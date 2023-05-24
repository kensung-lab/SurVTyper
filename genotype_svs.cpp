#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <chrono>

#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/tbx.h"
#include "utils.h"
#include "sam_utils.h"
#include "stat_tests.h"
#include "libs/cptl_stl.h"
#include "libs/ssw_cpp.h"
#include "libs/ssw.h"
#include "libs/seqan/store.h"
#include "libs/seqan/consensus.h"
#include "assemble.h"
#include "vcf_utils.h"

#include "genotype_deletions.h"


chr_seqs_map_t chr_seqs;
config_t config;
std::mutex mtx, consensus_log_mtx;
std::ofstream flog, consensus_flog;

std::string bam_fname, reference_fname;
bcf_hdr_t* vcf_header, * out_vcf_header;

std::vector<std::pair<uint64_t, uint32_t> > partial_sums;
std::vector<uint32_t> sampled_insert_sizes;
std::vector<uint32_t> depths, n_pairs_crossing;
std::vector<double> del_is_population;
// estimates how many pairs do we expect to see in a given inserted sequence
std::vector<std::vector<uint32_t> > dist_between_end_and_rnd;
std::vector<uint32_t> min_disc_pairs_for_ins_by_size, max_disc_pairs_for_ins_by_size;

std::mutex* mtx_contig;
std::vector<std::vector<std::string> > mate_seqs_buffers;


bool accept_aln_strict(StripedSmithWaterman::Alignment& aln) {
	return get_left_clip_size(aln) == 0 && get_right_clip_size(aln) == 0;
}


std::pair<std::string, double> calculate_small_dup_genotype(int alt_better_reads, int ref_better_reads) {
	if (alt_better_reads + ref_better_reads == 0) {
		return {"NO_GT", 0.0};
	}

	// in theory, we expect all alt_better for hom_alt calls, and a ratio of 1:1 alt_better:ref_better for het
	double ratio = (double) alt_better_reads/(alt_better_reads + ref_better_reads);
	if (ratio >= 0.75) {
		return {"HOM_ALT", ratio};
	} else if (ratio >= 0.25) {
		return {"HET", ratio};
	} else {
		return {"HOM_REF", ratio};
	}
}

int estimate_copy_number_large_dup(int alt_better_reads, int ref_better_reads) {
	if (alt_better_reads + ref_better_reads == 0) {
		return 0;
	}

	// in theory, we expect a ratio of 1:2 alt_better:ref_better for homalt. and 1:4 for het
	// so copy number is 2 + alt_better*4/ref_better
	return 2 + lround(alt_better_reads*4.0/ref_better_reads);
}

std::pair<std::string, double> calculate_small_ins_genotype(int alt_better_reads, int ref_better_reads) {
	if (alt_better_reads + ref_better_reads == 0) {
		return {"NO_GT", 0.0};
	}

	// in theory, we expect all alt_better for hom_alt calls, and a ratio of ~1:1 alt_better:ref_better for het
	double ratio = (double) alt_better_reads/(alt_better_reads + ref_better_reads);
	if (ratio >= 0.75) {
		return {"HOM_ALT", ratio};
	} else if (ratio >= 0.25) {
		return {"HET", ratio};
	} else {
		return {"HOM_REF", ratio};
	}
}

std::pair<std::string, double> calculate_large_ins_genotype(int alt_better_reads, int ref_better_reads) {
	if (alt_better_reads + ref_better_reads == 0) {
		return {"NO_GT", 0.0};
	}

	// in theory, we expect all alt_better for hom_alt calls, and a ratio of ~1:1 alt_better:ref_better for het
	double ratio = (double) alt_better_reads/(alt_better_reads + ref_better_reads);
	if (ratio >= 0.75) {
		return {"HOM_ALT", ratio};
	} else if (ratio >= 0.25) {
		return {"HET", ratio};
	} else {
		return {"HOM_REF", ratio};
	}
}



std::string build_consensus_contig(std::string contig, std::vector<std::string>& reads, std::vector<int>& reads_start, std::vector<int>& reads_end,
		std::string id = "") {
    seqan::FragmentStore<> store;
    int shift_by = 0;
    if (contig.empty()) {
    	shift_by = *std::min_element(reads_start.begin(), reads_start.end());
    }
    for (int i = 0; i < reads.size(); i++) {
        seqan::appendRead(store, reads[i]);
        seqan::appendAlignedRead(store, i, 0, reads_start[i]-shift_by, reads_end[i]-shift_by);
    }
    if (!contig.empty()) {
		seqan::appendRead(store, contig);
		seqan::appendAlignedRead(store, reads.size(), 0, 0, (int) contig.length());
    }

    seqan::ConsensusAlignmentOptions options;
    seqan::consensusAlignment(store, options);

    std::string gapped_consensus, consensus;
    getGappedConsensus(store, gapped_consensus, 0);
    consensus = gapped_consensus;
    consensus.erase(std::remove(consensus.begin(), consensus.end(), '-'), consensus.end());

	if (!id.empty()) {
		seqan::AlignedReadLayout layout;
		layoutAlignment(layout, store);
		consensus_log_mtx.lock();
		consensus_flog << id << std::endl;
		seqan::printAlignment(consensus_flog, layout, store, /*contigID=*/ 0, /*beginPos=*/ 0,
				/*endPos=*/ (int) gapped_consensus.length(), 0, (int) reads.size()+1);
		consensus_flog << std::endl;
		consensus_log_mtx.unlock();
	}

    return consensus;
}


void genotype_small_dup(std::string& contig_name, duplication_t* sv, open_samFile_t* bam_file,
                  StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner, stats_t stats) {

	std::stringstream log_ss;

	int dup_start = sv->start, dup_end = sv->end;
	int contig_len = chr_seqs.get_len(contig_name);

	int extend = config.read_len + 20;
	int dup_len = sv->len();

	// See comments for relative code in genotype_del
	dup_start++; dup_end++;

	int ref_start = std::max(0, dup_start-extend), ref_end = std::min(dup_end+extend, contig_len);
	int ref_len = ref_end - ref_start;
	char* ref_seq = new char[ref_len + 1];
	strncpy(ref_seq, chr_seqs.get_seq(contig_name)+ref_start, ref_len);
	ref_seq[ref_len] = 0;

	int left_extend = dup_start - ref_start, right_extend = ref_end - dup_end;

	std::vector<char*> alt_seqs;
	std::vector<int> alt_seqs_len;
	for (int copies = 1; copies*dup_len < config.read_len; copies++) {
		int alt_len = ref_len + copies*dup_len + 1;
		alt_seqs_len.push_back(alt_len);

		char* alt_seq = new char[ref_len + copies*dup_len + 1];
		int pos = 0;
		strncpy(alt_seq, chr_seqs.get_seq(contig_name)+ref_start, dup_end-ref_start);
		pos += dup_end - ref_start;
		for (int i = 0; i < copies; i++) {
			strncpy(alt_seq+pos, chr_seqs.get_seq(contig_name)+dup_start, dup_len);
			pos += dup_len;
		}
		strncpy(alt_seq+pos, chr_seqs.get_seq(contig_name)+dup_end, ref_end-dup_end);
		pos += ref_end - dup_end;
		alt_seq[pos] = 0;
		alt_seqs.push_back(alt_seq);
	}

	if (strchr(ref_seq, 'N') != NULL) return; // do not genotype if N region
	for (char* alt_seq : alt_seqs) {
		if (strchr(alt_seq, 'N') != NULL) return;
	}

	char region[100];
	sprintf(region, "%s:%d-%d", contig_name.c_str(), ref_start, ref_end);
	hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, region);
	bam1_t* read = bam_init1();

	std::vector<std::string> read_seqs, qnames;
	std::vector<bool> good_read;
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read)) continue;

		// do not consider reads that do not intersect the breakpoints
		if (get_unclipped_end(read) <= dup_start+config.min_clip_size) continue;
		if (get_unclipped_start(read) >= dup_start-config.min_clip_size && get_unclipped_end(read) <= dup_end+config.min_clip_size) continue;
		if (get_unclipped_start(read) >= dup_end-config.min_clip_size) continue;

		read_seqs.push_back(get_sequence(read));
		qnames.push_back(bam_get_qname(read));
		good_read.push_back(is_samechr(read) && !is_samestr(read) && !is_mate_unmapped(read));
	}

	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment alt_aln, ref_aln;
	// one vector per read, that has as many positins as alternative sequences.
	// Position [i][j] refers to read i, alternative sequence j
	std::vector<std::vector<bool> > alt_seqs_support(read_seqs.size()), same_support_seqs(read_seqs.size()); // whether read i supports alt seq j
	std::vector<std::vector<int> > alt_seqs_starts(read_seqs.size()); // start in alt seq j of read i
	std::vector<std::vector<int> > alt_seqs_score_diff(read_seqs.size());
	std::vector<bool> ref_accepted;
	for (int i = 0; i < read_seqs.size(); i++) {
		std::string seq = std::string(config.clip_penalty, 'N') + read_seqs[i] + std::string(config.clip_penalty, 'N');

		// align to REF (two breakpoints)
		aligner.Align(seq.c_str(), ref_seq, ref_len, filter, &ref_aln, 0);

		std::vector<int> alt_scores(alt_seqs.size());
		int alt_best_score = 0;
		for (int j = 0; j < alt_seqs.size(); j++) {
			aligner.Align(seq.c_str(), alt_seqs[j], alt_seqs_len[j], filter, &alt_aln, 0);
			if (accept_aln_strict(alt_aln)) {
				alt_scores[j] = alt_aln.sw_score;
				if (alt_best_score < alt_aln.sw_score) {
					alt_best_score = alt_scores[j];
				}
			}
			alt_seqs_starts[i].push_back(alt_aln.ref_begin + config.clip_penalty);
			alt_seqs_score_diff[i].push_back(alt_aln.sw_score-ref_aln.sw_score);
		}

		// initially this read does not support any alt seq
		alt_seqs_support[i].resize(alt_seqs.size());
		same_support_seqs[i].resize(alt_seqs.size());

		// TODO: !accept_aln_strict(ref_aln) - transfer this logic to deletions
		ref_accepted.push_back(accept_aln_strict(ref_aln));

		if (alt_best_score > ref_aln.sw_score) {
			for (int j = 0; j < alt_seqs.size(); j++) { // this read supports all alt sequences with maximal score
				alt_seqs_support[i][j] = (alt_best_score == alt_scores[j]);
			}
		} else if (alt_best_score < ref_aln.sw_score && accept_aln_strict(ref_aln)) {
			sv->ref_better++;
		} else if (alt_best_score == ref_aln.sw_score && accept_aln_strict(ref_aln) && accept_aln_strict(alt_aln)) {
			sv->same++;
			for (int j = 0; j < alt_seqs.size(); j++) { // this reads supports all alt sequences with maximal score
				same_support_seqs[i][j] = (alt_best_score == alt_scores[j]);
			}
		}
	}

	std::vector<int> alt_seqs_support_count(alt_seqs.size());
	for (int i = 0; i < read_seqs.size(); i++) {
		for (int j = 0; j < alt_seqs.size(); j++) {
			alt_seqs_support_count[j] += alt_seqs_support[i][j];
		}
	}

	// select best alt sequence
	int chosen_alt_seq_index = 0;
	for (int i = 0; i < alt_seqs.size(); i++) {
		if (alt_seqs_support_count[chosen_alt_seq_index] < alt_seqs_support_count[i]) chosen_alt_seq_index = i;
	}
	sv->alt_better = alt_seqs_support_count[chosen_alt_seq_index];
	char* best_alt_seq = alt_seqs[chosen_alt_seq_index];
	int best_alt_seq_len = strlen(best_alt_seq);

	// gather reads supporting best alt seq
	std::vector<std::string> alt_better_seqs, alt_better_qnames;
	std::vector<int> alt_better_score_diff;
	std::vector<bool> alt_better_is_strong;

	// we divide the reads into two sets (not necessarily disjoint) of reads that span the left bp and those that span the right bp
	// we create a consensus for each
	std::vector<std::string> lbp_alt_better_seqs, rbp_alt_better_seqs;
	std::vector<int> lbp_alt_better_starts, lbp_alt_better_ends, rbp_alt_better_starts, rbp_alt_better_ends;
	std::vector<bool> alt_better_supp_lbp, alt_better_supp_rbp;
	int compl_inside_dup = 0;
	for (int i = 0; i < read_seqs.size(); i++) {
		if (alt_seqs_support[i][chosen_alt_seq_index]) {
			alt_better_seqs.push_back(read_seqs[i]);
			alt_better_qnames.push_back(qnames[i]);
			int score_diff = alt_seqs_score_diff[i][chosen_alt_seq_index];
			alt_better_score_diff.push_back(score_diff);
			alt_better_is_strong.push_back(score_diff >= config.min_score_diff || !ref_accepted[i]);

			// we use a read for consensus generation only if it is not entirely in the duplication and
			// its mate is outside the duplication
			int aln_start = alt_seqs_starts[i][chosen_alt_seq_index];
			int aln_end = aln_start + config.read_len;
			if (aln_start < left_extend-config.min_clip_size) { // will build left bp consensus
				lbp_alt_better_seqs.push_back(read_seqs[i]);
				lbp_alt_better_starts.push_back(aln_start);
				lbp_alt_better_ends.push_back(aln_end);
				alt_better_supp_lbp.push_back(true);
			} else alt_better_supp_lbp.push_back(false);
			if (aln_end > best_alt_seq_len-right_extend+config.min_clip_size) { // supports right bp
				rbp_alt_better_seqs.push_back(read_seqs[i]);
				rbp_alt_better_starts.push_back(aln_start);
				rbp_alt_better_ends.push_back(aln_end);
				alt_better_supp_rbp.push_back(true);
			} else alt_better_supp_rbp.push_back(false);

			if (left_extend < aln_start && aln_end < best_alt_seq_len-right_extend) {
				compl_inside_dup++;
			}
		}
	}
	sv->compl_inside_dup = compl_inside_dup;

	bool not_enough_lbp_reads = lbp_alt_better_seqs.size() < 3, not_enough_rbp_reads = rbp_alt_better_seqs.size() < 3;

	for (int i = 0; i < read_seqs.size(); i++) {
		if (same_support_seqs[i][chosen_alt_seq_index]) {
			int aln_start = alt_seqs_starts[i][chosen_alt_seq_index];
			int aln_end = aln_start + config.read_len;
			if (not_enough_lbp_reads && aln_start < dup_start-ref_start-config.min_clip_size) { // supports left bp
				lbp_alt_better_seqs.push_back(read_seqs[i]);
				lbp_alt_better_starts.push_back(aln_start);
				lbp_alt_better_ends.push_back(aln_end);
			}
			if (not_enough_rbp_reads && aln_end > dup_end-ref_start+config.min_clip_size) { // supports right bp
				rbp_alt_better_seqs.push_back(read_seqs[i]);
				rbp_alt_better_starts.push_back(aln_start);
				rbp_alt_better_ends.push_back(aln_end);
			}
		}
	}

	if (alt_better_seqs.size() > 2*stats.max_depth || lbp_alt_better_seqs.size() > stats.max_depth || rbp_alt_better_seqs.size() > stats.max_depth) {
		auto gt = calculate_del_genotype(0, 0);
		sv->before_read_qc_gt = sv->called_gt = gt.first;

		for (char* alt_seq : alt_seqs) {
			delete[] alt_seq;
		}
		delete[] ref_seq;

		bam_destroy1(read);
		sam_itr_destroy(iter);
		return;
	}

	// calculate GT
	auto gt_before_qc = calculate_small_dup_genotype(sv->alt_better, sv->ref_better);
	sv->before_read_qc_gt = gt_before_qc.first;

	std::string lbp_alt_consensus, rbp_alt_consensus;
	log_ss << "ALT READS: " << lbp_alt_better_seqs.size() << "," << rbp_alt_better_seqs.size() << std::endl;
	if (gt_is_positive(sv->before_read_qc_gt) && lbp_alt_better_seqs.size() >= 3 && rbp_alt_better_seqs.size() >= 3) {
		lbp_alt_consensus = build_consensus_contig("", lbp_alt_better_seqs, lbp_alt_better_starts, lbp_alt_better_ends);
		rbp_alt_consensus = build_consensus_contig("", rbp_alt_better_seqs, rbp_alt_better_starts, rbp_alt_better_ends);
		log_ss << lbp_alt_consensus << " " << rbp_alt_consensus << std::endl;

		suffix_prefix_aln_t spa = aln_suffix_prefix(lbp_alt_consensus, rbp_alt_consensus, 1, -4, config.max_seq_error, config.min_clip_size);
		log_ss << lbp_alt_consensus.size() << " " << rbp_alt_consensus.size() << " " << spa.overlap << " " << spa.score << " " << spa.mismatches << std::endl;
		if (spa.overlap >= config.min_clip_size) {
			std::string junction_seq = lbp_alt_consensus.substr(0, lbp_alt_consensus.size()-spa.overlap) + rbp_alt_consensus;
			log_ss << "JUNCTION SEQ: " << junction_seq << std::endl;
			StripedSmithWaterman::Alignment ref_jun_aln;
			aligner.Align(junction_seq.c_str(), ref_seq, ref_len, filter, &ref_jun_aln, 0);
			log_ss << ref_len << " " << junction_seq.length() << std::endl;
			log_ss << ref_jun_aln.query_begin << " " << ref_jun_aln.query_end << " " << ref_jun_aln.cigar_string << std::endl;

			std::string left_anchor = junction_seq.substr(0, 50);
			std::string right_anchor = junction_seq.substr(junction_seq.length()-50);
			StripedSmithWaterman::Alignment la_aln, ra_aln;
			aligner.Align(left_anchor.c_str(), ref_seq, ref_len, filter, &la_aln, 0);
			aligner.Align(right_anchor.c_str(), ref_seq, ref_len, filter, &ra_aln, 0);
			sv->predicted_insertion_len = junction_seq.length() - (ra_aln.ref_end-la_aln.ref_begin+get_left_clip_size(la_aln)+get_right_clip_size(ra_aln));

			std::string padded_junction_seq = std::string(config.clip_penalty, 'N') + junction_seq + std::string(config.clip_penalty, 'N');

			char ref_cstr[100000];
			for (int i = 0; i < ref_len; i++) {
				ref_cstr[i] = toupper(ref_seq[i]);
			}
			ref_cstr[ref_len] = '\0';
			int* prefix_scores = smith_waterman_gotoh(ref_cstr, ref_len, padded_junction_seq.c_str(), padded_junction_seq.length(), 1, -4, -6, -1);

			for (int i = 0; i < ref_len/2; i++) {
				std::swap(ref_cstr[i], ref_cstr[ref_len-i-1]);
			}
			std::string padded_junction_seq_rev = std::string(padded_junction_seq.rbegin(), padded_junction_seq.rend());
			int* suffix_scores = smith_waterman_gotoh(ref_cstr, ref_len, padded_junction_seq_rev.c_str(), padded_junction_seq_rev.length(), 1, -4, -6, -1);

			int max_score = 0, split_i = 0;
			for (int i = 0; i < junction_seq.length(); i++) {

				int prefix_score = prefix_scores[config.clip_penalty+i-1], suffix_score = suffix_scores[junction_seq.length()-i-1+config.clip_penalty];
				if (prefix_score + suffix_score > max_score) {
					max_score = prefix_score + suffix_score;
					split_i = i;
				}
			}

			// ref_cstr was previously reverted, reset
			for (int i = 0; i < ref_len; i++) {
				ref_cstr[i] = toupper(ref_seq[i]);
			}
			ref_cstr[ref_len] = '\0';

			StripedSmithWaterman::Alignment lh_aln, rh_aln;
			char* padded_junction_seq_cstr = new char[padded_junction_seq.length()+1];
			strcpy(padded_junction_seq_cstr, padded_junction_seq.c_str());
			char terminator = '\0';
			std::swap(terminator, padded_junction_seq_cstr[config.clip_penalty+split_i]);
			aligner.Align(padded_junction_seq_cstr, ref_cstr, ref_len, filter, &lh_aln, 0);
			std::swap(terminator, padded_junction_seq_cstr[config.clip_penalty+split_i]);

			aligner.Align(padded_junction_seq_cstr+config.clip_penalty+split_i, ref_cstr, ref_len, filter, &rh_aln, 0);

			log_ss << lh_aln.cigar_string << " " << lh_aln.ref_begin << " " << lh_aln.ref_end << std::endl;
			log_ss << rh_aln.cigar_string << " " << rh_aln.ref_begin << " " << rh_aln.ref_end << std::endl;

			if (get_left_clip_size(lh_aln) > 0 || get_right_clip_size(rh_aln) > 0) sv->lowq_alt_allele = true;

			sv->remapped_start = ref_start + rh_aln.ref_begin-1;
			sv->remapped_end = ref_start + lh_aln.ref_end;
			sv->was_remapped = true;
		}

		int i = 0;
		StripedSmithWaterman::Alignment consensus_alt_aln;
		for (int i = 0; i < alt_better_seqs.size(); i++) {
			bool good_qual = false, strong = false;
//			if (ref_better_qnames.count(alt_better_qnames[i])) continue;
			std::string seq = alt_better_seqs[i];

			std::string padded_seq = std::string(config.clip_penalty, 'N') + seq + std::string(config.clip_penalty, 'N');
			if (alt_better_supp_lbp[i]) {
				harsh_aligner.Align(padded_seq.c_str(), lbp_alt_consensus.c_str(), lbp_alt_consensus.length(), filter, &consensus_alt_aln, 0);
				if (get_left_clip_size(consensus_alt_aln) <= config.clip_penalty && get_right_clip_size(consensus_alt_aln) <= config.clip_penalty) {
					// see deletions code for explanation on computation of n of mismatches
					int mismatches = consensus_alt_aln.mismatches - 2*config.clip_penalty + get_left_clip_size(consensus_alt_aln) +
							get_right_clip_size(consensus_alt_aln);

					if (mismatches <= seq.length()*config.max_seq_error) {
						good_qual = true;
						if (alt_better_is_strong[i]) strong = true;
					}
				}
			}
			if (alt_better_supp_rbp[i]) {
				harsh_aligner.Align(padded_seq.c_str(), rbp_alt_consensus.c_str(), rbp_alt_consensus.length(), filter, &consensus_alt_aln, 0);
				if (get_left_clip_size(consensus_alt_aln) <= config.clip_penalty && get_right_clip_size(consensus_alt_aln) <= config.clip_penalty) {
					// see deletions code for explanation on computation of n of mismatches
					int mismatches = consensus_alt_aln.mismatches - 2*config.clip_penalty + get_left_clip_size(consensus_alt_aln) +
							get_right_clip_size(consensus_alt_aln);

					if (mismatches <= seq.length()*config.max_seq_error) {
						good_qual = true;
						if (alt_better_is_strong[i]) strong = true;
					}
				}
			}

			if (good_qual) sv->alt_better_goodqual++;
			if (strong) sv->alt_better_strong++;
		}
	}

	// recalculate GT after excluding poor quality reads
	auto gt = calculate_small_dup_genotype(sv->alt_better_goodqual, sv->ref_better);
	sv->called_gt = gt.first;
	sv->ratio = gt.second;
	sv->copy_number = 0; // 2 + chosen_alt_seq_index+1;

	log_ss << sv->id << " " << contig_name << " " << sv->start << " " << sv->end << std::endl;
	log_ss << "REF: " << ref_seq << std::endl;
	for (int i = 0; i < alt_seqs.size(); i++) {
		log_ss << "ALT (" << i+1 << " extra copies) " << alt_seqs[i] << " " << alt_seqs_support_count[i] << "/" << sv->alt_better << std::endl;
	}
	log_ss << "CHOSEN ALT: " << alt_seqs_support_count[chosen_alt_seq_index] << "/" << sv->alt_better << std::endl;
	log_ss << "CONSENSUS READS: " << alt_better_seqs.size() << std::endl;
	log_ss << sv->alt_better << " " << sv->alt_better_goodqual << " " << sv->alt_better_strong << " " << sv->ref_better << " " << sv->same << std::endl;
	log_ss << sv->before_read_qc_gt << " " << sv->called_gt << std::endl;

	for (char* alt_seq : alt_seqs) {
		delete[] alt_seq;
	}
	delete[] ref_seq;

	bam_destroy1(read);
	sam_itr_destroy(iter);

//	mtx.lock();
//	flog << log_ss.str() << std::endl;
//	mtx.unlock();
}

void genotype_large_dup(std::string& contig_name, duplication_t* sv, open_samFile_t* bam_file,
                  StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner, stats_t stats) {

	std::stringstream log_ss;

	int dup_start = sv->start, dup_end = sv->end;
	int contig_len = chr_seqs.get_len(contig_name);

	int extend = config.read_len + 20;

	// See comments for relative code in genotype_del
	dup_start++; dup_end++;

	// all ranges will be start-inclusive and end-exclusive, i.e. [a,b)

	int ref_bp1_start = std::max(0, dup_start-extend), ref_bp1_end = std::min(dup_start+extend, contig_len);
	int ref_bp1_len = ref_bp1_end - ref_bp1_start;
	char* ref_bp1_seq = new char[ref_bp1_len + 1];
	strncpy(ref_bp1_seq, chr_seqs.get_seq(contig_name)+ref_bp1_start, ref_bp1_len);
	ref_bp1_seq[ref_bp1_len] = 0;

	int ref_bp2_start = std::max(0, dup_end-extend), ref_bp2_end = std::min(dup_end+extend, contig_len);
	int ref_bp2_len = ref_bp2_end - ref_bp2_start;
	char* ref_bp2_seq = new char[ref_bp2_len + 1];
	strncpy(ref_bp2_seq, chr_seqs.get_seq(contig_name)+ref_bp2_start, ref_bp2_len);
	ref_bp2_seq[ref_bp2_len] = 0;

	// build alt allele
	int alt_lh_len = dup_end - ref_bp2_start;
	int alt_rh_len = ref_bp1_end - dup_start;
	int alt_len = alt_lh_len + alt_rh_len;
	char* alt_seq = new char[alt_len + 1];
	strncpy(alt_seq, chr_seqs.get_seq(contig_name)+ref_bp2_start, alt_lh_len);
	strncpy(alt_seq+alt_lh_len, chr_seqs.get_seq(contig_name)+dup_start, alt_rh_len);
	alt_seq[alt_len] = 0;

	if (strchr(alt_seq, 'N') != NULL || strchr(ref_bp1_seq, 'N') != NULL || strchr(ref_bp2_seq, 'N') != NULL) { // do not genotype if N region
		return;
	}

	char l_region[100], r_region[100];
	sprintf(l_region, "%s:%d-%d", contig_name.c_str(), ref_bp1_start, ref_bp1_end);
	sprintf(r_region, "%s:%d-%d", contig_name.c_str(), ref_bp2_start, ref_bp2_end);
	char* regions[] = {l_region, r_region};

	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);
	bam1_t* read = bam_init1();

	std::vector<std::string> read_seqs, qnames;
	std::vector<bool> good_read;
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read)) continue;

		// do not consider reads that do not intersect the breakpoints
		if (get_unclipped_end(read) <= dup_start) continue;
		if (dup_start <= get_unclipped_start(read) && get_unclipped_end(read) <= dup_end) continue;
		if (get_unclipped_start(read) >= dup_end) continue;

		read_seqs.push_back(get_sequence(read));
		qnames.push_back(bam_get_qname(read));
		good_read.push_back(is_samechr(read) && !is_samestr(read) && !is_mate_unmapped(read));
	}

	char region[100];
	sprintf(region, "%s:%d-%d", contig_name.c_str(), dup_end-stats.max_is, dup_end);
	iter = sam_itr_querys(bam_file->idx, bam_file->header, region);
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read) || !is_samechr(read) || is_samestr(read) || bam_is_rev(read)) continue;

		hts_pos_t mate_end  = get_mate_end(read);

		if (read->core.mpos < dup_start || bam_endpos(read) > dup_end) continue;

		if (read->core.isize < stats.min_is && mate_end > dup_start && mate_end < dup_start+config.max_is) {
			log_ss << "OW " << bam_get_qname(read) << std::endl;
			sv->ow_pairs++;
		}
	}

	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment alt_aln, ref1_aln, ref2_aln;
	std::vector<std::string> alt_better_seqs, alt_better_qnames;
	std::vector<int> alt_better_starts, alt_better_ends, alt_better_scores;
	std::vector<bool> alt_is_strong;
	for (int i = 0; i < read_seqs.size(); i++) {
		std::string seq = std::string(config.clip_penalty, 'N') + read_seqs[i] + std::string(config.clip_penalty, 'N');

		// align to ALT
		aligner.Align(seq.c_str(), alt_seq, alt_len, filter, &alt_aln, 0);

		// align to REF (two breakpoints)
		aligner.Align(seq.c_str(), ref_bp1_seq, ref_bp1_len, filter, &ref1_aln, 0);
		aligner.Align(seq.c_str(), ref_bp2_seq, ref_bp2_len, filter, &ref2_aln, 0);

		StripedSmithWaterman::Alignment& ref_aln = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln : ref2_aln;
		if (alt_aln.sw_score > ref_aln.sw_score && accept_aln_strict(alt_aln)) {
			sv->alt_better++;
			alt_is_strong.push_back(/*alt_aln.sw_score-ref_aln.sw_score >= config.min_score_diff || */
					good_read[i] && !accept_aln_strict(ref_aln));
			alt_better_qnames.push_back(qnames[i]);
			alt_better_seqs.push_back(read_seqs[i]);
			alt_better_starts.push_back(alt_aln.ref_begin + config.clip_penalty);
			alt_better_ends.push_back(alt_aln.ref_begin + config.clip_penalty + config.read_len);
			alt_better_scores.push_back(alt_aln.sw_score);
		}
		else if (alt_aln.sw_score < ref_aln.sw_score && accept_aln_strict(ref_aln)) {
			sv->ref_better++;
		}
		else if (alt_aln.sw_score == ref_aln.sw_score && accept_aln_strict(alt_aln) && accept_aln_strict(ref_aln)) sv->same++;

		if (sv->ref_better > 4*stats.max_depth) {
			sv->gt_skipped = true;
			return;
		}
	}

	sv->before_qc_copy_number = estimate_copy_number_large_dup(sv->alt_better, sv->ref_better);
	if (sv->before_qc_copy_number == 2) {
		sv->before_read_qc_gt = "HOM_REF";
	} else if (sv->before_qc_copy_number >= 3) {
		sv->before_read_qc_gt = "SV_EXISTS";
	}
	sv->called_gt = sv->before_read_qc_gt;

	std::string consensus_alt_seq;
	if (gt_is_positive(sv->called_gt)) {
		std::vector<std::pair<int, int> > score_w_index;
		for (int  i = 0; i < alt_better_seqs.size(); i++) {
			score_w_index.push_back({alt_better_scores[i], i});
		}
		std::random_shuffle(score_w_index.begin(), score_w_index.end());
		std::sort(score_w_index.begin(), score_w_index.end(), [](const std::pair<int,int>& p1, const std::pair<int,int>& p2) {
			return p1.first > p2.first;
		});

		std::vector<std::string> subsampled_alt_better_seqs;
		std::vector<int> subsampled_alt_better_starts, subsampled_alt_better_ends;
		for (int i = 0; i < score_w_index.size() && i <= stats.max_depth; i++) {
			subsampled_alt_better_seqs.push_back(alt_better_seqs[score_w_index[i].second]);
			subsampled_alt_better_starts.push_back(alt_better_starts[score_w_index[i].second]);
			subsampled_alt_better_ends.push_back(alt_better_ends[score_w_index[i].second]);
		}
		consensus_alt_seq = build_consensus_contig(alt_seq, subsampled_alt_better_seqs, subsampled_alt_better_starts, subsampled_alt_better_ends);

		StripedSmithWaterman::Alignment consensus_alt_aln;
		for (int i = 0; i < alt_better_seqs.size(); i++) {
			std::string seq = alt_better_seqs[i];
			std::string padded_seq = std::string(config.clip_penalty, 'N') + seq + std::string(config.clip_penalty, 'N');
			harsh_aligner.Align(padded_seq.c_str(), consensus_alt_seq.c_str(), consensus_alt_seq.length(), filter, &consensus_alt_aln, 0);

			bool accepted = false, strong = false;
			if (get_left_clip_size(consensus_alt_aln) <= config.clip_penalty && get_right_clip_size(consensus_alt_aln) <= config.clip_penalty) {
				int mismatches = consensus_alt_aln.mismatches - 2*config.clip_penalty + get_left_clip_size(consensus_alt_aln) +
						get_right_clip_size(consensus_alt_aln);

				if (mismatches <= seq.length()*config.max_seq_error) {
					sv->alt_better_goodqual++;
					accepted = true;
					if (alt_is_strong[i]) {
						sv->alt_better_strong++;
						strong = true;
					}
				}
			}

			log_ss << alt_better_qnames[i] << " " << consensus_alt_aln.cigar_string << " " << consensus_alt_aln.ref_begin << " " <<
					consensus_alt_aln.ref_end << " " << accepted << " " << strong << std::endl;
		}
	}

	sv->copy_number = estimate_copy_number_large_dup(sv->alt_better_goodqual, sv->ref_better);
	if (sv->copy_number == 2) {
		sv->called_gt = "HOM_REF";
	} else if (sv->copy_number >= 3) {
		sv->called_gt = "SV_EXISTS";
	}

	log_ss << sv->id << " " << contig_name << " " << sv->start << " " << sv->end << std::endl;
	log_ss << "REF1: " << ref_bp1_seq << std::endl;
	log_ss << "REF2: " << ref_bp2_seq << std::endl;
	log_ss << "ALT: " << alt_seq << std::endl;
	log_ss << "CONSENSUS ALT: " << consensus_alt_seq << std::endl;
	log_ss << sv->alt_better << " " << sv->alt_better_goodqual << " " << sv->ref_better << " " << sv->same << std::endl;

//	mtx.lock();
//	flog << log_ss.str() << std::endl;
//	mtx.unlock();

	delete[] alt_seq;
	delete[] ref_bp1_seq;
	delete[] ref_bp2_seq;

	bam_destroy1(read);
	sam_itr_destroy(iter);
}

duplication_t* vcf_record_to_duplication(bcf1_t* vcf_dup) {
	bcf_unpack(vcf_dup, BCF_UN_INFO);
    duplication_t* dup = new duplication_t(vcf_dup->d.id, vcf_dup->pos, get_sv_end(vcf_dup, vcf_header));
    dup->bcf_entry = vcf_dup;
    return dup;
}

void genotype_dups(int id, std::string contig_name, std::vector<bcf1_t*> vcf_dups, stats_t stats) {

	mtx.lock();
	std::cerr << "Genotyping duplications in " << contig_name << std::endl;
	mtx.unlock();

	StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, true);
	StripedSmithWaterman::Aligner harsh_aligner(1, 4, 100, 1, true);
	open_samFile_t* bam_file = open_samFile(bam_fname);
	hts_set_fai_filename(bam_file->file, reference_fname.c_str());

	std::vector<duplication_t*> dups;
	for (bcf1_t* sv : vcf_dups) {
		dups.push_back(vcf_record_to_duplication(sv));
	}

	std::chrono::steady_clock::time_point begin, end;
	begin = std::chrono::steady_clock::now();
	depth_filter_dup(contig_name, chr_seqs.get_len(contig_name), dups, bam_file);
	end = std::chrono::steady_clock::now();
	std::cout << "Dups depths for " << contig_name << " calculated in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " [s]" << std::endl;

	begin = std::chrono::steady_clock::now();
	for (duplication_t* dup : dups) {
		if (dup->left_flanking_cov > stats.max_depth || dup->right_flanking_cov > stats.max_depth) { // skip very high depth regions, because they are very slow and will be filtered anyway
			dup->gt_skipped = true;
			continue;
		}
		if (dup->is_short(config)) {
			genotype_small_dup(contig_name, dup, bam_file, aligner, harsh_aligner, stats);
		} else {
			genotype_large_dup(contig_name, dup, bam_file, aligner, harsh_aligner, stats);
		}
	}
	end = std::chrono::steady_clock::now();
	std::cout << "Genotyped " << contig_name << " in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " [s]" << std::endl;

	mtx.lock();
	for (duplication_t* sv : dups) {
		std::string filter;
		if (!sv->gt_skipped && !sv->alt_better && !sv->ref_better) filter += "NO_USEFUL_READS;";
		int dup_len = sv->len();
		if (sv->left_flanking_cov > stats.max_depth || sv->right_flanking_cov > stats.max_depth ||
			sv->left_flanking_cov < stats.min_depth || sv->right_flanking_cov < stats.min_depth) filter += "ANOMALOUS_FLANKING_DEPTH;";
		if (sv->indel_left_cov < stats.min_depth || sv->indel_right_cov < stats.min_depth) filter += "ANOMALOUS_DUP_DEPTH;";
		if (gt_is_positive(sv->called_gt) && sv->end-sv->start >= config.min_size_for_depth_filtering &&
		(sv->mleft_flanking_cov*1.26>=sv->mindel_left_cov || sv->mindel_right_cov<=sv->mright_flanking_cov*1.26)) {
			filter += "DEPTH_FILTER;";
		}
		if (gt_is_positive(sv->called_gt)) {
			if (sv->alt_better_strong < 3) filter += "LOW_ALT_STRONG_READS;";
			if (sv->remapped_end-sv->remapped_start < config.min_sv_size) filter += "TOO_SHORT_AFTER_REMAPPING;";
			if (sv->lowq_alt_allele) filter += "POOR_JUNCTION_SEQ;";
			else if (!sv->was_remapped && sv->is_short(config) && sv->compl_inside_dup < 3) filter += "POOR_JUNCTION_SEQ;"; // this removes quite a few hard TPs
		}
		if (filter.empty()) {
			if (sv->called_gt != "NO_GT") filter = "PASS";
			else filter = "SKIPPED";
		}

		bcf_subset(out_vcf_header, sv->bcf_entry, 0, {});
		int gt_data[2];
		strgt2datagt(sv->called_gt, gt_data);
		bcf_update_genotypes(out_vcf_header, sv->bcf_entry, gt_data, 2);

		const char* ft_val[1];
		ft_val[0] = filter.c_str();
		int remapped_start_1based = sv->remapped_start+1, remapped_end_1based = sv->remapped_end+1;
		bcf_update_format_string(out_vcf_header, sv->bcf_entry, "FT", ft_val, 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "AR", &(sv->alt_better_goodqual), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "SR", &(sv->alt_better_strong), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "RR", &(sv->ref_better), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "CN", &(sv->copy_number), 1);
//		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "IL", &(sv->predicted_insertion_len), 1);
        bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "RS", &remapped_start_1based, 1);
        bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "RE", &remapped_end_1based, 1);

        int was_remapped = sv->was_remapped;
        bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "WR", &was_remapped, 1);
        bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "DFL", &(sv->left_flanking_cov), 1);
        bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "DDL", &(sv->indel_left_cov), 1);
        bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "DDR", &(sv->indel_right_cov), 1);
        bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "DFR", &(sv->right_flanking_cov), 1);

        bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "MDFL", &(sv->mleft_flanking_cov), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "MDDL", &(sv->mindel_left_cov), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "MDDR", &(sv->mindel_right_cov), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "MDFR", &(sv->mright_flanking_cov), 1);

        bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "CID", &(sv->compl_inside_dup), 1);
        bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "OW", &(sv->ow_pairs), 1);
	}
	mtx.unlock();

	close_samFile(bam_file);
}

void genotype_small_ins(std::string& contig_name, insertion_t* sv, open_samFile_t* bam_file, StripedSmithWaterman::Aligner& aligner,
		StripedSmithWaterman::Aligner& harsh_aligner, std::unordered_map<std::string, std::string>& mateseqs_map, stats_t stats) {

	std::stringstream log_ss;

	int ins_start = sv->start, ins_end = sv->end;
	int contig_len = chr_seqs.get_len(contig_name);

	int extend = config.read_len + 20;

	// build alt allele
	/*
	 * POS in VCF is the base BEFORE the insertion
	 * END seems to be the base BEFORE the reference resumes - i.e., for a "clean" insertion (no deletion),POS == END, otherwise the last base deleted
	 * As usual, in order to make intervals [ ), we increase the coordinates by 1
	 */
	ins_start++; ins_end++;

	char* contig_seq = chr_seqs.get_seq(contig_name);

	int alt_start = std::max(0, ins_start-extend);
	int alt_end = std::min(ins_end+extend, contig_len);
	int alt_lf_len = ins_start-alt_start, alt_rf_len = alt_end-ins_end;
	int alt_len = alt_lf_len + sv->ins_seq.length() + alt_rf_len;
	char* alt_seq = new char[alt_len + 1];
	strncpy(alt_seq, chr_seqs.get_seq(contig_name)+alt_start, alt_lf_len);
	strncpy(alt_seq+alt_lf_len, sv->ins_seq.c_str(), sv->ins_seq.length());
	strncpy(alt_seq+alt_lf_len+sv->ins_seq.length(), contig_seq+ins_end, alt_rf_len);
    alt_seq[alt_len] = '\0';

    int ref_bp1_start = alt_start, ref_bp1_end = std::min(ins_start+extend, contig_len);
    int ref_bp1_len = ref_bp1_end - ref_bp1_start;
    int ref_bp2_start = std::max(0, ins_end-extend), ref_bp2_end = alt_end;
    int ref_bp2_len = ref_bp2_end - ref_bp2_start;

    int ref_len = ins_start-ref_bp1_start + ref_bp2_end-ins_end;
    char* ref_seq = new char[ref_len + 1];
	strncpy(ref_seq, contig_seq+ref_bp1_start, ins_start-ref_bp1_start);
	strncpy(ref_seq+(ins_start-ref_bp1_start), contig_seq+ins_end, ref_bp2_end-ins_end);
	ref_seq[ref_len] = '\0';

	char l_region[100], r_region[100];
	sprintf(l_region, "%s:%d-%d", contig_name.c_str(), alt_start, ref_bp1_end);
	sprintf(r_region, "%s:%d-%d", contig_name.c_str(), ref_bp2_start, alt_end);
	char* regions[] = {l_region, r_region};

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);
    bam1_t* read = bam_init1();

    std::vector<std::string> read_seqs, qnames;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read)) continue;

		if (!bam_is_rev(read) && read->core.pos < ins_start && is_dc_pair(read)) {
			std::string mate_qname = bam_get_qname(read);
			if (read->core.tid == read->core.mtid) {
				mate_qname += (read->core.flag & BAM_FREAD1 ? "_2" : "_1");
			}
			std::string mate_seq = mateseqs_map[mate_qname];
			rc(mate_seq);

			read_seqs.push_back(mate_seq);
			qnames.push_back(mate_qname + "_MATE");
		} else if (bam_is_rev(read) && bam_endpos(read) > ins_end && is_dc_pair(read)) {
			std::string mate_qname = bam_get_qname(read);
			if (read->core.tid == read->core.mtid) {
				mate_qname += (read->core.flag & BAM_FREAD1 ? "_2" : "_1");
			}
			std::string mate_seq = mateseqs_map[mate_qname];

			read_seqs.push_back(mate_seq);
			qnames.push_back(mate_qname + "_MATE");
		}

		// do not consider reads that do not intersect the insertion
		if (get_unclipped_end(read) <= ins_start) continue;
		if (get_unclipped_start(read) >= ins_end) continue;

		read_seqs.push_back(get_sequence(read));
		qnames.push_back(bam_get_qname(read));
    }

    StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment alt_aln, ref1_aln, ref2_aln;
	std::vector<std::string> alt_better_seqs;
	std::vector<bool> alt_better_seqs_lc, alt_better_seqs_rc; // TODO: populate with actual values
	std::vector<bool> alt_strong;
	log_ss << sv->id << std::endl;
	int alt_better = 0, ref_better = 0, same = 0;
	int i = 0;
	for (std::string& seq : read_seqs) {
		std::string padded_seq = std::string(config.clip_penalty, 'N') + seq + std::string(config.clip_penalty, 'N');

		// align to ALT
		aligner.Align(padded_seq.c_str(), alt_seq, alt_len, filter, &alt_aln, 0);

		// align to REF (two breakpoints)
		aligner.Align(padded_seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_len, filter, &ref1_aln, 0);
		aligner.Align(padded_seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_len, filter, &ref2_aln, 0);

		StripedSmithWaterman::Alignment& ref_aln = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln : ref2_aln;

		if (alt_aln.sw_score > ref_aln.sw_score /*&& accept_aln_strict(alt_aln)*/) {
			alt_better++;
			alt_better_seqs.push_back(seq);
			alt_better_seqs_lc.push_back(false);
			alt_better_seqs_rc.push_back(false);
			alt_strong.push_back(alt_aln.sw_score-ref_aln.sw_score >= config.min_score_diff || !accept_aln_strict(ref_aln));
		} else if (alt_aln.sw_score < ref_aln.sw_score && accept_aln_strict(ref_aln)) {
			ref_better++;
		} else if (alt_aln.sw_score == ref_aln.sw_score && accept_aln_strict(alt_aln) && accept_aln_strict(ref_aln)) same++;

		log_ss << seq << " " << alt_aln.cigar_string << " " << ref_aln.cigar_string << " " << alt_aln.sw_score << " " << ref_aln.sw_score << std::endl;
	}
	sv->alt_better = alt_better; sv->ref_better = ref_better;

	auto gt1 = calculate_small_ins_genotype(alt_better, ref_better);
	sv->before_read_qc_gt = gt1.first;

	std::string alt_consensus;
	StripedSmithWaterman::Alignment consensus_alt_aln, diff_consensus_alt_aln;
	if (alt_better > 0 && alt_better <= 4*stats.max_depth && ref_better <= 2*stats.max_depth) {
		std::vector<StripedSmithWaterman::Alignment> consensus_contigs_alns;
		std::string consensus_log;
		alt_consensus = generate_consensus(alt_better_seqs, alt_better_seqs_lc, alt_better_seqs_rc,
				aligner, harsh_aligner, consensus_contigs_alns, consensus_log);

		log_ss << "ALT CONSENSUS: " << alt_consensus << std::endl;

		int i = 0;
		int alt_better = 0, diff_alt_better = 0;
		for (std::string& seq : alt_better_seqs) {
			std::string padded_seq = std::string(config.clip_penalty, 'N') + seq + std::string(config.clip_penalty, 'N');
			harsh_aligner.Align(padded_seq.c_str(), alt_consensus.c_str(), alt_consensus.length(), filter, &consensus_alt_aln, 0);


			if (get_left_clip_size(consensus_alt_aln) <= config.clip_penalty && get_right_clip_size(consensus_alt_aln) <= config.clip_penalty) {
				int mismatches = consensus_alt_aln.mismatches - 2*config.clip_penalty + get_left_clip_size(consensus_alt_aln) +
						get_right_clip_size(consensus_alt_aln);

				if (mismatches <= seq.length()*config.max_seq_error) {
					sv->alt_better_goodqual++;
					if (alt_strong[i]) sv->alt_better_strong++;
				}
			}
			i++;
		}

		std::string padded_alt_consensus = std::string(config.clip_penalty, 'N') + alt_consensus + std::string(config.clip_penalty, 'N');

		int contig_len = chr_seqs.get_len(contig_name);

		char ref_cstr[100000];
		int ref_realn_start = ins_start - 2*alt_consensus.length(), ref_realn_end = ins_end + 2*alt_consensus.length();
		if (ref_realn_start < 0) ref_realn_start = 0;
		if (ref_realn_end >= contig_len) ref_realn_end = contig_len-1;
		int ref_realn_len = ref_realn_end - ref_realn_start;
		strncpy(ref_cstr, chr_seqs.get_seq(contig_name)+ref_realn_start, ref_realn_len);
		for (int i = 0; i < ref_realn_len; i++) {
			ref_cstr[i] = toupper(ref_cstr[i]);
		}

		remap_consensus_to_reference(alt_consensus, ref_cstr, ref_realn_len, aligner, sv->remapped_cigar, sv->remapped_ins_len, sv->lowq_alt_allele);
	}
	delete[] ref_seq;

	auto gt2 = calculate_small_ins_genotype(sv->alt_better_goodqual, sv->ref_better);
	sv->called_gt = gt2.first, sv->ratio = gt2.second;

//	mtx.lock();
//	flog << log_ss.str() << std::endl;
//	mtx.unlock();
}

struct read_pos_w_strand {
	hts_pos_t start, end;
	bool is_fwd;

	read_pos_w_strand(hts_pos_t start, hts_pos_t end, bool is_fwd) : start(start), end(end), is_fwd(is_fwd) {}
};
void genotype_large_ins(std::string& contig_name, insertion_t* sv, open_samFile_t* bam_file, StripedSmithWaterman::Aligner& aligner,
		StripedSmithWaterman::Aligner& harsh_aligner, std::unordered_map<std::string, std::string>& mateseqs_map, stats_t stats) {

	std::stringstream log_ss;

	int ins_start = sv->start, ins_end = sv->end;
	int contig_len = chr_seqs.get_len(contig_name);

	int extend = config.read_len + 20;

	// build alt allele
	/*
	 * POS in VCF is the base BEFORE the insertion
	 * END seems to be the base BEFORE the reference resumes - i.e., for a "clean" insertion (no deletion), POS == END,
	 * otherwise the last base deleted
	 * As usual, in order to make intervals [ ), we increase the coordinates by 1
	 */
	ins_start++; ins_end++;

	char* contig_seq = chr_seqs.get_seq(contig_name);

	int alt_start = std::max(0, ins_start-extend);
	int alt_end = std::min(ins_end+extend, contig_len);
	int alt_lf_len = ins_start-alt_start, alt_rf_len = alt_end-ins_end;
	int ins_lh_len = std::min(extend, (int) sv->ins_seq.length());
	int ins_rh_len = ins_lh_len;

	int alt_bp1_len = alt_lf_len + ins_lh_len;
	char* alt_bp1_seq = new char[alt_bp1_len+1];
	strncpy(alt_bp1_seq, contig_seq+alt_start, alt_lf_len);
	strncpy(alt_bp1_seq+alt_lf_len, sv->ins_seq.c_str(), ins_lh_len);
	alt_bp1_seq[alt_bp1_len] = 0;

	int alt_bp2_len = ins_rh_len + alt_rf_len;
	char* alt_bp2_seq = new char[alt_bp2_len+1];
	strncpy(alt_bp2_seq, sv->ins_seq.c_str()+(sv->ins_seq.length()-ins_rh_len), ins_rh_len);
	strncpy(alt_bp2_seq+ins_rh_len, contig_seq+ins_end, alt_rf_len);
	alt_bp2_seq[alt_bp2_len] = 0;

	int ref_bp1_start = alt_start, ref_bp1_end = std::min(ins_start+extend, contig_len);
	int ref_bp1_len = ref_bp1_end - ref_bp1_start;
	int ref_bp2_start = std::max(0, ins_end-extend), ref_bp2_end = alt_end;
	int ref_bp2_len = ref_bp2_end - ref_bp2_start;

	//    if (strchr(alt_seq, 'N') != NULL || strchr(ref_bp1_seq, 'N') != NULL || strchr(ref_bp2_seq, 'N') != NULL) { // do not genotype if N region
	//		return;
	//	}

	char l_region[100], r_region[100];
	sprintf(l_region, "%s:%d-%d", contig_name.c_str(), ins_start-stats.max_is, ref_bp1_end);
	sprintf(r_region, "%s:%d-%d", contig_name.c_str(), ref_bp2_start, ins_end+stats.max_is);
	char* regions[] = {l_region, r_region};

	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);
	bam1_t* read = bam_init1();

	std::vector<std::string> read_seqs, qnames, mate_seqs, stable_seqs;
	std::vector<bool> stable_seqs_lc, stable_seqs_rc;
	std::vector<read_pos_w_strand> mate_seqs_stable_pos;
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read)) continue;

		if (is_dc_pair(read)) {
			std::string mate_qname = bam_get_qname(read);
			if (read->core.tid == read->core.mtid) {
				mate_qname += (read->core.flag & BAM_FREAD1 ? "_2" : "_1");
			}
			std::string mate_seq = mateseqs_map[mate_qname];

			if (!bam_is_rev(read)) rc(mate_seq);

			mate_seqs.push_back(mate_seq);
			mate_seqs_stable_pos.push_back({read->core.pos, bam_endpos(read), !bam_is_rev(read)});
			stable_seqs.push_back(get_sequence(read));
			stable_seqs_lc.push_back(get_left_clip_size(read) > 0);
			stable_seqs_rc.push_back(get_right_clip_size(read) > 0);

			read_seqs.push_back(mate_seq);
			qnames.push_back(mate_qname + "_MATE");
		}

		// do not consider reads that do not intersect the insertion
		if (get_unclipped_end(read) <= ins_start) continue;
		if (get_unclipped_start(read) >= ins_end) continue;

		read_seqs.push_back(get_sequence(read));
		qnames.push_back(bam_get_qname(read));
	}

	// find semi-mapped reads
	char region[100];
	sprintf(region, "%s:%d-%d", contig_name.c_str(), ins_start-stats.max_is, ins_end+stats.max_is);
	iter = sam_itr_querys(bam_file->idx, bam_file->header, region);
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read)) continue;
		bool is_clipped = is_left_clipped(read, config.min_clip_size) || is_right_clipped(read, config.min_clip_size);
		bool mate_stable =  bam_is_mrev(read) ?
				(get_mate_end(read)-ins_end > 0 && get_mate_end(read)-ins_end <= stats.max_is) :
				(ins_start-read->core.mpos > 0 && ins_start-read->core.mpos <= stats.max_is);
		if (is_clipped && mate_stable) {
			std::string seq = get_sequence(read);
			if (bam_is_rev(read)) rc(seq); // get it in fastq orientation
			if (!bam_is_mrev(read)) rc(seq);
			mate_seqs.push_back(seq);
			mate_seqs_stable_pos.push_back({read->core.mpos, get_mate_end(read), !bam_is_mrev(read)});
		}

		if (is_mate_unmapped(read) || is_clipped || is_mate_clipped(read) || !is_samechr(read) || bam_is_rev(read) || !bam_is_mrev(read))
			continue;

		int midpoint = (ins_start + ins_end)/2;
		int start = read->core.pos + config.read_len, end = read->core.mpos + config.read_len;
		if (start <= midpoint && midpoint <= end) sv->conc_pairs++;
	}

	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment alt1_aln, alt2_aln, ref1_aln, ref2_aln;
	log_ss << sv->id << std::endl;
	int alt_better = 0, ref_better = 0, same = 0;
	std::vector<std::string> alt_bp1_better_seqs, alt_bp2_better_seqs;
	std::vector<bool> alt_bp1_strong, alt_bp2_strong;
	int i = 0;
	for (std::string& seq : read_seqs) {
		// TODO: right now we follow what we did for deletions - just count alt better and ref better, without distinguishing which bp they support
		// if needed, in the future we may be stricter and compute individual GTs for the two breakpoints, and check whether they are the same

		// align to ALT
		aligner.Align(seq.c_str(), alt_bp1_seq, alt_bp1_len, filter, &alt1_aln, 0);
		aligner.Align(seq.c_str(), alt_bp2_seq, alt_bp2_len, filter, &alt2_aln, 0);
		StripedSmithWaterman::Alignment& alt_aln = alt1_aln.sw_score >= alt2_aln.sw_score ? alt1_aln : alt2_aln;

		// align to REF (two breakpoints)
		aligner.Align(seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_len, filter, &ref1_aln, 0);
		aligner.Align(seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_len, filter, &ref2_aln, 0);
		StripedSmithWaterman::Alignment& ref_aln = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln : ref2_aln;

		if (alt_aln.sw_score > ref_aln.sw_score && accept_aln(alt_aln, config.min_clip_size) && !accept_aln(ref_aln, config.min_clip_size)) {
			alt_better++;
			if (alt1_aln.sw_score >= alt2_aln.sw_score) {
				alt_bp1_better_seqs.push_back(seq);
				alt_bp1_strong.push_back(alt1_aln.sw_score-ref_aln.sw_score >= config.min_score_diff || !accept_aln(ref_aln, config.min_clip_size));
			}
			if (alt2_aln.sw_score >= alt1_aln.sw_score) {
				alt_bp2_better_seqs.push_back(seq);
				alt_bp2_strong.push_back(alt2_aln.sw_score-ref_aln.sw_score >= config.min_score_diff || !accept_aln(ref_aln, config.min_clip_size));
			}
			log_ss << qnames[i] << " " << alt_aln.cigar_string << " " << ref_aln.cigar_string << " " << alt_aln.sw_score-ref_aln.sw_score << std::endl;
		} else if (alt_aln.sw_score < ref_aln.sw_score && accept_aln(ref_aln, config.min_clip_size)) {
			ref_better++;
		} else if (alt_aln.sw_score == ref_aln.sw_score && accept_aln(alt_aln, config.min_clip_size) && accept_aln(ref_aln, config.min_clip_size)) same++;

		i++;
	}
	sv->alt_better = alt_better; sv->ref_better = ref_better;

	// fwd stable mates should only be mapped to ins_seq_lh, and rev stable mates only to ins_seq_rh
	std::string ins_seq_lh = sv->ins_seq.substr(0, std::min((int) sv->ins_seq.length(), config.max_is));
	std::string ins_seq_rh = sv->ins_seq.substr(std::max(0, (int) sv->ins_seq.length()-config.max_is));

	StripedSmithWaterman::Alignment aln;
	std::vector<std::string> accepted_mates_fwd_stable, accepted_mates_rev_stable, accepted_fwd_stable, accepted_rev_stable;
	std::vector<int> accepted_mates_fwd_start, accepted_mates_fwd_end, accepted_mates_rev_start, accepted_mates_rev_end;
	std::vector<bool> accepted_fwd_stable_lc, accepted_rev_stable_rc;
	int dc_pairs = 0;
	for (int i = 0; i < mate_seqs.size(); i++) {
		read_pos_w_strand& rpos = mate_seqs_stable_pos[i];
		std::string padded_seq = std::string(config.clip_penalty, 'N') + mate_seqs[i] + std::string(config.clip_penalty, 'N');
		bool accept;
		if (rpos.is_fwd) {
			aligner.Align(padded_seq.c_str(), ins_seq_lh.c_str(), ins_seq_lh.length(), filter, &aln, 0);
			accept = get_left_clip_size(aln) <= config.clip_penalty || aln.ref_begin <= 5;
			accept &= get_right_clip_size(aln) <= config.clip_penalty || aln.ref_end >= ins_seq_lh.length()-1-5;
		} else {
			aligner.Align(padded_seq.c_str(), ins_seq_rh.c_str(), ins_seq_rh.length(), filter, &aln, 0);
			accept = get_left_clip_size(aln) <= config.clip_penalty || aln.ref_begin <= 5;
			accept &= get_right_clip_size(aln) <= config.clip_penalty || aln.ref_end >= ins_seq_rh.length()-1-5;
		}

		if (aln.query_end-aln.query_begin+1 >= padded_seq.length()/2) {
			if (accept) {
				if (rpos.is_fwd) {
					sv->fwd_anchor_start = std::min(sv->fwd_anchor_start, rpos.start);
					sv->fwd_anchor_end = std::max(sv->fwd_anchor_end, rpos.end);
					accepted_mates_fwd_stable.push_back(mate_seqs[i]);
					if (i < stable_seqs.size()) {
						accepted_fwd_stable.push_back(stable_seqs[i]);
						accepted_fwd_stable_lc.push_back(stable_seqs_lc[i]);
					}
				} else {
					sv->rev_anchor_start = std::min(sv->rev_anchor_start, rpos.start);
					sv->rev_anchor_end = std::max(sv->rev_anchor_end, rpos.end);
					accepted_mates_rev_stable.push_back(mate_seqs[i]);
					if (i < stable_seqs.size()) {
						accepted_rev_stable.push_back(stable_seqs[i]);
						accepted_rev_stable_rc.push_back(stable_seqs_rc[i]);
					}
				}
			}
		}
	}

	auto gt1 = calculate_large_ins_genotype(alt_better, ref_better);
	sv->before_read_qc_gt = gt1.first;
	log_ss << alt_better << " " << ref_better << std::endl;
	log_ss << gt1.first << " " << gt1.second << std::endl;

	std::string alt_bp1_consensus, alt_bp2_consensus;
	StripedSmithWaterman::Alignment consensus_alt_aln;
	if (alt_better > 4*stats.max_depth || ref_better > 2*stats.max_depth ||
			accepted_fwd_stable.size() > stats.max_pairs_crossing || accepted_rev_stable.size() > stats.max_pairs_crossing) {
		sv->gt_skipped = true;
		return;
	}


	// see if we can build a full junction sequence
	std::vector<StripedSmithWaterman::Alignment> consensus_contigs_alns;
	std::string consensus_log;
	std::string junction_seq;

	std::vector<std::string> all_alt_seqs;
	all_alt_seqs.insert(all_alt_seqs.end(), accepted_fwd_stable.begin(), accepted_fwd_stable.end());
	all_alt_seqs.insert(all_alt_seqs.end(), alt_bp1_better_seqs.begin(), alt_bp1_better_seqs.end());
	all_alt_seqs.insert(all_alt_seqs.end(), accepted_mates_fwd_stable.begin(), accepted_mates_fwd_stable.end());
	all_alt_seqs.insert(all_alt_seqs.end(), accepted_mates_rev_stable.begin(), accepted_mates_rev_stable.end());
	all_alt_seqs.insert(all_alt_seqs.end(), alt_bp2_better_seqs.begin(), alt_bp2_better_seqs.end());
	all_alt_seqs.insert(all_alt_seqs.end(), accepted_rev_stable.begin(), accepted_rev_stable.end());
	std::vector<bool> all_alt_seqs_lc(all_alt_seqs.size(), false), all_alt_seqs_rc(all_alt_seqs.size(), false);
	junction_seq = generate_consensus(all_alt_seqs, all_alt_seqs_lc, all_alt_seqs_rc, aligner, harsh_aligner, consensus_contigs_alns, consensus_log);
	if (!junction_seq.empty()) {
		char ref_cstr[100000];
		int ref_realn_start = ins_start - 2*junction_seq.length(), ref_realn_end = ins_end + 2*junction_seq.length();
		if (ref_realn_start < 0) ref_realn_start = 0;
		if (ref_realn_end >= contig_len) ref_realn_end = contig_len-1;
		int ref_realn_len = ref_realn_end - ref_realn_start;
		strncpy(ref_cstr, chr_seqs.get_seq(contig_name)+ref_realn_start, ref_realn_len);
		for (int i = 0; i < ref_realn_len; i++) {
			ref_cstr[i] = toupper(ref_cstr[i]);
		}

		std::string remapped_cigar;
		int remapped_ins_len;
		bool lowq_alt_allele;
		remap_consensus_to_reference(junction_seq, ref_cstr, ref_realn_len, aligner, remapped_cigar, remapped_ins_len, lowq_alt_allele);

		if (!lowq_alt_allele) {
			sv->remapped_cigar = remapped_cigar;
			sv->remapped_ins_len = remapped_ins_len;
		}
	}


	// see if we can merge the two consensus into a full junction sequence

	consensus_contigs_alns.clear();
	std::vector<std::string> alt_bp1_building_seqs;
	std::vector<bool> alt_bp1_building_seqs_lc(alt_bp1_better_seqs.size(), false);
	std::vector<bool> alt_bp1_building_seqs_rc(alt_bp1_better_seqs.size()+accepted_fwd_stable.size(), false);
	alt_bp1_building_seqs.insert(alt_bp1_building_seqs.end(), alt_bp1_better_seqs.begin(), alt_bp1_better_seqs.end());
	alt_bp1_building_seqs.insert(alt_bp1_building_seqs.end(), accepted_fwd_stable.begin(), accepted_fwd_stable.end());
	alt_bp1_building_seqs_lc.insert(alt_bp1_building_seqs_lc.end(), accepted_fwd_stable_lc.begin(), accepted_fwd_stable_lc.end());
	alt_bp1_consensus = generate_consensus(alt_bp1_building_seqs, alt_bp1_building_seqs_lc, alt_bp1_building_seqs_rc,
			aligner, harsh_aligner, consensus_contigs_alns, consensus_log);

	i = 0;
	for (std::string& seq : alt_bp1_better_seqs) {
		std::string padded_seq = std::string(config.clip_penalty, 'N') + seq + std::string(config.clip_penalty, 'N');
		harsh_aligner.Align(padded_seq.c_str(), alt_bp1_consensus.c_str(), alt_bp1_consensus.length(), filter, &consensus_alt_aln, 0);

		if (get_left_clip_size(consensus_alt_aln) <= config.clip_penalty && get_right_clip_size(consensus_alt_aln) <= config.clip_penalty) {
			int mismatches = consensus_alt_aln.mismatches - 2*config.clip_penalty + get_left_clip_size(consensus_alt_aln) +
					get_right_clip_size(consensus_alt_aln);

			if (mismatches <= seq.length()*config.max_seq_error) {
				sv->alt_better_goodqual++;
				sv->alt_bp1_better_goodqual++;
				if (alt_bp1_strong[i]) {
					sv->alt_better_strong++;
					sv->alt_bp1_better_strong++;
				}
			}
		}
		i++;
	}

	consensus_contigs_alns.clear();
	std::vector<std::string> alt_bp2_building_seqs;
	std::vector<bool> alt_bp2_building_seqs_lc(alt_bp2_better_seqs.size()+accepted_rev_stable.size(), false);
	std::vector<bool> alt_bp2_building_seqs_rc(alt_bp2_better_seqs.size(), false);
	alt_bp2_building_seqs.insert(alt_bp2_building_seqs.end(), alt_bp2_better_seqs.begin(), alt_bp2_better_seqs.end());
	alt_bp2_building_seqs.insert(alt_bp2_building_seqs.end(), accepted_rev_stable.begin(), accepted_rev_stable.end());
	alt_bp2_building_seqs_rc.insert(alt_bp2_building_seqs_rc.end(), accepted_rev_stable_rc.begin(), accepted_rev_stable_rc.end());
	alt_bp2_consensus = generate_consensus(alt_bp2_building_seqs, alt_bp2_building_seqs_lc, alt_bp2_building_seqs_rc,
			aligner, harsh_aligner, consensus_contigs_alns, consensus_log);

	i = 0;
	for (std::string& seq : alt_bp2_better_seqs) {
		std::string padded_seq = std::string(config.clip_penalty, 'N') + seq + std::string(config.clip_penalty, 'N');
		harsh_aligner.Align(padded_seq.c_str(), alt_bp2_consensus.c_str(), alt_bp2_consensus.length(), filter, &consensus_alt_aln, 0);

		if (get_left_clip_size(consensus_alt_aln) <= config.clip_penalty && get_right_clip_size(consensus_alt_aln) <= config.clip_penalty) {
			int mismatches = consensus_alt_aln.mismatches - 2*config.clip_penalty + get_left_clip_size(consensus_alt_aln) +
					get_right_clip_size(consensus_alt_aln);

			if (mismatches <= seq.length()*config.max_seq_error) {
				sv->alt_better_goodqual++;
				sv->alt_bp2_better_goodqual++;
				if (alt_bp2_strong[i]) {
					sv->alt_better_strong++;
					sv->alt_bp2_better_strong++;
				}
			}
		}
		i++;
	}

	if (sv->remapped_cigar.empty()) {
		junction_seq = "";

		// check if one of the consensuses is completely contained in the other
		std::string ref = alt_bp1_consensus;
		std::string query = alt_bp2_consensus;
		if (ref.length() < query.length()) std::swap(ref, query);
		if (!query.empty() && is_contained(ref, query, query.length()*config.max_seq_error)) {
			junction_seq = ref;
		}

		if (junction_seq.empty()) {
			suffix_prefix_aln_t spa = aln_suffix_prefix(alt_bp1_consensus, alt_bp2_consensus, 1, -4, config.max_seq_error, config.min_clip_size);
			if (spa.overlap) {
				junction_seq = alt_bp1_consensus + alt_bp2_consensus.substr(spa.overlap);
			}
		}

		if (!junction_seq.empty()) {
			char ref_cstr[100000];
			int ref_realn_start = ins_start - 2*junction_seq.length(), ref_realn_end = ins_end + 2*junction_seq.length();
			if (ref_realn_start < 0) ref_realn_start = 0;
			if (ref_realn_end >= contig_len) ref_realn_end = contig_len-1;
			int ref_realn_len = ref_realn_end - ref_realn_start;
			strncpy(ref_cstr, chr_seqs.get_seq(contig_name)+ref_realn_start, ref_realn_len);
			for (int i = 0; i < ref_realn_len; i++) {
				ref_cstr[i] = toupper(ref_cstr[i]);
			}

			remap_consensus_to_reference(junction_seq, ref_cstr, ref_realn_len, aligner, sv->remapped_cigar, sv->remapped_ins_len, sv->lowq_alt_allele);
		}
	}

	if (accepted_mates_fwd_stable.size() + accepted_mates_rev_stable.size() < 1000) {
		StripedSmithWaterman::Alignment aln;
		std::vector<StripedSmithWaterman::Alignment> consensus_contigs_alns1, consensus_contigs_alns2;
		std::string consensus_log;

		std::string fwd_stable_mates_consensus = generate_reference_guided_consensus(sv->ins_seq, accepted_mates_fwd_stable, aligner, harsh_aligner, consensus_contigs_alns1, consensus_log);
		for (std::string& seq : accepted_mates_fwd_stable) {
			std::string padded_seq = std::string(config.clip_penalty, 'N') + seq + std::string(config.clip_penalty, 'N');
			harsh_aligner.Align(padded_seq.c_str(), fwd_stable_mates_consensus.c_str(), fwd_stable_mates_consensus.length(), filter, &aln, 0);

			bool accept = get_left_clip_size(aln) <= config.clip_penalty || aln.ref_begin <= 10;
			accept &= get_right_clip_size(aln) <= config.clip_penalty || aln.ref_end >= fwd_stable_mates_consensus.length()-1-10;
			int mismatches = aln.mismatches - 2*config.clip_penalty + std::min(config.clip_penalty, get_left_clip_size(aln)) +
					std::min(config.clip_penalty, get_right_clip_size(aln));
			if (accept && mismatches <= seq.length()*config.max_seq_error) {
				sv->disc_pairs_fwd_stable++;
			}

			log_ss << "DISC " << seq << " " << aln.cigar_string << " " << aln.ref_begin << " " << aln.ref_end << " ";
			log_ss << "STABLE FWD " << (aln.query_end-aln.query_begin+1 >= padded_seq.length()/2) << " ";
			log_ss << (accept && mismatches <= seq.length()*config.max_seq_error) << std::endl;
		}

		std::string rev_stable_mates_consensus = generate_reference_guided_consensus(sv->ins_seq, accepted_mates_rev_stable, aligner, harsh_aligner, consensus_contigs_alns2, consensus_log);
		for (std::string& seq : accepted_mates_rev_stable) {
			std::string padded_seq = std::string(config.clip_penalty, 'N') + seq + std::string(config.clip_penalty, 'N');
			harsh_aligner.Align(padded_seq.c_str(), rev_stable_mates_consensus.c_str(), rev_stable_mates_consensus.length(), filter, &aln, 0);

			bool accept = get_left_clip_size(aln) <= config.clip_penalty || aln.ref_begin <= 10;
			accept &= get_right_clip_size(aln) <= config.clip_penalty || aln.ref_end >= rev_stable_mates_consensus.length()-1-10;
			int mismatches = aln.mismatches - 2*config.clip_penalty + std::min(config.clip_penalty, get_left_clip_size(aln)) +
					std::min(config.clip_penalty, get_right_clip_size(aln));
			if (accept && mismatches <= seq.length()*config.max_seq_error) {
				sv->disc_pairs_rev_stable++;
			}

			log_ss << "DISC " << seq << " " << aln.cigar_string << " " << aln.ref_begin << " " << aln.ref_end << " ";
			log_ss << "STABLE REV " << (aln.query_end-aln.query_begin+1 >= padded_seq.length()/2) << " ";
			log_ss << (accept && mismatches <= seq.length()*config.max_seq_error) << std::endl;
		}
	}

	auto gt2 = calculate_large_ins_genotype(sv->alt_better_goodqual, sv->ref_better);
	sv->called_gt = gt2.first, sv->ratio = gt2.second;
	log_ss << alt_better << " " << ref_better << std::endl;
	log_ss << gt2.first << " " << gt2.second << std::endl;
	log_ss << "ALT BP1 CONSENSUS: " << alt_bp1_consensus << std::endl;
	log_ss << "ALT BP2 CONSENSUS: " << alt_bp2_consensus << std::endl;
	log_ss << "JUNCTION SEQ: " << junction_seq << std::endl;

	mtx.lock();
	flog << log_ss.str() << std::endl;
	mtx.unlock();
}


insertion_t* vcf_record_to_insertion(bcf1_t* vcf_ins) {
	bcf_unpack(vcf_ins, BCF_UN_INFO);
    insertion_t* ins = new insertion_t(vcf_ins->d.id, vcf_ins->pos, get_sv_end(vcf_ins, vcf_header), get_ins_seq(vcf_ins, vcf_header));
    ins->bcf_entry = vcf_ins;
    return ins;
}

void genotype_inss(int id, std::string contig_name, std::vector<bcf1_t*> vcf_inss, stats_t stats, std::string workdir) {

	mtx.lock();
	std::cerr << "Genotyping insertions in " << contig_name << std::endl;
	mtx.unlock();

	// read mateseqs
	std::unordered_map<std::string, std::string> mateseqs;
	std::ifstream mateseqs_fin(workdir + "/mateseqs/" + contig_name + ".txt");
	std::string qname, seq;
	while (mateseqs_fin >> qname >> seq) {
		mateseqs[qname] = seq;
	}

	auto is_small = [](insertion_t* ins) { return ins->ins_seq.length() <= config.read_len-2*config.min_clip_size; };

	StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, true);
	StripedSmithWaterman::Aligner harsh_aligner(1, 4, 100, 1, true);
	open_samFile_t* bam_file = open_samFile(bam_fname);
	hts_set_fai_filename(bam_file->file, reference_fname.c_str());

	std::vector<insertion_t*> insertions, shorter_inss;
	for (bcf1_t* vcf_ins : vcf_inss) {
		insertion_t* ins = vcf_record_to_insertion(vcf_ins);
		if (is_small(ins)) {
			genotype_small_ins(contig_name, ins, bam_file, aligner, harsh_aligner, mateseqs, stats);
			shorter_inss.push_back(ins);
		} else {
			genotype_large_ins(contig_name, ins, bam_file, aligner, harsh_aligner, mateseqs, stats);
		}
		insertions.push_back(ins);
	}

	depth_filter_ins(contig_name, chr_seqs.get_len(contig_name), insertions, bam_file, config.max_is);
	calculate_confidence_interval_size(contig_name, shorter_inss, bam_file, config, stats);

	close_samFile(bam_file);

	for (insertion_t* sv : insertions) {
		int disc_pairs = sv->disc_pairs_fwd_stable + sv->disc_pairs_rev_stable;
		int min_dp = min_disc_pairs_for_ins_by_size[std::min(stats.max_is, (int) sv->ins_seq.length())];
		int max_dp = 2*max_disc_pairs_for_ins_by_size[std::min(stats.max_is, (int) sv->ins_seq.length())];
		std::string evidence = "SR";

		std::string filter;
		if (sv->gt_skipped) {
			sv->called_gt = "NO_GT";
			filter = "SKIPPED";
		} else {
			if (!sv->alt_better_goodqual && !sv->ref_better) filter += "NO_USEFUL_READS;";

			if (sv->left_flanking_cov > stats.max_depth || sv->right_flanking_cov > stats.max_depth ||
				sv->left_flanking_cov < stats.min_depth || sv->right_flanking_cov < stats.min_depth) filter += "ANOMALOUS_FLANKING_DEPTH;";
			if (gt_is_positive(sv->called_gt) && sv->alt_better_strong < 3) filter += "LOW_ALT_STRONG_READS;";

			// rescue high-DP insertions
			if (!gt_is_positive(sv->called_gt) && filter.empty() && disc_pairs/double(sv->conc_pairs) > 0.25) {
				sv->called_gt = calculate_large_ins_genotype(disc_pairs, sv->conc_pairs).first;
				evidence = "DP_ONLY";
			}

			if (gt_is_positive(sv->called_gt)) {
				if (is_small(sv)) {
					int ins_len = sv->ins_seq.length();
					if (sv->called_gt == "HOM_ALT" && ins_len > -sv->min_conf_size) filter += "SIZE_FILTER;";
					else if (sv->called_gt == "HET" && ins_len > -sv->min_conf_size*2) filter += "SIZE_FILTER;";
				} else {
					if (disc_pairs < min_dp) filter += "NOT_ENOUGH_DISC_PAIRS;";
				}
				if (sv->lowq_alt_allele) filter += "POOR_JUNCTION_SEQ;";
				if (sv->remapped_ins_len < 50) filter += "REMAPPED_INS_LEN_LT50BP;";
				if (sv->fwd_anchor_end > 0 && sv->rev_anchor_end > 0 && sv->fwd_anchor_end-sv->rev_anchor_start > config.min_mh_len) filter += "MH_TOO_LONG;";
			}

			if (filter.empty()) {
				filter = "PASS";
			}
		}

		bcf_subset(out_vcf_header, sv->bcf_entry, 0, {});
		int gt_data[2];
		strgt2datagt(sv->called_gt, gt_data);
		bcf_update_genotypes(out_vcf_header, sv->bcf_entry, gt_data, 2);

		const char* ft_val[1];
		ft_val[0] = filter.c_str();
		bcf_update_format_string(out_vcf_header, sv->bcf_entry, "FT", ft_val, 1);

		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "AR", &(sv->alt_better_goodqual), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "SR", &(sv->alt_better_strong), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "RR", &(sv->ref_better), 1);

//		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "AR1", &(sv->alt_bp1_better_goodqual), 1);
//		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "SR1", &(sv->alt_bp1_better_strong), 1);
//		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "AR2", &(sv->alt_bp2_better_goodqual), 1);
//		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "SR2", &(sv->alt_bp2_better_strong), 1);

		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "MDFL", &(sv->mleft_flanking_cov), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "MDFR", &(sv->mright_flanking_cov), 1);

		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "DP", &(disc_pairs), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "DPFS", &(sv->disc_pairs_fwd_stable), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "DPRS", &(sv->disc_pairs_rev_stable), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "mDP", &(min_dp), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "MDP", &(max_dp), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "CP", &(sv->conc_pairs), 1);

		if (sv->fwd_anchor_end > 0) {
			const char* fan_val[1];
			std::string fan_val_str = std::to_string(sv->fwd_anchor_start) + "-" + std::to_string(sv->fwd_anchor_end);
			fan_val[0] = fan_val_str.c_str();
			bcf_update_format_string(out_vcf_header, sv->bcf_entry, "FAN", fan_val, 1);
		}
		if (sv->rev_anchor_end > 0) {
			const char* ran_val[1];
			std::string ran_val_str = std::to_string(sv->rev_anchor_start) + "-" + std::to_string(sv->rev_anchor_end);
			ran_val[0] = ran_val_str.c_str();
			bcf_update_format_string(out_vcf_header, sv->bcf_entry, "RAN", ran_val, 1);
		}

		if (is_small(sv)) {
			bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "MCS", &(sv->min_conf_size), 1);
		}
		if (sv->remapped_ins_len != INT32_MAX) {
			bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "RIL", &(sv->remapped_ins_len), 1);
			const char* rcg_val[1];
			rcg_val[0] = sv->remapped_cigar.c_str();
			bcf_update_format_string(out_vcf_header, sv->bcf_entry, "RCG", rcg_val, 1);
		}

		const char* ed_val[1];
		ed_val[0] = evidence.c_str();
		bcf_update_format_string(out_vcf_header, sv->bcf_entry, "ED", ed_val, 1);
	}
}

int max_size_to_test() {
    return config.max_is;
}

void estimate_stats(int id, std::string contig_name, std::vector<int> rnd_positions) {

    sort(rnd_positions.begin(), rnd_positions.end());

    open_samFile_t* bam_file = open_samFile(bam_fname);
    if (hts_set_fai_filename(bam_file->file, fai_path(reference_fname.c_str())) != 0) {
        throw "Failed to read reference " + reference_fname;
    }

    std::vector<char*> regions;
    for (int rnd_pos : rnd_positions) {
        std::stringstream ss;
        ss << contig_name << ":" << rnd_pos-config.max_is << "-" << rnd_pos+10;
        char* region = new char[1000];
        strcpy(region, ss.str().c_str());
        regions.push_back(region);
    }

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
    bam1_t* read = bam_init1();

    uint64_t sum_is = 0;
    uint32_t n_is = 0;

    int curr_pos = 0;
    std::vector<uint32_t> rnd_positions_depth(rnd_positions.size());
    std::vector<uint32_t> rnd_positions_n_pairs_crossing(rnd_positions.size());
    std::vector<std::vector<uint32_t> > rnd_positions_dist_between_end_and_rnd(rnd_positions.size());
    std::vector<uint32_t> del_pop, temp_insert_sizes;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read)) continue;

        while (curr_pos < rnd_positions.size() && read->core.pos > rnd_positions[curr_pos]) curr_pos++;
        if (curr_pos >= rnd_positions.size()) break;

        int rlen = read->core.l_qseq;
        if ((read->core.flag & BAM_FPROPER_PAIR) && !bam_is_rev(read) && read->core.isize > 0 && read->core.isize <= config.max_is) {
            if (read->core.pos+rlen/2 <= rnd_positions[curr_pos] && rnd_positions[curr_pos] <= read->core.pos+read->core.isize-rlen/2) {
            	temp_insert_sizes.push_back(read->core.isize);
                sum_is += read->core.isize;
                n_is++;
            }
        }

        if (is_samechr(read) && !is_samestr(read) && read->core.isize > 0 && read->core.isize < config.max_is + max_size_to_test()) {
            int start = read->core.pos + read->core.l_qseq / 2;
            int end = read->core.pos + read->core.isize - read->core.l_qseq / 2;
            if (start >= end) continue;

            for (int i = curr_pos; i < rnd_positions.size() && rnd_positions[i] < end; i++) {
                if (start <= rnd_positions[i] && rnd_positions[i] <= end) {
                    del_pop.push_back(read->core.isize);
                    rnd_positions_n_pairs_crossing[i]++;
                    rnd_positions_dist_between_end_and_rnd[i].push_back(end-rnd_positions[i]);
                }
            }
        }

        if (read->core.tid != read->core.mtid || read->core.qual < 20) continue;

        hts_pos_t endpos = bam_endpos(read);
        for (int i = curr_pos; i < rnd_positions.size() && rnd_positions[i] < endpos; i++) {
            if (read->core.pos <= rnd_positions[i] && bam_endpos(read) >= rnd_positions[i]) {
                rnd_positions_depth[i]++;
            }
        }
    }

    mtx.lock();
    partial_sums.push_back({sum_is, n_is});
    for (uint32_t d : rnd_positions_depth) {
        if (d > 0) depths.push_back(d);
    }
    del_is_population.insert(del_is_population.end(), del_pop.begin(), del_pop.end());
    for (uint32_t n : rnd_positions_n_pairs_crossing) {
    	if (n > 0) n_pairs_crossing.push_back(n);
    }
    sampled_insert_sizes.insert(sampled_insert_sizes.end(), temp_insert_sizes.begin(), temp_insert_sizes.end());
    for (auto& v : rnd_positions_dist_between_end_and_rnd) {
    	if (v.empty()) continue;
    	dist_between_end_and_rnd.push_back(v);
    }
    mtx.unlock();

    for (char* region : regions) {
        delete[] region;
    }
    hts_itr_destroy(iter);
    bam_destroy1(read);

    close_samFile(bam_file);
}

void extract_dc_reads(int id, std::string contig_name, std::string workdir) {

	mtx.lock();
	std::cout << "Categorizing " << contig_name << std::endl;
	mtx.unlock();

    open_samFile_t* bam_file = open_samFile(bam_fname);
    if (hts_set_fai_filename(bam_file->file, fai_path(reference_fname.c_str())) != 0) {
        throw "Failed to read reference " + reference_fname;
    }

    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, contig_name.c_str());
    if (iter == NULL) { // no reads
    	close_samFile(bam_file);
    	return;
    }

    bam1_t* read = bam_init1();
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_mate_unmapped(read) || !is_primary(read) || !is_dc_pair(read)) continue;

		std::string qname = bam_get_qname(read), read_seq = get_sequence(read);
		if (bam_is_rev(read)) {
			rc(read_seq);
		}
		if (is_samechr(read)) {
			if (read->core.flag & BAM_FREAD1) qname += "_1";
			else qname += "_2";
		}
		mtx_contig[read->core.mtid].lock();
		mate_seqs_buffers[read->core.mtid].push_back(qname + " " + read_seq);
		mtx_contig[read->core.mtid].unlock();
    }

    bam_destroy1(read);
	hts_itr_destroy(iter);

	close_samFile(bam_file);
}

int main(int argc, char* argv[]) {

    std::string sv_fname = argv[1];
    std::string workdir = argv[2];

    flog.open(workdir + "/log");
    consensus_flog.open(workdir + "/consensus_log");

    config.parse(workdir + "/config.txt");
    bam_fname = config.bam_fname;
    reference_fname = config.reference_fname;

    std::string contig_name;
    std::ifstream rnd_pos_fin(workdir + "/random_pos.txt");
    int pos;
    std::unordered_map<std::string, std::vector<int> > rnd_pos_map;
    while (rnd_pos_fin >> contig_name >> pos) {
        rnd_pos_map[contig_name].push_back(pos);
    }

    open_samFile_t* bam_file = open_samFile(bam_fname);
	if (hts_set_fai_filename(bam_file->file, fai_path(config.reference_fname.c_str())) != 0) {
		throw "Failed to read reference " + config.reference_fname;
	}

    std::cout << "Estimating BAM statistics." << std::endl;
    ctpl::thread_pool thread_pool1(config.threads);
    std::vector<std::future<void> > futures;
    for (auto& e : rnd_pos_map) {
        std::future<void> future = thread_pool1.push(estimate_stats, e.first, e.second);
        futures.push_back(std::move(future));
    }
    thread_pool1.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        try {
            futures[i].get();
        } catch (char const* s) {
            std::cerr << s << std::endl;
        }
    }
    futures.clear();

    uint64_t sum_is = 0;
    uint32_t n_is = 0;
    for (std::pair<uint64_t, uint32_t>& p : partial_sums) {
        sum_is += p.first;
        n_is += p.second;
    }
    int pop_avg_crossing_is = sum_is / n_is;

    uint64_t tot = 0;
    for (uint32_t is : sampled_insert_sizes) {
    	tot += is;
    }
    pop_avg_crossing_is = tot/sampled_insert_sizes.size();

    uint64_t tot_lt = 0, tot_gt = 0;
    int n_lt = 0, n_gt = 0;
    for (uint32_t is : sampled_insert_sizes) {
    	if (is <= pop_avg_crossing_is) {
    		tot_lt += (is - pop_avg_crossing_is)*(is - pop_avg_crossing_is);
    		n_lt++;
    	}
    	if (is >= pop_avg_crossing_is) {
    		tot_gt += (is - pop_avg_crossing_is)*(is - pop_avg_crossing_is);
    		n_gt++;
    	}
    }
    int stddev_lt = std::sqrt(tot_lt/n_lt), stddev_gt = std::sqrt(tot_gt/n_gt);
    int max_is = pop_avg_crossing_is+3*stddev_gt;

    std::sort(depths.begin(), depths.end());
    int min_depth = depths[depths.size()/100];
    int median_depth = depths[depths.size()/2];
    int max_depth = depths[depths.size()-depths.size()/100];

    std::sort(n_pairs_crossing.begin(), n_pairs_crossing.end());
    int min_pairs_crossing = n_pairs_crossing[n_pairs_crossing.size()/100];
    int max_pairs_crossing = n_pairs_crossing[n_pairs_crossing.size()-1-(n_pairs_crossing.size()/100)];
    std::cout << "MIN PAIRS CROSSING " << min_pairs_crossing << std::endl;

    std::vector<std::vector<uint32_t>> pairs_crossing_dists(max_is+1);
    for (int i = 0; i < dist_between_end_and_rnd.size(); i++) {
    	std::vector<uint32_t> dist(max_is+1);
    	for (uint32_t val : dist_between_end_and_rnd[i]) {
    		if (val <= max_is) dist[val]++;
    	}
    	pairs_crossing_dists[0].push_back(dist[0]);
    	for (int j = 1; j <= max_is; j++) {
    		dist[j] += dist[j-1];
    		pairs_crossing_dists[j].push_back(dist[j]);
    	}
    }

    for (int i = 0; i <= max_is; i++) {
    	std::sort(pairs_crossing_dists[i].begin(), pairs_crossing_dists[i].end());
    	min_disc_pairs_for_ins_by_size.push_back(pairs_crossing_dists[i][pairs_crossing_dists[i].size()/100]);
    	max_disc_pairs_for_ins_by_size.push_back(pairs_crossing_dists[i][pairs_crossing_dists[i].size()-1-(pairs_crossing_dists[i].size()/100)]);
    }

    stats_t stats(pop_avg_crossing_is, pop_avg_crossing_is-3*stddev_lt, max_is, min_depth, median_depth, max_depth, min_pairs_crossing, max_pairs_crossing);

    std::ofstream stats_fout(workdir + "/stats.txt");
    stats_fout << "pop_avg_crossing_is " << stats.pop_avg_crossing_is << std::endl;
    stats_fout << "min_is " << stats.min_is << std::endl;
    stats_fout << "max_is " << stats.max_is << std::endl;
    stats_fout << "min_depth " << stats.min_depth << std::endl;
    stats_fout << "median_depth " << stats.median_depth << std::endl;
    stats_fout << "max_depth " << stats.max_depth << std::endl;
    stats_fout << "min_pairs_crossing " << stats.min_pairs_crossing << std::endl;
    stats_fout << "max_pairs_crossing " << stats.max_pairs_crossing << std::endl;

    std::random_shuffle(del_is_population.begin(), del_is_population.end());
    del_is_population.resize(config.max_pop_size);

    chr_seqs.read_fasta_into_map(reference_fname);

    htsFile* sv_vcf_file = bcf_open(sv_fname.c_str(), "r");
    if (sv_vcf_file == NULL) {
        throw std::runtime_error("Unable to open file " + sv_fname + ".");
    }

    vcf_header = bcf_hdr_read(sv_vcf_file);
    if (vcf_header == NULL) {
        throw std::runtime_error("Failed to read the VCF header.");
    }

    bcf1_t* vcf_record = bcf_init();
    std::unordered_map<std::string, std::vector<bcf1_t*> > dels_by_chr, dups_by_chr, inss_by_chr;
    bool has_insertions = false;
    while (bcf_read(sv_vcf_file, vcf_header, vcf_record) == 0) {
        std::string contig_name = bcf_seqname(vcf_header, vcf_record);
        std::string sv_type = get_sv_type(vcf_record, vcf_header);
        if (sv_type == "DEL") {
            dels_by_chr[contig_name].push_back(bcf_dup(vcf_record));
        } else if (sv_type == "DUP") {
            dups_by_chr[contig_name].push_back(bcf_dup(vcf_record));
        } else if (sv_type == "INS") {
        	inss_by_chr[contig_name].push_back(bcf_dup(vcf_record));
        	has_insertions = true;
        }
    }

    if (has_insertions) {
    	mtx_contig = new std::mutex[bam_file->header->n_targets];
    	mate_seqs_buffers.resize(bam_file->header->n_targets);

    	ctpl::thread_pool thread_pool2(config.threads);
		for (int i = 0; i < bam_file->header->n_targets; i++) {
			std::string contig_name = bam_file->header->target_name[i];
			std::future<void> future = thread_pool2.push(extract_dc_reads, contig_name, workdir);
			futures.push_back(std::move(future));
		}
		thread_pool2.stop(true);
		for (int i = 0; i < futures.size(); i++) {
			try {
				futures[i].get();
			} catch (char const* s) {
				std::cerr << s << std::endl;
			}
		}
		futures.clear();

		for (int i = 0; i < bam_file->header->n_targets; i++) {
			std::string contig_name = bam_file->header->target_name[i];
			std::ofstream mate_seqs_fout(workdir + "/mateseqs/" + contig_name + ".txt");
			for (std::string& seq : mate_seqs_buffers[i]) {
				mate_seqs_fout << seq << std::endl;
			}
		}

		delete[] mtx_contig;
    }

    std::string out_vcf_fname = workdir + "/genotyped.vcf.gz";
    htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
    out_vcf_header = generate_out_hdr(config.sample_name, vcf_header);
    if (bcf_hdr_write(out_vcf_file, out_vcf_header) != 0) {
    	throw std::runtime_error("Failed to read the VCF header.");
    }

    // genotype chrs in descending order of svs
    std::vector<std::pair<int, std::string> > dels_by_chr_nums, dups_by_chr_nums, inss_by_chr_nums;
    for (auto& e : dels_by_chr) {
    	dels_by_chr_nums.push_back({e.second.size(), e.first});
    }
    std::sort(dels_by_chr_nums.begin(), dels_by_chr_nums.end(), std::greater<>());
    for (auto& e : dups_by_chr) {
		dups_by_chr_nums.push_back({e.second.size(), e.first});
	}
    std::sort(dups_by_chr_nums.begin(), dups_by_chr_nums.end(), std::greater<>());
    for (auto& e : inss_by_chr) {
		inss_by_chr_nums.push_back({e.second.size(), e.first});
	}
	std::sort(inss_by_chr_nums.begin(), inss_by_chr_nums.end(), std::greater<>());

    ctpl::thread_pool thread_pool3(config.threads);
    for (auto& p : dels_by_chr_nums) {
    	std::string contig_name = p.second;
        std::future<void> future = thread_pool3.push(genotype_dels, contig_name, chr_seqs.get_seq(contig_name),
        		chr_seqs.get_len(contig_name), dels_by_chr[contig_name], vcf_header, out_vcf_header, stats, config);
        futures.push_back(std::move(future));
    }
    for (auto& p : dups_by_chr_nums) {
    	std::string contig_name = p.second;
    	std::future<void> future = thread_pool3.push(genotype_dups, contig_name, dups_by_chr[contig_name], stats);
    	futures.push_back(std::move(future));
    }
    for (auto& p : inss_by_chr_nums) {
		std::string contig_name = p.second;
		std::future<void> future = thread_pool3.push(genotype_inss, contig_name, inss_by_chr[contig_name], stats, workdir);
		futures.push_back(std::move(future));
	}

    thread_pool3.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        try {
            futures[i].get();
        } catch (char const* s) {
            std::cerr << s << std::endl;
        }
    }
    futures.clear();

    // print contigs in vcf order
    int n_seqs;
    const char** seqnames = bcf_hdr_seqnames(vcf_header, &n_seqs);
    for (int i = 0; i < n_seqs; i++) {
    	std::string contig_name = seqnames[i];
    	std::vector<bcf1_t*> contig_svs;
    	if (dels_by_chr.count(contig_name) > 0) contig_svs.insert(contig_svs.end(), dels_by_chr[contig_name].begin(), dels_by_chr[contig_name].end());
    	if (dups_by_chr.count(contig_name) > 0) contig_svs.insert(contig_svs.end(), dups_by_chr[contig_name].begin(), dups_by_chr[contig_name].end());
    	if (inss_by_chr.count(contig_name) > 0) contig_svs.insert(contig_svs.end(), inss_by_chr[contig_name].begin(), inss_by_chr[contig_name].end());
    	std::sort(contig_svs.begin(), contig_svs.end(), [](const bcf1_t* b1, const bcf1_t* b2) {return b1->pos < b2->pos;});

		for (auto& vcf_record : contig_svs) {
			bcf_update_info_int32(out_vcf_header, vcf_record, "AC", NULL, 0);
			bcf_update_info_int32(out_vcf_header, vcf_record, "AN", NULL, 0);
			if (bcf_write(out_vcf_file, out_vcf_header, vcf_record) != 0) {
				throw std::runtime_error("Failed to write VCF record to " + out_vcf_fname);
			}
		}
    }
    delete[] seqnames;

    bcf_destroy(vcf_record);
    bcf_hdr_destroy(vcf_header);
    hts_close(sv_vcf_file);

    hts_close(out_vcf_file);

    tbx_index_build(out_vcf_fname.c_str(), 0, &tbx_conf_vcf);

}
