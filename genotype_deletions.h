#ifndef GENOTYPE_DELETIONS_H_
#define GENOTYPE_DELETIONS_H_

extern config_t config;

#include "sam_utils.h"
#include "utils.h"

std::pair<std::string, double> calculate_del_genotype(int alt_better_reads, int ref_better_reads) {
	if (alt_better_reads + ref_better_reads == 0) {
		return {"NO_GT", 0.0};
	}

	// in theory, we expect all alt_better for hom_alt calls, and a ratio of 1:2 alt_better:ref_better for het
	double ratio = (double) alt_better_reads/(alt_better_reads + ref_better_reads);
	if (ratio >= 0.75) {
		return {"HOM_ALT", ratio};
	} else if (ratio >= 0.17) {
		return {"HET", ratio};
	} else {
		return {"HOM_REF", ratio};
	}
}

bool overlap_read_consensus_lh(bam1_t* read, int lh_remap_start, std::string& alt_consensus) {
	if (is_left_clipped(read, 1) || lh_remap_start - read->core.pos < config.min_clip_size) return false;

	std::string read_seq = get_sequence(read);
	suffix_prefix_aln_t spa = aln_suffix_prefix(read_seq, alt_consensus, 1, -4, config.max_seq_error, config.min_clip_size);
	return spa.overlap;
}

bool overlap_read_consensus_rh(bam1_t* read, int rh_remap_start, std::string& alt_consensus) {
	if (is_right_clipped(read, 1) || bam_endpos(read) - rh_remap_start < config.min_clip_size) return false;

	std::string read_seq = get_sequence(read);
	suffix_prefix_aln_t spa = aln_suffix_prefix(alt_consensus, read_seq, 1, -4, config.max_seq_error, config.min_clip_size);
	return spa.overlap;
}

void genotype_del(deletion_t* sv, std::string& contig_name, char* contig_seq, int contig_len, open_samFile_t* bam_file,
		StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner, stats_t stats, config_t config,
		std::string& log_string) {

	std::stringstream log_ss;

    int del_start = sv->start, del_end = sv->end;

    int extend = config.read_len + 20;

    // build alt allele
    /* POS in VCF is the base BEFORE the deletion - i.e., the first deleted base is POS+1.
     * Therefore, we want the ALT allele to *include* base POS
     * (note that POS is 1-based in the VCF file, but htslib kindly returns the 0-based coordinate here).
     * As for the END coordinate, my current understanding (which may change) is that it represents the last base deleted.
     * Therefore, the ALT allele should NOT include base END, i.e. it should start at END+1.
     * Here we shift both coordinates by 1, to make them the base immediately AFTER the breakpoints, which is a bit more intuitive for me. */
    del_start++; del_end++;

    // all ranges will be start-inclusive and end-exclusive, i.e. [a,b)
    int alt_start = std::max(0, del_start-extend);
    int alt_end = std::min(del_end+extend, contig_len);
    int alt_lh_len = del_start-alt_start, alt_rh_len = alt_end-del_end;
    int alt_len = alt_lh_len + sv->ins_seq.length() + alt_rh_len;
    char* alt_seq = new char[alt_len + 1];
    strncpy(alt_seq, contig_seq+alt_start, alt_lh_len);
    strncpy(alt_seq+alt_lh_len, sv->ins_seq.c_str(), sv->ins_seq.length());
    strncpy(alt_seq+alt_lh_len+sv->ins_seq.length(), contig_seq+del_end, alt_rh_len);
    alt_seq[alt_len] = 0;

    // extract ref breakpoint allele(s) - will be useful for consensus generation
    int ref_bp1_start = alt_start, ref_bp1_end = std::min(del_start+extend, contig_len);
    int ref_bp1_len = ref_bp1_end - ref_bp1_start;
    char* ref_bp1_seq = new char[ref_bp1_len + 1];
    strncpy(ref_bp1_seq, contig_seq+ref_bp1_start, ref_bp1_len);
    ref_bp1_seq[ref_bp1_len] = 0;

    int ref_bp2_start = std::max(0, del_end-extend), ref_bp2_end = alt_end;
    int ref_bp2_len = ref_bp2_end - ref_bp2_start;
    char* ref_bp2_seq = new char[ref_bp2_len + 1];
    strncpy(ref_bp2_seq, contig_seq+ref_bp2_start, ref_bp2_len);
    ref_bp2_seq[ref_bp2_len] = 0;

    if (strchr(alt_seq, 'N') != NULL || strchr(ref_bp1_seq, 'N') != NULL || strchr(ref_bp2_seq, 'N') != NULL) { // do not genotype if N region
        return;
    }

    char l_region[100], r_region[100];
    sprintf(l_region, "%s:%d-%d", contig_name.c_str(), alt_start, ref_bp1_end);
    sprintf(r_region, "%s:%d-%d", contig_name.c_str(), ref_bp2_start, alt_end);
    char* regions[] = {l_region, r_region};

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);

    bam1_t* read = bam_init1();

    /* Gather reads to remap */
    std::vector<bam1_t*> reads_to_remap;
    std::vector<bool> good_read, read_fwd;
    std::unordered_set<std::string> ref_better_qnames;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;

        // do not consider reads that do not intersect the breakpoints
        if (get_unclipped_end(read) <= del_start) continue;
        if (del_start <= get_unclipped_start(read) && get_unclipped_end(read) <= del_end) {
        	ref_better_qnames.insert(bam_get_qname(read));
        	continue;
        }
        if (get_unclipped_start(read) >= del_end) continue;

        reads_to_remap.push_back(bam_dup1(read));
        good_read.push_back(is_samechr(read) && !is_samestr(read) && !is_mate_unmapped(read));
        read_fwd.push_back(!bam_is_rev(read));
    }

    /* Find reads that align better to the alternative allele */
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alt_aln, ref1_aln, ref2_aln;
    int ref_better = 0, alt_better = 0, same = 0;
    std::vector<std::string> alt_better_seqs, alt_better_qnames;
    std::vector<std::string> alt_better_original_alt_aln, alt_better_original_ref_aln;
    std::vector<int> alt_better_score_diff, alt_better_mapq;
    std::vector<bool> alt_better_good_read, alt_better_is_fwd;
    std::vector<int> alt_better_seqs_start, alt_better_seqs_end;
    std::vector<bam1_t*> ref_better_reads, alt_better_reads;
    for (int i = 0; i < reads_to_remap.size(); i++) {
    	std::string seq = get_sequence(reads_to_remap[i]);

        // align to ALT
        aligner.Align(seq.c_str(), alt_seq, alt_len, filter, &alt_aln, 0);

        // align to REF (two breakpoints)
        aligner.Align(seq.c_str(), ref_bp1_seq, ref_bp1_len, filter, &ref1_aln, 0);
        aligner.Align(seq.c_str(), ref_bp2_seq, ref_bp2_len, filter, &ref2_aln, 0);

        StripedSmithWaterman::Alignment& ref_aln = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln : ref2_aln;
        if (alt_aln.sw_score > ref_aln.sw_score && accept_aln(alt_aln, config.min_clip_size)) {
            alt_better++;
            alt_better_is_fwd.push_back(read_fwd[i]);
            alt_better_seqs.push_back(seq);
            alt_better_qnames.push_back(bam_get_qname(reads_to_remap[i]));
            alt_better_seqs_start.push_back(alt_aln.ref_begin - get_left_clip_size(alt_aln));
            alt_better_seqs_end.push_back(alt_aln.ref_end + get_right_clip_size(alt_aln));

            alt_better_original_alt_aln.push_back(alt_aln.cigar_string);
            alt_better_original_ref_aln.push_back(ref_aln.cigar_string);
            alt_better_score_diff.push_back(alt_aln.sw_score-ref_aln.sw_score);
            alt_better_mapq.push_back(reads_to_remap[i]->core.qual);
            alt_better_good_read.push_back(good_read[i]);
            alt_better_reads.push_back(reads_to_remap[i]);
        }
        else if (alt_aln.sw_score < ref_aln.sw_score && accept_aln(ref_aln, config.min_clip_size)) {
        	ref_better++;
        	ref_better_reads.push_back(reads_to_remap[i]);
        	ref_better_qnames.insert(bam_get_qname(reads_to_remap[i]));
        }
        else if (alt_aln.sw_score == ref_aln.sw_score && accept_aln(alt_aln, config.min_clip_size) && accept_aln(ref_aln, config.min_clip_size)) same++;

        i++;

        if (alt_better + ref_better + same > 4*stats.max_depth) {
        	auto gt = calculate_del_genotype(0, 0);

        	log_ss << "SKIPPING " << sv->id << std::endl;

        	delete[] alt_seq;
			delete[] ref_bp1_seq;
			delete[] ref_bp2_seq;

			bam_destroy1(read);
			sam_itr_destroy(iter);
        	return;
        }
    }

    log_ss << sv->id << std::endl;

    /* Build a consensus alternative allele and count reads that support it */
    StripedSmithWaterman::Alignment consensus_alt_aln;
    std::string alt_consensus;
	std::vector<StripedSmithWaterman::Alignment> consensus_contigs_alns;
	std::string consensus_log;
	alt_consensus = generate_reference_guided_consensus(alt_seq, alt_better_seqs, aligner, harsh_aligner, consensus_contigs_alns, consensus_log, false);

	std::vector<bam1_t*> alt_better_accepted_reads;
	std::vector<StripedSmithWaterman::Alignment> alt_better_accepted_alns;
	for (int i = 0; i < alt_better_reads.size(); i++) {
		std::string seq = alt_better_seqs[i];

		if (ref_better_qnames.count(alt_better_qnames[i])) continue;

		std::string padded_seq = std::string(config.clip_penalty, 'N') + seq + std::string(config.clip_penalty, 'N');
		harsh_aligner.Align(padded_seq.c_str(), alt_consensus.c_str(), alt_consensus.length(), filter, &consensus_alt_aln, 0);

		// we force N to match with any base. However, they will still be counted in aln.mismatches, so we remove them
		if (get_left_clip_size(consensus_alt_aln) <= config.clip_penalty && get_right_clip_size(consensus_alt_aln) <= config.clip_penalty) {
			// (twice, once for each tail)
			// sometimes the N paddings are clipped, (only case I can think of is boundaries of the alt_consensus reached)
			// in that case we removed too many mismatches and we put them back
			int mismatches = consensus_alt_aln.mismatches - 2*config.clip_penalty + get_left_clip_size(consensus_alt_aln) +
					get_right_clip_size(consensus_alt_aln);

			if (mismatches <= seq.length()*config.max_seq_error) {
				sv->alt_better++;
				if (alt_better_is_fwd[i]) sv->alt_better_fwd++;
				else sv->alt_better_rev++;
				if (alt_better_score_diff[i] >= config.min_score_diff && alt_better_good_read[i]) {
					log_ss << alt_better_qnames[i] << " " << seq << "\n" << alt_better_original_alt_aln[i] << " " << alt_better_original_ref_aln[i] << " ";
					log_ss << alt_better_score_diff[i] << " " << mismatches << "/" << (seq.length()*config.max_seq_error) << " ";
					log_ss << consensus_alt_aln.cigar_string << " " << alt_better_mapq[i] << std::endl;
					sv->alt_better_strong++;
					sv->alt_better_max_mq = std::max(sv->alt_better_max_mq, alt_better_mapq[i]);
				}

				// realign without padding for saving as evidence
				harsh_aligner.Align(seq.c_str(), alt_consensus.c_str(), alt_consensus.length(), filter, &consensus_alt_aln, 0);
				alt_better_accepted_reads.push_back(alt_better_reads[i]);
				alt_better_accepted_alns.push_back(consensus_alt_aln);
			}
		}
	}

    sv->ref_better = ref_better;
	sv->same = same;

    log_ss << "REF1: " << ref_bp1_seq << std::endl;
    log_ss << "REF2: " << ref_bp2_seq << std::endl;
    log_ss << "ALT: " << alt_seq << std::endl;
    log_ss << "CONS_ALT: " << alt_consensus << std::endl;
    log_ss << sv->alt_better << " " << sv->ref_better << " " << sv->same << std::endl;

    auto gt = calculate_del_genotype(sv->alt_better, sv->ref_better);
    sv->called_gt = gt.first;
    sv->ratio = gt.second;

	// realign consensus ALT
    int split_i = 0;
    if (gt_is_positive(sv->called_gt)) {
		StripedSmithWaterman::Alignment alt_cons_realn;
		std::string padded_alt_consensus = std::string(config.clip_penalty, 'N') + alt_consensus + std::string(config.clip_penalty, 'N');
		int ref_realn_start = del_start-padded_alt_consensus.length(), ref_realn_end = del_end+padded_alt_consensus.length();
		if (ref_realn_start < 0) ref_realn_start = 0;
		if (ref_realn_end >= contig_len) ref_realn_end = contig_len-1;
		int ref_realn_len = ref_realn_end - ref_realn_start;
		if (sv->len() < config.min_size_for_depth_filtering) {
			aligner.Align(padded_alt_consensus.c_str(), contig_seq+ref_realn_start, ref_realn_len, filter, &alt_cons_realn, 0);
			if (alt_cons_realn.query_begin > 0 && alt_cons_realn.query_end < alt_consensus.length()-1) sv->lowq_alt_allele = true;
		}

		char ref_cstr[100000];
		if (sv->len() < config.min_size_for_depth_filtering) {
			strncpy(ref_cstr, contig_seq+ref_realn_start, ref_realn_len);
		} else { // too long to realing to the whole deletion
			int ref_bp1_realn_start = ref_realn_start;
			int ref_bp1_realn_end = std::min(int(del_start+padded_alt_consensus.length()), contig_len);
			strncpy(ref_cstr, contig_seq+ref_bp1_realn_start, ref_bp1_realn_end-ref_bp1_realn_start);
			int ref_bp2_realn_start = std::max(0, int(del_end-padded_alt_consensus.length()));
			int ref_bp2_realn_end = std::min(int(del_end+padded_alt_consensus.length()), contig_len);
			strncpy(ref_cstr+ref_bp1_realn_end-ref_bp1_realn_start, contig_seq+ref_bp2_realn_start, ref_bp2_realn_end-ref_bp2_realn_start);
			ref_realn_len = strlen(ref_cstr);
		}
		for (int i = 0; i < ref_realn_len; i++) {
			ref_cstr[i] = toupper(ref_cstr[i]);
		}
		int* prefix_scores = smith_waterman_gotoh(ref_cstr, ref_realn_len, padded_alt_consensus.c_str(), padded_alt_consensus.length(), 1, -4, -6, -1);

		for (int i = 0; i < ref_realn_len/2; i++) {
			std::swap(ref_cstr[i], ref_cstr[ref_realn_len-i-1]);
		}
		std::string padded_alt_consensus_rev = std::string(padded_alt_consensus.rbegin(), padded_alt_consensus.rend());
		int* suffix_scores = smith_waterman_gotoh(ref_cstr, ref_realn_len, padded_alt_consensus_rev.c_str(), padded_alt_consensus_rev.length(), 1, -4, -6, -1);

		int max_score = 0;
		for (int i = 0; i < alt_consensus.length(); i++) {
			int prefix_score = prefix_scores[config.clip_penalty+i-1], suffix_score = suffix_scores[alt_consensus.length()-i-1+config.clip_penalty];
			if (prefix_score + suffix_score > max_score) {
				max_score = prefix_score + suffix_score;
				split_i = i;
			}
		}

		if (sv->len() < config.min_size_for_depth_filtering) {
			StripedSmithWaterman::Alignment lh_aln, rh_aln;
			char* padded_alt_consensus_cstr = new char[padded_alt_consensus.length()+1];
			strcpy(padded_alt_consensus_cstr, padded_alt_consensus.c_str());
			char terminator = '\0';
			std::swap(terminator, padded_alt_consensus_cstr[config.clip_penalty+split_i]);
			aligner.Align(padded_alt_consensus_cstr, contig_seq+ref_realn_start, ref_realn_len, filter, &lh_aln, 0);
			std::swap(terminator, padded_alt_consensus_cstr[config.clip_penalty+split_i]);

			aligner.Align(padded_alt_consensus_cstr+config.clip_penalty+split_i, contig_seq+ref_realn_start, ref_realn_len, filter, &rh_aln, 0);

			sv->was_remapped = true;
			sv->remapped_start = ref_realn_start + lh_aln.ref_end;
			sv->remapped_end = ref_realn_start + rh_aln.ref_begin-1;
			sv->rs_cigar = lh_aln.cigar_string;
			sv->re_cigar = rh_aln.cigar_string;
			sv->full_cigar = alt_cons_realn.cigar_string;

			// Find reads vouching for the remapping location of the left half of the consensus
			char region[100];
			int lh_remap_start = ref_realn_start + lh_aln.ref_begin + config.clip_penalty;
			sprintf(region, "%s:%d-%d", contig_name.c_str(), lh_remap_start, lh_remap_start);
			iter = sam_itr_querys(bam_file->idx, bam_file->header, region);
			while (sam_itr_next(bam_file->file, iter, read) >= 0) {
				if (!is_primary(read) || is_unmapped(read)) continue;
				sv->vouching_for_lh_consensus_remapping += overlap_read_consensus_lh(read, lh_remap_start, alt_consensus);
			}

			int rh_remap_end = ref_realn_start + rh_aln.ref_end - config.clip_penalty;
			sprintf(region, "%s:%d-%d", contig_name.c_str(), rh_remap_end, rh_remap_end);
			iter = sam_itr_querys(bam_file->idx, bam_file->header, region);
			while (sam_itr_next(bam_file->file, iter, read) >= 0) {
				if (!is_primary(read) || is_unmapped(read)) continue;
				sv->vouching_for_rh_consensus_remapping += overlap_read_consensus_rh(read, rh_remap_end, alt_consensus);
			}

			log_ss << "LEFT HALF: " << lh_aln.cigar_string << " " << lh_aln.sw_score << std::endl;
			log_ss << padded_alt_consensus.substr(0, config.clip_penalty+split_i) << " " << config.clip_penalty+split_i << std::endl;
			log_ss << "RIGHT HALF: " << rh_aln.cigar_string  << " " << rh_aln.sw_score << std::endl;
			log_ss << padded_alt_consensus.substr(config.clip_penalty+split_i) << " " << padded_alt_consensus.length()-(config.clip_penalty+split_i) << std::endl;
			log_ss << "FULL: " << alt_cons_realn.cigar_string << " " << alt_cons_realn.sw_score << std::endl;

			int lr_mismatches = lh_aln.mismatches + rh_aln.mismatches - 2*config.clip_penalty -
					count_inserted_bps(lh_aln) - count_inserted_bps(rh_aln);
			int lr_indels = count_indels(lh_aln) + count_indels(rh_aln) + 1;
			int full_mismatches = alt_cons_realn.mismatches - 2*config.clip_penalty - count_inserted_bps(alt_cons_realn); // not correct if clipped, but it does not matter (see (b) below)
			int full_indels = count_indels(alt_cons_realn);

			// TODO: we need a model here to pit the split alignment vs the full alignment.
			// Currently, let's (very arbitrarily) assume that 1 indel = 2.5 mismatches
			double INDEL_TO_MISMATCH_CONVERSION = 2.5;
			double lr_delta = lr_mismatches + INDEL_TO_MISMATCH_CONVERSION*lr_indels;
			int full_delta = full_mismatches + INDEL_TO_MISMATCH_CONVERSION*full_indels;

			// either (a) the split alignment is better than the full alignment (need a better model here), or
			// (b) full alignment is clipped, or
			// (c) the full alignment actually include the deletion (often the full alignment is actually identical to the split alignment)
			sv->split_aln_accepted = false;
			if (lr_delta < full_delta ||
				get_left_clip_size(alt_cons_realn) > 0 || get_right_clip_size(alt_cons_realn) > 0 ||
				largest_del(alt_cons_realn) >= config.min_sv_size) {
				sv->split_aln_accepted = true;
			}
			log_ss << lr_mismatches << " " << lr_indels << " " << lr_delta << " ";
			log_ss << full_mismatches << " " << full_indels << " " << full_delta << std::endl;
		}

		delete[] prefix_scores;
		delete[] suffix_scores;
    }

    if (config.save_evidence) {
		std::string evidence_fname = config.workdir + "/evidence/" + sv->id + ".bam";
		bam_hdr_t* evidence_hdr = bam_hdr_dup(bam_file->header);
		std::string alt_consensus_len_str = std::to_string(alt_consensus.length());
		sam_hdr_add_line(evidence_hdr, "SQ", "SN", "ALT", "LN", alt_consensus_len_str.c_str(), NULL);
		samFile* evidence_bam = get_writer(evidence_fname, evidence_hdr);
		for (bam1_t* read : ref_better_reads) {
			if (sam_write1(evidence_bam, evidence_hdr, read) < 0) {
				throw "Failed to write evidence to " + evidence_fname;
			}
		}

		std::vector<bam1_t*> alt_better_accepted_reads_reset;
		for (int i = 0; i < alt_better_accepted_reads.size(); i++) {
			bam1_t* read = alt_better_accepted_reads[i];
			bam1_t* new_read = bam_init1();
			StripedSmithWaterman::Alignment& aln = alt_better_accepted_alns[i];
			uint32_t* cigar = NULL;
			size_t n_cigar = 0;
			sam_parse_cigar(aln.cigar_string.c_str(), NULL, &cigar, &n_cigar);

			bam_set1(new_read, strlen(bam_get_qname(read)), bam_get_qname(read), read->core.flag, sam_hdr_name2tid(evidence_hdr, "ALT"),
					aln.ref_begin, read->core.qual, n_cigar, cigar, read->core.mtid, read->core.mpos, 0, read->core.l_qseq,
					get_sequence(read).c_str(), (const char*) bam_get_qual(read), 0);
			bam_aux_update_int(new_read, "AS", aln.sw_score);
			bam_aux_update_int(new_read, "NM", aln.mismatches);
			alt_better_accepted_reads_reset.push_back(new_read);
		}

		if (split_i) {
			int ref_lh_realn_start = std::max(0, int(del_start-alt_consensus.length()));
			int ref_lh_realn_end, ref_rh_realn_start;
			int ref_rh_realn_end = std::min(int(del_end+alt_consensus.length()), contig_len);
			if (sv->len() < config.min_size_for_depth_filtering) {
				ref_lh_realn_end = ref_rh_realn_end;
				ref_rh_realn_start = ref_lh_realn_start;
			} else {
				ref_lh_realn_end = std::min(int(del_start+alt_consensus.length()), contig_len);
				ref_rh_realn_start = std::max(0, int(del_end-alt_consensus.length()));
			}

			StripedSmithWaterman::Alignment lh_aln, rh_aln;
			std::string lh_alt_consensus = alt_consensus.substr(0, split_i);
			std::string rh_alt_consensus = alt_consensus.substr(split_i);
			aligner.Align(lh_alt_consensus.c_str(), contig_seq+ref_lh_realn_start, ref_lh_realn_end-ref_lh_realn_start, filter, &lh_aln, 0);
			aligner.Align(rh_alt_consensus.c_str(), contig_seq+ref_rh_realn_start, ref_rh_realn_end-ref_rh_realn_start, filter, &rh_aln, 0);

			bam1_t* new_read = bam_init1();
			uint32_t* cigar = NULL;
			size_t n_cigar = 0;
			sam_parse_cigar(lh_aln.cigar_string.c_str(), NULL, &cigar, &n_cigar);
			bam_set1(new_read, 10, "ALT_PREFIX", 66, sam_hdr_name2tid(evidence_hdr, contig_name.c_str()),
					ref_lh_realn_start+lh_aln.ref_begin, 60, n_cigar, cigar, sam_hdr_name2tid(evidence_hdr, contig_name.c_str()),
					ref_rh_realn_start+rh_aln.ref_begin, rh_aln.ref_end-lh_aln.ref_begin, lh_alt_consensus.length(),
					lh_alt_consensus.c_str(), NULL, 0);
			alt_better_accepted_reads_reset.push_back(new_read);

			new_read = bam_init1();
			cigar = NULL;
			n_cigar = 0;
			sam_parse_cigar(rh_aln.cigar_string.c_str(), NULL, &cigar, &n_cigar);
			bam_set1(new_read, 10, "ALT_SUFFIX", 130, sam_hdr_name2tid(evidence_hdr, contig_name.c_str()),
					ref_rh_realn_start+rh_aln.ref_begin, 60, n_cigar, cigar, sam_hdr_name2tid(evidence_hdr, contig_name.c_str()),
					ref_lh_realn_start+lh_aln.ref_begin, lh_aln.ref_begin-rh_aln.ref_end, rh_alt_consensus.length(),
					rh_alt_consensus.c_str(), NULL, 0);
			alt_better_accepted_reads_reset.push_back(new_read);
		}

		std::sort(alt_better_accepted_reads_reset.begin(), alt_better_accepted_reads_reset.end(),
				[](bam1_t* b1, bam1_t* b2) {return b1->core.pos < b2->core.pos;});
		for (bam1_t* r : alt_better_accepted_reads_reset) {
			if (sam_write1(evidence_bam, evidence_hdr, r) < 0) {
				throw "Failed to write evidence to " + evidence_fname;
			}
		}
		sam_close(evidence_bam);

		std::string alt_consensus_fname = config.workdir + "/evidence/" + sv->id + ".fa";
		std::ofstream alt_consensus_fout(alt_consensus_fname);
		alt_consensus_fout << ">" << sv->id << std::endl;
		alt_consensus_fout << alt_consensus << std::endl;
		alt_consensus_fout.close();
	}

    delete[] alt_seq;
    delete[] ref_bp1_seq;
    delete[] ref_bp2_seq;

    bam_destroy1(read);
    sam_itr_destroy(iter);

    log_string = log_ss.str();
}

void genotype_dels(int id, std::string contig_name, char* contig_seq, int contig_len, std::vector<bcf1_t*> vcf_dels,
		bcf_hdr_t* in_vcf_header, bcf_hdr_t* out_vcf_header, stats_t stats, config_t config) {

    std::cout << "Genotyping deletions in " << contig_name << std::endl;

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, true);
    StripedSmithWaterman::Aligner harsh_aligner(1, 4, 100, 1, true);
    open_samFile_t* bam_file = open_samFile(config.bam_fname);
    hts_set_fai_filename(bam_file->file, config.reference_fname.c_str());

    std::vector<deletion_t*> dels, shorter_dels, longer_dels;
    for (bcf1_t* sv : vcf_dels) {
		deletion_t* del = vcf_record_to_deletion(sv, in_vcf_header);
		std::string log_string;
		genotype_del(del, contig_name, contig_seq, contig_len, bam_file, aligner, harsh_aligner, stats, config, log_string);
//		mtx.lock();
//		flog << log_string << std::endl << std::endl;
//		mtx.unlock();
		dels.push_back(del);
		if (del->is_short(config)) shorter_dels.push_back(del);
		else longer_dels.push_back(del);
    }

    calculate_confidence_interval_size(contig_name, shorter_dels, bam_file, config, stats);
	depth_filter_del(contig_name, contig_len, dels, bam_file);

	int min_disc_pairs = std::max(3, int(stats.median_depth+5)/10);

    for (deletion_t* sv : dels) {
        std::string filter;

        int del_len = sv->remapped_end-sv->remapped_start - sv->ins_seq.length();

        // ANOMALOUS COVERAGE - do not trust no matter what
        if (sv->left_flanking_cov > stats.max_depth || sv->right_flanking_cov > stats.max_depth ||
			sv->left_flanking_cov < stats.min_depth || sv->right_flanking_cov < stats.min_depth) {
        	filter += "ANOMALOUS_FLANKING_DEPTH;";
        	sv->anomalous_cov = true;
        }
		if (sv->indel_left_cov > stats.max_depth || sv->indel_right_cov > stats.max_depth) {
			filter += "ANOMALOUS_DEL_DEPTH;";
			sv->anomalous_cov = true;
		}

		// If positive SR-based GT, try SR-based filters
		if (gt_is_positive(sv->called_gt)) {
			if (sv->is_short(config)) { // size filter
				if (sv->called_gt == "HOM_ALT" && del_len > sv->max_conf_size) filter += "SIZE_FILTER;";
				else if (sv->called_gt == "HET" && del_len > sv->max_conf_size*2) filter += "SIZE_FILTER;";
			}
			if (sv->alt_better_strong < 3) filter += "LOW_ALT_STRONG_READS;";
			if (sv->remapped_end-sv->remapped_start < config.min_sv_size) filter += "TOO_SHORT_AFTER_REMAPPING;";
			if (!sv->split_aln_accepted) filter += "WEAK_SUPPORT_BY_ALT_CONTIG;";
			int size_diff = abs((sv->end-sv->start) - (sv->remapped_end-sv->remapped_start));
			if (size_diff >= config.max_size_diff) filter += "REMAPPED_SIZE_DIFF;";
			if (sv->was_remapped && (sv->vouching_for_lh_consensus_remapping < 3 || sv->vouching_for_rh_consensus_remapping < 3))
				filter += "NO_READS_SUPPORTING_REMAPPING;";
		}
        if (sv->alt_better > 2*stats.max_depth) filter += "TOO_MANY_ALT_READS;";
        sv->filter = filter;
    }

	calculate_ptn_ratio(contig_name, dels, bam_file, config, stats);

	for (deletion_t* sv : dels) {
		std::string filter = sv->filter;
        if (gt_is_positive(sv->called_gt) && sv->disc_pairs < 3 && !sv->is_short(config)) {
        	filter += "NOT_ENOUGH_DISC_PAIRS;";
        }

        // RESCUE DELETIONS THAT DO NOT HAVE SPLIT READS BUT ARE HIGHLY SUPPORTED BY DISCORDANT PAIRS
        bool dp_only = false;
        if (!sv->anomalous_cov && (!gt_is_positive(sv->called_gt) || !filter.empty()) && sv->disc_pairs >= min_disc_pairs) {
        	if (sv->mleft_flanking_cov*0.25 >= sv->mindel_left_cov && sv->mright_flanking_cov*0.25 >= sv->mindel_right_cov) {
        		sv->called_gt = "HOM_ALT";
        		filter = "";
        	} else {
        		sv->called_gt = "HET";
        		filter = "";
        	}
        	dp_only = true;
        	if (gt_is_positive(sv->called_gt) && sv->disc_pairs+sv->conc_pairs < stats.min_pairs_crossing) filter += "NOT_ENOUGH_PAIRS;";
        }

        if (!gt_is_positive(sv->called_gt) && !sv->alt_better && !sv->ref_better) filter += "NO_USEFUL_READS;";

        // force depth filter for dp_only variants
        if (gt_is_positive(sv->called_gt) && (sv->end-sv->start >= config.min_size_for_depth_filtering || dp_only)) {
        	if (sv->mleft_flanking_cov*0.75<sv->mindel_left_cov || sv->mright_flanking_cov*0.75<sv->mindel_right_cov) {
				filter += "DEPTH_FILTER;";
        	}
		}

        if (filter.empty()) filter = "PASS";

        bcf_subset(out_vcf_header, sv->bcf_entry, 0, {});
        int gt_data[2];
        strgt2datagt(sv->called_gt, gt_data);
        bcf_update_genotypes(out_vcf_header, sv->bcf_entry, gt_data, 2);

        const char* ft_val[1];
        ft_val[0] = filter.c_str();

        const char* ed_val[1];
        ed_val[0] = dp_only ? "DP_ONLY" : "SR";

        int remapped_start_1based = sv->remapped_start+1, remapped_end_1based = sv->remapped_end+1;
        bcf_update_format_string(out_vcf_header, sv->bcf_entry, "FT", ft_val, 1);
        bcf_update_format_string(out_vcf_header, sv->bcf_entry, "ED", ed_val, 1);
        bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "AR", &(sv->alt_better), 1);
        bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "ARF", &(sv->alt_better_fwd), 1);
        bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "ARR", &(sv->alt_better_rev), 1);
        bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "RR", &(sv->ref_better), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "SR", &(sv->alt_better_strong), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "ARMQ", &(sv->alt_better_max_mq), 1);
		if (sv->was_remapped) {
			bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "RS", &remapped_start_1based, 1);
			bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "RE", &remapped_end_1based, 1);
			const char* rsc_val[1];
			rsc_val[0] = sv->rs_cigar.c_str();
			bcf_update_format_string(out_vcf_header, sv->bcf_entry, "RSC", rsc_val, 1);
			const char* rec_val[1];
			rec_val[0] = sv->re_cigar.c_str();
			bcf_update_format_string(out_vcf_header, sv->bcf_entry, "REC", rec_val, 1);
			const char* rfc_val[1];
			rfc_val[0] = sv->full_cigar.c_str();
			bcf_update_format_string(out_vcf_header, sv->bcf_entry, "RFC", rfc_val, 1);
			bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "VL", &(sv->vouching_for_lh_consensus_remapping), 1);
			bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "VR", &(sv->vouching_for_rh_consensus_remapping), 1);
		}

        bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "DFL", &(sv->left_flanking_cov), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "DDL", &(sv->indel_left_cov), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "DDR", &(sv->indel_right_cov), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "DFR", &(sv->right_flanking_cov), 1);

		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "MDFL", &(sv->mleft_flanking_cov), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "MDDL", &(sv->mindel_left_cov), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "MDDR", &(sv->mindel_right_cov), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "MDFR", &(sv->mright_flanking_cov), 1);

		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "MCS", &(sv->max_conf_size), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "DP", &(sv->disc_pairs), 1);
		bcf_update_format_int32(out_vcf_header, sv->bcf_entry, "CP", &(sv->conc_pairs), 1);
    }

    close_samFile(bam_file);
}


#endif /* GENOTYPE_DELETIONS_H_ */
