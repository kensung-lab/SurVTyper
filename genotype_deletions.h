#ifndef GENOTYPE_DELETIONS_H_
#define GENOTYPE_DELETIONS_H_

extern config_t config;

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

void genotype_del(std::string& contig_name, deletion_t* sv, open_samFile_t* bam_file, StripedSmithWaterman::Aligner& aligner,
		StripedSmithWaterman::Aligner& harsh_aligner, stats_t stats, chr_seqs_map_t& chr_seqs, std::string& log_string) {

	std::stringstream log_ss;

    int del_start = sv->start, del_end = sv->end;
    int contig_len = chr_seqs.get_len(contig_name);

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
    strncpy(alt_seq, chr_seqs.get_seq(contig_name)+alt_start, alt_lh_len);
    strncpy(alt_seq+alt_lh_len, sv->ins_seq.c_str(), sv->ins_seq.length());
    strncpy(alt_seq+alt_lh_len+sv->ins_seq.length(), chr_seqs.get_seq(contig_name)+del_end, alt_rh_len);
    alt_seq[alt_len] = 0;

    // extract ref breakpoint allele(s) - will be useful for consensus generation
    int ref_bp1_start = alt_start, ref_bp1_end = std::min(del_start+extend, contig_len);
    int ref_bp1_len = ref_bp1_end - ref_bp1_start;
    char* ref_bp1_seq = new char[ref_bp1_len + 1];
    strncpy(ref_bp1_seq, chr_seqs.get_seq(contig_name)+ref_bp1_start, ref_bp1_len);
    ref_bp1_seq[ref_bp1_len] = 0;

    int ref_bp2_start = std::max(0, del_end-extend), ref_bp2_end = alt_end;
    int ref_bp2_len = ref_bp2_end - ref_bp2_start;
    char* ref_bp2_seq = new char[ref_bp2_len + 1];
    strncpy(ref_bp2_seq, chr_seqs.get_seq(contig_name)+ref_bp2_start, ref_bp2_len);
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
    std::vector<std::string> read_seqs, qnames;
    std::vector<bool> good_read, read_fwd;
    std::unordered_set<std::string> ref_better_qnames;
    std::vector<int> mapq;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;

        // do not consider reads that do not intersect the breakpoints
        if (get_unclipped_end(read) <= del_start) continue;
        if (del_start <= get_unclipped_start(read) && get_unclipped_end(read) <= del_end) {
        	ref_better_qnames.insert(bam_get_qname(read));
        	continue;
        }
        if (get_unclipped_start(read) >= del_end) continue;

        read_seqs.push_back(get_sequence(read));
        qnames.push_back(bam_get_qname(read));
        good_read.push_back(is_samechr(read) && !is_samestr(read) && !is_mate_unmapped(read));
        read_fwd.push_back(!bam_is_rev(read));
        mapq.push_back(read->core.qual);
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
    int i = 0;
    for (std::string& seq : read_seqs) {
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
            alt_better_qnames.push_back(qnames[i]);
            alt_better_seqs_start.push_back(alt_aln.ref_begin - get_left_clip_size(alt_aln));
            alt_better_seqs_end.push_back(alt_aln.ref_end + get_right_clip_size(alt_aln));

            alt_better_original_alt_aln.push_back(alt_aln.cigar_string);
            alt_better_original_ref_aln.push_back(ref_aln.cigar_string);
            alt_better_score_diff.push_back(alt_aln.sw_score-ref_aln.sw_score);
            alt_better_mapq.push_back(mapq[i]);
            alt_better_good_read.push_back(good_read[i]);
        }
        else if (alt_aln.sw_score < ref_aln.sw_score && accept_aln(ref_aln, config.min_clip_size)) {
        	ref_better++;
        	ref_better_qnames.insert(qnames[i]);
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

	i = 0;
	for (std::string& seq : alt_better_seqs) {
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
			}
		}
		i++;
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
    if (gt_is_positive(sv->called_gt) && sv->len() < config.min_size_for_depth_filtering) {
		StripedSmithWaterman::Alignment alt_cons_realn;
		std::string padded_alt_consensus = std::string(config.clip_penalty, 'N') + alt_consensus + std::string(config.clip_penalty, 'N');
		int ref_realn_start = del_start-padded_alt_consensus.length(), ref_realn_end = del_end+padded_alt_consensus.length();
		if (ref_realn_start < 0) ref_realn_start = 0;
		if (ref_realn_end >= contig_len) ref_realn_end = contig_len-1;
		int ref_realn_len = ref_realn_end - ref_realn_start;
		aligner.Align(padded_alt_consensus.c_str(), chr_seqs.get_seq(contig_name)+ref_realn_start, ref_realn_len, filter, &alt_cons_realn, 0);
		if (alt_cons_realn.query_begin > 0 && alt_cons_realn.query_end < alt_consensus.length()-1) sv->lowq_alt_allele = true;

		char ref_cstr[100000];
		strncpy(ref_cstr, chr_seqs.get_seq(contig_name)+ref_realn_start, ref_realn_len);
		for (int i = 0; i < ref_realn_len; i++) {
			ref_cstr[i] = toupper(ref_cstr[i]);
		}
		int* prefix_scores = smith_waterman_gotoh(ref_cstr, ref_realn_len, padded_alt_consensus.c_str(), padded_alt_consensus.length(), 1, -4, -6, -1);

		for (int i = 0; i < ref_realn_len/2; i++) {
			std::swap(ref_cstr[i], ref_cstr[ref_realn_len-i-1]);
		}
		std::string padded_alt_consensus_rev = std::string(padded_alt_consensus.rbegin(), padded_alt_consensus.rend());
		int* suffix_scores = smith_waterman_gotoh(ref_cstr, ref_realn_len, padded_alt_consensus_rev.c_str(), padded_alt_consensus_rev.length(), 1, -4, -6, -1);

		int max_score = 0, split_i = 0;
		for (int i = 0; i < alt_consensus.length(); i++) {
			int prefix_score = prefix_scores[config.clip_penalty+i-1], suffix_score = suffix_scores[alt_consensus.length()-i-1+config.clip_penalty];
			if (prefix_score + suffix_score > max_score) {
				max_score = prefix_score + suffix_score;
				split_i = i;
			}
		}

		StripedSmithWaterman::Alignment lh_aln, rh_aln;
		char* padded_alt_consensus_cstr = new char[padded_alt_consensus.length()+1];
		strcpy(padded_alt_consensus_cstr, padded_alt_consensus.c_str());
		char terminator = '\0';
		std::swap(terminator, padded_alt_consensus_cstr[config.clip_penalty+split_i]);
		aligner.Align(padded_alt_consensus_cstr, chr_seqs.get_seq(contig_name)+ref_realn_start, ref_realn_len, filter, &lh_aln, 0);
		std::swap(terminator, padded_alt_consensus_cstr[config.clip_penalty+split_i]);

		aligner.Align(padded_alt_consensus_cstr+config.clip_penalty+split_i, chr_seqs.get_seq(contig_name)+ref_realn_start, ref_realn_len, filter, &rh_aln, 0);

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

		delete[] prefix_scores;
		delete[] suffix_scores;
    }

    delete[] alt_seq;
    delete[] ref_bp1_seq;
    delete[] ref_bp2_seq;

    bam_destroy1(read);
    sam_itr_destroy(iter);

    log_string = log_ss.str();
}


#endif /* GENOTYPE_DELETIONS_H_ */
