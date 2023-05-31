#ifndef GENOTYPE_INSERTIONS_H_
#define GENOTYPE_INSERTIONS_H_

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

void genotype_small_ins(std::string& contig_name, insertion_t* sv, char* contig_seq, int contig_len, open_samFile_t* bam_file,
		StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner,
		std::unordered_map<std::string, std::string>& mateseqs_map, stats_t stats) {

	std::stringstream log_ss;

	int ins_start = sv->start, ins_end = sv->end;

	int extend = config.read_len + 20;

	// build alt allele
	/*
	 * POS in VCF is the base BEFORE the insertion
	 * END seems to be the base BEFORE the reference resumes - i.e., for a "clean" insertion (no deletion),POS == END, otherwise the last base deleted
	 * As usual, in order to make intervals [ ), we increase the coordinates by 1
	 */
	ins_start++; ins_end++;

	int alt_start = std::max(1, ins_start-extend);
	int alt_end = std::min(ins_end+extend, contig_len);
	int alt_lf_len = ins_start-alt_start, alt_rf_len = alt_end-ins_end;
	int alt_len = alt_lf_len + sv->ins_seq.length() + alt_rf_len;
	char* alt_seq = new char[alt_len + 1];
	strncpy(alt_seq, contig_seq+alt_start, alt_lf_len);
	strncpy(alt_seq+alt_lf_len, sv->ins_seq.c_str(), sv->ins_seq.length());
	strncpy(alt_seq+alt_lf_len+sv->ins_seq.length(), contig_seq+ins_end, alt_rf_len);
    alt_seq[alt_len] = '\0';

    int ref_bp1_start = alt_start, ref_bp1_end = std::min(ins_start+extend, contig_len);
    int ref_bp1_len = ref_bp1_end - ref_bp1_start;
    int ref_bp2_start = std::max(0, ins_end-extend), ref_bp2_end = alt_end;
    int ref_bp2_len = ref_bp2_end - ref_bp2_start;

    int ref_len = ref_bp2_end - ref_bp1_start;
    char* ref_seq = new char[ref_len + 1];
	strncpy(ref_seq, contig_seq+ref_bp1_start, ref_bp2_end-ref_bp1_start);
	ref_seq[ref_len] = '\0';

	char l_region[100], r_region[100];
	sprintf(l_region, "%s:%d-%d", contig_name.c_str(), alt_start, ref_bp1_end);
	sprintf(r_region, "%s:%d-%d", contig_name.c_str(), ref_bp2_start, alt_end);
	char* regions[] = {l_region, r_region};

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);
    bam1_t* read = bam_init1();

    std::vector<bam1_t*> reads_to_remap;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read)) continue;

		if ((!bam_is_rev(read) && read->core.pos < ins_start && is_dc_pair(read)) ||
			 (bam_is_rev(read) && bam_endpos(read) > ins_end && is_dc_pair(read))) {
			std::string mate_qname = bam_get_qname(read);
			if (read->core.tid == read->core.mtid) {
				mate_qname += (read->core.flag & BAM_FREAD1 ? "_2" : "_1");
			}
			std::string mate_seq = mateseqs_map[mate_qname];
			if (!bam_is_rev(read)) rc(mate_seq);

			int flag = BAM_FUNMAP; // so that htslib stops complaining about no cigar
			int mapq = read->core.qual;
			std::string qual(mate_seq.length(), '?');

			bam1_t* new_read = bam_init1();
			bam_set1(new_read, mate_qname.length(), mate_qname.c_str(), flag, read->core.mtid,
					read->core.mpos, mapq, 0, NULL, read->core.tid, read->core.pos, 0, mate_seq.length(),
					mate_seq.c_str(), qual.c_str(), 0);
			reads_to_remap.push_back(new_read);
		}

		// do not consider reads that do not intersect the insertion
		if (get_unclipped_end(read) <= ins_start) continue;
		if (get_unclipped_start(read) >= ins_end) continue;

		reads_to_remap.push_back(bam_dup1(read));
    }

    StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment alt_aln, ref1_aln, ref2_aln;
	std::vector<std::string> alt_better_seqs;
    std::vector<bam1_t*> ref_better_reads, alt_better_reads;
	std::vector<bool> alt_better_seqs_lc, alt_better_seqs_rc; // TODO: populate with actual values
	std::vector<bool> alt_strong;
	log_ss << sv->id << std::endl;
	int alt_better = 0, ref_better = 0, same = 0;
	int i = 0;
	for (int i = 0; i < reads_to_remap.size(); i++) {
		std::string seq = get_sequence(reads_to_remap[i]);
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
			alt_better_reads.push_back(reads_to_remap[i]);
			alt_better_seqs_lc.push_back(false);
			alt_better_seqs_rc.push_back(false);
			alt_strong.push_back(alt_aln.sw_score-ref_aln.sw_score >= config.min_score_diff || !accept_aln_strict(ref_aln));
		} else if (alt_aln.sw_score < ref_aln.sw_score && accept_aln_strict(ref_aln)) {
			ref_better++;
			ref_better_reads.push_back(reads_to_remap[i]);
		} else if (alt_aln.sw_score == ref_aln.sw_score && accept_aln_strict(alt_aln) && accept_aln_strict(ref_aln)) same++;

		log_ss << seq << " " << alt_aln.cigar_string << " " << ref_aln.cigar_string << " " << alt_aln.sw_score << " " << ref_aln.sw_score << std::endl;
	}
	sv->ref_better = ref_better;

	auto gt1 = calculate_small_ins_genotype(alt_better, ref_better);
	sv->before_read_qc_gt = gt1.first;

	std::string alt_consensus;
	std::vector<bam1_t*> alt_better_accepted_reads;
	std::vector<StripedSmithWaterman::Alignment> alt_better_accepted_alns;
	StripedSmithWaterman::Alignment consensus_alt_aln, diff_consensus_alt_aln;
	if (alt_better > 0 && alt_better <= 4*stats.max_depth && ref_better <= 2*stats.max_depth) {
		std::vector<StripedSmithWaterman::Alignment> consensus_contigs_alns;
		std::string consensus_log;
		alt_consensus = generate_consensus(alt_better_seqs, alt_better_seqs_lc, alt_better_seqs_rc,
				aligner, harsh_aligner, consensus_contigs_alns, consensus_log);

		log_ss << "ALT CONSENSUS: " << alt_consensus << std::endl;

		int i = 0;
		int alt_better = 0, diff_alt_better = 0;
		for (int i = 0; i < alt_better_reads.size(); i++) {
			std::string seq = alt_better_seqs[i];
			std::string padded_seq = std::string(config.clip_penalty, 'N') + seq + std::string(config.clip_penalty, 'N');
			harsh_aligner.Align(padded_seq.c_str(), alt_consensus.c_str(), alt_consensus.length(), filter, &consensus_alt_aln, 0);

			if (get_left_clip_size(consensus_alt_aln) <= config.clip_penalty && get_right_clip_size(consensus_alt_aln) <= config.clip_penalty) {
				int mismatches = consensus_alt_aln.mismatches - 2*config.clip_penalty + get_left_clip_size(consensus_alt_aln) +
						get_right_clip_size(consensus_alt_aln);

				if (mismatches <= seq.length()*config.max_seq_error) {
					sv->alt_better++;
					if (alt_strong[i]) sv->alt_better_strong++;

					harsh_aligner.Align(seq.c_str(), alt_consensus.c_str(), alt_consensus.length(), filter, &consensus_alt_aln, 0);
					alt_better_accepted_reads.push_back(alt_better_reads[i]);
					alt_better_accepted_alns.push_back(consensus_alt_aln);
				}
			}
		}

		char ref_cstr[100000];
		int ref_realn_start = ins_start - 2*alt_consensus.length(), ref_realn_end = ins_end + 2*alt_consensus.length();
		if (ref_realn_start < 0) ref_realn_start = 0;
		if (ref_realn_end >= contig_len) ref_realn_end = contig_len-1;
		int ref_realn_len = ref_realn_end - ref_realn_start;
		strncpy(ref_cstr, contig_seq+ref_realn_start, ref_realn_len);
		for (int i = 0; i < ref_realn_len; i++) {
			ref_cstr[i] = toupper(ref_cstr[i]);
		}

		remap_consensus_to_reference(alt_consensus, ref_cstr, ref_realn_len, aligner, sv->remapped_cigar, sv->remapped_ins_len, sv->lowq_alt_allele);
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
			StripedSmithWaterman::Alignment& aln = alt_better_accepted_alns[i];
			uint32_t* cigar = NULL;
			size_t n_cigar = 0;
			sam_parse_cigar(aln.cigar_string.c_str(), NULL, &cigar, &n_cigar);

			bam1_t* new_read = bam_init1();
			bam_set1(new_read, strlen(bam_get_qname(read)), bam_get_qname(read), 0, sam_hdr_name2tid(evidence_hdr, "ALT"),
					aln.ref_begin, read->core.qual, n_cigar, cigar, read->core.mtid, read->core.mpos, 0, read->core.l_qseq,
					get_sequence(read).c_str(), (const char*) bam_get_qual(read), 0);
			bam_aux_update_int(read, "AS", aln.sw_score);
			bam_aux_update_int(read, "NM", aln.mismatches);
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

		for (bam1_t* r : alt_better_accepted_reads_reset) {
			bam_destroy1(r);
		}

		std::string alt_consensus_fname = config.workdir + "/evidence/" + sv->id + ".fa";
		std::ofstream alt_consensus_fout(alt_consensus_fname);
		alt_consensus_fout << ">" << sv->id << std::endl;
		alt_consensus_fout << alt_consensus << std::endl;
		alt_consensus_fout.close();
	}

	delete[] ref_seq;

	auto gt2 = calculate_small_ins_genotype(sv->alt_better, sv->ref_better);
	sv->called_gt = gt2.first, sv->ratio = gt2.second;

//	mtx.lock();
//	flog << log_ss.str() << std::endl;
//	mtx.unlock();
}


#endif /* GENOTYPE_INSERTIONS_H_ */
