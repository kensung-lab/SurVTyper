#ifndef VCF_UTILS_H_
#define VCF_UTILS_H_


deletion_t* vcf_record_to_deletion(bcf1_t* vcf_del, bcf_hdr_t* vcf_header) {
	bcf_unpack(vcf_del, BCF_UN_INFO);
    deletion_t* del = new deletion_t(vcf_del->d.id, vcf_del->pos, get_sv_end(vcf_del, vcf_header), get_ins_seq(vcf_del, vcf_header), "");
    del->bcf_entry = vcf_del;
    return del;
}

bcf_hdr_t* generate_out_hdr(std::string& sample_name, bcf_hdr_t* vcf_header) {
	// remove existing FORMAT lines - we are going to replace them with ours
	bcf_hdr_t* new_header = bcf_hdr_subset(vcf_header, 0, {}, {});
	bcf_hdr_add_sample(new_header, sample_name.c_str());

	bcf_hdr_remove(new_header, BCF_HL_FMT, NULL);
	int len;

	const char* gt_tag = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, gt_tag, &len));

	const char* ft_tag = "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Filter. PASS indicates a reliable call.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, ft_tag, &len));

	const char* ed_tag = "##FORMAT=<ID=ED,Number=1,Type=String,Description=\"Evidence used to call the deletion. "
			"Most reliable calls require SR (split reads). DP_ONLY means that the call was made using discordant pairs only.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, ed_tag, &len));

	const char* ar_tag = "##FORMAT=<ID=AR,Number=1,Type=Integer,Description=\"Number of good quality reads that support the alternative allele.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, ar_tag, &len));

	const char* armq_tag = "##FORMAT=<ID=ARMQ,Number=1,Type=Integer,Description=\"Maximum MAPQ among strong reads.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, armq_tag, &len));

	const char* arf_tag = "##FORMAT=<ID=ARF,Number=1,Type=Integer,Description=\"Number of good quality forward reads that support the alternative allele.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, arf_tag, &len));

	const char* arr_tag = "##FORMAT=<ID=ARR,Number=1,Type=Integer,Description=\"Number of good quality reverse reads that support the alternative allele.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, arr_tag, &len));

	const char* sr_tag = "##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Number of reads that strongly support the alternative allele.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, sr_tag, &len));

	const char* ar1_tag = "##FORMAT=<ID=AR1,Number=1,Type=Integer,Description=\"Number of good quality reads that support the alternative allele.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, ar1_tag, &len));

	const char* sr1_tag = "##FORMAT=<ID=SR1,Number=1,Type=Integer,Description=\"Number of reads that strongly support the alternative allele.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, sr1_tag, &len));

	const char* ar2_tag = "##FORMAT=<ID=AR2,Number=1,Type=Integer,Description=\"Number of good quality reads that support the alternative allele.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, ar2_tag, &len));

	const char* sr2_tag = "##FORMAT=<ID=SR2,Number=1,Type=Integer,Description=\"Number of reads that strongly support the alternative allele.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, sr2_tag, &len));

	const char* rr_tag = "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"Number of reads that support the reference allele.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, rr_tag, &len));

	const char* rs_tag = "##FORMAT=<ID=RS,Number=1,Type=Integer,Description=\"Start of deletion after remapping.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, rs_tag, &len));

	const char* re_tag = "##FORMAT=<ID=RE,Number=1,Type=Integer,Description=\"End of deletion after remapping.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, re_tag, &len));

	const char* rsc_tag = "##FORMAT=<ID=RSC,Number=1,Type=String,Description=\"\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, rsc_tag, &len));

	const char* rec_tag = "##FORMAT=<ID=REC,Number=1,Type=String,Description=\"\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, rec_tag, &len));

	const char* rfc_tag = "##FORMAT=<ID=RFC,Number=1,Type=String,Description=\"\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, rfc_tag, &len));

	const char* cn_tag = "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, cn_tag, &len));

	const char* ril_tag = "##FORMAT=<ID=RIL,Number=1,Type=Integer,Description=\"Remapped insertion length.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, ril_tag, &len));

	const char* rcg_tag = "##FORMAT=<ID=RCG,Number=.,Type=String,Description=\"Remapped cigar.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, rcg_tag, &len));

	// will remove later
	const char* wr_tag = "##FORMAT=<ID=WR,Number=1,Type=Integer,Description=\"Was remapped.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, wr_tag, &len));

	const char* vl_tag = "##FORMAT=<ID=VL,Number=1,Type=Integer,Description=\"Number of reads vouching for the remapping of the left half of the consensus alternative allele.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, vl_tag, &len));

	const char* vr_tag = "##FORMAT=<ID=VR,Number=1,Type=Integer,Description=\"Number of reads vouching for the remapping of the right half of the consensus alternative allele.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, vr_tag, &len));

	const char* dfl_tag = "##FORMAT=<ID=DFL,Number=1,Type=Integer,Description=\"Average depth of left flanking region.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, dfl_tag, &len));

	const char* ddl_tag = "##FORMAT=<ID=DDL,Number=1,Type=Integer,Description=\"Average depth of left portion of the duplication.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, ddl_tag, &len));

	const char* ddr_tag = "##FORMAT=<ID=DDR,Number=1,Type=Integer,Description=\"Average depth of right portion of the duplication.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, ddr_tag, &len));

	const char* dfr_tag = "##FORMAT=<ID=DFR,Number=1,Type=Integer,Description=\"Average depth of right flanking region.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, dfr_tag, &len));

	const char* mdfl_tag = "##FORMAT=<ID=MDFL,Number=1,Type=Integer,Description=\"Median depth of left flanking region.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, mdfl_tag, &len));

	const char* mddl_tag = "##FORMAT=<ID=MDDL,Number=1,Type=Integer,Description=\"Median depth of left portion of the duplication.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, mddl_tag, &len));

	const char* mddr_tag = "##FORMAT=<ID=MDDR,Number=1,Type=Integer,Description=\"Median depth of right portion of the duplication.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, mddr_tag, &len));

	const char* mdfr_tag = "##FORMAT=<ID=MDFR,Number=1,Type=Integer,Description=\"Median depth of right flanking region.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, mdfr_tag, &len));

	const char* cid_tag = "##FORMAT=<ID=CID,Number=1,Type=Integer,Description=\"Completely inside duplications.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, cid_tag, &len));

	const char* ow_tag = "##FORMAT=<ID=OW,Number=1,Type=Integer,Description=\"Outward pairs.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, ow_tag, &len));

	const char* mcs_tag = "##FORMAT=<ID=MCS,Number=1,Type=Integer,Description=\"Maximum size.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, mcs_tag, &len));

	const char* cp_tag = "##FORMAT=<ID=CP,Number=1,Type=Integer,Description=\"Concordant pairs.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, cp_tag, &len));

	const char* dp_tag = "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Discordant pairs.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, dp_tag, &len));

	const char* dpfs_tag = "##FORMAT=<ID=DPFS,Number=1,Type=Integer,Description=\"Discordant pairs with stable mate on forward strand.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, dpfs_tag, &len));

	const char* dprs_tag = "##FORMAT=<ID=DPRS,Number=1,Type=Integer,Description=\"Discordant pairs with stable mate on reverse strand.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, dprs_tag, &len));

	const char* mindp_tag = "##FORMAT=<ID=mDP,Number=1,Type=Integer,Description=\"Minimum number of discordant pairs.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, mindp_tag, &len));

	const char* maxdp_tag = "##FORMAT=<ID=MDP,Number=1,Type=Integer,Description=\"Maximum number of discordant pairs.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, maxdp_tag, &len));

	const char* sdp_tag = "##FORMAT=<ID=SDP,Number=1,Type=Integer,Description=\"Strong discordant pairs.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, sdp_tag, &len));

	const char* fwd_anchor_tag = "##FORMAT=<ID=FAN,Number=1,Type=String,Description=\"Forward anchor.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, fwd_anchor_tag, &len));

	const char* rev_anchor_tag = "##FORMAT=<ID=RAN,Number=1,Type=String,Description=\"Reverse anchor.\">";
	bcf_hdr_add_hrec(new_header, bcf_hdr_parse_line(new_header, rev_anchor_tag, &len));

	return new_header;
}

#endif /* VCF_UTILS_H_ */
