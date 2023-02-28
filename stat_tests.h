#ifndef SURVINDEL2_STAT_TESTS_H
#define SURVINDEL2_STAT_TESTS_H

#include <set>
#include <cmath>

#include "htslib/sam.h"
#include "sam_utils.h"
#include "ks-test.h"

//struct read_w_cached_info_t {
//    bam1_t* read;
//    hts_pos_t start, end, isize;
//    int64_t as;
//    int aln_len;
//    bool left_clipped, right_clipped, is_rev, is_mrev;
//    int references = 0;
//
//    read_w_cached_info_t() {}
//    read_w_cached_info_t(bam1_t* read) : read(bam_dup1(read)), as(get_AS_tag(read)), start(read->core.pos), end(bam_endpos(read)),
//                                         aln_len(get_aligned_portion_len(read)), left_clipped(is_left_clipped(read, 0)),
//                                         right_clipped(is_right_clipped(read, 0)), is_rev(bam_is_rev(read)), is_mrev(bam_is_mrev(read)),
//                                         isize(read->core.isize) {}
//};
//
//
//int cleanup_bin_size(int read_len) {
//    return 200;
//}
//
//std::vector<read_w_cached_info_t*> cleanup_reads(std::vector<read_w_cached_info_t*>& reads_w_as, hts_pos_t region_start,
//                                                hts_pos_t region_end, hts_pos_t original_dup_left_bp, hts_pos_t original_dup_right_bp,
//                                                stats_t stats, config_t config, std::vector<double>& max_allowed_frac_normalized) {
//    int BIN_SIZE = cleanup_bin_size(config.read_len);
//    hts_pos_t region_len = region_end - region_start;
//    int n_bins = (region_len-1)/BIN_SIZE + 1;
//    double desired_depth = stats.min_depth;
//
//    // compute coverage (in term of bases) that we want to achieve for each bin
//    std::vector<uint32_t> bin_desired_cov_bases(n_bins);
//    for (int i = 0; i < n_bins-1; i++) {
//        bin_desired_cov_bases[i] = desired_depth*BIN_SIZE;
//    }
//    bin_desired_cov_bases[n_bins-1] = desired_depth*(region_len%BIN_SIZE == 0 ? BIN_SIZE : region_len%BIN_SIZE);
//
//    std::random_shuffle(reads_w_as.begin(), reads_w_as.end());
//    std::sort(reads_w_as.begin(), reads_w_as.end(), [](read_w_cached_info_t* r1, read_w_cached_info_t* r2) {
//        return r1->aln_len-r1->as < r2->aln_len-r2->as;
//    });
//
//    auto choose_best_bin = [&region_len, &n_bins](hts_pos_t rs_in_region, hts_pos_t re_in_region, int bin_size) {
//        int rs_bin = rs_in_region/bin_size, re_bin = re_in_region/bin_size;
//        int best_overlap = 0, best_bin;
//        for (int i = rs_bin; i <= re_bin && i < n_bins; i++) {
//            int bin_overlap = overlap(rs_in_region, re_in_region, i*bin_size,
//                                      std::min(hts_pos_t(i+1)*bin_size, region_len));
//            if (bin_overlap > best_overlap) {
//                best_overlap = bin_overlap;
//                best_bin = i;
//            }
//        }
//        return best_bin;
//    };
//
//    std::vector<read_w_cached_info_t*> solid_reads_w_as;
//    std::vector<int> region_coverage(region_len+1);
//    for (read_w_cached_info_t* read_w_as : reads_w_as) {
//        hts_pos_t rs = read_w_as->start, re = read_w_as->end;
//        if (read_w_as->left_clipped && abs(rs - original_dup_left_bp) >= 5) continue;
//        if (read_w_as->right_clipped && abs(re - original_dup_right_bp) >= 5) continue;
//
//        hts_pos_t rs_in_region = std::max(hts_pos_t(0), rs-region_start), re_in_region = std::min(re-region_start, region_len);
//        if (re_in_region <= 0 || rs_in_region >= region_len) continue;
//
//        bool is_solid = false;
//        for (int i = rs_in_region; i <= re_in_region; i++) {
//            if (region_coverage[i] < desired_depth) {
//                is_solid = true;
//                break;
//            }
//        }
//
//        if (is_solid) {
//            for (int i = rs_in_region; i <= re_in_region; i++) {
//                region_coverage[i]++;
//            }
//            solid_reads_w_as.push_back(read_w_as);
//        }
//    }
//
//    std::vector<std::vector<int64_t>> solid_as_diff_dists(n_bins);
//    std::vector<std::vector<read_w_cached_info_t*> > binned_reads(n_bins);
//    for (read_w_cached_info_t* read_w_as : solid_reads_w_as) {
//        hts_pos_t rs = read_w_as->start, re = read_w_as->end;
//        hts_pos_t rs_in_region = std::max(hts_pos_t(0), rs-region_start), re_in_region = std::min(re-region_start, region_len);
//        int best_bin = choose_best_bin(rs_in_region, re_in_region, BIN_SIZE);
//        solid_as_diff_dists[best_bin].push_back(read_w_as->aln_len-read_w_as->as);
//    }
//    for (read_w_cached_info_t* read_w_as : reads_w_as) {
//        hts_pos_t rs = read_w_as->start, re = read_w_as->end;
//        hts_pos_t rs_in_region = std::max(hts_pos_t(0), rs-region_start), re_in_region = std::min(re-region_start, region_len);
//        if (re_in_region <= 0 || rs_in_region >= region_len) continue;
//        int best_bin = choose_best_bin(rs_in_region, re_in_region, BIN_SIZE);
//        binned_reads[best_bin].push_back(read_w_as);
//    }
//
//    std::vector<read_w_cached_info_t*> kept_reads;
//    for (int i = 0; i < n_bins; i++) {
//        std::vector<int64_t>& bin_as_dists = solid_as_diff_dists[i];
//        int n = bin_as_dists.size();
//        if (n < 3) continue;
//
//        std::sort(bin_as_dists.begin(), bin_as_dists.end());
//        int64_t median_as_diff = bin_as_dists[n/2];
//
//        int median_or_better = 0;
//        for (read_w_cached_info_t* read_w_as : binned_reads[i]) {
//            if (read_w_as->aln_len-read_w_as->as <= median_as_diff) {
//                median_or_better++;
//            }
//        }
//
//        std::vector<int> as_diff_count(max_allowed_frac_normalized.size());
//        for (read_w_cached_info_t* read_w_as : binned_reads[i]) {
//            if (read_w_as->left_clipped && abs(read_w_as->start - original_dup_left_bp) >= 5) continue;
//            if (read_w_as->right_clipped && abs(read_w_as->end - original_dup_right_bp) >= 5) continue;
//
//            int64_t as_diff = read_w_as->aln_len - read_w_as->as;
//            int64_t shifted_as_diff = std::max(int64_t(0), as_diff-median_as_diff);
//            if (shifted_as_diff >= max_allowed_frac_normalized.size()) continue;
//            if (as_diff_count[shifted_as_diff]+1 > max_allowed_frac_normalized[shifted_as_diff] * median_or_better) continue;
//
//            kept_reads.push_back(read_w_as);
//            as_diff_count[shifted_as_diff]++;
//        }
//    }
//
//    return kept_reads;
//}
//
//void compute_cleaned_up_depth(duplication_t* dup, open_samFile_t* bam_file,
//       std::vector<read_w_cached_info_t*>& lf_reads, std::vector<read_w_cached_info_t*>& ldup_reads,
//       std::vector<read_w_cached_info_t*>& rdup_reads, std::vector<read_w_cached_info_t*>& rf_reads,
//       stats_t& stats, config_t& config, std::vector<double>& max_allowed_frac_normalized, std::string workdir) {
//    hts_pos_t left_flanking_start = dup->start - 5000, left_flanking_end = dup->start;
//    hts_pos_t dup_left_start = dup->start, dup_left_end = std::min(dup->end, dup->start+5000);
//    hts_pos_t dup_right_start = std::max(dup->start, dup->end-5000), dup_right_end = dup->end;
//    hts_pos_t right_flanking_start = dup->end, right_flanking_end = dup->end+5000;
//
//    // cleanup reads
//    std::vector<read_w_cached_info_t*> left_flanking_reads =
//            cleanup_reads(lf_reads, left_flanking_start, left_flanking_end, dup->original_start, dup->original_end,
//                          stats, config, max_allowed_frac_normalized);
//    std::vector<read_w_cached_info_t*> dup_left_reads =
//            cleanup_reads(ldup_reads, dup_left_start, dup_left_end, dup->original_start, dup->original_end, stats, config,
//                          max_allowed_frac_normalized);
//    std::vector<read_w_cached_info_t*> dup_right_reads =
//            cleanup_reads(rdup_reads, dup_right_start, dup_right_end, dup->original_start, dup->original_end, stats, config,
//                          max_allowed_frac_normalized);
//    std::vector<read_w_cached_info_t*> right_flanking_reads =
//            cleanup_reads(rf_reads, right_flanking_start, right_flanking_end, dup->original_start, dup->original_end,
//                          stats, config,max_allowed_frac_normalized);
//
//    std::set<read_w_cached_info_t*> kept_reads;
//    kept_reads.insert(left_flanking_reads.begin(), left_flanking_reads.end());
//    kept_reads.insert(dup_left_reads.begin(), dup_left_reads.end());
//    kept_reads.insert(dup_right_reads.begin(), dup_right_reads.end());
//    kept_reads.insert(right_flanking_reads.begin(), right_flanking_reads.end());
//
//    uint32_t left_flanking_coverage[5000], dup_left_coverage[5000], dup_right_coverage[5000], right_flanking_coverage[5000];
//    memset(left_flanking_coverage, 0, 5000 * sizeof(left_flanking_coverage[0]));
//    memset(dup_left_coverage, 0, 5000 * sizeof(dup_left_coverage[0]));
//    memset(dup_right_coverage, 0, 5000 * sizeof(dup_right_coverage[0]));
//    memset(right_flanking_coverage, 0, 5000 * sizeof(right_flanking_coverage[0]));
//    for (read_w_cached_info_t* read_w_as : kept_reads) {
//        hts_pos_t rs = read_w_as->start, re = read_w_as->end;
//
//        // note: this may overestimate the coverage, for example if the read has deletions. But it is much simpler and faster than considering cigar
//        hts_pos_t b = std::max(hts_pos_t(0), rs-left_flanking_start), e = std::min(hts_pos_t(5000), re-left_flanking_start);
//        for (hts_pos_t i = b; i < e; i++) {
//            left_flanking_coverage[i]++;
//        }
//
//        b = std::max(hts_pos_t(0), rs-dup_left_start), e = std::min(hts_pos_t(5000), re-dup_left_start);
//        for (hts_pos_t i = b; i < e; i++) {
//            dup_left_coverage[i]++;
//        }
//
//        b = std::max(hts_pos_t(0), rs-dup_right_start), e = std::min(hts_pos_t(5000), re-dup_right_start);
//        for (hts_pos_t i = b; i < e; i++) {
//            dup_right_coverage[i]++;
//        }
//
//        b = std::max(hts_pos_t(0), rs-right_flanking_start), e = std::min(hts_pos_t(5000), re-right_flanking_start);
//        for (hts_pos_t i = b; i < e; i++) {
//            right_flanking_coverage[i]++;
//        }
//    }
//
//    // trying median
//    hts_pos_t dup_tested_len = std::min(dup->end-dup->start, hts_pos_t(5000));
//    std::sort(left_flanking_coverage, left_flanking_coverage+5000, std::greater<uint32_t>());
//    int end = std::find(left_flanking_coverage, left_flanking_coverage+5000, 0)-left_flanking_coverage;
//    dup->left_flanking_cov = left_flanking_coverage[end/2];
//    std::sort(dup_left_coverage, dup_left_coverage+dup_tested_len, std::greater<uint32_t>());
//    dup->indel_left_cov = dup_left_coverage[dup_tested_len/2];
//    std::sort(dup_right_coverage, dup_right_coverage+dup_tested_len, std::greater<uint32_t>());
//    dup->indel_right_cov = dup_right_coverage[dup_tested_len/2];
//    std::sort(right_flanking_coverage, right_flanking_coverage+5000, std::greater<uint32_t>());
//    end = std::find(right_flanking_coverage, right_flanking_coverage+5000, 0)-right_flanking_coverage;
//    dup->right_flanking_cov = right_flanking_coverage[end/2];
//
//    int ow_pairs = 0;
//    for (read_w_cached_info_t* r : dup_left_reads) {
//        hts_pos_t mate_end = r->start + r->isize;
//        if (r->start >= dup->original_start && r->end-dup->original_start <= config.max_is && r->is_rev && !r->is_mrev
//        && mate_end <= dup->original_end && dup->original_end-mate_end <= config.max_is && r->isize > 0) {
//            ow_pairs++;
//        }
//    }
//    dup->ow_pairs = ow_pairs;
//
////    if (dup->left_flanking_cov*1.26<dup->indel_left_cov || dup->indel_right_cov>dup->right_flanking_cov*1.26) {
////        std::string fname = workdir + "/workspace/dup-bams/" + dup->id + ".bam";
////        samFile* writer = get_writer(fname, bam_file->header);
////        for (read_w_cached_info_t* r : kept_reads) {
////            sam_write1(writer, bam_file->header, r->read);
////        }
////        sam_close(writer);
////    }
//}
//
//void clear_rcis(std::vector<read_w_cached_info_t*>& testable_dups_lf_reads, std::vector<read_w_cached_info_t*>& testable_dups_ldup_reads,
//                std::vector<read_w_cached_info_t*>& testable_dups_rdup_reads, std::vector<read_w_cached_info_t*>& testable_dups_rf_reads) {
//    for (read_w_cached_info_t* r : testable_dups_lf_reads) {
//        r->references--;
//        if (r->references == 0) delete r;
//    }
//    for (read_w_cached_info_t* r : testable_dups_ldup_reads) {
//        r->references--;
//        if (r->references == 0) delete r;
//    }
//    for (read_w_cached_info_t* r : testable_dups_rdup_reads) {
//        r->references--;
//        if (r->references == 0) delete r;
//    }
//    for (read_w_cached_info_t* r : testable_dups_rf_reads) {
//        r->references--;
//        if (r->references == 0) delete r;
//    }
//    std::vector<read_w_cached_info_t*>().swap(testable_dups_lf_reads);
//    std::vector<read_w_cached_info_t*>().swap(testable_dups_ldup_reads);
//    std::vector<read_w_cached_info_t*>().swap(testable_dups_rdup_reads);
//    std::vector<read_w_cached_info_t*>().swap(testable_dups_rf_reads);
//}
//
//void depth_filter_dup_w_cleanup(std::string contig_name, std::vector<duplication_t*>& duplications,
//                                open_samFile_t* bam_file, stats_t stats, config_t config,
//                                std::vector<double>& max_allowed_frac_normalized, std::string workdir) {
//    // we have our population-wise distribution of AS differences
//    // any read with AS >= the median for its best bin is considered AS diff = 0
//    // from there, we compute how many reads we allow for each AS diff > 0, based on the population-wise distribution
//    std::vector<char*> regions;
//    for (duplication_t* dup : duplications) {
//        std::stringstream ss;
//        ss << contig_name << ":" << dup->start-5000 << "-" << std::min(dup->start+5000, dup->end);
//        char* region = new char[1000];
//        strcpy(region, ss.str().c_str());
//        regions.push_back(region);
//
//        ss.str(std::string());
//        ss << contig_name << ":" << std::max(dup->end-5000, dup->start) << "-" << dup->end+5000;
//        region = new char[1000];
//        strcpy(region, ss.str().c_str());
//        regions.push_back(region);
//    }
//    std::sort(duplications.begin(), duplications.end(), [](const duplication_t* dup1, const duplication_t* dup2) {
//        return dup1->start < dup2->start;
//    });
//
//    std::vector<std::pair<duplication_t*, int> > duplications_endsorted(duplications.size());
//    for (int i = 0; i < duplications.size(); i++) {
//        duplications_endsorted[i] = std::make_pair(duplications[i], i);
//    }
//    std::sort(duplications_endsorted.begin(), duplications_endsorted.end(),
//              [] (std::pair<duplication_t*, int>& p1, std::pair<duplication_t*, int>& p2)
//              { return p1.first->end < p2.first->end; });
//
//    std::vector<std::vector<read_w_cached_info_t*> > dups_lf_reads(duplications.size()), dups_ldup_reads(duplications.size());
//    std::vector<std::vector<read_w_cached_info_t*> > dups_rdup_reads(duplications.size()), dups_rf_reads(duplications.size());
//    int curr_pos = 0;
//    int processed_dups = 0;
//
//    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
//    bam1_t* read = bam_init1();
//    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
//        if (is_unmapped(read) || is_mate_unmapped(read)) continue;
//        if (read->core.tid != read->core.mtid || read->core.qual < 20) continue; // do not count low qual reads and inter-chr pairs
//
//        hts_pos_t rs = read->core.pos, re = bam_endpos(read);
//        while (curr_pos < duplications.size() && rs > duplications[curr_pos]->end+5000) curr_pos++;
//
//        while (processed_dups < duplications.size() && rs > duplications_endsorted[processed_dups].first->end+5000) {
//            int idx = duplications_endsorted[processed_dups].second;
//            compute_cleaned_up_depth(duplications[idx], bam_file, dups_lf_reads[idx],
//                                     dups_ldup_reads[idx], dups_rdup_reads[idx],
//                                     dups_rf_reads[idx], stats, config, max_allowed_frac_normalized, workdir);
//            clear_rcis(dups_lf_reads[idx],dups_ldup_reads[idx],dups_rdup_reads[idx],dups_rf_reads[idx]);
//            processed_dups++;
//        }
//
//        read_w_cached_info_t* rci = new read_w_cached_info_t(read);
//        for (int i = curr_pos; i < duplications.size() && re > duplications[i]->start-5000; i++) {
//            duplication_t* dup = duplications[i];
//            hts_pos_t left_flanking_start = dup->start - 5000, left_flanking_end = dup->start;
//            hts_pos_t dup_left_start = dup->start, dup_left_end = std::min(dup->end, dup->start+5000);
//            hts_pos_t dup_right_start = std::max(dup->start, dup->end-5000), dup_right_end = dup->end;
//            hts_pos_t right_flanking_start = dup->end, right_flanking_end = dup->end+5000;
//
//            bool overlaps_lf = overlap(rs, re, left_flanking_start, left_flanking_end) > 0;
//            bool overlaps_ldup = overlap(rs, re, dup_left_start, dup_left_end) > 0;
//            bool overlaps_rdup = overlap(rs, re, dup_right_start, dup_right_end) > 0;
//            bool overlaps_rf = overlap(rs, re, right_flanking_start, right_flanking_end) > 0;
//
//            if (overlaps_lf || overlaps_ldup || overlaps_rdup || overlaps_rf) {
//                if (overlaps_lf) {
//                    dups_lf_reads[i].push_back(rci);
//                }
//                if (overlaps_ldup) {
//                    dups_ldup_reads[i].push_back(rci);
//                }
//                if (overlaps_rdup) {
//                    dups_rdup_reads[i].push_back(rci);
//                }
//                if (overlaps_rf) {
//                    dups_rf_reads[i].push_back(rci);
//                }
//                rci->references += overlaps_lf + overlaps_ldup + overlaps_rdup + overlaps_rf;
//            }
//        }
//    }
//    while (processed_dups < duplications.size()) {
//        int idx = duplications_endsorted[processed_dups].second;
//        compute_cleaned_up_depth(duplications[idx], bam_file, dups_lf_reads[idx],
//                                 dups_ldup_reads[idx], dups_rdup_reads[idx],
//                dups_rf_reads[idx], stats, config, max_allowed_frac_normalized, workdir);
//        clear_rcis(dups_lf_reads[idx],dups_ldup_reads[idx],dups_rdup_reads[idx],dups_rf_reads[idx]);
//        processed_dups++;
//    }
//
//    for (char* region : regions) {
//        delete[] region;
//    }
//    hts_itr_destroy(iter);
//    bam_destroy1(read);
//}

void depth_filter_indel(std::string contig_name, hts_pos_t contig_len, std::vector<indel_t*>& indels, open_samFile_t* bam_file, int flanking_size = 5000) {
    std::vector<char*> regions;
    for (indel_t* indel : indels) {
        std::stringstream ss;
        ss << contig_name << ":" << indel->start - flanking_size << "-" << std::min(indel->start + flanking_size, indel->end);
        char* region = new char[1000];
        strcpy(region, ss.str().c_str());
        regions.push_back(region);

        ss.str(std::string());
        ss << contig_name << ":" << std::max(indel->end - flanking_size, indel->start) << "-" << indel->end + flanking_size;
        region = new char[1000];
        strcpy(region, ss.str().c_str());
        regions.push_back(region);
    }
    std::sort(indels.begin(), indels.end(), [](const indel_t* indel1, const indel_t* indel2) {
        return indel1->start < indel2->start;
    });

    std::vector<int64_t> flanking_left_cov(indels.size()), del_left_cov(indels.size());
    std::vector<int64_t> del_right_cov(indels.size()), flanking_right_cov(indels.size());

    std::vector<std::vector<int64_t> > _flanking_left_cov(indels.size(), std::vector<int64_t>(flanking_size)), _indel_left_cov(indels.size(), std::vector<int64_t>(flanking_size));
	std::vector<std::vector<int64_t> > _indel_right_cov(indels.size(), std::vector<int64_t>(flanking_size)), _flanking_right_cov(indels.size(), std::vector<int64_t>(flanking_size));

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
    bam1_t* read = bam_init1();
    int curr_pos = 0;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read)) continue;
        if (read->core.tid != read->core.mtid || read->core.qual < 20) continue; // do not count low qual reads and inter-chr pairs

        hts_pos_t rs = read->core.pos, re = bam_endpos(read);
        while (curr_pos < indels.size() && indels[curr_pos]->end + flanking_size < rs) curr_pos++;

        if (curr_pos == indels.size()) break;

        for (int i = curr_pos; i < indels.size() && indels[i]->start - flanking_size < re; i++) {
            indel_t* indel = indels[i];
            hts_pos_t indel_tested_len = std::min(indel->end-indel->start, hts_pos_t(flanking_size));

            hts_pos_t left_flanking_start = std::max(hts_pos_t(0), indel->start-flanking_size);
            hts_pos_t left_flanking_len = indel->start - left_flanking_start;
			hts_pos_t indel_lcov_start = indel->start;
			hts_pos_t indel_rcov_start = std::max(indel->start, indel->end-flanking_size);
			hts_pos_t right_flanking_start = indel->end;

            hts_pos_t b = std::max(hts_pos_t(0), rs-left_flanking_start), e = std::min(left_flanking_len, re-left_flanking_start);
            for (hts_pos_t j = b; j < e; j++) {
            	_flanking_left_cov[i][j]++;
			}

            b = std::max(hts_pos_t(0), rs-indel_lcov_start), e = std::min(indel_tested_len, re-indel_lcov_start);
			for (hts_pos_t j = b; j < e; j++) {
				_indel_left_cov[i][j]++;
			}

			b = std::max(hts_pos_t(0), rs-indel_rcov_start), e = std::min(indel_tested_len, re-indel_rcov_start);
			for (hts_pos_t j = b; j < e; j++) {
				_indel_right_cov[i][j]++;
			}

			b = std::max(hts_pos_t(0), rs-right_flanking_start), e = std::min(hts_pos_t(flanking_size), re-right_flanking_start);
			for (hts_pos_t j = b; j < e; j++) {
				_flanking_right_cov[i][j]++;
			}

            flanking_left_cov[i] += overlap(indel->start - flanking_size, indel->start, rs, re);
            del_left_cov[i] += overlap(indel->start, std::min(indel->start + flanking_size, indel->end), rs, re);
            del_right_cov[i] += overlap(std::max(indel->end - flanking_size, indel->start), indel->end, rs, re);
            flanking_right_cov[i] += overlap(indel->end, indel->end + flanking_size, rs, re);
        }
    }

    for (int i = 0; i < indels.size(); i++) {
        indel_t* indel = indels[i];
        hts_pos_t left_flanking_region = std::max(std::min(indel->start, hts_pos_t(flanking_size)), hts_pos_t(1));
        indel->left_flanking_cov = flanking_left_cov[i] / left_flanking_region;
        hts_pos_t indel_len = std::max(std::min(indel->len(), hts_pos_t(flanking_size)), hts_pos_t(1));
        indel->indel_left_cov = del_left_cov[i] / indel_len;
        indel->indel_right_cov = del_right_cov[i] / indel_len;
        hts_pos_t right_flanking_region = std::max(std::min(contig_len-indel->end, hts_pos_t(flanking_size)), hts_pos_t(1));
        indel->right_flanking_cov = flanking_right_cov[i] / right_flanking_region;

        // trying median
		hts_pos_t indel_tested_len = std::min(indel->end-indel->start, hts_pos_t(flanking_size));
		std::sort(_flanking_left_cov[i].begin(), _flanking_left_cov[i].end(), std::greater<int64_t>());
		int end = std::find(_flanking_left_cov[i].begin(), _flanking_left_cov[i].end(), 0)-_flanking_left_cov[i].begin();
		indel->mleft_flanking_cov = _flanking_left_cov[i][end/2];
		std::sort(_indel_left_cov[i].begin(), _indel_left_cov[i].end(), std::greater<int64_t>());
		indel->mindel_left_cov = _indel_left_cov[i][indel_tested_len/2];
		std::sort(_indel_right_cov[i].begin(), _indel_right_cov[i].end(), std::greater<int64_t>());
		indel->mindel_right_cov = _indel_right_cov[i][indel_tested_len/2];
		std::sort(_flanking_right_cov[i].begin(), _flanking_right_cov[i].end(), std::greater<int64_t>());
		end = std::find(_flanking_right_cov[i].begin(), _flanking_right_cov[i].end(), 0)-_flanking_right_cov[i].begin();
		indel->mright_flanking_cov = _flanking_right_cov[i][end/2];
    }

    for (char* region : regions) {
        delete[] region;
    }
    hts_itr_destroy(iter);
    bam_destroy1(read);
}

void depth_filter_del(std::string contig_name, hts_pos_t contig_len, std::vector<deletion_t*>& deletions, open_samFile_t* bam_file, int flanking_size = 5000) {
    std::vector<indel_t*> testable_indels(deletions.begin(), deletions.end());
    depth_filter_indel(contig_name, contig_len, testable_indels, bam_file, flanking_size);
}
void depth_filter_dup(std::string contig_name, hts_pos_t contig_len, std::vector<duplication_t*>& duplications, open_samFile_t* bam_file, int flanking_size = 5000) {
    std::vector<indel_t*> testable_indels(duplications.begin(), duplications.end());
    depth_filter_indel(contig_name, contig_len, testable_indels, bam_file, flanking_size);
}
void depth_filter_ins(std::string contig_name, hts_pos_t contig_len, std::vector<insertion_t*>& insertions, open_samFile_t* bam_file, int flanking_size = 5000) {
    std::vector<indel_t*> testable_indels(insertions.begin(), insertions.end());
    depth_filter_indel(contig_name, contig_len, testable_indels, bam_file, flanking_size);
}

void calculate_ptn_ratio(std::string contig_name, std::vector<deletion_t*>& deletions, open_samFile_t* bam_file, config_t config, stats_t stats) {

	std::sort(deletions.begin(), deletions.end(), [](const deletion_t* d1, const deletion_t* d2) {
		return (d1->start+d1->end)/2 < (d2->start+d2->end)/2;
	});

	std::vector<hts_pos_t> midpoints, sizes;
	std::vector<char*> regions;
	for (deletion_t* deletion : deletions) {
		hts_pos_t midpoint = (deletion->start+deletion->end)/2, size = deletion->end-deletion->start;
		midpoints.push_back(midpoint);
		sizes.push_back(size);

		std::stringstream ss;
		ss << contig_name << ":" << deletion->start-config.max_is << "-" << deletion->start;
		char* region = new char[ss.str().length()+1];
		strcpy(region, ss.str().c_str());
		regions.push_back(region);

		ss.str("");
		ss << contig_name << ":" << midpoint-config.max_is << "-" << midpoint;
		region = new char[ss.str().length()+1];
		strcpy(region, ss.str().c_str());
		regions.push_back(region);
	}
	if (regions.empty()) return;

	std::vector<int> tmp_disc_pairs_count(deletions.size());
	std::vector<std::vector<int> > disc_pairs_to_del;

	int curr_pos = 0;
	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
	bam1_t* read = bam_init1();
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		while (curr_pos < deletions.size() && midpoints[curr_pos] < read->core.pos) curr_pos++;

		if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read)) continue;
		if (!is_samechr(read) || is_samestr(read) || bam_is_rev(read) || read->core.isize <= 0) continue;
//		if (read->core.qual < 20) continue;

		hts_pos_t start = read->core.pos + read->core.l_qseq/2;
		hts_pos_t end = read->core.pos + read->core.isize - read->core.l_qseq/2;

		std::vector<int> _disc_pairs_to_del;
		for (int i = curr_pos; i < deletions.size() && midpoints[i] < read->core.pos+read->core.isize; i++) {
			// normally a discordant pair should not be bigger than max_is + del size
			// however, we are a bit permissive (max_is tolerance) because some imprecise deletions may be underestimated
//			if (read->core.isize > 2*config.max_is+sizes[i]) continue;
			if (read->core.isize > 1*config.max_is+sizes[i]) continue;
			if (start < deletions[i]->start && deletions[i]->end < end && read->core.isize > stats.max_is) {
				tmp_disc_pairs_count[i]++;
				_disc_pairs_to_del.push_back(i);
			}
			else if (start <= midpoints[i] && midpoints[i] <= end && read->core.isize <= stats.max_is) deletions[i]->conc_pairs++;
		}

		if (_disc_pairs_to_del.size() == 1) {
			deletions[_disc_pairs_to_del[0]]->disc_pairs++;
		} else if (_disc_pairs_to_del.size() > 1) {
			disc_pairs_to_del.push_back(_disc_pairs_to_del);
		}
	}

	// for DP that support multiple deletions, choose the one(s) with maximal support
	// If a deletion is accepted by SR, force it to have maximal support
	for (int i = 0; i < deletions.size(); i++) {
		if (gt_is_positive(deletions[i]->called_gt) && deletions[i]->filter.empty()) {
			tmp_disc_pairs_count[i] = INT32_MAX;
		}
	}

	for (auto& v : disc_pairs_to_del) {
		int max_dp = 0;
		for (int i : v) {
			if (tmp_disc_pairs_count[i] > max_dp) max_dp = tmp_disc_pairs_count[i];
		}
		for (int i : v) {
			if (tmp_disc_pairs_count[i] == max_dp) deletions[i]->disc_pairs++;
		}
	}

	for (char* region : regions) {
		delete[] region;
	}
}

void calculate_confidence_interval_size(std::string contig_name, std::vector<indel_t*>& indels, open_samFile_t* bam_file,
                                        config_t config, stats_t stats, std::vector<double>& population, bool do_ks_test = false) {

	std::sort(indels.begin(), indels.end(), [](const indel_t* d1, const indel_t* d2) {
        return (d1->start+d1->end)/2 < (d2->start+d2->end)/2;
    });

    std::vector<hts_pos_t> midpoints, sizes;
    std::vector<uint64_t> sums(indels.size()), sq_sums(indels.size());
    std::vector<uint32_t> ns(indels.size());
    std::vector<std::vector<double> > is_dists(indels.size());
    std::vector<char*> regions;
    for (indel_t* indel : indels) {
        hts_pos_t midpoint = (indel->start+indel->end)/2, size = indel->end-indel->start;
        midpoints.push_back(midpoint);
        sizes.push_back(size);

        std::stringstream ss;
        ss << contig_name << ":" << indel->start-(config.max_is+size) << "-" << indel->start;
        char* region = new char[1000];
        strcpy(region, ss.str().c_str());
        regions.push_back(region);
    }
    if (regions.empty()) return;

    int curr_pos = 0;
    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
    bam1_t* read = bam_init1();
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        while (curr_pos < indels.size() && midpoints[curr_pos] < read->core.pos) curr_pos++;

        if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read)) continue;
        if (!is_samechr(read) || is_samestr(read) || bam_is_rev(read) || read->core.isize <= 0) continue;

        hts_pos_t start = read->core.pos + read->core.l_qseq/2;
        hts_pos_t end = read->core.pos + read->core.isize - read->core.l_qseq/2;

//        for (int i = curr_pos; i < indels.size() && midpoints[i] < end; i++) {
//            if (start <= midpoints[i] && midpoints[i] <= end) {
//                is_dists[i].push_back(read->core.isize);
//            }
//        }

        if (read->core.qual < 20) continue;

        for (int i = curr_pos; i < indels.size() && midpoints[i] < read->core.pos+read->core.isize; i++) {
            if (start <= midpoints[i] && midpoints[i] <= end && read->core.isize <= config.max_is+sizes[i]) {
                sums[i] += read->core.isize;
                sq_sums[i] += read->core.isize*read->core.isize;
                ns[i]++;
            }
        }
    }

    for (int i = 0; i < indels.size(); i++) {
        indel_t* indel = indels[i];
        uint32_t n = ns[i];
        uint64_t sum = sums[i], sq_sum = sq_sums[i];
        if (n >= 4) {
            int avg_is = sum/n;
            int var_is = (sq_sum - sum*sum/n)/(n-1);
            int confidence_ival = 2.576 * sqrt(var_is/n);
            int predicted_size = avg_is - stats.pop_avg_crossing_is;
            indel->max_conf_size = predicted_size + confidence_ival;
            indel->min_conf_size = predicted_size - confidence_ival;
        }
    }

    for (char* region : regions) {
        delete[] region;
    }
    hts_itr_destroy(iter);
    bam_destroy1(read);
}

void calculate_confidence_interval_size(std::string contig_name, std::vector<deletion_t*>& indels, open_samFile_t* bam_file,
                                        config_t config, stats_t stats, std::vector<double>& population, bool do_ks_test = false) {
	std::vector<indel_t*> deletions(indels.begin(), indels.end());
	calculate_confidence_interval_size(contig_name, deletions, bam_file, config, stats, population, do_ks_test);
}
void calculate_confidence_interval_size(std::string contig_name, std::vector<insertion_t*>& indels, open_samFile_t* bam_file,
                                        config_t config, stats_t stats, std::vector<double>& population, bool do_ks_test = false) {
	std::vector<indel_t*> insertions(indels.begin(), indels.end());
	calculate_confidence_interval_size(contig_name, insertions, bam_file, config, stats, population, do_ks_test);
}

#endif //SURVINDEL2_STAT_TESTS_H
