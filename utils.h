#ifndef SURVTYPER_UTILS_H
#define SURVTYPER_UTILS_H

#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unistd.h>
#include <cmath>

#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "libs/ssw_cpp.h"
#include "libs/ssw.h"
KSEQ_INIT(int, read)

struct config_t {
	std::string bam_fname, reference_fname, workdir, sample_name;
    int threads, read_len;
    int max_is;
    double max_seq_error = 0.04;
    int clip_penalty = 7;
    int min_size_for_depth_filtering = 1000;
    int max_pop_size = 10000;
    int min_score_diff = 15;
    int min_sv_size = 50;
    int min_clip_size = 15;
    int min_mh_len = 100;
    int max_size_diff = 200;
    bool save_evidence = false;

    void parse(std::string config_file) {
        std::unordered_map<std::string, std::string> config_params;
        std::ifstream fin(config_file);
        std::string name, value;
        while (fin >> name >> value) {
            config_params[name] = value;
        }
        fin.close();

        bam_fname = config_params["bam_file"];
        workdir = config_params["workdir"];
        reference_fname = config_params["reference"];
        sample_name = config_params["sample_name"];
        threads = stoi(config_params["threads"]);
        read_len = stoi(config_params["read_len"]);
        max_is = stoi(config_params["max_is"]);
        max_size_diff = stoi(config_params["max_size_diff"]);
        save_evidence = stoi(config_params["save_evidence"]);
    };
};

struct stats_t {
    int pop_avg_crossing_is;
    int min_is, max_is;
    int min_depth, max_depth, median_depth;
    int min_pairs_crossing, max_pairs_crossing;

    stats_t(int pop_avg_crossing_is, int min_is, int max_is, int min_depth, int median_depth, int max_depth, int min_pairs_crossing, int max_pairs_crossing) :
    	pop_avg_crossing_is(pop_avg_crossing_is), min_is(min_is), max_is(max_is), min_depth(min_depth), median_depth(median_depth),
		max_depth(max_depth), min_pairs_crossing(min_pairs_crossing), max_pairs_crossing(max_pairs_crossing) {}
};


struct chr_seq_t {
    char* seq;
    hts_pos_t len;

    chr_seq_t(char* seq, hts_pos_t len) : seq(seq), len(len) {}
    ~chr_seq_t() {delete[] seq;}
};
struct chr_seqs_map_t {
    std::unordered_map<std::string, chr_seq_t*> seqs;

    void read_fasta_into_map(std::string& reference_fname) {
        FILE* fasta = fopen(reference_fname.c_str(), "r");
        kseq_t* seq = kseq_init(fileno(fasta));
        while (kseq_read(seq) >= 0) {
            std::string seq_name = seq->name.s;
            char* chr_seq = new char[seq->seq.l + 1];
            strcpy(chr_seq, seq->seq.s);
            seqs[seq_name] = new chr_seq_t(chr_seq, seq->seq.l);
        }
        kseq_destroy(seq);
        fclose(fasta);
    }

    char* get_seq(std::string seq_name) {
        return seqs[seq_name]->seq;
    }

    hts_pos_t get_len(std::string seq_name) {
        return seqs[seq_name]->len;
    }

    void clear() {
        for (auto& e : seqs) {
            delete e.second;
            e.second = NULL;
        }
    }

    ~chr_seqs_map_t() {
        clear();
    }
};


std::string get_sv_type(bcf1_t* sv, bcf_hdr_t* hdr) {
    char* data = NULL;
    int len = 0;
    if (bcf_get_info_string(hdr, sv, "SVTYPE", &data, &len) < 0) {
    	std::cout << "SVTYPE not found" << std::endl;
        throw std::runtime_error("Failed to determine SVTYPE for sv " + std::string(sv->d.id));
    }
    std::string svtype = data;
    delete[] data;
    return svtype;
}

int get_sv_end(bcf1_t* sv, bcf_hdr_t* hdr) {
    int* data = NULL;
    int size = 0;
    bcf_get_info_int32(hdr, sv, "END", &data, &size);
    if (size > 0) {
        int end = data[0];
        delete[] data;
        return end-1; // return 0-based
    }

    bcf_get_info_int32(hdr, sv, "SVLEN", &data, &size);
    if (size > 0) {
        int svlen = data[0];
        delete[] data;
        return sv->pos + abs(svlen);
    }

    throw std::runtime_error("SV " + std::string(sv->d.id) + "has no END or SVLEN annotation.");
}

std::string get_ins_seq(bcf1_t* sv, bcf_hdr_t* hdr) {
	// priority to the ALT allele, if it is not symbolic and longer than just the padding base
	char c = toupper(sv->d.allele[1][0]);
	if ((c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N') && strlen(sv->d.allele[1]) > 1) {
		return sv->d.allele[1];
	}

	// otherwise, look for SVINSSEQ (compliant with Manta)
	char* data = NULL;
	int size = 0;
	bcf_get_info_string(hdr, sv, "SVINSSEQ", (void**) &data, &size);
	if (data) return data;

	return "";
}

int get_left_clip_size(StripedSmithWaterman::Alignment& aln) {
    uint32_t l = aln.cigar[0];
    return cigar_int_to_op(l) == 'S' ? cigar_int_to_len(l) : 0;
}
int get_right_clip_size(StripedSmithWaterman::Alignment& aln) {
    uint32_t r = aln.cigar[aln.cigar.size()-1];
    return cigar_int_to_op(r) == 'S' ? cigar_int_to_len(r) : 0;
}

bool gt_is_positive(std::string& gt) {
	return gt == "HET" || gt == "HOM_ALT" || gt == "SV_EXISTS";
}

struct indel_t {
    std::string id;
    hts_pos_t start, end, remapped_start, remapped_end;
    int alt_better = 0, ref_better = 0, same = 0,  alt_better_strong = 0, alt_better_goodqual = 0;
    int alt_better_fwd = 0, alt_better_rev = 0;
    int vouching_for_lh_consensus_remapping = 0, vouching_for_rh_consensus_remapping = 0;

    int alt_better_max_mq = 0;
    std::string rs_cigar, re_cigar, full_cigar;
    std::string filter;
    bool anomalous_cov = false;

    double mm_rate = 0.0;
    int left_flanking_cov = 0, indel_left_cov = 0, indel_right_cov = 0, right_flanking_cov = 0;
    int mleft_flanking_cov = 0, mindel_left_cov = 0, mindel_right_cov = 0, mright_flanking_cov = 0; // median version
    int min_conf_size = 0, max_conf_size = 0;

    bool lowq_alt_allele = false;
    bool was_remapped = false;

    std::string called_gt = "NO_GT", before_read_qc_gt = "NO_GT";
    double ratio = 0.0;
    bool gt_skipped = false;

    bcf1_t* bcf_entry = NULL;

    indel_t(std::string id, hts_pos_t start, hts_pos_t end)
            : id(id), start(start), end(end), remapped_start(start), remapped_end(end) {}

	hts_pos_t len() { return end-start; }
	virtual std::string indel_type() { return ""; };
};

struct deletion_t : indel_t {
	std::string ins_seq;

	bool split_aln_accepted = true;
    double p_val = -1.0, p_val2 = 1.0;
    int n_is_kstest = 0;

    int disc_pairs = 0, conc_pairs = 0;
    std::vector<std::string> disc_pairs_qnames;

    deletion_t(std::string id, hts_pos_t start, hts_pos_t end, std::string ins_seq, std::string benchmark_gt) :
               indel_t(id, start, end), ins_seq(ins_seq) {};

    std::string indel_type() override { return "DEL"; }

    bool is_short(config_t& config) {
		return len() <= config.max_is;
	}
};

struct duplication_t : indel_t {
	int copy_number = 0, before_qc_copy_number = 0;
	int compl_inside_dup = 0, ow_pairs = 0;
	int predicted_insertion_len = 0;

	duplication_t(std::string id, hts_pos_t start, hts_pos_t end) : indel_t(id, start, end) {}

	bool is_short(config_t& config) {
		return len() <= config.read_len-2*config.min_clip_size;
	}
};

struct insertion_t : indel_t {
	int alt_bp1_better_strong = 0, alt_bp1_better_goodqual = 0, alt_bp2_better_strong = 0, alt_bp2_better_goodqual = 0;
	int disc_pairs_fwd_stable = 0, disc_pairs_rev_stable = 0, conc_pairs = 0;
	std::string ins_seq;
	hts_pos_t fwd_anchor_start = INT32_MAX, fwd_anchor_end = 0, rev_anchor_start = INT32_MAX, rev_anchor_end = 0;

	std::string remapped_cigar;
	int remapped_ins_len = INT32_MAX;

	insertion_t(std::string id, hts_pos_t start, hts_pos_t end, std::string ins_seq) : indel_t(id, start, end), ins_seq(ins_seq) {}
};

int64_t overlap(hts_pos_t s1, hts_pos_t e1, hts_pos_t s2, hts_pos_t e2) {
    int64_t overlap = std::min(e1, e2) - std::max(s1, s2);
    return std::max(int64_t(0), overlap);
}

void strgt2datagt(std::string gt, int data[2]) {
	if (gt == "NO_GT") {
		data[0] = data[1] = bcf_gt_missing;
	} else if (gt == "HOM_REF") {
		data[0] = data[1] = bcf_gt_unphased(0);
	} else if (gt == "HET") {
		data[0] = bcf_gt_unphased(0);
		data[1] = bcf_gt_unphased(1);
	} else if (gt == "HOM_ALT") {
		data[0] = data[1] = bcf_gt_unphased(1);
	} else if (gt == "SV_EXISTS") {
		data[0] = bcf_gt_unphased(1);
		data[1] = bcf_gt_missing;
	}
}

struct suffix_prefix_aln_t {
    int overlap, score, mismatches;

    suffix_prefix_aln_t(int overlap, int score, int mismatches) : overlap(overlap), score(score), mismatches(mismatches) {}
};

// Finds the best alignment between a suffix of s1 and a prefix of s2
// Disallows gaps
suffix_prefix_aln_t aln_suffix_prefix(std::string& s1, std::string& s2, int match_score, int mismatch_score, double max_seq_error,
                                      int min_overlap = 1, int max_overlap = INT32_MAX, int max_mismatches = INT32_MAX) {
    int best_score = 0, best_aln_mismatches = 0;
    int overlap = 0;

    for (int i = std::max(0, (int) s1.length()-max_overlap); i < s1.length()-min_overlap+1; i++) {
        if (i+s2.length() < s1.length()) continue;

        int sp_len = s1.length()-i;
        if (best_score >= sp_len*match_score) break; // current best score is unbeatable

        int mismatches = 0;
        const char* s1_suffix = s1.data()+i;
        const char* s2_prefix = s2.data();
        while (*s1_suffix) {
            if (*s1_suffix != *s2_prefix) mismatches++;
            s1_suffix++; s2_prefix++;
        }

        int score = (sp_len-mismatches)*match_score + mismatches*mismatch_score;

        int max_acceptable_mm = max_seq_error == 0.0 ? 0 : std::max(1.0, sp_len*max_seq_error);
        if (best_score < score && mismatches <= max_acceptable_mm && mismatches <= max_mismatches) {
            best_score = score;
            best_aln_mismatches = mismatches;
            overlap = sp_len;
        }
    }
    return suffix_prefix_aln_t(overlap, best_score, best_aln_mismatches);
}

bool is_contained(std::string& ref, std::string& query, int max_mismatches) {
	for (int i = 0; i < ref.length()-query.length(); i++) {
		int mm = 0;
		const char* r = ref.data()+i, * q = query.data();
		for (int j = 0; j < query.length(); j++) {
			if (*r != *q) mm++;
			r++; q++;
			if (mm > max_mismatches) break;
		}
		if (mm <= max_mismatches) return true;
	}
	return false;
}

template<typename T>
inline T max(T a, T b, T c) { return std::max(std::max(a,b), c); }

template<typename T>
inline T max(T a, T b, T c, T d) { return std::max(std::max(a,b), std::max(c,d)); }

int score(char a, char b, int match_score, int mismatch_penalty) {
	return (toupper(a) == toupper(b) || a == 'N' || b == 'N') ? match_score : mismatch_penalty;
}
int* smith_waterman_gotoh(const char* ref, int ref_len, const char* read, int read_len, int match_score, int mismatch_penalty, int gap_open, int gap_extend) {
	const int INF = 1000000;

	int** dab = new int*[ref_len+1];
	int** dag = new int*[ref_len+1];
	int** dgb = new int*[ref_len+1];
	for (int i = 0; i <= ref_len; i++) {
		dab[i] = new int[read_len+1];
		dag[i] = new int[read_len+1];
		dgb[i] = new int[read_len+1];
		std::fill(dab[i], dab[i]+read_len+1, 0);
		std::fill(dag[i], dag[i]+read_len+1, 0);
		std::fill(dgb[i], dgb[i]+read_len+1, 0);
	}

	for (int i = 1; i <= ref_len; i++) {
		dab[i][0] = -INF;
		dag[i][0] = -INF;
		dgb[i][0] = gap_open + (i-1)*gap_extend;
	}
	for (int i = 1; i <= read_len; i++) {
		dab[0][i] = -INF;
		dag[0][i] = gap_open + (i-1)*gap_extend;
		dgb[0][i] = -INF;
	}

	for (int i = 1; i <= ref_len; i++) {
		for (int j = 1; j <= read_len; j++) {
			dab[i][j] = score(ref[i-1], read[j-1], match_score, mismatch_penalty) + max(dab[i-1][j-1], dag[i-1][j-1], dgb[i-1][j-1], 0);
			dag[i][j] = max(gap_open + dab[i][j-1], gap_extend + dag[i][j-1], gap_open + dgb[i][j-1]);
			dgb[i][j] = max(gap_open + dab[i-1][j], gap_open + dag[i-1][j], gap_extend + dgb[i-1][j]);
		}
	}

	int* prefix_scores = new int[read_len];
	std::fill(prefix_scores, prefix_scores+read_len, 0);
	for (int i = 1; i <= ref_len; i++) {
		for (int j = 1; j <= read_len; j++) {
			prefix_scores[j-1] = std::max(prefix_scores[j-1], dab[i][j]);
		}
	}

	for (int i = 0; i <= ref_len; i++) {
		delete[] dab[i];
		delete[] dag[i];
		delete[] dgb[i];
	}
	delete[] dab;
	delete[] dag;
	delete[] dgb;

	for (int i = 1; i < read_len; i++) {
		prefix_scores[i] = std::max(prefix_scores[i], prefix_scores[i-1]);
	}

	return prefix_scores;
}

bool accept_aln(StripedSmithWaterman::Alignment& aln, int min_clip_size) {
    return get_left_clip_size(aln) < min_clip_size && get_right_clip_size(aln) < min_clip_size;
}
bool accept_aln_strict(StripedSmithWaterman::Alignment& aln) {
	return get_left_clip_size(aln) == 0 && get_right_clip_size(aln) == 0;
}

int count_indels(StripedSmithWaterman::Alignment& aln) {
	int indels = 0;
	for (uint32_t c : aln.cigar) {
		char op = bam_cigar_opchr(c);
		if (bam_cigar_opchr(c)) {
			if (op == 'S' || op == 'I' || op == 'D') indels++;
		}
	}
	return indels;
}

int count_inserted_bps(StripedSmithWaterman::Alignment& aln) {
	int indels = 0;
	for (uint32_t c : aln.cigar) {
		char op = bam_cigar_opchr(c);
		if (bam_cigar_opchr(c)) {
			if (op == 'I') indels += bam_cigar_oplen(c);
		}
	}
	return indels;
}

int largest_del(StripedSmithWaterman::Alignment& aln) {
	int ldel = 0;
	for (uint32_t c : aln.cigar) {
		char op = bam_cigar_opchr(c);
		if (bam_cigar_opchr(c)) {
			if (op == 'D') ldel = std::max(ldel, (int) bam_cigar_oplen(c));
		}
	}
	return ldel;
}

#endif //SURVTYPER_UTILS_H
