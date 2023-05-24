#ifndef SURVTYPER_SAM_UTILS_H
#define SURVTYPER_SAM_UTILS_H

#include <sstream>

bool is_unmapped(bam1_t* r) {
    return r->core.flag & BAM_FUNMAP;
}
bool is_mate_unmapped(bam1_t* r) {
    return r->core.flag & BAM_FMUNMAP;
}

bool is_primary(bam1_t* r) {
    return !(r->core.flag & BAM_FSECONDARY) && !(r->core.flag & BAM_FSUPPLEMENTARY);
}

bool is_samechr(bam1_t* r) {
    return r->core.tid == r->core.mtid;
}
bool is_samestr(bam1_t* r) {
    return (r->core.flag & BAM_FREVERSE) == (r->core.flag & BAM_FMREVERSE);
}
bool is_dc_pair(bam1_t* r) {
    return !is_samechr(r) || std::abs(r->core.isize) > 100000 || is_unmapped(r) != is_mate_unmapped(r);
}

int get_left_clip_size(bam1_t* r) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[0]) == 'S' ? bam_cigar_oplen(cigar[0]): 0;
}
int get_right_clip_size(bam1_t* r) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[r->core.n_cigar-1]) == 'S' ? bam_cigar_oplen(cigar[r->core.n_cigar-1]): 0;
}

bool is_left_clipped(bam1_t* r, int min_clip_len) {
    return get_left_clip_size(r) >= min_clip_len;
}
bool is_right_clipped(bam1_t* r, int min_clip_len) {
    return get_right_clip_size(r) >= min_clip_len;
}

hts_pos_t get_unclipped_start(bam1_t* r) {
    return r->core.pos - get_left_clip_size(r);
}
hts_pos_t get_unclipped_end(bam1_t* r) {
    return bam_endpos(r) + get_right_clip_size(r);
}

char* get_mc_tag(bam1_t* r) {
	const uint8_t* mc_tag = bam_aux_get(r, "MC");
	if (mc_tag == NULL) {
		throw "Read " + std::string(bam_get_qname(r)) + " does not have the MC tag.";
	}
	return bam_aux2Z(mc_tag);
}
bool is_mate_left_clipped(bam1_t* r) {
    char* mc_tag_str = get_mc_tag(r);
    int i = 0;
    while (mc_tag_str[i] >= '0' && mc_tag_str[i] <= '9') i++;
    return mc_tag_str[i] == 'S';
}
bool is_mate_right_clipped(bam1_t* r) {
	char* mc_tag_str = get_mc_tag(r);
    int i = strlen(mc_tag_str)-1;
    return mc_tag_str[i] == 'S';
}
bool is_mate_clipped(bam1_t* r) {
	return is_mate_left_clipped(r) && is_mate_right_clipped(r);
}

char get_base(const uint8_t* seq, int i) {
    char nucl2chr[16];
    nucl2chr[1] = 'A'; nucl2chr[2] = 'C'; nucl2chr[4] = 'G'; nucl2chr[8] = 'T'; nucl2chr[15] = 'N';
    return nucl2chr[bam_seqi(seq, i)];
}

hts_pos_t get_mate_end(bam1_t* r) {
	uint8_t *mcs = bam_aux_get(r, "MC");
	if (mcs == NULL) return r->core.mpos; // if no MC, return mpos

	char* mc = bam_aux2Z(mcs);
	uint32_t* cigar = NULL;
	size_t size = 0;
	sam_parse_cigar(mc, NULL, &cigar, &size);

	return r->core.mpos + bam_cigar2rlen(size, cigar);
}

std::string get_sequence(bam1_t* r) {
    char seq[100000];
    const uint8_t* bam_seq = bam_get_seq(r);
    for (int i = 0; i < r->core.l_qseq; i++) {
        seq[i] = get_base(bam_seq, i);
    }
    seq[r->core.l_qseq] = '\0';
    return std::string(seq);
}

std::string get_quals(bam1_t* r) {
	char qual[100000];
	const uint8_t* bam_qual = bam_get_qual(r);
	for (int i = 0; i < r->core.l_qseq; i++) {
		qual[i] = bam_qual[i]+33;
	}
	qual[r->core.l_qseq] = '\0';
	return std::string(qual);
}

std::string get_cigar_code(bam1_t* r) {
    const uint32_t* cigar = bam_get_cigar(r);
    std::stringstream ss;
    for (int i = 0; i < r->core.n_cigar; i++) {
        ss << bam_cigar_oplen(cigar[i]) << bam_cigar_opchr(cigar[i]);
    }
    return ss.str();
}

std::pair<int, int> get_mismatches_prefix(bam1_t* read, int prefix_len) {
    int indels = 0, mismatches = 0;
    int original_prefix_len = prefix_len;

    uint32_t* cigar = bam_get_cigar(read);
    for (int j = 0; j < read->core.n_cigar && prefix_len > 0; j++) {
        uint32_t op = cigar[j];
        char opchr = bam_cigar_opchr(op);
        int oplen = bam_cigar_oplen(op);
        if (opchr == 'I' || opchr == 'D') {
            indels += std::min(oplen, prefix_len);
        }
        if (opchr != 'D') {
            prefix_len -= oplen;
        }
    }

    prefix_len = original_prefix_len;
    if (bam_cigar_opchr(cigar[0]) == 'S') prefix_len -= bam_cigar_oplen(cigar[0]) ;

    char* md_tag = bam_aux2Z(bam_aux_get(read, "MD"));
    int md_tag_len = strlen(md_tag);
    int n = 0;
    for (int i = 0; i < md_tag_len && prefix_len > 0; i++) {
        if (md_tag[i] >= '0' && md_tag[i] <= '9') {
            n = n*10 + md_tag[i]-'0';
        } else {
            prefix_len -= n;
            if (prefix_len <= 0) break;
            n = 0;

            if (md_tag[i] == '^') {
                while (md_tag[i] < '0' || md_tag[i] > '9') i++;
            } else {
                mismatches++;
                prefix_len--;
            }
        }
    }

    return {mismatches, indels};
}

std::pair<int, int> get_mismatches_suffix(bam1_t* read, int prefix_len) {
    int indels = 0, mismatches = 0;
    int original_prefix_len = prefix_len;

    uint32_t* cigar = bam_get_cigar(read);
    for (int j = read->core.n_cigar-1; j >= 0 && prefix_len > 0; j--) {
        uint32_t op = cigar[j];
        char opchr = bam_cigar_opchr(op);
        int oplen = bam_cigar_oplen(op);
        if (opchr == 'I' || opchr == 'D') {
            indels += std::min(oplen, prefix_len);
        }
        if (opchr != 'D') {
            prefix_len -= oplen;
        }

    }

    prefix_len = original_prefix_len;
    if (bam_cigar_opchr(cigar[0]) == 'S') prefix_len -= bam_cigar_oplen(cigar[0]) ;

    char* md_tag = bam_aux2Z(bam_aux_get(read, "MD"));
    int md_tag_len = strlen(md_tag);
    int n = 0, m = 1;
    for (int i = md_tag_len-1; i >= 0 && prefix_len > 0; i--) {
        if (md_tag[i] >= '0' && md_tag[i] <= '9') {
            n += m * (md_tag[i]-'0');
            m *= 10;
        } else {
            prefix_len -= n;
            if (prefix_len <= 0) break;
            n = 0, m = 1;

            if (i == 0 || (md_tag[i] != '^' && md_tag[i-1] >= '0' && md_tag[i-1] <= '9')) {
                mismatches++;
            }
        }
    }

    return {mismatches, indels};
}

int64_t get_AS_tag(bam1_t* read) {
    uint8_t* aux_get = bam_aux_get(read, "AS");
    return bam_aux2i(aux_get);
}

int get_aligned_portion_len(bam1_t* read) {
    return read->core.l_qseq - get_left_clip_size(read) - get_right_clip_size(read);
}

void rc(std::string& read) {
    int len = read.length();
    for (int i = 0; i < len/2; i++) {
        std::swap(read[i], read[len-i-1]);
    }
    for (int i = 0; i < len; i++) {
        if (read[i] == 'A') read[i] = 'T';
        else if (read[i] == 'C') read[i] = 'G';
        else if (read[i] == 'G') read[i] = 'C';
        else if (read[i] == 'T') read[i] = 'A';
        else read[i] = 'N';
    }
}

struct open_samFile_t {
    samFile* file;
    bam_hdr_t* header;
    hts_idx_t* idx;

    open_samFile_t() {}

    open_samFile_t(samFile* file, bam_hdr_t* header, hts_idx_t* idx) : file(file), header(header), idx(idx) {}
};

open_samFile_t* open_samFile(std::string fname_str, bool index_file = false) {
    const char* fname = fname_str.c_str();
    open_samFile_t* sam_file = new open_samFile_t;
    sam_file->file = sam_open(fname, "r");
    if (sam_file->file == NULL) {
        throw "Could not open " + std::string(fname);
    }

    if (index_file) {
        int code = sam_index_build(fname, 0);
        if (code != 0) {
            throw "Cannot index " + std::string(fname);
        }
    }

    sam_file->idx = sam_index_load(sam_file->file, sam_file->file->fn);
    if (sam_file->idx == NULL) {
        throw "Unable to open index for " + std::string(fname);
    }

    sam_file->header = sam_hdr_read(sam_file->file);
    if (sam_file->header == NULL) {
        throw "Unable to open header for " + std::string(fname);
    }

    return sam_file;
}

void close_samFile(open_samFile_t* f) {
    hts_idx_destroy(f->idx);
    bam_hdr_destroy(f->header);
    sam_close(f->file);
    delete f;
}

bool is_homopolymer(const char* seq, int len) {
	int a = 0, c = 0, g = 0, t = 0;
	for (int i = 0; i < len; i++) {
		char b = std::toupper(seq[i]);
		if (b == 'A') a++;
		else if (b == 'C') c++;
		else if (b == 'G') g++;
		else if (b == 'T') t++;
	}
	return max(a, c, g, t)/double(a+c+g+t) >= 0.8;
}

bool is_homopolymer(std::string& seq) {
	return is_homopolymer(seq.data(), seq.length());
}

samFile* get_writer(std::string path, bam_hdr_t* header, bool bam = true) {
    samFile* writer = sam_open(path.c_str(), bam ? "wb" : "w");
    if (sam_hdr_write(writer, header) != 0) {
        throw "Could not write file " + path;
    }
    return writer;
}

#endif //SURVTYPER_SAM_UTILS_H
