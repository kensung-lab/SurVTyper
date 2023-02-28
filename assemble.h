#include <vector>

#include "utils.h"

extern config_t config;

struct edge_t {
	int next, score, overlap;

	edge_t() : next(0), score(0), overlap(0) {}
	edge_t(int next, int score, int overlap) : next(next), score(score), overlap(overlap) {}
};

struct path_permission_t {
	bool can_start_path, can_end_path;
};
struct seq_w_pp_t {
	std::string seq;
	path_permission_t clip_pair;

	seq_w_pp_t() : seq(), clip_pair() {}
	seq_w_pp_t(std::string& seq, bool can_start_path, bool can_end_path) : seq(seq) {
		clip_pair.can_start_path = can_start_path;
		clip_pair.can_end_path = can_end_path;
	}
};

bool accept(StripedSmithWaterman::Alignment& aln, std::string& log, double max_seq_error = 1.0, std::string qual_ascii = "", int min_avg_base_qual = 0) {
	// the 'N' we insterted for padding are counted as mismatches, so they must be removed from the count
	int lc_size = get_left_clip_size(aln), rc_size = get_right_clip_size(aln);
	if (lc_size > config.clip_penalty && rc_size > config.clip_penalty) return false; // do not accept if both sides clipped

	int left_padding_matched = std::max(0, config.clip_penalty - lc_size);
	int right_padding_matched = std::max(0, config.clip_penalty - rc_size);
	int mismatches = aln.mismatches - left_padding_matched - right_padding_matched;
	double mismatch_rate = double(mismatches)/(aln.query_end-aln.query_begin-left_padding_matched-right_padding_matched);
	return lc_size <= config.clip_penalty && rc_size <= config.clip_penalty && mismatch_rate <= max_seq_error;
};
bool accept(StripedSmithWaterman::Alignment& aln, double max_seq_error = 1.0, std::string qual_ascii = "", int min_avg_base_qual = 0) {
	std::string temp;
	return accept(aln, temp, max_seq_error, qual_ascii, min_avg_base_qual);
}

void add_alignment(std::string& reference, std::string& query, std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> >& accepted_alns,
		std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> >& rejected_alns, StripedSmithWaterman::Aligner& aligner) {
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	std::string padded_query = std::string(config.clip_penalty, 'N') + query + std::string(config.clip_penalty, 'N');
	aligner.Align(padded_query.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
	aln.ref_begin += config.clip_penalty - get_left_clip_size(aln);
	aln.ref_end -= config.clip_penalty - get_right_clip_size(aln);
	if (accept(aln)) {
		accepted_alns.push_back({query, aln});
	} else {
		rejected_alns.push_back({query, aln});
	}
}

void remap_consensus_to_reference(std::string& junction, char* ref, int ref_len, StripedSmithWaterman::Aligner& aligner,
		std::string& remapped_cigar, int& remapped_ins_len, bool& low_qual_alt_allele) {
	StripedSmithWaterman::Alignment aln;
	StripedSmithWaterman::Filter filter;
	std::string padded_junction = std::string(config.clip_penalty, 'N') + junction + std::string(config.clip_penalty, 'N');
	aligner.Align(padded_junction.c_str(), ref, ref_len, filter, &aln, 0);

	remapped_ins_len = 0;
	remapped_cigar = "";
	low_qual_alt_allele = false;

	if (aln.query_begin > 0 && aln.query_end < padded_junction.length()-1) {
		low_qual_alt_allele = true;
	} else if (aln.query_begin == 0 && aln.query_end == padded_junction.length()-1) {
		remapped_ins_len = padded_junction.length()-(aln.ref_end-aln.ref_begin+1);
		remapped_cigar = aln.cigar_string;
	} else {
		if (aln.query_begin == 0) {
			StripedSmithWaterman::Alignment rh_aln;
			std::string right_half = padded_junction.substr(aln.query_end);
			aligner.Align(right_half.c_str(), ref+aln.ref_end, ref_len-aln.ref_end, filter, &rh_aln, 0);
			if (get_right_clip_size(rh_aln) > 0) low_qual_alt_allele = true; // bad anchor
			else if (rh_aln.query_end-rh_aln.query_begin < config.min_clip_size+config.clip_penalty) low_qual_alt_allele = true; // anchor is too small
			remapped_ins_len = std::max(0, get_left_clip_size(rh_aln) - rh_aln.ref_begin);
			remapped_cigar = aln.cigar_string + "," + rh_aln.cigar_string + "," + std::to_string(rh_aln.ref_begin);
		} else if (aln.query_end == padded_junction.length()-1) {
			StripedSmithWaterman::Alignment lh_aln;
			std::string left_half = padded_junction.substr(0, aln.query_begin);
			aligner.Align(left_half.c_str(), ref, aln.ref_begin, filter, &lh_aln, 0);
			if (get_left_clip_size(lh_aln) > 0) low_qual_alt_allele = true; // bad anchor
			else if (lh_aln.query_end-lh_aln.query_begin < config.min_clip_size+config.clip_penalty) low_qual_alt_allele = true; // anchor is too small
			remapped_ins_len = std::max(0, get_right_clip_size(lh_aln) - (aln.ref_begin-1-lh_aln.ref_end));
			remapped_cigar = lh_aln.cigar_string + "," + aln.cigar_string + "," + std::to_string(aln.ref_begin-1-lh_aln.ref_end);
		}
	}
}

void build_aln_guided_graph(std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> >& alns, std::vector<int>& out_edges,
		std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev) {
	std::sort(alns.begin(), alns.end(),
			[](const std::pair<std::string, StripedSmithWaterman::Alignment>& aln1, const std::pair<std::string, StripedSmithWaterman::Alignment>& aln2) {
		return aln1.second.ref_begin < aln2.second.ref_begin;
	});

	for (int i = 0; i < alns.size(); i++) {
		for (int j = i+1; j < alns.size() && alns[i].second.ref_end-alns[j].second.ref_begin >= config.min_clip_size; j++) {
			suffix_prefix_aln_t spa = aln_suffix_prefix(alns[i].first, alns[j].first, 1, -4, config.max_seq_error, config.min_clip_size);
			if (spa.overlap) {
				out_edges[i]++;
				l_adj[i].push_back({j, spa.score, spa.overlap});
				l_adj_rev[j].push_back({i, spa.score, spa.overlap});
			}
		}
	}

	// two major differences with the regular assembly:
	// 1 - by how the graph is defined, no cycle is possible here
	// 2 - in regular assembly, we only report contigs made of at last 2 reads. Here we report even single reads
}

std::vector<int> find_rev_topological_order(int n, std::vector<int>& out_edges, std::vector<std::vector<edge_t> >& l_adj_rev) {

	std::queue<int> sinks;
	for (int i = 0; i < n; i++) {
		if (!out_edges[i]) sinks.push(i);
	}

	std::vector<int> rev_topological_order;
	while (!sinks.empty()) {
		int s = sinks.front();
		sinks.pop();
		rev_topological_order.push_back(s);
		for (edge_t& e : l_adj_rev[s]) {
			out_edges[e.next]--;
			if (out_edges[e.next] == 0) sinks.push(e.next);
		}
	}
	return rev_topological_order;
}

void correct_contig(std::string& contig, std::vector<std::string>& reads, StripedSmithWaterman::Aligner& harsh_aligner) {
	std::vector<int> As(contig.length()), Cs(contig.length()), Gs(contig.length()), Ts(contig.length());
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	for (std::string& read : reads) {
		std::string padded_read = std::string(config.clip_penalty, 'N') + read + std::string(config.clip_penalty, 'N');
		harsh_aligner.Align(padded_read.c_str(), contig.c_str(), contig.length(), filter, &aln, 0);
		if (accept(aln)) {
			int left_padding_aligned = config.clip_penalty - get_left_clip_size(aln);
			for (int i = 0; i < read.length(); i++) {
				char c = read[i];
				if (c == 'A') As[i+aln.ref_begin+left_padding_aligned]++;
				else if (c == 'C') Cs[i+aln.ref_begin+left_padding_aligned]++;
				else if (c == 'G') Gs[i+aln.ref_begin+left_padding_aligned]++;
				else if (c == 'T') Ts[i+aln.ref_begin+left_padding_aligned]++;
			}
		}
	}

	for (int i = 0; i < contig.length(); i++) {
		int max_freq = max(As[i], Cs[i], Gs[i], Ts[i]);
		if (max_freq == 0) continue;
		if (max_freq == As[i]) contig[i] = 'A';
		else if (max_freq == Cs[i]) contig[i] = 'C';
		else if (max_freq == Gs[i]) contig[i] = 'G';
		else if (max_freq == Ts[i]) contig[i] = 'T';
	}
}

void build_graph(std::vector<std::string>& read_seqs, std::vector<int>& order, std::vector<int>& out_edges,
		std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev,
		int max_mismatches, int min_overlap, std::stringstream& ss_graph) {

	int n = read_seqs.size();

	for (int i = 0; i < n; i++) {
		for (int j = i+1; j < n; j++) {
			std::string& s1 = read_seqs[i];
			std::string& s2 = read_seqs[j];

			int max_overlap = std::min(s1.length(), s2.length())-1;
			suffix_prefix_aln_t spa1 = aln_suffix_prefix(s1, s2, 1, -4, 1.0, min_overlap, max_overlap, max_mismatches);
			bool spa1_homopolymer = is_homopolymer(s2.c_str(), spa1.overlap);
			suffix_prefix_aln_t spa2 = aln_suffix_prefix(s2, s1, 1, -4, 1.0, min_overlap, max_overlap, max_mismatches);
			bool spa2_homopolymer = is_homopolymer(s1.c_str(), spa2.overlap);
			if (spa1.overlap && spa2.overlap) {
				if (spa1.score >= spa2.score && order[i] <= order[j]) {
					ss_graph << i << " -> " << j << " " << spa1.score << " " << spa1.overlap << " " << spa1_homopolymer << std::endl;
					if (!spa1_homopolymer) {
						out_edges[i]++;
						l_adj[i].push_back({j, spa1.score, spa1.overlap});
						l_adj_rev[j].push_back({i, spa1.score, spa1.overlap});
					}
				} else if (spa1.score < spa2.score && order[j] <= order[i]) {
					ss_graph << j << " -> " << i << " " << spa2.score << " " << spa2.overlap << " " << spa2_homopolymer << std::endl;
					if (!spa2_homopolymer) {
						out_edges[j]++;
						l_adj[j].push_back({i, spa2.score, spa2.overlap});
						l_adj_rev[i].push_back({j, spa2.score, spa2.overlap});
					}
				}
			} else if (spa1.overlap && order[i] <= order[j]) {
				ss_graph << i << " -> " << j << " " << spa1.score << " " << spa1.overlap << " " << spa1_homopolymer << std::endl;
				if (!spa1_homopolymer) {
					out_edges[i]++;
					l_adj[i].push_back({j, spa1.score, spa1.overlap});
					l_adj_rev[j].push_back({i, spa1.score, spa1.overlap});
				}
			} else if (spa2.overlap && order[j] <= order[i]) {
				ss_graph << j << " -> " << i << " " << spa2.score << " " << spa2.overlap << " " << spa2_homopolymer << std::endl;
				if (!spa2_homopolymer) {
					out_edges[j]++;
					l_adj[j].push_back({i, spa2.score, spa2.overlap});
					l_adj_rev[i].push_back({j, spa2.score, spa2.overlap});
				}
			}
		}
	} ss_graph << std::endl;
}

std::vector<std::string> assemble_reads(std::vector<seq_w_pp_t>& left_stable_read_seqs, std::vector<seq_w_pp_t>& unstable_read_seqs,
		std::vector<seq_w_pp_t>& right_stable_read_seqs, StripedSmithWaterman::Aligner& harsh_aligner, config_t config, std::stringstream& ss_graph) {

	std::vector<std::string> read_seqs;
	std::vector<path_permission_t> path_permissions;
	std::vector<int> order;
	for (seq_w_pp_t& s : left_stable_read_seqs) {
		read_seqs.push_back(s.seq);
		path_permissions.push_back(s.clip_pair);
		order.push_back(1);
	}
	for (seq_w_pp_t& s : unstable_read_seqs) {
		read_seqs.push_back(s.seq);
		path_permissions.push_back(s.clip_pair);
		order.push_back(2);
	}
	for (seq_w_pp_t& s : right_stable_read_seqs) {
		read_seqs.push_back(s.seq);
		path_permissions.push_back(s.clip_pair);
		order.push_back(3);
	}

	int n = read_seqs.size();
	std::vector<int> out_edges(n);
	std::vector<std::vector<edge_t> > l_adj(n), l_adj_rev(n);

	for (int i = 0; i < n; i++) {
		ss_graph << i << " " << read_seqs[i] << " " << order[i] << std::endl;
	} ss_graph << std::endl;

	build_graph(read_seqs, order, out_edges, l_adj, l_adj_rev, 1, config.min_clip_size, ss_graph);

	std::vector<int> rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);

	if (rev_topological_order.size() < n) {
		build_graph(read_seqs, order, out_edges, l_adj, l_adj_rev, 0.0, config.min_clip_size, ss_graph);

		int min_overlap = config.min_clip_size;
		for (; min_overlap <= config.read_len/2; min_overlap += 10) {
			ss_graph << "HAS A CYCLE, TRYING NO MISMATCHES AND min_overlap = " << min_overlap << std::endl;

			for (int i = 0; i < n; i++) {
				l_adj[i].erase(std::remove_if(l_adj[i].begin(), l_adj[i].end(),
						[&min_overlap](edge_t& e) { return e.overlap < min_overlap; }), l_adj[i].end());
				l_adj_rev[i].erase(std::remove_if(l_adj_rev[i].begin(), l_adj_rev[i].end(),
						[&min_overlap](edge_t& e) { return e.overlap < min_overlap; }), l_adj_rev[i].end());
				out_edges[i] = l_adj[i].size();
			}

			rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);

			if (rev_topological_order.size() == n) {
				ss_graph << "CYCLE DISAPPEARED!" << std::endl;
				break;
			}
		}
		if (rev_topological_order.size() < n) {
			return {"HAS_CYCLE"};
		}
	}

	std::vector<std::string> assembled_sequences;
	std::vector<bool> used(n);
	while (true) {
		// compute longest paths
		std::vector<int> best_scores(n);
		std::vector<edge_t> best_edges(n);
		for (int i : rev_topological_order) {
			if (used[i]) continue;
			if (best_scores[i] == 0 && !path_permissions[i].can_end_path) continue; // sink and cannot end path => discard
			for (edge_t& e : l_adj_rev[i]) {
				if (best_scores[e.next] < e.score + best_scores[i]) {
					best_scores[e.next] = e.score + best_scores[i];
					best_edges[e.next] = {i, e.score, e.overlap};
				}
			}
		}

		int best_score = 0, curr_vertex = 0;
		for (int i = 0; i < best_scores.size(); i++) {
			if (path_permissions[i].can_start_path && best_score < best_scores[i]) {
				best_score = best_scores[i];
				curr_vertex = i;
			}
		}
		if (best_score == 0) break;

		std::string assembled_sequence = read_seqs[curr_vertex];
		std::vector<std::string> used_reads; // track reads used to build this contig, so that we can use them for correction
		used_reads.push_back(read_seqs[curr_vertex]);
		while (best_edges[curr_vertex].overlap) {
			ss_graph << curr_vertex << " -> " << best_edges[curr_vertex].next << " , " << best_edges[curr_vertex].overlap << " ";
			ss_graph << read_seqs[curr_vertex] << std::endl;
			used[curr_vertex] = true;
			int overlap = best_edges[curr_vertex].overlap;
			curr_vertex = best_edges[curr_vertex].next;
			assembled_sequence += read_seqs[curr_vertex].substr(overlap);
			used_reads.push_back(read_seqs[curr_vertex]);
		}
		used[curr_vertex] = true;
		correct_contig(assembled_sequence, used_reads, harsh_aligner);
		assembled_sequences.push_back(assembled_sequence);
		ss_graph << "-> " << read_seqs[curr_vertex] << std::endl;
		ss_graph << assembled_sequence << std::endl << std::endl;
	}

	for (int i = 0; i < read_seqs.size(); i++) {
		if (!used[i]) assembled_sequences.push_back(read_seqs[i]);
	}

	return assembled_sequences;
}

std::vector<std::string> generate_reference_guided_contigs(std::string reference, std::vector<std::string>& seqs,
		StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner,
		std::vector<StripedSmithWaterman::Alignment>& consensus_contigs_alns, std::string& consensus_log) {

	std::stringstream ss;
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;

	std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> > accepted_alns, rejected_alns;
	for (std::string& seq : seqs) {
		add_alignment(reference, seq, accepted_alns, rejected_alns, aligner);
	}

	ss << "ACCEPTED ALNS" << std::endl;
	for (auto& e : accepted_alns) {
		ss << e.first << " " << e.second.ref_begin << " " << e.second.cigar_string << std::endl;
	} ss << std::endl;

	int n = accepted_alns.size();
	std::vector<int> out_edges(n);
	std::vector<std::vector<edge_t> > l_adj(n), l_adj_rev(n);
	build_aln_guided_graph(accepted_alns, out_edges, l_adj, l_adj_rev);

	ss << "GRAPH" << std::endl;
	for (int i = 0; i < l_adj.size(); i++) {
		for (auto& e : l_adj[i]) {
			ss << i << " -> " << e.next << ", " << e.score << std::endl;
		}
	} ss << std::endl;

	std::vector<int> rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);

	std::vector<std::string> assembled_sequences;
	std::vector<bool> used(n);
	while (true) {
		// compute longest paths
		std::vector<int> best_scores(n);
		std::vector<edge_t> best_edges(n);
		for (int i : rev_topological_order) {
			if (used[i]) continue;
			for (edge_t& e : l_adj_rev[i]) {
				if (!used[e.next] && best_scores[e.next] < e.score + best_scores[i]) {
					best_scores[e.next] = e.score + best_scores[i];
					best_edges[e.next] = {i, e.score, e.overlap};
				}
			}
		}

		int best_score = 0, curr_vertex = 0;
		for (int i = 0; i < best_scores.size(); i++) {
			if (best_score < best_scores[i]) {
				best_score = best_scores[i];
				curr_vertex = i;
			}
		}
		if (best_score == 0) break;

		std::string assembled_sequence = accepted_alns[curr_vertex].first;
		std::vector<std::string> used_reads; // track reads used to build this contig, so that we can use them for correction
		used_reads.push_back(accepted_alns[curr_vertex].first);
		ss << curr_vertex << " ";
		while (best_edges[curr_vertex].overlap) {
			used[curr_vertex] = true;
			int overlap = best_edges[curr_vertex].overlap;
			curr_vertex = best_edges[curr_vertex].next;
			assembled_sequence += accepted_alns[curr_vertex].first.substr(overlap);
			used_reads.push_back(accepted_alns[curr_vertex].first);
			ss << curr_vertex << " ";
		}
		used[curr_vertex] = true;

		std::string corrected_assembled_sequence = assembled_sequence;
		correct_contig(corrected_assembled_sequence, used_reads, harsh_aligner);
		assembled_sequences.push_back(corrected_assembled_sequence);
		ss << corrected_assembled_sequence << std::endl;
	}
	for (int i = 0; i < n; i++) {
		if (!used[i]) {
			assembled_sequences.push_back(accepted_alns[i].first);
		}
	}

	// retain assembled sequences that align without clipping and do not overlap a higher rated sequence
	std::vector<std::string> retained_assembled_sequences;
	for (std::string& assembled_sequence : assembled_sequences) {
		aligner.Align(assembled_sequence.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		if (!accept(aln)) continue;

		bool overlaps = false;
		for (StripedSmithWaterman::Alignment& existing_aln : consensus_contigs_alns) {
			if (std::max(existing_aln.ref_begin, aln.ref_begin) <= std::min(existing_aln.ref_end, aln.ref_end)) {
				overlaps = true;
				break;
			}
		}
		if (!overlaps) {
			consensus_contigs_alns.push_back(aln);
			retained_assembled_sequences.push_back(assembled_sequence);
		}
	}

	if (retained_assembled_sequences.empty()) return {};

	// try scaffolding using rejected reads
	std::stringstream ss_graph;
	std::vector<seq_w_pp_t> rejected_alns_w_pp, temp1, temp2; // TODO: separate reads into stable left etc. For now we use no partial ordering
	for (auto& e : rejected_alns) rejected_alns_w_pp.push_back({e.first, true, true});
	std::vector<std::string> scaffolds = assemble_reads(temp1, rejected_alns_w_pp, temp2, harsh_aligner, config, ss_graph);

	ss << "Assembled " << n << " reads into " << assembled_sequences.size() << " sequences." << std::endl;
	for (std::string a : assembled_sequences) {
		aligner.Align(a.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		ss << a.length() << "," << aln.ref_begin << "-" << aln.ref_end << "," << accept(aln) << " ";
	}
	ss << std::endl;
	for (int i = 0; i < retained_assembled_sequences.size(); i++) {
		StripedSmithWaterman::Alignment& aln = consensus_contigs_alns[i];
		ss << retained_assembled_sequences[i].length() << "," << aln.ref_begin << "-" << aln.ref_end << "," << accept(aln) << "," << aln.cigar_string << " ";
	}
	ss << std::endl << std::endl;

	ss << "RETAINED CONTIGS: " << retained_assembled_sequences.size() <<  std::endl;
	for (std::string& a : retained_assembled_sequences) ss << a << std::endl;
	ss << std::endl;

	ss << "REJECTED READS: " << rejected_alns.size() << std::endl;
	for (auto& e : rejected_alns) {
		ss << e.first << " " << e.second.ref_begin << " " << e.second.cigar_string << std::endl;
	} ss << std::endl;

	ss << "SCAFFOLDS: " << scaffolds.size() << std::endl;
	for (std::string& a : scaffolds) ss << a << std::endl;
	ss << std::endl;


	std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> > contigs_sorted_by_pos;
	for (int i = 0; i < consensus_contigs_alns.size(); i++) {
		contigs_sorted_by_pos.push_back({retained_assembled_sequences[i], consensus_contigs_alns[i]});
	}
	std::sort(contigs_sorted_by_pos.begin(), contigs_sorted_by_pos.end(),
			[](std::pair<std::string, StripedSmithWaterman::Alignment>& p1, std::pair<std::string, StripedSmithWaterman::Alignment>& p2) {
		return p1.second.ref_begin < p2.second.ref_begin;
	});

	std::vector<int> linked(contigs_sorted_by_pos.size()-1, -1);
	std::vector<std::pair<int, int> > link_overlaps(contigs_sorted_by_pos.size()-1);
	for (int i = 0; i < scaffolds.size(); i++) {
		std::string& scaffold = scaffolds[i];
		int best_link = -1, best_link_w = 0;
		std::pair<int, int> best_link_overlap;
		for (int j = 0; j < contigs_sorted_by_pos.size()-1; j++) {
			if (linked[j] >= 0) continue;

			suffix_prefix_aln_t spa1 = aln_suffix_prefix(contigs_sorted_by_pos[j].first, scaffold, 1, -4, config.max_seq_error, config.min_clip_size);
			suffix_prefix_aln_t spa2 = aln_suffix_prefix(scaffold, contigs_sorted_by_pos[j+1].first, 1, -4, config.max_seq_error, config.min_clip_size);
			if (spa1.overlap && spa2.overlap && best_link_w < spa1.score+spa2.score) {
				best_link = j, best_link_w = spa1.score+spa2.score;
				best_link_overlap = {spa1.overlap, spa2.overlap};
			}
		}

		if (best_link >= 0) {
			linked[best_link] = i;
			link_overlaps[best_link] = best_link_overlap;
		}
	}

	std::vector<std::string> scaffolded_sequences;
	std::string curr_seq = contigs_sorted_by_pos[0].first;
	for (int i = 1; i < contigs_sorted_by_pos.size(); i++) {
		if (linked[i-1] == -1) {
			scaffolded_sequences.push_back(curr_seq);
			curr_seq = contigs_sorted_by_pos[i].first;
		} else {
			std::string link = scaffolds[linked[i-1]];
			auto& lo = link_overlaps[i-1];
			link = link.substr(lo.first);
			curr_seq += link + contigs_sorted_by_pos[i].first.substr(lo.second);
			ss << "Linking " << contigs_sorted_by_pos[i-1].first << " and " << contigs_sorted_by_pos[i].first << " with " << scaffolds[linked[i-1]] << std::endl;
			ss << lo.first << " " << link.length()-lo.first-lo.second << std::endl;
		}
	}
	scaffolded_sequences.push_back(curr_seq);

	// Use scaffolds to left-extend the first sequence
	int best_overlap = 0, best_overlap_i = -1;
	for (int i = 0; i < scaffolds.size(); i++) {
		suffix_prefix_aln_t spa = aln_suffix_prefix(scaffolds[i], scaffolded_sequences[0], 1, -4, config.max_seq_error, config.min_clip_size);
		if (best_overlap < spa.overlap) {
			best_overlap = spa.overlap;
			best_overlap_i = i;
		}
	}
	if (best_overlap_i >= 0) {
		scaffolded_sequences[0] = scaffolds[best_overlap_i].substr(0, scaffolds.size()-best_overlap) + scaffolded_sequences[0];
		scaffolds[best_overlap_i] = "";
	}

	// Use scaffolds to right-extend the last sequence
	best_overlap = 0, best_overlap_i = -1;
	for (int i = 0; i < scaffolds.size(); i++) {
		suffix_prefix_aln_t spa = aln_suffix_prefix(scaffolded_sequences[scaffolded_sequences.size()-1], scaffolds[i], 1, -4, config.max_seq_error, config.min_clip_size);
		if (best_overlap < spa.overlap) {
			best_overlap = spa.overlap;
			best_overlap_i = i;
		}
	}
	if (best_overlap_i >= 0) {
		scaffolded_sequences[scaffolded_sequences.size()-1] = scaffolded_sequences[scaffolded_sequences.size()-1] +
				scaffolds[best_overlap_i].substr(best_overlap);
		scaffolds[best_overlap_i] = "";
	}

	std::vector<StripedSmithWaterman::Alignment> scaffolded_seqs_alns;
	ss << "SCAFFOLDED SEQUENCES: " << scaffolded_sequences.size() << std::endl;
	bool scaffolding_failed = false;
	for (std::string s : scaffolded_sequences) {
		aligner.Align(s.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		if (!accept(aln)) {
			scaffolding_failed = true;
			break;
		}
		scaffolded_seqs_alns.push_back(aln);
		ss << s << " " << s.length() << "," << aln.ref_begin << "-" << aln.ref_end << "," << accept(aln) << std::endl;
	} ss << std::endl;

	if (!scaffolding_failed)
	for (int i = 0; i < scaffolded_seqs_alns.size(); i++) {
		for (int j = i+1; j < scaffolded_seqs_alns.size(); j++) {
			if (std::max(scaffolded_seqs_alns[i].ref_begin, scaffolded_seqs_alns[j].ref_begin) <= std::min(scaffolded_seqs_alns[i].ref_end, scaffolded_seqs_alns[j].ref_end)) {
				scaffolding_failed = true;
				break;
			}
		}
		if (scaffolding_failed) break;
	}

	ss << reference << std::endl;

	consensus_log = ss.str();

	if (!scaffolding_failed) scaffolded_seqs_alns.swap(consensus_contigs_alns);
	return scaffolding_failed ? retained_assembled_sequences : scaffolded_sequences;
}

std::string generate_reference_guided_consensus(std::string reference, std::vector<std::string>& seqs,
		StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner,
		std::vector<StripedSmithWaterman::Alignment>& consensus_contigs_alns, std::string& consensus_log, bool include_ref_padding = true) {

	if (seqs.empty()) return reference;

	std::vector<std::string> consensus_contigs = generate_reference_guided_contigs(reference, seqs, aligner, harsh_aligner, consensus_contigs_alns, consensus_log);

	if (consensus_contigs_alns.empty()) return reference;

	std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> > contigs_sorted_by_pos;
	for (int i = 0; i < consensus_contigs.size(); i++) {
		contigs_sorted_by_pos.push_back({consensus_contigs[i], consensus_contigs_alns[i]});
	}
	std::sort(contigs_sorted_by_pos.begin(), contigs_sorted_by_pos.end(),
			[](std::pair<std::string, StripedSmithWaterman::Alignment>& p1, std::pair<std::string, StripedSmithWaterman::Alignment>& p2) {
		return p1.second.ref_begin < p2.second.ref_begin;
	});

	std::string corrected_reference;
	if (include_ref_padding) corrected_reference = reference.substr(0, contigs_sorted_by_pos[0].second.ref_begin);
	corrected_reference += contigs_sorted_by_pos[0].first;
	for (int i = 1; i < contigs_sorted_by_pos.size(); i++) {
		corrected_reference += reference.substr(contigs_sorted_by_pos[i-1].second.ref_end+1, contigs_sorted_by_pos[i].second.ref_begin-contigs_sorted_by_pos[i-1].second.ref_end);
		corrected_reference += contigs_sorted_by_pos[i].first;
	}
	if (include_ref_padding) corrected_reference += reference.substr(contigs_sorted_by_pos[contigs_sorted_by_pos.size()-1].second.ref_end+1);
	consensus_log += "\nCONSENSUS REF: " + corrected_reference + "\n";

	return corrected_reference;
}

std::string generate_consensus(std::vector<std::string>& seqs, std::vector<bool>& seqs_lc, std::vector<bool>& seqs_rc,
		StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner,
		std::vector<StripedSmithWaterman::Alignment>& consensus_contigs_alns, std::string& consensus_log) {

	if (seqs.empty()) return "";

	std::vector<seq_w_pp_t> seqs_w_pp, temp1, temp2;
	for (int i = 0; i < seqs.size(); i++) seqs_w_pp.push_back({seqs[i], !seqs_lc[i], !seqs_rc[i]});
	std::stringstream ss_graph;
	std::vector<std::string> contigs = assemble_reads(temp1, seqs_w_pp, temp2, harsh_aligner, config, ss_graph);
	std::string primary_assembled_sequence = contigs[0];
	return primary_assembled_sequence;
}

