#ifndef ISO_EXP_MODEL_H
///Define this macro to prevent from including this header file more than once.
#define ISO_EXP_MODEL_H

#include "isoform_opt.h"
#include "math_utils.h"
#include "fasta.h"

class iso_exp_model_2 { // for expression_analysis
public:
	class _isoform_str {
	public:
		string name;
		int length;
		_isoform_str():name(""),length(0){};
	};
	class _gene_str {
	public:
		string name;
		int unique_count;
		double total_count;
		vector<_isoform_str> isoforms;
		vector<double> isoform_expression;
		map<string, int> isoform_map;
		map<pair<vector<int>, vector<int>>, int> collapsed_reads;
		vector<vector<double>> A;
		vector<double> N;
		vector<double> W;
		_gene_str():name(""),unique_count(0),total_count(0){
			isoforms.clear();
			isoform_map.clear();
			collapsed_reads.clear();
			A.clear();
			N.clear();
			W.clear();
		};
		void accumulate_read(const vector<int> &isoform_indices, const vector<int> &start_pos, const vector<int> &insert_lens) {
			__ASSERT(isoform_indices.size() == insert_lens.size(), "error: illegal arguments for _gene_str::accumulate_read()\n");
			vector<int> st_pos(isoforms.size(), 0);
			vector<int> ins_lens(isoforms.size(), 0);
			for (int i = 0; i < (int)isoform_indices.size(); i++) {
				st_pos[isoform_indices[i]] = start_pos[i];
				ins_lens[isoform_indices[i]] = insert_lens[i];
			}
			if (collapsed_reads.count(pair<vector<int>,vector<int>>(st_pos, ins_lens)) == 0) {
				collapsed_reads[pair<vector<int>,vector<int>>(st_pos, ins_lens)] = 1;
			} else {
				collapsed_reads[pair<vector<int>,vector<int>>(st_pos, ins_lens)]++;
			}
			unique_count++;
			total_count+=1;
		};
		void get_rates(int total_supporting_reads, map<int, double> &p_table, int read_length) {
			int I = (int)isoforms.size();
			A.clear();
			N.clear();
			W.clear();
			for (auto itr = collapsed_reads.begin(); itr != collapsed_reads.end(); itr++) {
				vector<double> rate(I, 0);
				bool positive_rate = false;
				for (int i = 0; i < I; i++) {
					int insert_len = itr->first.second[i];
					if (p_table.count(insert_len) > 0) {
						rate[i] = p_table[insert_len] * total_supporting_reads / 1000 / 1000000;
						positive_rate = true;
					} else rate[i] = 0;
				}
				if (positive_rate) {
					int index = -1;
					for (int j = 0; j < (int)A.size(); j++) {
						if (vector_parallel(A[j], rate)) {
							index = j;
							N[j] += itr->second;
							for (int i = 0; i < I; i++) A[j][i] += rate[i];
							break;
						}
					}
					if (index == -1) {
						A.push_back(rate);
						N.push_back(itr->second);
					}
				}
			}
			for (int i = 0; i < I; i++) {
				double temp_sum = 0;
				for (auto itr = p_table.begin(); itr != p_table.end(); itr++) {
					if (itr->first <= isoforms[i].length) {
						temp_sum += total_supporting_reads * itr->second * (isoforms[i].length - read_length + 1) / 1000 / 1000000;
					}
				}
				W.push_back(temp_sum);
			}
		};
		void compute() {
			if (N.empty()) {
				isoform_expression.resize(isoforms.size());
				for (int i = 0; i < (int)isoform_expression.size(); i++) {
					isoform_expression[i] = 0;
				}
			} else {
				solve_likelihood_t(A, N, isoform_expression, &W);
			}
		};
		void output_data(ofstream &ofs) {
			int I = (int)isoforms.size();
			int J = (int)N.size();
			ofs << name << "\t" << I << "\t" << J << endl;
			for (int i = 0; i < I; i++) {
				ofs << isoforms[i].name << "\t";
			}
			ofs << endl;
			for (int j = 0; j < J; j++) {
				for (int i = 0; i < I; i++) {
					ofs << A[j][i] << "\t";
				}
				ofs << N[j] << endl;
			}
			for (int i = 0; i < I; i++) {
				ofs << W[i] << "\t";
			}
			ofs << endl;
		};
		void output_exp(ofstream &ofs) {		
			int I = (int)isoforms.size();
			double gene_exp = 0;
			for (int i = 0; i < I; i++) {
				gene_exp += isoform_expression[i];
			}
			ofs << " " << name << "\t" << gene_exp << "\t" << I << "\t";
			for (int i = 0; i < I; i++) {
				ofs << isoforms[i].name;
				if (i < I-1) ofs << ",";
			}
			ofs << "\t";
			for (int i = 0; i < I; i++) {
				ofs << isoform_expression[i];
				if (i < I-1) ofs << ",";
			}
			ofs << endl;
		}
	};
	vector<_gene_str> genes;
	map<string, int> gene_map;
	iso_exp_model_2() {
		genes.clear();
		gene_map.clear();
	};
	void build_from_fasta(const fasta &myfasta) {
		for (int i = 0; i < (int)myfasta.tags.size(); i++) {
			size_t pos = myfasta.tags[i].find_first_of("$$");
			string gene_name = (pos == string::npos) ? myfasta.tags[i] : myfasta.tags[i].substr(0,pos);
			string isoform_name = (pos == string::npos) ? myfasta.tags[i] : myfasta.tags[i].substr(pos+2);
			int gene_index = -1;
			if (gene_map.count(gene_name) > 0) {
				gene_index = gene_map[gene_name];
			} else {
				_gene_str new_gene;
				new_gene.name = gene_name;
				genes.push_back(new_gene);
				gene_index = (int)genes.size()-1;
				gene_map[gene_name] = gene_index;
			}
			__ASSERT (genes[gene_index].isoform_map.count(isoform_name) == 0, string("error: duplicate isoform ") + isoform_name);
			_isoform_str new_isoform;
			new_isoform.name = isoform_name;
			new_isoform.length = (int)myfasta.sequences[i].length();
			genes[gene_index].isoforms.push_back(new_isoform);
			int isoform_index = (int)genes[gene_index].isoforms.size()-1;
			genes[gene_index].isoform_map[isoform_name] = isoform_index;
		}
	};
	void check_genefile(genefile &mygenefile) {
		for (int i = 0; i < (int)genes.size(); i++) {
			for (int j = 0; j < (int)genes[i].isoforms.size(); j++) {
				string chr;
				gene_struct gene_s;
				mygenefile.search_by_name(genes[i].isoforms[j].name, chr, gene_s);
				if (chr == "") {
					cout << "warning: isoform " << genes[i].isoforms[j].name << " not found in annotation file.\n";
				}
				if (gene_s.geneName != genes[i].name) {
					cout << "warning: gene name for isoform " << genes[i].isoforms[j].name << " mismatch: " << gene_s.geneName << " vs " << genes[i].name << ".\n";
				}
				if (genes[i].isoforms[j].length != round_double(interval_set(gene_s.exons).length())) {
					cout << "warning: length for isoform " << genes[i].isoforms[j].name << " mismatch: " << genes[i].isoforms[j].length << " vs " << round_double(interval_set(gene_s.exons).length()) << ".\n";
				}
			}
		}
	};
	void get_rates(map<int, double> &p_table, int read_length) {
		int total_supporting_reads = 0;
		for (int i = 0; i < (int)genes.size(); i++) {
			total_supporting_reads += genes[i].unique_count;
		}		
		for (int i = 0; i < (int)genes.size(); i++) {
			genes[i].get_rates(total_supporting_reads, p_table, read_length);
		}
	};
	void compute() {
		for (int i = 0; i < (int)genes.size(); i++) {
			genes[i].compute();
		}
	};
	void output_data(const string &file_name) {
		ofstream ofs(file_name);
		for (int i = 0; i < (int)genes.size(); i++) {
			genes[i].output_data(ofs);
		}
		ofs.close();
	};
	void output_exp(const string &file_name) {
		ofstream ofs(file_name);
		for (int i = 0; i < (int)genes.size(); i++) {
			genes[i].output_exp(ofs);
		}
		ofs.close();
	};
};

class iso_exp_model {
public:
	ej_gene *gene;
	int read_length;
	map<int, double> *p_table;
	vector<vector<int> > reads_start, reads_end;
	vector<vector<double> > sampling_rate;
	vector<vector<double> > collapsed_sampling_rate;
	vector<vector<double> > normalized_collapsed_sampling_rate;
	vector<int> collapsed_count;
	vector<vector<double> > possible_sampling_rate;
	vector<vector<double> > collapsed_possible_sampling_rate;
	vector<vector<double> > normalized_collapsed_possible_sampling_rate;
	vector<double> possible_category_counts;

	iso_exp_model(ej_gene &gene1, const int read_length1, map<int, double> &p_table1) : gene(&gene1), read_length(read_length1), p_table(&p_table1) {};
	inline void accumulate_read(const vector<int> &start, const vector<int> &end);
	inline void gen_possible_rates();
	inline void get_rates();
	inline void compute(const double mapped_reads);
	inline void collapse_rates();
	inline void collapse_possible_rates();
	inline void compute_possible_category_counts();
};

inline void iso_exp_model::accumulate_read(const vector<int> &start, const vector<int> &end) {
	__ASSERT(start.size() == end.size() && start.size() == gene->trans.size(), "internal error: illegal arguments for isoform_exp_model::accumulate_read");
	reads_start.push_back(start);
	reads_end.push_back(end);
}

inline void iso_exp_model::gen_possible_rates() {
	int min_insert_length = 1000000, max_insert_length = -1;
	for (map<int, double>::iterator mi = p_table->begin(); mi != p_table->end(); mi++) {
		min_insert_length = min(min_insert_length, mi->first);
		max_insert_length = max(max_insert_length, mi->first);
	}
	__ASSERT(0 < min_insert_length && min_insert_length < max_insert_length, "error insert length.\n");
	int m = (int)gene->trans.size(), n = (int)gene->subexon_sets.size();
	
	//if (m == 1 || n == 1) return;
	//if (sampling_rate.empty()) return;

	vector<vector<bool> > A;
	A.resize(m);
	for (int i = 0; i < m; i++) A[i] = gene->trans[i].subexon_set_vector;
	vector<int> L = gene->len_subexon_sets;
	int r = read_length;
	int sumL = sum(L);
	vector<int> L_A;
	L_A.resize(m);
	for (int i = 0; i < m; i++) {
		L_A[i] = 0;
		for (int j = 0; j < n; j++) {
			if (A[i][j]) L_A[i] += L[j];
		}
	}

//---------------------------------------------------------
// possible_reads[j][i] = position of read j in isoform i,
// -1 if read j is not in isoform i
//---------------------------------------------------------
	vector<vector<int> > possible_reads;
	possible_reads.resize(sumL * m);
	for (int i = 0; i < sumL * m; i++) {
		possible_reads[i].resize(m);
		for (int j = 0; j < m; j++) {
			possible_reads[i][j] = -1;
		}
	}
	vector<int> pos;
	pos.resize(m);
	for (int i = 0; i < m; i++) {
		pos[i] = -1;
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < L[i]; j++) {
			for (int k = 0; k < m; k++) {
				if (A[k][i]) pos[k] ++;
			}
			for (int k = 0; k < m; k++) {
				if (!A[k][i]) {
					continue;
				}
	            
				int rest_exon = i;
				int rest_len = L[i] - j;
				while (rest_exon < n - 1 && rest_len < r) {
					rest_exon = rest_exon + 1;
					if (A[k][rest_exon]) {
						rest_len = rest_len + L[rest_exon];
					}
				}
				if (rest_len < r) {
					continue;
				}
	            
				vector<int> read_vector;
				read_vector.resize(m);
				for (int l = 0; l < m; l++) {
					read_vector[l] = -1;
				}
				for (int l = 0; l < m; l++) {
					bool has_read = true;
					for (int ll = i; ll <= rest_exon; ll++) {
						if (A[l][ll] != A[k][ll]) {
							has_read = false;
							break;
						}
					}
					if (has_read) {						
						read_vector[l] = pos[l];
					}
				}
				possible_reads[k * sumL + pos[k]] = read_vector;
			}
		}
	}

	vector<pair<vector<int>, vector<double> > > possible_paired_end_reads;
	possible_paired_end_reads.clear();
	for (int ii = 0; ii < m; ii++) {
		for (int pos = 0; pos <= L_A[ii] - min_insert_length; pos++) {
			int i = ii * sumL + pos;
			for (int len = min_insert_length; len <= max_insert_length; len++) {
                int j = i + len - r;
				if (j <= i || j > ii * sumL + L_A[ii] - r + 1) {
                    continue;
				}
                vector<int> first_end = possible_reads[i];
                vector<int> second_end = possible_reads[j];
				bool has_read = false;
				for (int k = 0; k < m; k++) {
					if (first_end[k] == -1 || second_end[k] == -1) {
						first_end[k] = second_end[k] = -1;
					} else {
						has_read = true;
					}
				}
				if (!has_read) {
                    continue;
				}

				bool valid_read = true;
				for (int k = 0; k < m; k++) {
					if (first_end[k] > second_end[k]) {
						valid_read = false;
						break;
					}
				}
				if (valid_read) {
                    vector<double> p;
					p.resize(m);
					for (int k = 0; k < m; k++) {
						p[k] = 0;
					}
					bool has_read = false;
					for (int k = 0; k < m; k++) {
						if (first_end[k] > -1 && second_end[k] >= first_end[k] && second_end[k] + r - first_end[k] >= min_insert_length && second_end[k] + r - first_end[k] <= max_insert_length && p_table->count(second_end[k] + r - first_end[k]) > 0) {
                            p[k] = (*p_table)[second_end[k] + r - first_end[k]];
							has_read = true;
						}
					}
					if (!has_read) {
                        continue;
					}
					vector<int> position;
					position.resize(2*m);
					for (int k = 0; k < m; k++) {
						position[k] = first_end[k];
						position[k+m] = second_end[k];
					}
					possible_paired_end_reads.push_back(pair<vector<int>, vector<double> >(position, p));
				}
			}
		}
	}

	sort(possible_paired_end_reads.begin(), possible_paired_end_reads.end());
	possible_sampling_rate.clear();
	for (int i = 0; i < (int)possible_paired_end_reads.size(); i++) {
		if (i == 0 || possible_paired_end_reads[i - 1].first != possible_paired_end_reads[i].first) {
			possible_sampling_rate.push_back(possible_paired_end_reads[i].second);
		}
	}

	collapse_possible_rates();
}

inline void iso_exp_model::compute(const double mapped_reads) {

	int sumL = sum(gene->len_subexon_sets);
	double effective_length = 0;
	for (map<int, double>::iterator mi = p_table->begin(); mi != p_table->end(); mi++) {
		if (sumL - mi->first + 1 <= 0) {
			break;
		}
		effective_length += mi->second * (sumL - mi->first + 1);
	}
	if (effective_length <= epsilon) {
		if (sampling_rate.size() == 0) {
			cout << "warning: sampling_rate.size() != 0 for gene: " << gene->name << endl;
		}
		effective_length = 1;
	}
	gene->gene_expression = (double)sampling_rate.size() / effective_length * 1000 / mapped_reads * 1000000;

	if (gene->iem->sampling_rate.size() <= 20 || gene->subexon_sets.size() == 1) {
		return;
	}

	__ASSERT(gene->trans.size() > 1, string("internal error: more than one subexon_sets for single isoform gene: ") + gene->name);

	gen_possible_rates();
	compute_possible_category_counts();
	gene->tran_expression.resize(gene->trans.size());

	vector<vector<double> > A;
	A.resize(gene->trans.size());
	for (int i = 0; i < (int)A.size(); i++) {
		A[i].resize(collapsed_possible_sampling_rate.size());
		for (int j = 0; j < (int)A[i].size(); j++) {
			A[i][j] = collapsed_possible_sampling_rate[j][i] / 1000 * (double)mapped_reads / 1000000;
		}
	}
	solve_likelihood(A, possible_category_counts, gene->tran_expression);
//	gene->gene_expression = sum(gene->tran_expression);
}

inline void iso_exp_model::get_rates() {
	sampling_rate.clear();
	for (int i = 0; i < (int)reads_start.size(); i++) {
		vector<double> rate;
		rate.resize(gene->trans.size());
		bool positive_rate = false;
		for (int j = 0; j < (int)gene->trans.size(); j++) {
			rate[j] = 0;
			if (reads_start[i][j] >= 0 && reads_end[i][j] >= 0) {
				int insert_len = reads_end[i][j] - reads_start[i][j] + read_length;
				if (p_table->count(insert_len) > 0)
					rate[j] = (*p_table)[insert_len]; //or p_table->operator[](insert_len)
				else 
					rate[j] = 0;
			}
			if (rate[j] > 0) positive_rate = true;
		}
		if (positive_rate) sampling_rate.push_back(rate);
	}

	collapse_rates();
}

inline void iso_exp_model::collapse_rates() {
	collapsed_sampling_rate.clear();
	normalized_collapsed_sampling_rate.clear();
	collapsed_count.clear();
	for (int i = 0; i < (int)sampling_rate.size(); i++) {
		vector<double> rate = sampling_rate[i];
		L1_normalize(rate);
		int j = 0;
		for (j = 0; j < (int)normalized_collapsed_sampling_rate.size(); j++) {
			if (nearly_equal(rate, normalized_collapsed_sampling_rate[j])) {
				for (int k = 0; k < (int)sampling_rate[i].size(); k++) {
					collapsed_sampling_rate[j][k] += sampling_rate[i][k];
				}
				collapsed_count[j] += 1;
				break;
			}
		}
		if (j == (int)normalized_collapsed_sampling_rate.size()) {
			normalized_collapsed_sampling_rate.push_back(rate);
			collapsed_sampling_rate.push_back(sampling_rate[i]);
			collapsed_count.push_back(1);
		}
	}
}

inline void iso_exp_model::collapse_possible_rates() {
	collapsed_possible_sampling_rate.clear();
	normalized_collapsed_possible_sampling_rate.clear();
	for (int i = 0; i < (int)possible_sampling_rate.size(); i++) {
		vector<double> rate = possible_sampling_rate[i];
		L1_normalize(rate);
		int j = 0;
		for (j = 0; j < (int)normalized_collapsed_possible_sampling_rate.size(); j++) {
			if (nearly_equal(rate, normalized_collapsed_possible_sampling_rate[j])) {
				for (int k = 0; k < (int)possible_sampling_rate[i].size(); k++) {
					collapsed_possible_sampling_rate[j][k] += possible_sampling_rate[i][k];
				}
				break;
			}
		}
		if (j == (int)normalized_collapsed_possible_sampling_rate.size()) {
			normalized_collapsed_possible_sampling_rate.push_back(rate);
			collapsed_possible_sampling_rate.push_back(possible_sampling_rate[i]);
		}
	}
}

inline void iso_exp_model::compute_possible_category_counts() {
	possible_category_counts.resize(collapsed_possible_sampling_rate.size());
	for (int i = 0; i < (int)possible_category_counts.size(); i++) {
		possible_category_counts[i] = 0;
	}
	for (int i = 0; i < (int)normalized_collapsed_sampling_rate.size(); i++) {
		int j = 0;
		for (j = 0; j < (int)normalized_collapsed_possible_sampling_rate.size(); j++) {
			if (nearly_equal(normalized_collapsed_sampling_rate[i], normalized_collapsed_possible_sampling_rate[j])) {
				possible_category_counts[j] += collapsed_count[i];
				break;
			}
		}
//		__ASSERT(j < (int)normalized_collapsed_possible_sampling_rate.size(), "internal error: category not found.\n");
		if (j == (int)normalized_collapsed_possible_sampling_rate.size()) {
			printf("warning: category not found, probably mapping error (mapped to N), gene: %s, read count: %d\n", gene->name.c_str(), collapsed_count[i]);
		}
	}
}

#endif // ISOFORM_OPT_MODEL_H
