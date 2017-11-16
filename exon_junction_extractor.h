/*
exon_junction_extractor.h - Header file for exon and junction sequence extraction
written by JIANG Hui, 
Institute for Computational and Mathematical Engineering, Stanford University
May, 2008 -
*/

#ifndef EXON_JUNCTION_EXTRACTOR_H
///Define this macro to prevent from including this header file more than once.
#define EXON_JUNCTION_EXTRACTOR_H

#include "genefile.h"
#include "math_utils.h"
#include "isoform_opt.h"

class ej_trans {
public:
	string name;
	string chr;
	bool strand;
	vector<pair<int, int> > exons;
	vector<pair<int, int> > exon_coords;
	vector<bool> subexon_vector;
	vector<bool> subexon_set_vector;
};

class iso_exp_model;

class ej_gene {
public:
	iso_exp_model *iem; //remember to release this pointers

	string name;
	vector<ej_trans> trans;

	vector<pair<int, int> > subexons;
	vector<bool> subexons_has_break;

	vector<pair<int, int> > subexon_sets;
	vector<int> len_subexon_sets;
	vector<int> len_subexon_set_junc;

	pair<int, int> pos_seq;
	vector<pair<int, int> > pos_subexon_sets;
	vector<pair<int, int> > pos_junc_first;
	vector<vector<pair<int, int> > > pos_junctions;

	int unique_count;
	int total_count;
	vector<int> subexon_set_unique_count;
	vector<vector<int> > subexon_set_junc_unique_count;

	double gene_expression;
	vector<double> tran_expression;
	vector<vector<double> > covar;
	vector<vector<double> > obs_covar;
	double gene_var;

//for differential testing and Baysian inference
	vector<vector<double> > rates;
	vector<double> counts;
};

class exon_junction_extractor {
public:
	int read_len;
	int hanging_len;
	vector<ej_gene> genes;
	map<string, int> genes_map;

	inline exon_junction_extractor(){};

	inline bool check_data() {
		printf("checking data.\n");

		//check duplicate transcripts
		set<string> trans_names;
		for (int i = 0; i < (int)genes.size(); i++) {
			for (int j = 0; j < (int)genes[i].trans.size(); j++) {
				if (trans_names.count(genes[i].trans[j].name) == 0) {
					trans_names.insert(genes[i].trans[j].name);
				} else {
					printf("duplicate transcript: %s, removed.\n", genes[i].trans[j].name.c_str());
					genes[i].trans.erase(genes[i].trans.begin() + j);
					j--;
				}
			}
			if (genes[i].trans.size() == 0) {
				genes.erase(genes.begin() + i);
				i--;
			}
		}

		//print statistics
		int total_trans = 0, max_num_trans = 0, unique_gene = 0;
		for (int i = 0; i < (int)genes.size(); i++) {
			__ASSERT(genes[i].trans.size() >= 1, "internal error: genes[i].trans.size() == 0");
			total_trans += (int)genes[i].trans.size();
			max_num_trans = max(max_num_trans, (int)genes[i].trans.size());
			unique_gene += ((int)genes[i].trans.size() == 1);
		}
		printf("total %d genes, %d transcripts, %d genes have unique transcript, max num transcript a gene is %d.\n", (int)genes.size(), total_trans, unique_gene, max_num_trans);

		//check transcript consistency
		for (int i = 0; i < (int)genes.size(); i++) {
			string chr = "";
			bool strand = true;
			for (int j = 0; j < (int)genes[i].trans.size(); j++) {
				if (j == 0) {
					chr = genes[i].trans[j].chr;
					strand = genes[i].trans[j].strand;
				} else if (chr != genes[i].trans[j].chr) {
					printf("warning: chromosome inconsistency for transcript %s, removed.\n", genes[i].trans[j].name.c_str());
					genes[i].trans.erase(genes[i].trans.begin() + j);
					j--;
				} else if (strand != genes[i].trans[j].strand) {
					printf("warning: strand inconsistency for transcript %s, removed.\n", genes[i].trans[j].name.c_str());
					genes[i].trans.erase(genes[i].trans.begin() + j);
					j--;
				}
			}
			if (genes[i].trans.size() == 0) {
				genes.erase(genes.begin() + i);
				i--;
			}
		}

		//check transcript coordinates
		for (int i = 0; i < (int)genes.size(); i++) {
			for (int j = 0; j < (int)genes[i].trans.size(); j++) {
				bool good_exon = true;
				for (int k = 0; k < (int)genes[i].trans[j].exons.size(); k++) {
					if (!(genes[i].trans[j].exons[k].first >= 0 && genes[i].trans[j].exons[k].second >= genes[i].trans[j].exons[k].first && (k == 0 || genes[i].trans[j].exons[k].first >= genes[i].trans[j].exons[k-1].second))) {
						good_exon = false;
						break;
					}
					if (genes[i].trans[j].exons[k].second == genes[i].trans[j].exons[k].first) {
						printf("warning: exon length is zero for transcript %s, exon removed.\n", genes[i].trans[j].name.c_str());
						genes[i].trans[j].exons.erase(genes[i].trans[j].exons.begin() + k);
						k--;
						continue;
					}
					if (k > 0 && genes[i].trans[j].exons[k].first == genes[i].trans[j].exons[k-1].second) {
						printf("warning: intron length is zero for transcript %s, exon merged.\n", genes[i].trans[j].name.c_str());
						genes[i].trans[j].exons[k-1].second = genes[i].trans[j].exons[k].second; 
						genes[i].trans[j].exons.erase(genes[i].trans[j].exons.begin() + k);
						k--;
						continue;
					}
				}
				if (!good_exon) {
					printf("warning: transcript coordinates inconsistency for transcript %s, removed.\n", genes[i].trans[j].name.c_str());
					genes[i].trans.erase(genes[i].trans.begin() + j);
					j--;
					continue;
				}
				if (genes[i].trans[j].exons.size() == 0) {
					genes[i].trans.erase(genes[i].trans.begin() + j);
					j--;
					continue;
				}
			}
			if (genes[i].trans.size() == 0) {
				genes.erase(genes.begin() + i);
				i--;
				continue;
			}
		}

		//check transcript overlap
		int num_overlap_gene = 0;
		for (int i = 0; i < (int)genes.size(); i++) {
			interval_set intersect_is;
			interval_set union_is;
			for (int j = 0; j < (int)genes[i].trans.size(); j++) {
				interval_set temp_is;
				for (int k = 0; k < (int)genes[i].trans[j].exons.size(); k++) {
					temp_is.union_with(interval_set(genes[i].trans[j].exons[k].first, genes[i].trans[j].exons[k].second));
				}
				union_is.union_with(temp_is);
				if (j == 0) {
					intersect_is.union_with(temp_is);
				} else {
					intersect_is.intersect_with(temp_is);
				}
			}
			__ASSERT(intersect_is.check_valid() && union_is.check_valid(), string("internal error: !interval_set::check_valid(), gene ") + genes[i].name + "\n");
			__ASSERT(union_is.length() > 0, string("interval error: union_is.length() == 0, gene ") + genes[i].name + "\n");
			if (intersect_is.length() / union_is.length() < 0.1) {
				num_overlap_gene++;
				if (num_overlap_gene < 10) {
					printf("warning: transcripts of gene %s overlap too small.\n", genes[i].name.c_str());
				}
			}
		}
		printf ("%d genes have too small overlap.\n", num_overlap_gene);

		return true;
	};

	inline bool generate_result() {
		printf("generating results.\n");

		printf("generating subexons.\n");
		for (int i = 0; i < (int)genes.size(); i++) {
			interval_set is;
			for (int j = 0; j < (int)genes[i].trans.size(); j++) {
				interval_set temp_is(genes[i].trans[j].exons);				
				__ASSERT(temp_is.check_valid(), string("internal error, invalid exon positions, gene ") + genes[i].name + "\n");
				is.break_union_with(temp_is);
			}
			__ASSERT(is.check_int(), string("internal error: !interval_set::check_int(), gene ") + genes[i].name + "\n");
			if (!is.convert_to_int_pairs(genes[i].subexons))
				panic(string("internal error: failed converting to int pair, gene ") + genes[i].name + "\n");
		}

		printf("checking subexons.\n");
		for (int i = 0; i < (int)genes.size(); i++) {
			for (int j = 0; j < (int)genes[i].subexons.size(); j++) {
				__ASSERT(genes[i].subexons[j].first >= 0 && genes[i].subexons[j].second > genes[i].subexons[j].first && (j == 0 || genes[i].subexons[j].first >= genes[i].subexons[j-1].second), string("internal error, wrong subexons coordinates for gene ") + genes[i].name.c_str() + "\n");
			}
			interval_set union_is;
			for (int j = 0; j < (int)genes[i].trans.size(); j++) {
				union_is.union_with(interval_set(genes[i].trans[j].exons));
				interval_set temp_is(genes[i].trans[j].exons);
				temp_is.intersect_with(interval_set(genes[i].subexons));
				__ASSERT(interval_set(genes[i].trans[j].exons).length() == temp_is.length(), string("internal error, wrong subexons size for transcripts ") + genes[i].trans[j].name.c_str() + "\n");
			}
			__ASSERT(union_is.length() == interval_set(genes[i].subexons).length(), string("internal error, wrong subexons size for gene ") + genes[i].name.c_str() + "\n");
		}

		printf("generating subexon vectors.\n");
		for (int i = 0; i < (int)genes.size(); i++) {
			for (int j = 0; j < (int)genes[i].trans.size(); j++) {
				genes[i].trans[j].subexon_vector.clear();
				for (int k = 0; k < (int)genes[i].subexons.size(); k++) {
					interval_set temp_is(genes[i].trans[j].exons);
					temp_is.intersect_with(interval_set(genes[i].subexons[k].first, genes[i].subexons[k].second));
					if (temp_is.length() == 0) {
						genes[i].trans[j].subexon_vector.push_back(false);
					} else if (temp_is.length() == genes[i].subexons[k].second - genes[i].subexons[k].first) {
						genes[i].trans[j].subexon_vector.push_back(true);
					} else {
						panic(string("internal error: failed generating subexon vector, transcript ") + genes[i].trans[j].name + "\n");
					}
				}
			}
		}

		printf("checking subexon vectors.\n");
		for (int i = 0; i < (int)genes.size(); i++) {
			for (int j = 0; j < (int)genes[i].trans.size(); j++) {
				__ASSERT(genes[i].trans[j].subexon_vector.size() == genes[i].subexons.size(), string("internal error, wrong subexon_vector size, transcript ") + genes[i].trans[j].name + "\n");
			}
			for (int k = 0; k < (int)genes[i].subexons.size(); k++) {
				bool has_one = false;
				for (int j = 0; j < (int)genes[i].trans.size(); j++) {
					if (genes[i].trans[j].subexon_vector[k]) {
						has_one = true;
						break;
					}
				}
				__ASSERT(has_one, string("internal error, wrong !has_one, gene ") + genes[i].name + "\n");
				if (k > 0 && genes[i].subexons[k].first == genes[i].subexons[k-1].second) {
					bool has_break = false;
					for (int j = 0; j < (int)genes[i].trans.size(); j++) {
						if (genes[i].trans[j].subexon_vector[k] != genes[i].trans[j].subexon_vector[k-1]) {
							has_break = true;
							break;
						}
					}					
					__ASSERT(has_break, string("internal error, wrong !has_break, gene ") + genes[i].name + "\n");
				}
			}
		}
		return true;
	}

	inline bool load_refFlat_file(const string file_name) {
		printf("loading refFlat file: %s\n", file_name.c_str());
		genefile gf;
		gf.read_from_text_file(file_name);
		int num_dup_trans = 0;
		for (int i = 0; i < (int)gf.num_seq; i++) {
			for (int j = 0; j < (int)gf.sequences[i].num_genes; j++) {
				int k = 0;
				for (k = 0; k < (int) genes.size(); k++) {
					if (genes[k].name == gf.sequences[i].genes[j].geneName) break;
				}
				if (k == (int)genes.size()) {
					ej_gene temp_gene;
					genes.push_back(temp_gene);
					genes[k].name = gf.sequences[i].genes[j].geneName;
				}
				int l = 0;
				for (l = 0; l < (int)genes[k].trans.size(); l++) {
					if (genes[k].trans[l].name == gf.sequences[i].genes[j].name) break;
				}
				if (l == (int)genes[k].trans.size()) {
					ej_trans temp_trans;
					genes[k].trans.push_back(temp_trans);
					genes[k].trans[l].chr = gf.sequences[i].genes[j].chrom;
					genes[k].trans[l].name = gf.sequences[i].genes[j].name;
					genes[k].trans[l].strand = gf.sequences[i].genes[j].strand;
					genes[k].trans[l].exons = gf.sequences[i].genes[j].exons;
				} else {
					num_dup_trans++;
					if (num_dup_trans < 10) {
						printf("warning: duplicate transcript %s, skipped\n", gf.sequences[i].genes[j].name.c_str());
					}
				}
			}
		}
		printf("%d duplicate transcripts.\n", num_dup_trans);
		return check_data();
	};

	inline bool write_to_text_file(const string subexon_file_name) {
		if (!generate_result()) return false;
		printf("writing subexon file: %s\n", subexon_file_name.c_str());
		ofstream subexon_ofs(subexon_file_name.c_str());
		for (int i = 0; i < (int)genes.size(); i++) {
			subexon_ofs << genes[i].name << "\t" <<  genes[i].trans[0].chr << "\t" << (genes[i].trans[0].strand?"+":"-") << "\t" << (int)genes[i].trans.size() << "\t" << (int)genes[i].subexons.size() << "\t";
			for (int j = 0; j < (int)genes[i].subexons.size(); j++) {
				subexon_ofs << genes[i].subexons[j].first << ",";
			}
			subexon_ofs << "\t";
			for (int j = 0; j < (int)genes[i].subexons.size(); j++) {
				subexon_ofs << genes[i].subexons[j].second << ",";
			}
			subexon_ofs << endl;
			for (int j = 0; j < (int)genes[i].trans.size(); j++) {
				subexon_ofs << genes[i].trans[j].name << "\t" << genes[i].trans[j].chr << "\t" << (genes[i].trans[j].strand?"+":"-") << "\t" << j+1 << "\t" << (int)genes[i].trans[j].exons.size() << "\t";
				for (int k = 0; k < (int)genes[i].trans[j].subexon_vector.size(); k++) {
					subexon_ofs << (genes[i].trans[j].subexon_vector[k]?"1":"0") << ",";
				}
				subexon_ofs << endl;
			}
		}
		subexon_ofs.close();
		return true;
	};

	inline bool load_from_file(const string subexon_file_name) {
		genes.clear();
		printf("loading subexon file: %s\n", subexon_file_name.c_str());
		ifstream ifs(subexon_file_name.c_str());
		string readline;
		int num_trans = 0;
		while(getline(ifs, readline)) {
			vector<string> tokens = string_tokenize(readline, " \t");
			if (tokens.size() < 7) {
				printf("error reading subexon file");
				return false;
			}
			ej_gene gene;
			gene.name = tokens[0];
			ej_trans trans;
			trans.chr = tokens[1];
			if (tokens[2] != "+" && tokens[2] != "-") {
				printf("error reading subexon file");
				return false;
			}
			trans.strand = (tokens[2] == "+");
			if (!is_int(tokens[3]) || str2int(tokens[3]) <= 0) {
				printf("error reading subexon file");
				return false;
			}
			gene.trans.resize(str2int(tokens[3]));
			if (!is_int(tokens[4]) || str2int(tokens[4]) <= 0) {
				printf("error reading subexon file");
				return false;
			}
			gene.subexons.resize(str2int(tokens[4]));
			vector<int> firsts = str2int_vec(tokens[5]), seconds = str2int_vec(tokens[6]);
			if (firsts.size() < gene.subexons.size() || seconds.size() < gene.subexons.size()) {
				printf("error reading subexon file");
				return false;
			}
			for (int i = 0; i < (int)gene.subexons.size(); i++) {
				gene.subexons[i].first = firsts[i];
				gene.subexons[i].second = seconds[i];
			}
			for (int i = 0; i < (int)gene.trans.size(); i++) {
				getline(ifs, readline);
				tokens = string_tokenize(readline, " \t");
				gene.trans[i] = trans;
				gene.trans[i].name = tokens[0];
				if (trans.chr != tokens[1] || (trans.strand?"+":"-") != tokens[2] || str2int(tokens[3]) != i + 1) {
					printf("error reading subexon file");
					return false;
				}
				vector<int> subexon_vector = str2int_vec(tokens[5]);
				if (subexon_vector.size() < gene.subexons.size()) {
					printf("error reading subexon file");
					return false;
				}
				gene.trans[i].subexon_vector.resize(gene.subexons.size());
				for (int j = 0; j < (int)gene.subexons.size(); j++) {
					gene.trans[i].subexon_vector[j] = (subexon_vector[j] == 1);
				}
//regenerate trans.exons
				vector<pair<int, int> > trans_subexons;
				for (int j = 0; j < (int)gene.subexons.size(); j++) {
					if (gene.trans[i].subexon_vector[j]) trans_subexons.push_back(gene.subexons[j]);
				}
				interval_set is(trans_subexons);
				is.remove_empty_intervals();
				is.remove_redundant_inner_points();
				gene.trans[i].exons.clear();
				if (!is.convert_to_int_pairs(gene.trans[i].exons)) {
					panic("!is.convert_to_int_pairs(gene.trans[i].exons)");
				}
			}
			genes.push_back(gene);
			num_trans += (int)gene.trans.size();
		}
		ifs.close();
		printf("%d genes loaded, total %d transcripts.\n", (int)genes.size(), num_trans);
		return true;
	}

	inline void prepare_map() {
		printf("prepare map...");
		genes_map.clear();
		for (int i = 0; i < (int) genes.size(); i++) {
			genes_map[genes[i].name] = i;
		}
		printf("done.\n");
	}

	inline void prepare_acceptor(const int index, const int hanging_length = 20, const int separate_N = 5, const bool include_subexons = true, const bool include_junctions = true, const bool combine_exons = true) { 
		ej_gene &gene = genes[index];
		int n = (int)gene.subexons.size();
		gene.subexons_has_break.resize(n);

		for (int j = 0; j < n; j++) {
			gene.subexons_has_break[j] = false;
			if (j < n - 1) {
				for (int k = 0; k < (int)gene.trans.size(); k++) {
					if (gene.trans[k].subexon_vector[j] != gene.trans[k].subexon_vector[j+1]) {
						gene.subexons_has_break[j] = true;
						break;
					}
				}
			}
		}

		gene.subexon_sets.clear();
		gene.len_subexon_sets.clear();
		if (combine_exons) {
			int start = 0, len = 0;
			for (int i = 0; i < n; i++) {
				len += gene.subexons[i].second - gene.subexons[i].first;
				if (gene.subexons_has_break[i] || i == n - 1) {
					gene.subexon_sets.push_back(pair<int, int>(start, i));
					gene.len_subexon_sets.push_back(len);
					len = 0;
					start = i + 1;
				}
			}
		} else {
			for (int i = 0; i < n; i++) {
				gene.subexon_sets.push_back(pair<int, int>(i, i));
				gene.len_subexon_sets.push_back(gene.subexons[i].second - gene.subexons[i].first);
			}
		}

		int m = (int)gene.subexon_sets.size();

		for (int j = 0; j < (int)gene.trans.size(); j++) {
			gene.trans[j].subexon_set_vector.resize(m);
			for (int i = 0; i < m; i++) {
				gene.trans[j].subexon_set_vector[i] = gene.trans[j].subexon_vector[gene.subexon_sets[i].first];
			}
		}

//compute trans.exon_coords
		for (int j = 0; j < (int)gene.trans.size(); j++) {
			gene.trans[j].exon_coords.clear();
			int coord = 0;
			for (int i = 0; i < (int)gene.subexons.size(); i++) {
				int len = gene.subexons[i].second - gene.subexons[i].first;
				if (gene.trans[j].subexon_vector[i]) {
					gene.trans[j].exon_coords.push_back(pair<int, int>(coord, coord + len));
				}
				coord += len;
			}
			interval_set is(gene.trans[j].exon_coords);
			is.remove_empty_intervals();
			is.remove_redundant_inner_points();
			gene.trans[j].exon_coords.clear();
			if (!is.convert_to_int_pairs(gene.trans[j].exon_coords)) {
				panic("internal error: !is.convert_to_int_pairs(gene.trans[j].exon_coords)");
			}
			__ASSERT((int)interval_set(gene.trans[j].exon_coords).length() == (int)interval_set(gene.trans[j].exons).length(), "internal error: wrong exon_coords length");
		}

		gene.subexon_set_unique_count.resize(m);
		gene.subexon_set_junc_unique_count.resize(m);
		for (int j = 0; j < m; j++) {
			gene.subexon_set_junc_unique_count[j].resize(m);
			gene.subexon_set_unique_count[j] = 0;
			for (int k = 0; k < m; k++) {
				gene.subexon_set_junc_unique_count[j][k] = 0;
			}
		}

		gene.len_subexon_set_junc.resize(m);
		gene.pos_subexon_sets.resize(m);
		gene.pos_junc_first.resize(m);
		gene.pos_junctions.resize(m);
		for (int i = 0; i < m; i++) {
			gene.pos_junctions[i].resize(m);
		}

		int pos = 0;
		gene.pos_seq.first = pos;
		for (int j = 0; j < m; j++) {
			if (j > 0 && include_subexons) pos += separate_N;
			gene.pos_subexon_sets[j].first = pos;
			if (include_subexons) pos += gene.len_subexon_sets[j];
			gene.pos_subexon_sets[j].second = pos;
		}
		gene.pos_seq.second = pos;

		for (int j = 0; j < m; j++) {
			gene.len_subexon_set_junc[j] = min(hanging_length, gene.len_subexon_sets[j]);
		}

		for (int j = 0; j < m - 1; j++) {
			gene.pos_junc_first[j].first = pos;
			for (int k = j + 1; k < m; k++) {
				if (include_junctions && (include_subexons || j > 0 || k > 1)) pos += separate_N;
				gene.pos_junctions[j][k].first = pos;
				if (include_junctions) pos += gene.len_subexon_set_junc[j] + gene.len_subexon_set_junc[k];
				gene.pos_junctions[j][k].second = pos;
			}
			gene.pos_junc_first[j].second = pos;
		}
	}

	inline void get_sequence(string &sequence, const int index, const fasta &fa, const int hanging_length = 20, const int separate_N = 5, const bool include_subexons = true, const bool include_junctions = true, const bool combine_exons = true) {
		prepare_acceptor(index, hanging_length, separate_N, include_subexons, include_junctions, combine_exons);
		sequence = "";
		string seperator = "";
		seperator.insert(seperator.begin(), separate_N, 'N');
		for (int i = 0; i < (int)fa.tags.size(); i++) {
			if (genes[index].trans[0].chr == fa.tags[i]) {
				vector<string> seq_subexon_sets;
				seq_subexon_sets.resize(genes[index].subexon_sets.size());
				for (int j = 0; j < (int)genes[index].subexon_sets.size(); j++) {
					if (j > 0 && include_subexons) sequence += seperator;
					seq_subexon_sets[j] = "";
					for (int k = genes[index].subexon_sets[j].first; k <= genes[index].subexon_sets[j].second; k++) {
						seq_subexon_sets[j] += fa.sequences[i].substr(genes[index].subexons[k].first, genes[index].subexons[k].second - genes[index].subexons[k].first);
					}
					if (include_subexons) sequence += seq_subexon_sets[j];
				}

				if (include_junctions) {
					vector<string> lefts, rights;
					for (int j = 0; j < (int)genes[index].subexon_sets.size(); j++) {
						lefts.push_back(seq_subexon_sets[j].substr(genes[index].len_subexon_sets[j] - genes[index].len_subexon_set_junc[j], genes[index].len_subexon_set_junc[j]));
						rights.push_back(seq_subexon_sets[j].substr(0, genes[index].len_subexon_set_junc[j]));
					}

					for (int j = 0; j < (int)genes[index].subexon_sets.size() - 1; j++) {
						for (int k = j + 1; k < (int)genes[index].subexon_sets.size(); k++) {
							if (include_subexons || j > 0 || k > 1) sequence += seperator;
							sequence += lefts[j];
							sequence += rights[k];
						}
					}
				}

				break;
			}
		}
	}

	inline void get_ready(const string subexon_file_name, const int read_length = 25, const int hanging_length = 20, const int separate_N = 5, const bool include_subexons = true, const bool include_junctions = true, const bool combine_exons = true) { 
		read_len = read_length;
		hanging_len = hanging_length;
		load_from_file(subexon_file_name);
		printf("preparing acceptors...");
		for (int i = 0; i < (int)genes.size(); i++) {
			prepare_acceptor(i, hanging_length, separate_N, include_subexons, include_junctions, combine_exons);
		}
		printf("done.\n");
		prepare_map();
	}

	inline bool locate_coord(const string &gene_name, const int coord, ej_gene* &p_gene, int &subexon_set_id1, int &subexon_set_id2, int &gene_pos, int &chr_pos, vector<int> &tran_pos) {
//		__ASSERT(genes_map.count(gene_name) == 1, string("gene not found: ") + gene_name);
		__ASSERT(genes_map.count(gene_name) <= 1, string("multiple genes found: ") + gene_name);

		if (genes_map.count(gene_name) == 0) {
			p_gene = NULL;
			return false;
		}

		ej_gene &gene = genes[genes_map[gene_name]];
		p_gene = &gene;

		int m = (int)gene.subexon_sets.size();
		int i, j, k, l, ll;
		if (coord >= gene.pos_seq.first && coord < gene.pos_seq.second) {
			for (i = 0; i < m; i++) {
				if (coord >= gene.pos_subexon_sets[i].first && coord < gene.pos_subexon_sets[i].second) {
					int pos = coord - gene.pos_subexon_sets[i].first;
					for (j = gene.subexon_sets[i].first; j <= gene.subexon_sets[i].second; j++) {
						if (pos < gene.subexons[j].second - gene.subexons[j].first) {
							subexon_set_id1 = i;
							subexon_set_id2 = -1;
							gene_pos = 0;
							for (k = 0; k < j; k++) gene_pos += gene.subexons[k].second - gene.subexons[k].first;
							gene_pos += pos;
							chr_pos = gene.subexons[j].first + pos;
							tran_pos.clear();
							tran_pos.resize(gene.trans.size());
							for (k = 0; k < (int)gene.trans.size(); k++) {
								if (gene.trans[k].subexon_set_vector[i]) {
									tran_pos[k] = gene_pos;
									for (l = 0; l < gene.subexon_sets[i].first; l++) {
										if (!gene.trans[k].subexon_vector[l]) {
											tran_pos[k] -= gene.subexons[l].second - gene.subexons[l].first;
										}
									}
								} else {
									tran_pos[k] = -1;
								}
							}
							return true;
						} else {
							pos -= gene.subexons[j].second - gene.subexons[j].first;
						}
					}
					//__ASSERT(j <= gene.subexon_sets[i].second, string("internal error, pos not found, gene: ") + gene_name);
					if (j > gene.subexon_sets[i].second) {
						cout << string("warning: pos not found, gene: ") + gene_name << endl;
						return false;
					}
					return false;
				}
			}
			//__ASSERT(i < m, string("internal error, pos not found, gene: ") + gene_name);
			if (i >= m) {
				cout << string("warning: pos not found, gene: ") + gene_name << endl;
				return false;
			}
			return false;
		} else {
			for (i = 0; i < m - 1; i++) {
				if (coord >= gene.pos_junc_first[i].first && coord < gene.pos_junc_first[i].second) {
					for (j = i+1; j < m; j++) {
						if (coord >= gene.pos_junctions[i][j].first && coord < gene.pos_junctions[i][j].second) {
							if (coord < gene.pos_junctions[i][j].first + gene.len_subexon_set_junc[i]) {
								int pos = gene.pos_junctions[i][j].first + gene.len_subexon_set_junc[i] - coord;
								for (k = gene.subexon_sets[i].second; k >= gene.subexon_sets[i].first; k--) {
									if (pos <= gene.subexons[k].second - gene.subexons[k].first) {
										subexon_set_id1 = i;
										subexon_set_id2 = j;
										gene_pos = 0;
										for (l = 0; l <= k; l++) gene_pos += gene.subexons[l].second - gene.subexons[l].first;
										gene_pos -= pos;
										chr_pos = gene.subexons[k].second - pos;
										tran_pos.clear();
										tran_pos.resize(gene.trans.size());
										for (l = 0; l < (int)gene.trans.size(); l++) {
											bool f_has_r = gene.trans[l].subexon_set_vector[i] && gene.trans[l].subexon_set_vector[j];
											for (ll = i + 1; ll < j; ll++) {
												f_has_r &= (!gene.trans[l].subexon_set_vector[ll]);
											}
											if (f_has_r) {
												tran_pos[l] = gene_pos;
												for (ll = 0; ll < gene.subexon_sets[i].first; ll++) {
													if (!gene.trans[l].subexon_vector[ll]) {
														tran_pos[l] -= gene.subexons[ll].second - gene.subexons[ll].first;
													}
												}
											} else {
												tran_pos[l] = -1;
											}
										}
										return true;
									} else {
										pos -= gene.subexons[k].second - gene.subexons[k].first;
									}
								}
								//__ASSERT(k >= gene.subexon_sets[i].first, string("internal error, pos not found, gene: ") + gene_name);
								if (k < gene.subexon_sets[i].first) {
									cout << string("warning: pos not found, gene: ") + gene_name << endl;
									return false;
								}
								return false;
							} else { // coord >= gene.pos_junctions[i][j].second - gene.len_subexon_set_junc[j]) {
								__ASSERT(coord >= gene.pos_junctions[i][j].second - gene.len_subexon_set_junc[j] && gene.pos_junctions[i][j].second - gene.pos_junctions[i][j].first == gene.len_subexon_set_junc[i] + gene.len_subexon_set_junc[j], string("internal error, wrong junction length, gene: ") + gene_name);
								int pos = coord - (gene.pos_junctions[i][j].second - gene.len_subexon_set_junc[j]);
								for (k = gene.subexon_sets[j].first; k <= gene.subexon_sets[j].second; k++) {
									if (pos < gene.subexons[k].second - gene.subexons[k].first) {
										subexon_set_id1 = i;
										subexon_set_id2 = j;
										gene_pos = 0;
										for (l = 0; l < k; l++) gene_pos += gene.subexons[l].second - gene.subexons[l].first;
										gene_pos += pos;
										chr_pos = gene.subexons[k].first + pos;
										tran_pos.clear();
										tran_pos.resize(gene.trans.size());
										for (l = 0; l < (int)gene.trans.size(); l++) {
											bool f_has_r = gene.trans[l].subexon_set_vector[i] && gene.trans[l].subexon_set_vector[j];
											for (ll = i + 1; ll < j; ll++) {
												f_has_r &= (!gene.trans[l].subexon_set_vector[ll]);
											}
											if (f_has_r) {
												tran_pos[l] = gene_pos;
												for (ll = 0; ll < gene.subexon_sets[j].first; ll++) {
													if (!gene.trans[l].subexon_vector[ll]) {
														tran_pos[l] -= gene.subexons[ll].second - gene.subexons[ll].first;
													}
												}
											} else {
												tran_pos[l] = -1;
											}
										}
										return true;
									} else {
										pos -= gene.subexons[k].second - gene.subexons[k].first;
									}
								}
								//__ASSERT(k <= gene.subexon_sets[j].second, string("internal error, pos not found, gene: ") + gene_name);
								if (k > gene.subexon_sets[j].second) {
									cout << string("warning: pos not found, gene: ") + gene_name << endl;
									return false;
								}
								return false;
							}
						}
					}
					//__ASSERT(j < m, string("internal error, pos not found, gene: ") + gene_name);
					if (j >= m) {
						cout << string("warning: pos not found, gene: ") + gene_name << endl;
						return false;
					}
					return false;
				}
			}
			//__ASSERT(i < m - 1, string("internal error, pos not found, gene: ") + gene_name);
			if (i >= m - 1) {
				cout << string("warning: pos not found, gene: ") + gene_name << endl;
				return false;
			}
		}
		return false;
	}

/* waited to replace old ones later

	inline void handler(const string &gene_name, const int coord, const bool unique_read) {
		ej_gene *p_gene;
		int subexon_set_id1, subexon_set_id2, gene_pos, chr_pos;
		if (!locate_coord(gene_name, coord, p_gene, subexon_set_id1, subexon_set_id2, gene_pos, chr_pos)) return;
		p_gene->total_count += 1;
		if (!unique_read) return;
		p_gene->unique_count += 1;
		__ASSERT(subexon_set_id1 >= 0 && subexon_set_id1 < (int)p_gene->subexon_sets.size(), "internal error: bad subexon_set_id1");
		if (subexon_set_id2 == -1) {
			p_gene->subexon_set_unique_count[subexon_set_id1]++;
		} else {
			__ASSERT(subexon_set_id2 >= 0 && subexon_set_id2 < (int)p_gene->subexon_sets.size(), "internal error: bad subexon_set_id2");
			p_gene->subexon_set_junc_unique_count[subexon_set_id1][subexon_set_id2]++;
		}
	}

	inline pair<string, int> map_coord_to_chr(const string &gene_name, const int coord) {
		ej_gene *p_gene;
		int subexon_set_id1, subexon_set_id2, gene_pos, chr_pos;
		pair<string, int> result;
		result.first = "";
		result.second = -1;
		if (!locate_coord(gene_name, coord, p_gene, subexon_set_id1, subexon_set_id2, gene_pos, chr_pos)) return result;
		__ASSERT(subexon_set_id1 >= 0 && subexon_set_id1 < (int)p_gene->subexon_sets.size(), "internal error: bad subexon_set_id1");
		if (subexon_set_id2 == -1) {
			result.first = p_gene->trans[0].chr;
			result.second = chr_pos;
		} else {
			__ASSERT(subexon_set_id2 >= 0 && subexon_set_id2 < (int)p_gene->subexon_sets.size(), "internal error: bad subexon_set_id2");
		}
		return result;
	}
*/

	inline void handler(const string &gene_name, const int coord, const bool unique_read) {
//		__ASSERT(genes_map.count(gene_name) == 1, string("gene not found: ") + gene_name);
		__ASSERT(genes_map.count(gene_name) <= 1, string("multiple genes found: ") + gene_name);
		
		if (genes_map.count(gene_name) == 0) {
			return;
		}

		ej_gene &gene = genes[genes_map[gene_name]];
		gene.total_count += 1;

		if (!unique_read) {
			return;
		}

		gene.unique_count += 1;

		int m = (int)gene.subexon_sets.size();
		int i, j;
		if (coord >= gene.pos_seq.first && coord <= gene.pos_seq.second) {
			for (i = 0; i < m; i++) {
				if (coord >= gene.pos_subexon_sets[i].first && coord <= gene.pos_subexon_sets[i].second) {
					gene.subexon_set_unique_count[i]++;
					break;
				}
			}
			//__ASSERT(i < m, string("internal error, pos not found, gene: ") + gene_name);
			if (i >= m) {
				cout << string("warning: pos not found, gene: ") + gene_name << endl;
				return;
			}
		} else {
			for (i = 0; i < m - 1; i++) {
				if (coord >= gene.pos_junc_first[i].first && coord <= gene.pos_junc_first[i].second) {
					for (j = i+1; j < m; j++) {
						if (coord >= gene.pos_junctions[i][j].first && coord <= gene.pos_junctions[i][j].second) {
							gene.subexon_set_junc_unique_count[i][j]++;
							break;
						}
					}
					//__ASSERT(j < m, string("internal error, pos not found, gene: ") + gene_name);
					if (j >= m) {
						cout << string("warning: pos not found, gene: ") + gene_name << endl;
						return;
					}
					break;
				}
			}
			//__ASSERT(i < m, string("internal error, pos not found, gene: ") + gene_name);
			if (i >= m) {
				cout << string("warning: pos not found, gene: ") + gene_name << endl;
				return;
			}
		}
	}

	inline pair<string, int> map_coord_to_chr(const string &gene_name, const int coord) {
		//__ASSERT(genes_map.count(gene_name) == 1, string("gene not found: ") + gene_name);
		__ASSERT(genes_map.count(gene_name) <= 1, string("multiple genes found: ") + gene_name);

		if (genes_map.count(gene_name) == 0) {
			return pair<string, int>("", -1);
		}

		ej_gene &gene = genes[genes_map[gene_name]];
		int m = (int)gene.subexon_sets.size();
		int i, j;
		pair<string, int> result;
		result.first = "";
		result.second = -1;
		if (coord >= gene.pos_seq.first && coord <= gene.pos_seq.second) {
			for (i = 0; i < m; i++) {
				if (coord >= gene.pos_subexon_sets[i].first && coord <= gene.pos_subexon_sets[i].second) {
					int pos = coord - gene.pos_subexon_sets[i].first;
					for (j = gene.subexon_sets[i].first; j <= gene.subexon_sets[i].second; j++) {
						if (pos < gene.subexons[j].second - gene.subexons[j].first) {
							result.first = gene.trans[0].chr;
							result.second = gene.subexons[j].first + pos;
							break;
						} else {
							pos -= gene.subexons[j].second - gene.subexons[j].first;
						}
					}
					//__ASSERT(j <= gene.subexon_sets[i].second, string("internal error, pos not found, gene: ") + gene_name);
					if (j > gene.subexon_sets[i].second) {
						cout << string("warning: pos not found, gene: ") + gene_name << endl;
						return pair<string, int>("", -1);
					}
					break;
				}
			}
			//__ASSERT(i < m, string("internal error, pos not found, gene: ") + gene_name);
			if (i >= m) {
				cout << string("warning: pos not found, gene: ") + gene_name << endl;
				return pair<string, int>("", -1);
			}
		} else {
			for (i = 0; i < m - 1; i++) {
				if (coord >= gene.pos_junc_first[i].first && coord <= gene.pos_junc_first[i].second) {
					for (j = i+1; j < m; j++) {
						if (coord >= gene.pos_junctions[i][j].first && coord <= gene.pos_junctions[i][j].second) {
							break;
						}
					}
					//__ASSERT(j < m, string("internal error, pos not found, gene: ") + gene_name);
					if (j >= m) {
						cout << string("warning: pos not found, gene: ") + gene_name << endl;
						return pair<string, int>("", -1);
					}
					break;
				}
			}
			//__ASSERT(i < m, string("internal error, pos not found, gene: ") + gene_name);
			if (i >= m) {
				cout << string("warning: pos not found, gene: ") + gene_name << endl;
				return pair<string, int>("", -1);
			}
		}
		return result;
	}
	
	inline void compute_expression(const int index, const int mapped_reads) {
		ej_gene &gene = genes[index];
		gene.tran_expression.resize(gene.trans.size());
		if (gene.subexon_sets.size() == 1 && gene.trans.size() > 1) cout << "warning: duplicated isoforms for gene: " << gene.name << endl;
//		if (gene.subexon_sets.size() > 1 && gene.trans.size() == 1) panic(string("internal error: more than one subexon_sets for single isoform gene: ") + gene.name);
		if (gene.subexon_sets.size() > 1 && gene.trans.size() == 1) cout << "warning: more than one subexon_sets for single isoform gene: " << gene.name << endl;
		if (gene.subexon_sets.size() == 1 && gene.trans.size() == 1) {
			int effective_length = gene.len_subexon_sets[0] - read_len + 1;
			if (effective_length <= 0) {
				__ASSERT(gene.subexon_set_unique_count[0] == 0, string("internal error: gene.subexon_set_unique_count[0] != 0 for gene: " + gene.name));
				effective_length = 1;
			}
			gene.gene_expression = gene.tran_expression[0] = (double)gene.subexon_set_unique_count[0] / mapped_reads / effective_length * 1000 * 1000000;
			double p = (double)gene.subexon_set_unique_count[0] / mapped_reads;
			gene.gene_var = 1e18 * p * (1-p) / mapped_reads / effective_length / effective_length;
		} else {
			int m = (int)gene.trans.size();
			int n = (int)gene.subexon_set_unique_count.size();
			vector<int> N;
			vector<int> L;
			vector<vector<int> > A;
			N.clear();
			L.clear();
			A.clear();
			for (int i = 0; i < n; i++) {
				int effective_length = gene.len_subexon_sets[i] - read_len + 1;
				if (effective_length <= 0) {
					if (gene.subexon_set_unique_count[i] > 0) {
						cout << "warning: gene.subexon_set_unique_count[i] != 0 for gene: " << gene.name << endl;
					}
				} else {
					N.push_back(gene.subexon_set_unique_count[i]);
					L.push_back(effective_length);
					vector<int> temp_A;
					for (int j = 0; j < m; j++) temp_A.push_back((int)gene.trans[j].subexon_set_vector[i]);
					A.push_back(temp_A);
				}
			}
			for (int i = 0; i < n; i++) {
				for (int j = i + 1; j < n; j++) {
					int effective_length = hanging_len * 2 - read_len + 1;
					if (effective_length <= 0) {
						__ASSERT(gene.subexon_set_junc_unique_count[i][j] == 0, string("internal error: gene.subexon_set_junc_unique_count[i][j] != 0 for gene: " + gene.name));
					} else {
						N.push_back(gene.subexon_set_junc_unique_count[i][j]);
						L.push_back(effective_length);
						vector<int> temp_A;
						for (int k = 0; k < m; k++) {
							int temp = (int)gene.trans[k].subexon_set_vector[i] * (int)gene.trans[k].subexon_set_vector[j];
							for (int l = i + 1; l < j; l++) {
								temp *= (1 - (int)gene.trans[k].subexon_set_vector[l]);
							}
							temp_A.push_back(temp);
						}
						A.push_back(temp_A);
					}
				}
			}
//collapse A
			for (int i = (int)A.size() - 1; i >= 0; i--) {
//remove zero
				bool is_zero = true;
				for (int j = 0; j < m; j++) {
					if (A[i][j] != 0) {
						is_zero = false;
						break;
					}
				}
				if (is_zero) {
					A.erase(A.begin() + i);
					N.erase(N.begin() + i);
					L.erase(L.begin() + i);
					continue;
				}
//merge duplicates
				int duplicate = -1;
				for (int j = 0; j < i; j++) {
					if (A[i] == A[j]) {
						duplicate = j;
						break;
					}
				}
				if (duplicate > 0) {
					N[duplicate] += N[i];
					L[duplicate] += L[i];
					A.erase(A.begin() + i);
					N.erase(N.begin() + i);
					L.erase(L.begin() + i);
					continue;
				}
			}

			vector<double> double_N;
			vector<vector<double> > double_A_t;
			double_A_t.resize(m);

			for (int i = 0; i < (int)N.size(); i++) double_N.push_back((double)N[i]);
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < (int)A.size(); j++) {
					double_A_t[i].push_back((double)A[j][i] * (double)L[j] / 1000 * (double)mapped_reads / 1000000);
				}
			}

			gene.rates = transpose(double_A_t);
			gene.counts = double_N;

			//iso_opt.verbose = true;
			//cout << "compute expression for gene: " << gene.name << endl;
			solve_likelihood(double_A_t, double_N, gene.tran_expression);
			gene.gene_expression = sum(gene.tran_expression);
			//gene.covar = iso_opt.inv_fisher(gene.tran_expression);
//			gene.obs_covar = iso_opt.inv_obs_fisher(gene.tran_expression);
			//gene.gene_var = sum_matrix(gene.covar);
		}
	}
};

#endif //#ifdef EXON_JUNCTION_EXTRACTOR_H
