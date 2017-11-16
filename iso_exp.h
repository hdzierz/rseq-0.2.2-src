/*
iso_exp.h - Header file for isoform expression analysis
written by JIANG Hui, 
Institute for Computational and Mathematical Engineering, Stanford University
May, 2008 -
*/

#ifndef ISO_EXP_H
///Define this macro to prevent from including this header file more than once.
#define ISO_EXP_H

#include "stl.h"
#include "type.h"
#include "string_operation.h"
#include "file_operation.h"
#include "bar.h"
#include "genefile.h"
#include "fasta.h"
#include "exon_junction_extractor.h"
#include "map_result.h"
#include "probe.h"
#include "math_utils.h"
#include "iso_exp_model.h"
#include "probe_match.h"
#include "chr_region.h"
#include "table.h"

class task_set {
public:
    set<string> tasks;
    task_set() {tasks.clear();}
    void operator = (string task) {tasks.insert(task);}
    bool operator == (string task) const {return tasks.count(task) > 0;}
    string to_string() const {
        string result = "";
        for (auto ti = tasks.begin(); ti != tasks.end(); ti++) {
            result += (*ti) + " ";
        }
        return result;
    }
    friend ostream& operator<< (std::ostream& stream, const task_set& ts) {return stream << ts.to_string();}
};

class gene_info_struct{
public:
	int length;
	int reduce_len;
	int unique_count;
	double total_count;
	gene_info_struct(){
		length = 0;
		reduce_len = 0;
		unique_count = 0;
		total_count = 0;
	}
};

class paired_target_struct {
public:
	string read_id;
	string gene_name;
	string isoform_name;
	bool reverse_strand;
	int isoform_index;
	int pos1, pos2, read_len1, read_len2;
	int insert_len;
};

class local_data{
public:
	fasta myfasta; //for all tasks
	genefile mygenefile; // for all tasks
	string output_prefix; // can not be specified using command-line options in multiple thread mode
	vector<string> other_args;
	vector<vector<int> > nt_counts; // for scan_reads
	ofstream ofs;
	map<pair<int, pair<int, int> >, int> exon_counts; // for exon_usage
	map<int, int> insert_len_hist; // for paired_end and expression_analysis
	iso_exp_model_2 iem2; // for expression_analysis
	map<int, double> p_table; // for paired_end and expression_analysis
	map<string, int> gene_counts; // for two_sample_diff_gene_exp
	map<string, vector<pair<set<int>, set<int> > > > exon_usage_sets; // for two_sample_diff_exon_usage
	map<string, vector<pair<int, int> > > exon_usage_set_lengths; // for two_sample_diff_exon_usage
	map<string, vector<pair<int, int> > > exon_usage_set_counts; // for two_sample_diff_exon_usage
	int total_read_counts; // for sequential
	local_data() {
		total_read_counts = 0;
	}
};

#define MAX_READ_LEN 1024 // for scan_reads and map_stat
#define MAX_INS_LEN 10000 // for paired_end and expression_analysis

///////////////////////////////////////////////////////
// command-line parameters
int read_length = 25;
int num_mismatch = 3;
string direction = "b";
string reduce_len_file = "";
string annotation_file = "";
string reference_file = "";
string ins_len_p_file = "";
string exon_file = "";
set<string> gene_name;
int version = 1;
bool t_test = false;
bool multiple = false;
bool quick = false;
bool do_stop = true;
int report_progress = 0;
int num_itr = 1000;
bool reset_srand = false;
bool old = false;
int num_max_isoform = 3;
int num_random_chr = 1;
int chr_random_size = 1000000;
int num_random_gene = 1;
int num_random_iso = 2;
int num_random_reads = 10000;
task_set task;
///////////////////////////////////////////////////////

exon_junction_extractor eje; // for comp_exp, paired_end
map<string, gene_info_struct> genes; // for comp_exp
//map<string, gene_info_struct> transcripts; // for comp_exp
ofstream ofs_F, ofs_R, ofs_M; // for extract, convert_coord, paired_end
mutex samples_lock; // for do_sequential_tasks
int active_samples = 0; // for do_sequential_tasks
int waiting_samples = 0; // for do_sequential_tasks
bool start_output = false; // for do_sequential_tasks
vector<local_data*> samples; // for do_sequential_tasks
condition_variable samples_cond_var; // for do_sequential_tasks
map<string, bool> diff_gene_rejected; // for two_sample_diff_gene_exp
map<string, vector<int> > diff_exon_usage_results; // for two_sample_diff_exon_usage, 0: uncertain, 1: diff_usage, 2: no_diff_usage, 3: < 1RPKM

thread_local int sample_id; // for do_sequential_tasks
thread_local local_data* local = NULL; // for do_sequential_tasks
thread_local int min_read_len = MAX_READ_LEN+1, max_read_len = -1, num_total_read = 0, num_good_read = 0; // for scan_reads
thread_local int64 total_nt = 0, total_gc = 0; // for scan_reads
thread_local int total_reads = 0, unmapped_reads = 0, mapped_reads = 0, unique_reads = 0, total_targets = 0; // for handler, single reads
thread_local int half_mapped_reads = 0, diff_chr_reads = 0, same_dir_reads = 0, neg_dist_reads = 0, large_dist_reads = 0, paired_reads = 0, other_reads = 0; // for handler, paired reads
thread_local bool first_phase = true; // for comp_exp [and potentiall expression_analysis]
thread_local int intergene_reads = 0, intron_reads = 0, UTR5_reads = 0, UTR3_reads = 0, exon_reads = 0; // for exon_usage
thread_local int real_read_length = 0; // for map_stat
thread_local int map_stat_counts[MAX_READ_LEN][5][5]; // for map_stat

bool gtf_convert(int output_format) { //1: refFlat, 2: BED
	genefile mygenefile;
	if (!mygenefile.read_from_text_file(local->other_args[0])) return false;
	if (local->output_prefix == "") local->output_prefix = local->other_args[0];
	if (output_format == 1) {
		mygenefile.version = 1.0f;
		mygenefile.write_to_text_file(local->output_prefix + ".refFlat.txt");
	} else if (output_format == 2) {
		mygenefile.version = 2.0f;
		mygenefile.write_to_text_file(local->output_prefix + ".BED");
	}
	return true;
}

void annotate_transcripts() {
	fasta myfasta;
	myfasta.read_from_file(local->other_args[0]);
	myfasta.build_tag_map();
	table mytable;
	mytable.read_from_file(local->other_args[1]);
	mytable.build_content_map();
	int num_multiple_matches = 0;
	for (int i = 0; i < (int)myfasta.tags.size(); i++) {
		vector<string> tokens = string_tokenize(myfasta.tags[i]);
		string id = "";
		if (tokens.size() == 0) {
			printf("warning: sequence tag is empty.\n");
		} else {
			id = tokens[0];
		}
		string gene_name = id;
		if (mytable.content_map.count(id) == 0) {
			num_multiple_matches++;
			if (num_multiple_matches <= 5) printf("warning: transcript %s not found.\n", id.c_str());
		} else {
			vector<pair<int, int> > targets = mytable.content_map[id];
			__ASSERT(targets.size() > 0, "internal error: target is empty");
			set<string> gene_names;
			for (int j = 0; j < (int)targets.size(); j++) {
				int x = targets[j].first;
				__ASSERT(x < (int)mytable.data.size(), "internal error: target is out of table.\n");
				if (mytable.data[x].size() < 5 || mytable.data[x][4] == "") printf("warning: gene name does not exist for transcript %s.\n", id.c_str());
				gene_names.insert(mytable.data[x][4]);
			}
			if (gene_names.size() > 1) {
				printf("warning: multiple matches for transcript %s, first one is taken.\n", id.c_str());
			}
			gene_name = *(gene_names.begin());
		}
		myfasta.tags[i] = gene_name + "$$" + id;
		myfasta.sequences[i] = toupper(myfasta.sequences[i]);
	}
	if (num_multiple_matches > 5) printf("%d warnings suppressed.\n", num_multiple_matches - 5);
	if (local->output_prefix == "") local->output_prefix = local->other_args[0];
	myfasta.write_to_file(local->output_prefix + ".new.fa");
}

void random_genome() {
	if (local->output_prefix == "") local->output_prefix = local->other_args[0];
	ofstream ofs((local->output_prefix + ".chrfilelist.txt").c_str());
	for (int i = 0; i < num_random_chr; i++) {
		fasta myfasta;
		string sequence = "";
		for (int j = 0; j < chr_random_size; j++) {
			sequence += random_NT();
		}
		myfasta.tags.push_back("chr" + int2str(i+1));
		myfasta.sequences.push_back(sequence);
		myfasta.write_to_file(local->output_prefix + ".chr" + int2str(i+1) + ".fa");
		ofs << local->other_args[0] + ".chr" + int2str(i+1) + ".fa" << endl;
	}
	ofs.close();
	printf("%d chromsomes of size %d generated.\n", num_random_chr, chr_random_size);
}

void random_annotation() {
	genefile mygenefile;
	mygenefile.num_seq = num_random_chr;
	mygenefile.sequences.resize(num_random_chr);
	for (int i = 0; i < num_random_chr; i++) {
		mygenefile.sequences[i].name = "chr" + int2str(i+1);
	}
	for (int i = 0; i < num_random_gene; i++) {
		while(true) {
			int chr = rand_int(0, num_random_chr-1);
			bool strand = (rand_int(2) == 1);
			vector<pair<int, int> > exons;
			int temp = round_double(rand_normal(10, 3));
			int num_exon = max(1, temp);
			exons.resize(num_exon);
			int start = rand_int(chr_random_size) - 1;
			for (int j = 0; j < num_exon; j++) {
				int temp = round_double(rand_normal(400, 200));
				exons[j].first = (j==0)? start : exons[j-1].second + max(100, temp);
				temp = round_double(rand_normal(80, 20));
				exons[j].second = exons[j].first + max(30, temp);
			}
			if (exons[num_exon-1].second > chr_random_size) continue;
			for (int j = 0; j < num_random_iso; j++) {
				gene_struct gene;
				gene.geneName = "gene" + int2str(i+1);
				gene.name = gene.geneName + "_isoform" + int2str(j+1);
				gene.chrom = mygenefile.sequences[chr].name;
				gene.strand = strand;
				for (int k = 0; k < num_exon; k++) {
					if (rand_int(2) > 1) gene.exons.push_back(exons[k]);
				}
				if (gene.exons.empty()) gene.exons.push_back(exons[0]);
				gene.txStart = gene.exons[0].first;
				gene.txEnd = gene.exons.back().second;
				gene.cdsStart = gene.txStart;
				gene.cdsEnd = gene.txEnd;
				mygenefile.sequences[chr].genes.push_back(gene);
				mygenefile.sequences[chr].num_genes++;
			}
			break;
		}
	}
	if (local->output_prefix == "") local->output_prefix = local->other_args[0];
	mygenefile.write_to_text_file(local->output_prefix);
	printf("%d genes generated.\n", num_random_gene);
}

void random_expression() {
	genefile mygenefile;
	mygenefile.read_from_text_file(local->other_args[0]);
	if (local->output_prefix == "") local->output_prefix = local->other_args[0];
	ofstream ofs(local->output_prefix + ".exp");
	int num_trans = 0;
	for (uint i = 0; i < mygenefile.num_seq; i++) {
		for (uint j = 0; j < mygenefile.sequences[i].num_genes; j++) {
			double exp = rand_exponential(10);
			ofs << mygenefile.sequences[i].genes[j].name << "\t" << exp << endl;
			num_trans++;
		}
	}
	ofs.close();
	printf("%d transcript expression values generated.\n", num_trans);
}

void random_reads() {
	fasta myfasta;
	myfasta.read_from_file(local->other_args[0]);
	for (auto itr = myfasta.tags.begin(); itr != myfasta.tags.end(); itr++) {
			size_t pos1 = itr->find("$$");
			if (pos1 != string::npos) {
				size_t pos2 = itr->find("$$", pos1+1);
				if (pos2 == string::npos) {
					//panic(string("bad format in sequence tag : ") + *itr);
					*itr = itr->substr(pos1 + 2);
				} else {
					*itr = itr->substr(pos1 + 2, pos2 - pos1 - 2);
				}
			}
	}
	myfasta.build_tag_map();
	
	ifstream ifs(local->other_args[1]);
	string readline;
	vector<double> prob;
	vector<int> index;
	while (getline(ifs, readline)) {
		vector<string> tokens = string_tokenize(readline);
		if (tokens.size() != 2 || !is_num(tokens[1])) panic(string("bad format in gene expression : ") + readline);
		if (myfasta.tag_map.count(tokens[0]) == 0) panic(string("gene sequence not found : ") + tokens[0]);
		index.push_back(myfasta.tag_map[tokens[0]]);
		prob.push_back(str2double(tokens[1]) * (myfasta.sequences[myfasta.tag_map[tokens[0]]].length() - read_length + 1));
	}
	ifs.close();

	vector<int> read_count = mnrnd(num_random_reads, prob);
	fasta output_fasta;
	int n = 1;
	for (uint i = 0; i < index.size(); i++) {
		for (int j = 0; j < read_count[i]; j++) {
			int pos = rand_int(0, (int)myfasta.sequences[index[i]].length() - read_length);
			string read = myfasta.sequences[index[i]].substr(pos, read_length);
			bool strand = (rand_double() > 0.5);
			if (!strand) read = get_reverse2(read);
			output_fasta.tags.push_back(myfasta.tags[index[i]] + "_" + int2str(n) + "_" + int2str(pos) + (strand?"F":"R"));
			output_fasta.sequences.push_back(read);
			n++;
		}
	}
	if (local->output_prefix == "") local->output_prefix = local->other_args[1];
	output_fasta.write_to_file(local->output_prefix + ".reads.fa");
}

void denovo() {	
	if (!is_region(local->other_args[0])) {
		panic(string("wrong genomic region: ") + local->other_args[0]);
	}
	chr_region region(local->other_args[0]);

	genefile mygenefile;
	if (!mygenefile.read_from_text_file(local->other_args[1])) {
		panic(string("failed reading read file: ") + local->other_args[1]);
	}

	vector<gene_struct> reads;
	if (!mygenefile.search_by_region(region.chr, region.start, region.end, reads)) {
		panic(string("failed retrieval reads"));
	}
	printf("%d reads retrived.\n", (int)reads.size());

	int len = region.end - region.start;
	vector<int> coverage(len, 0);
	map<int, map<int, int> > splices;
	for (int i = 0; i < (int)reads.size(); i++) {
		for (int j = 0; j < (int)reads[i].exons.size(); j++) {
			for (int k = reads[i].exons[j].first; k < reads[i].exons[j].second - 1; k++) {
				coverage[k - region.start]++;
			}
			if (j > 0) {
				map<int, int> targets;
				if (splices.count(reads[i].exons[j-1].second - region.start) > 0) {
					targets = splices[reads[i].exons[j-1].second - region.start];
				}
				if (targets.count(reads[i].exons[j].first - region.start) > 0) {
					targets[reads[i].exons[j].first - region.start] = targets[reads[i].exons[j].first - region.start] + 1;
				} else {
					targets[reads[i].exons[j].first - region.start] = 1;
				}
				splices[reads[i].exons[j-1].second - region.start] = targets;
			}
		}
	}

	double average_coverage;
	while (1) {
		int nonzero_len = 0;
		int sum_coverage = 0;
		for (int i = 0; i < (int)coverage.size(); i++) {
			if (coverage[i] > 0) {
				nonzero_len++;
				sum_coverage += coverage[i];
			}
		}
		if (nonzero_len == 0) {
			panic("no read");
		}
		average_coverage = (double)sum_coverage / nonzero_len;
		printf("covered length = %d, average coverage = %lf.\n", nonzero_len, average_coverage);
		bool changed = false;
		for (int i = 0; i < (int)coverage.size(); i++) {
			if (coverage[i] > 0 && coverage[i] <= max(2, (int)(average_coverage / 20))) {
				coverage[i] = 0;
				changed = true;
			}
		}

		nonzero_len = 0;
		for (int i = 0; i < (int)coverage.size(); i++) {
			if (coverage[i] > 0) {
				nonzero_len++;
			} else {
				if (nonzero_len > 0 && nonzero_len < 10) {
					for (int j = i - nonzero_len; j < i; j++) {
						coverage[j] = 0;
						changed = true;
					}
				}
				nonzero_len = 0;
			}
		}
		if (!changed) break;
	}

	if (local->output_prefix == "") local->output_prefix = local->other_args[1];

	ofstream ofs1((local->output_prefix + ".coverage.txt").c_str());
	int zero_len = 0;
	int nonzero_len = 0;
	int sum_coverage = 0;
	for (int i = 0; i < (int)coverage.size(); i++) {
		if (coverage[i] > 0) {
			if (zero_len > 0) {
				cout << i - zero_len + region.start << " -> " << i + 1 + region.start << " : " << 0 << " x " << zero_len<< endl;
			}
			zero_len = 0;
			//cout << coverage[i] << ",";
			nonzero_len++;
			sum_coverage += coverage[i];
		} else {
			if (nonzero_len > 0) {
				cout << i - nonzero_len + region.start << " -> " << i + 1 + region.start << " : " << (double)sum_coverage / nonzero_len << " x " << nonzero_len << endl;
				ofs1 << region.chr << "\t" << i - nonzero_len + region.start << "\t" << i + 1 + region.start << "\t" << (double)sum_coverage / nonzero_len << endl;
			}
			nonzero_len = 0;
			sum_coverage = 0;
			zero_len++;
		}
	}
	cout << endl;
	ofs1.close();

	for (map<int, map<int, int> >::iterator mi = splices.begin(); mi != splices.end(); mi++) {
		for (map<int, int>::iterator mi1 = mi->second.begin(); mi1 != mi->second.end(); mi1++) {
			if (mi1->second <= max(2, (int)(average_coverage / 20)) || !(mi1->first - mi->first > 30 && mi1->first - mi->first < 400000)) {
				mi1->second = 0;
			}
		}
	}

	ofstream ofs2((local->output_prefix + ".splices.refFlat.txt").c_str());
	for (map<int, map<int, int> >::iterator mi = splices.begin(); mi != splices.end(); mi++) {
		for (map<int, int>::iterator mi1 = mi->second.begin(); mi1 != mi->second.end(); mi1++) {
			if (mi1->second > 0) {
				cout << mi->first + region.start << " -> " << mi1->first + region.start << " x " << mi1->second << endl;
				ofs2 << mi1->second << "\t" << region.chr << "\t+\t" << mi->first + region.start - 1 << "\t" << mi1->first + region.start + 1 << "\t";
				ofs2 << mi->first + region.start - 1 << "\t" << mi1->first + region.start + 1 << "\t2\t";
				ofs2 << mi->first + region.start - 1 << "," << mi1->first + region.start << ",\t" << mi->first + region.start << "," << mi1->first + region.start + 1 << ",\n";
			}
		}
	}
	ofs2.close();
}

void mapability() {
	vector<uint64> mapability_hash_table;
	map<int64, int> map_count;
	int method = 1; //1: hash, 2: map
	vector<int64> reads;
	fasta myfasta;
	ifstream ifs;
	string readline;
	vector<string> tags;
	vector<int> lengths;

	cout << "reading inputs\n";
	reads.clear();
	ifs.open(local->other_args[0].c_str());
	while (getline(ifs, readline)) {
		myfasta.clear();
		myfasta.read_from_file(readline);
		for (int i = 0; i < (int)myfasta.sequences.size(); i++) {
			tags.push_back(myfasta.tags[i]);
			lengths.push_back((int)myfasta.sequences[i].length());
			reads.reserve(reads.size() + myfasta.sequences[i].length());
			int64 index = 0;
			for (int j = 0; j < (int)myfasta.sequences[i].length(); j++) {
				if (j > 0 && j % 10000000 == 0) cout << ".";
				flush(cout);
				char ch = (char)toupper(myfasta.sequences[i][j]);
				index <<= 2;
				index &= (((int64)(1) << (read_length * 2)) - 1);
				if (ch == 'C') index += 1;
				else if (ch == 'G') index += 2;
				else if (ch == 'T') index += 3;
				if (j >= read_length - 1) reads.push_back(index);
			}
		}
		cout << endl;
	}
	ifs.close();
	ifs.clear();
	cout << "done.\n";
	cout << (int64)reads.size() << " valid reads.\n";

	cout << "building index...\n";
	size_t nn = 0;
	if (method == 1) {
		int logNN = (int)ceil(log((double)reads.size() * 2.0) / log(double(2.0)));
		if (logNN >= 32) nn = (uint)(-1);
		else nn = 1 << logNN;
		mapability_hash_table.resize(nn);
		for (size_t i = 0; i < nn; i++) {
			mapability_hash_table[i] = uint64(-1);
		}
	} else if (method == 2) {
		map_count.clear();
	}
	for (uint i = 0; i < reads.size(); i++) {
		int64 index = reads[i];
		if (i > 0 && i % 10000000 == 0) cout << ".";
		flush(cout);
		if (method == 1) {
			size_t key = joaat_hash((uchar*)(&index), 8) % nn;
			size_t origin = key;
			while (true) {
				if (mapability_hash_table[key] == (uint64)(-1)) {
					mapability_hash_table[key] = (((uint64)(i)) << 32) + 1;
					break;
				} else if (reads[(uint)(mapability_hash_table[key] >> 32)] == index) {
					mapability_hash_table[key]++;
					break;
				} else {
					key ++;
					if (key == nn) key = 0;
					__ASSERT (origin != key, "internal error, hash full.\n");
				}
			}
		} else if (method == 2) {
			if (map_count.count(index) == 0) map_count[index] = 1;
			else map_count[index]++;
		}
	}
	cout << "done.\n";

	ofstream ofs;
	if (local->output_prefix == "") local->output_prefix = local->other_args[0];
	ofs.open((local->output_prefix + ".map").c_str());
	cout << "writing outputs...\n";
	uint pos = 0;
	for (int i = 0; i < (int)tags.size(); i++) {
		FILE *file = fopen((tags[i]+".cs").c_str(), "w");
		for (int j = 0; j <= lengths[i] - read_length; j++) {
			if (pos > 0 && pos % 10000000 == 0) cout << ".";
			flush(cout);
			int64 index = reads[pos++];
			if (method == 1) {
				uint key = joaat_hash((uchar*)(&index), 8) % nn;
				uint origin = key;
				while (true) {
					__ASSERT(mapability_hash_table[key] != (uint64)(-1), "internal error, key not found.\n");
					if (reads[(uint)(mapability_hash_table[key] >> 32)] == index) {
						if (((mapability_hash_table[key] << 32) >> 32) > 1) {
							ofs << tags[i] << "\t" << j << "\t" << ((mapability_hash_table[key] << 32) >> 32) - 1 << endl;
						}
						uchar score = (uchar) (255.0 / ((mapability_hash_table[key] << 32) >> 32));
						fwrite(&score, 1, 1, file);
						break;
					} else {
						key ++;
						if (key == nn) key = 0;
						__ASSERT (origin != key, "internal error, hash full.\n");
					}
				}
			} else if (method == 2) {
				if (map_count[index] > 1) {
					ofs << myfasta.tags[i] << "\t" << j << "\t" << map_count[index] - 1 << endl;
				}
				uchar score = (uchar) (255.0 / (map_count[index]));
				fwrite(&score, 1, 1, file);
			}
		}
		fclose(file);
	}
	__ASSERT(pos == reads.size(), "internal error, wrong pos.\n");
	ofs.close();
	ofs.clear();
	cout << "done.\n";
	if (method == 1) {
		mapability_hash_table.clear();
	} else if (method == 2) map_count.clear();
	reads.clear();
}

void differential() {
	ofstream ofs;

	if (local->output_prefix == "") local->output_prefix = local->other_args[0];

	if (gene_name.size() != 1) {
		ofs.open((local->output_prefix + ".diff.xls").c_str());
//		ofs << "gene-name\tstat-T\tp-value-T\tstat-G\tp-value-G\tnum-isoforms\tiso-exp-pooled\tgene-exp-pooled\tiso-exp-1\tgene-exp-1\tiso-exp-2\tgene-exp-2" << endl;
		ofs << "gene-name\tstat-T\tp-value-T\tnum-isoforms\tiso-exp-1\tgene-exp-1\tiso-exp-2\tgene-exp-2" << endl;	
	}

	ifstream ifs[2];
	string readline;
	vector<string> tokens;

//mapped reads
	int NN[2], read_len[2];

	for(int j = 0; j < 2; j++) {
		ifs[j].open(local->other_args[j].c_str());
		if (old) {
			getline(ifs[j], readline);
			__ASSERT (2 == sscanf(readline.c_str(), " %d %d\n", &NN[j], &read_len[j]), "error reading file head.\n");
		}
	}

	while (true) {
//gene name
		string name[2];
//#exons
		int m[2];
//#isoforms
		int n[2];
//#categories
		int k[2];
//#possible categories
		int kk[2];
//isoform names
		vector<string> trans[2];
//exon lengths
		vector<int> L[2];
//RPKM
		double rpkm[2];
//MLE
		vector<double> mle[2];
//SUM MLE
		double sum_mle[2];
//sampling rates
		vector<vector<double> > rates[2];
//read counts
		vector<double> counts[2];
//possible sampling rates
		vector<vector<double> > possible_rates[2];
		vector<vector<double> > possible_rates_t[2];
//idicator matrix
		vector<vector<int> > D[2];

		for (int j = 0; j < 2; j++) {
			if (!getline(ifs[j], readline)) {
				goto finished;
			}
			char buf[1025];

			if (old) {
				__ASSERT (5 == sscanf(readline.c_str(), " %s %d %d %d %d\n", buf, &m[j], &n[j], &k[j], &kk[j]), "error reading gene information.\n");
			} else {
				__ASSERT (3 == sscanf(readline.c_str(), " %s %d %d\n", buf, &n[j], &k[j]), "error reading gene information.\n");
			}

			name[j] = buf;
			getline(ifs[j], readline);
			trans[j] = string_tokenize(readline);

			if (old) {
				getline(ifs[j], readline);
				tokens = string_tokenize(readline);
				__ASSERT(m[j] == (int)tokens.size(), "error reading exon lengths.\n");
				for (int i = 0; i < m[j]; i++) {
					__ASSERT(is_int(tokens[i]), "error: wrong exon length.\n");
					L[j].push_back(str2int(tokens[i]));			
				}
				getline(ifs[j], readline);
				__ASSERT (1 == sscanf(readline.c_str(), "RPKM %lf\n", &rpkm[j]), "error reading RPKM.\n");
				if (kk[j] > 0) {
					getline(ifs[j], readline);
					tokens = string_tokenize(readline);
					__ASSERT((int)tokens.size() == n[j] + 1 && tokens[0] == "MLE", "error reading MLE.\n");
					for (int i = 0; i < n[j]; i++) {
						mle[j].push_back(str2double(tokens[i+1]));
					}
					getline(ifs[j], readline);
					__ASSERT (1 == sscanf(readline.c_str(), "SUM %lf\n", &sum_mle[j]), "error reading SUM MLE.\n");
				}
			}

			for (int i = 0; i < k[j]; i++){
				getline(ifs[j], readline);
				tokens = string_tokenize(readline);
				__ASSERT((int)tokens.size() == n[j] + 1, "error reading sampling rates.\n");
				vector<double> rate_vector;
				for (int l = 0; l < n[j]; l++) {
					rate_vector.push_back(str2double(tokens[l]));
				}
				counts[j].push_back(str2double(tokens[n[j]]));
				rates[j].push_back(rate_vector);
			}

			if (old) {
				for (int i = 0; i < kk[j]; i++){
					getline(ifs[j], readline);
					tokens = string_tokenize(readline);
					__ASSERT((int)tokens.size() == n[j], "error reading sampling rates.\n");
					vector<double> rate_vector;
					for (int l = 0; l < n[j]; l++) {
						rate_vector.push_back(str2double(tokens[l]));
					}
					possible_rates[j].push_back(rate_vector);
				}
				for (int i = 0; i < n[j]; i++){
					getline(ifs[j], readline);
					tokens = string_tokenize(readline);
					__ASSERT((int)tokens.size() == m[j], "error reading sampling rates.\n");
					vector<int> D_vector;
					for (int l = 0; l < m[j]; l++) {
						D_vector.push_back(str2int(tokens[l]));
					}
					D[j].push_back(D_vector);
				}
			}
		}

		if (old) {
			__ASSERT (name[0] == name[1] && m[0] == m[1] && n[0] == n[1] && trans[0] == trans[1] && L[0] == L[1] && D[0] == D[1], "error: inconsistent gene information.\n");
		} else {
			__ASSERT (name[0] == name[1] && n[0] == n[1] && trans[0] == trans[1], "error: inconsistent gene information.\n");
		}

		if (gene_name.empty() || gene_name.count(tolower(name[0])) > 0) {
			if (old) {
				if (m[0] > 1 && n[0] == 1) {
					panic(string("error: m > 1 && n == 1 for gene ") + name[0]);
				}
				if (m[0] == 1 && n[0] > 1) {
					cout << "warning: m == 1 && n > 1 for gene " << name[0] << endl;
				}
			} else {
				if (k[0] > 1 && n[0] == 1) {
					panic(string("error: k > 1 && n == 1 for gene ") + name[0]);
				}
				if (k[0] == 1 && n[0] > 1) {
					cout << "warning: k == 1 && n > 1 for gene " << name[0] << endl;
				}
			}

			if (gene_name.size() == 1) {
				cout << name[0] << endl;
			}

			if (old) {
	//			if (kk[0] > 0 && kk[1] > 0 && n[0] <= 5) {
	//			if (kk[0] > 0 && kk[1] > 0) {
				if (kk[0] == 0 || kk[1] == 0 || (quick && ((n[0] > num_max_isoform && num_max_isoform > 0) || sum(counts[0]) < 20 || sum(counts[1]) < 20))) {
					if (gene_name.size() != 1) {
						ofs << " " << name[0] << "\tNA\tNA\t" << n[0] << "\n";
					}
					continue;
				}
			} else {
				if (k[0] == 0 || k[1] == 0 || (quick && ((n[0] > num_max_isoform && num_max_isoform > 0) || sum(counts[0]) < 20 || sum(counts[1]) < 20))) {
					if (gene_name.size() != 1) {
						ofs << " " << name[0] << "\tNA\tNA\t" << n[0] << "\n";
					}
					continue;
				}
			}

			vector<double> N[2];
			//int num_reads[2];
			vector<double> X_opt[2];
			for (int j = 0; j < 2; j++) {
				
				if (old) {
					vector<vector<double> > normalized_rates = rates[j], normalized_possible_rates = possible_rates[j];
					for (int i = 0; i < k[j]; i++) {
						L1_normalize(normalized_rates[i]);
					}
					for (int i = 0; i < kk[j]; i++) {
						L1_normalize(normalized_possible_rates[i]);
					}
					N[j].resize(kk[j]);
					for (int i = 0; i < kk[j]; i++) {
						N[j][i] = 0;
					}
					for (int i = 0; i < k[j]; i++) {
						vector<pair<double, int> > distances;
						for (int l = 0; l < kk[j]; l++) {
							distances.push_back(pair<double, int>(L1_dist(normalized_rates[i], normalized_possible_rates[l]), l));
						}
						sort(distances.begin(), distances.end());
						if (kk[j] > 0 && distances[0].first < 1e-2) {
							N[j][distances[0].second] += counts[j][i];
						} else {
							printf("warning: category not found, gene: %s\n", name[0].c_str());
						}

	/*					bool found = false;
						for (int l = 0; l < kk[j]; l++) {
							if (nearly_equal(normalized_rates[i], normalized_possible_rates[l], 1e-2)) {
								found = true;
								N[j][l] += counts[j][i];
								break;
							}
						}
						if (!found) {
							printf("warning: category not found, gene: %s\n", name[0].c_str());
						}*/
					}

	/*				if (quick) { //remove shared parts of isoforms //remove small categories with zero reads
						int i = 0;
						while (i < kk[j]) {
							bool shared_category = true;
							for (int l = 1; l < n[j]; l++) {
								if (fabs(possible_rates[j][i][l] - possible_rates[j][i][l-1]) > epsilon) {
									shared_category = false;
									break;
								}
							}
							if (shared_category || N[j][i] == 0 && L1_norm(possible_rates[j][i]) < 10) {
								possible_rates[j].erase(possible_rates[j].begin() + i);
								N[j].erase(N[j].begin() + i);
								kk[j] = kk[j] - 1;
							} else {
								i++;
							}
						}
					}*/

					for (int i = 0; i < kk[j]; i++) {
						for (int l = 0; l < n[j]; l++) {
							possible_rates[j][i][l] = possible_rates[j][i][l] / 1e3 * NN[j] / 1e6;
						}
					}
				} else {
					possible_rates[j] = rates[j];
					kk[j] = k[j];
					N[j] = counts[j];
				}

				//num_reads[j] = (int)sum(N[j]);
				possible_rates_t[j] = transpose(possible_rates[j]);
				solve_likelihood(possible_rates_t[j], N[j], X_opt[j]);
				
				for (int i = 0; i < n[j]; i++) {
					if (old && !quick && fabs(X_opt[j][i] - mle[j][i]) > 1) {
						printf("warning: inaccurate caculation, gene: %s\n", name[0].c_str());
					}
				}

				if (gene_name.size() == 1) {
					cout << "X_opt[" << j << "]:\t";
					for (int i = 0; i < n[j]; i++) {
						cout << X_opt[j][i] << "\t";
					}
					cout << endl;
					
					if (old) {
						cout << "iso_opt[" << j << "]:\t";
						for (int i = 0; i < n[j]; i++) {
							cout << mle[j][i] << "\t";
						}
						cout << endl;
					}
				}
			}

			if (quick && (sum(N[0]) < 10 || sum(N[1]) < 10)) {
				if (gene_name.size() != 1) {
					ofs << " " << name[0] << "\tNA\tNA\t" << n[0] << "\n";
				}
				continue;
			}

			double T = 0;
			double sum1 = sum(X_opt[0]), sum2 = sum(X_opt[1]);

			if (quick && !(sum1 > 0 && sum2 > 0)) {
				panic(string("internal error: gene not expressed, zero reads. gene: ") + name[0]);
			}

			for (int i = 0; i < n[0]; i++) {
				T += fabs(X_opt[0][i] / sum1 - X_opt[1][i] / sum2);
			}
//				double G = fabs(sum1 - sum2);
			if (gene_name.size() == 1) {
				cout << "T: " << T << endl;
//					cout << "G: " << G << endl;
			}

			if (T < epsilon || (quick && T < .5)) {
				if (gene_name.size() != 1) {
					ofs << " " << name[0] << "\t" << T << "\tNA\t" << n[0] << "\t";
				}
				for (int j = 0; j < 2; j++) {
					for (int i = 0; i < n[j]; i++) {
						ofs << X_opt[j][i] << ",";
					}
					ofs << "\t" << sum(X_opt[j]) << "\t";
				}
				ofs << endl;
				continue;
			}

			vector<double> combined_N = N[0];
			combined_N.insert(combined_N.end(), N[1].begin(), N[1].end());
			vector<vector<double> > combined_possible_rates = possible_rates[0];
			combined_possible_rates.insert(combined_possible_rates.end(), possible_rates[1].begin(), possible_rates[1].end());
			vector<double> combined_X_opt;
			solve_likelihood_t(combined_possible_rates, combined_N, combined_X_opt);
			if (gene_name.size() == 1) {
				cout << "combined_rates:\n";
				for (int i = 0; i < n[0]; i++) {
					for (int j = 0; j < kk[0]+kk[1]; j++) {
						cout << combined_possible_rates[j][i] << "\t";
					}
					cout << endl;
				}
				cout << "combined_N:\t";
				for (int i = 0; i < kk[0]+kk[1]; i++) {
					cout << combined_N[i] << "\t";
				}
				cout << endl;
				cout << "combined_X_opt:\t";
				for (int i = 0; i < n[0]; i++) {
					cout << combined_X_opt[i] << "\t";
				}
				cout << endl;
			}
			
/*				vector<double> p_rates[2];
			for (int j = 0; j < 2; j++) {
				p_rates[j].resize(kk[j]);
				for (int i = 0; i < kk[j]; i++) {
					p_rates[j][i] = 0;
					for (int l = 0; l < n[j]; l++) {
						p_rates[j][i] += possible_rates_t[j][l][i] * combined_X_opt[l];
					}
				}
			}
			vector<double> combined_p_rates = p_rates[0];
			combined_p_rates.insert(combined_p_rates.end(), p_rates[1].begin(), p_rates[1].end());*/

			vector<double> T_null, G_null;
			const double p = 0.05;
			const double c = 2 * p / (1 - p);
			const double a = 10;
			int itr;
			double p_value_T;
//				double p_value_G;

/*			vector<double> read_count_rates;
			read_count_rates.push_back((double)num_reads[0]/(num_reads[0] + num_reads[1]));
			read_count_rates.push_back((double)num_reads[1]/(num_reads[0] + num_reads[1]));*/
//			read_count_rates.push_back((double)NN[0]/(NN[0]+NN[1]));
//			read_count_rates.push_back((double)NN[1]/(NN[0]+NN[1]));
/*			if (gene_name.size() == 1) {
				cout << "read_count_rates:\t" << read_count_rates[0] << "\t" << read_count_rates[1] << endl;
			}
			vector<vector<double> > separate_combined_possible_rates_t[2];
			for (int j = 0; j < 2; j++) {
				separate_combined_possible_rates_t[j].resize(n[j]);
				for (int i = 0; i < n[j]; i++) {
					separate_combined_possible_rates_t[j][i].resize(kk[0]+kk[1]);
					for (int l = 0; l < kk[0]+kk[1]; l++) {
						separate_combined_possible_rates_t[j][i][l] = combined_possible_rates[l][i] * (read_count_rates[0] * NN[0] + read_count_rates[1] * NN[1]) / (NN[0] + NN[1]);
					}
				}
			}*/

			if (reset_srand) {
				srand(8729343);
			}

			for (itr = 0; itr < num_itr; itr++) {
				vector<double> temp_X_opt[2];
				vector<double> categories[2];
				for (int j = 0; j < 2; j++) {
//					vector<int> iso = mnrnd(num_reads[j], combined_X_opt);
					categories[j].resize(kk[j]);
					for (int l = 0; l < kk[j]; l++) {
						if (N[j][l] == 0) {
							categories[j][l] = 0;
						} else {
							double pois_rate = inner_prod(possible_rates[j][l], combined_X_opt);
							pois_rate *= sum(X_opt[j]) / sum(combined_X_opt);
							pois_rate *= N[j][l] / inner_prod(possible_rates[j][l], X_opt[j]);
							categories[j][l] = rand_poisson(pois_rate);
						}
					}
                    
/*					for (int i = 0; i < n[j]; i++) {
					vector<int> temp = mnrnd(iso[i], possible_rates_t[j][i]);
					for (int l = 0; l < kk[j]; l++) {
						category[l] += temp[l];
					}
				}*/
//						vector<int> counts = mnrnd(num_reads[j], p_rates[j]);
//						vector<double> category(kk[j], 0);
//						for (int l = 0; l < kk[j]; l++) {
//							category[l] = counts[l];
//						}
					solve_likelihood(possible_rates_t[j], categories[j], temp_X_opt[j]);
				}

/*					vector<int> counts = mnrnd(num_reads[0] + num_reads[1], combined_p_rates);
				vector<int> counts = mnrnd(num_reads[0] + num_reads[1], combined_p_rates);
				vector<double> category;
				category.resize(kk[0]);
				for (int l = 0; l < kk[0]; l++) {
					category[l] = counts[l];
				}
				solve_likelihood(possible_rates_t[0], category, temp_X_opt[0]);
				category.resize(kk[1]);
				for (int l = 0; l < kk[1]; l++) {
					category[l] = counts[l + kk[0]];
				}
				solve_likelihood(possible_rates_t[1], category, temp_X_opt[1]);*/

/*
				for (int j = 0; j < 2; j++) {
					for (int l = 0; l < kk[j]; l++) {
						vector<int> counts = mnrnd((int)N[j][l], read_count_rates);
						categories[0].push_back(counts[0]);
						categories[1].push_back(counts[1]);
					}
				}

				for (int j = 0; j < 2; j++) {
					solve_likelihood(separate_combined_possible_rates_t[j], categories[j], temp_X_opt[j]);
				}*/
				
				double sum1 = sum(temp_X_opt[0]), sum2 = sum(temp_X_opt[1]);
				if (!(sum1 > 0 && sum2 > 0)) {
					printf("warning: gene not expressed, zero read sampled. gene: %s\n", name[0].c_str());
					sum1 += epsilon;
					sum2 += epsilon;
				}
				double temp_T = 0;
				for (int i = 0; i < n[0]; i++) {
					temp_T += fabs(temp_X_opt[0][i] / sum1 - temp_X_opt[1][i] / sum2);
				}
//					double temp_G = fabs(sum1 - sum2);

				if (gene_name.size() == 1 && itr < 3) {
					for (int j = 0; j < 2; j++) {
/*						cout << "rates[" << j << "]\n";
						for (int i = 0; i < n[0]; i++) {
							for (int l = 0; l < kk[0]+kk[1]; l++) {
								cout << separate_combined_possible_rates_t[j][i][l] << "\t";
							}
							cout << endl;
						}*/
						cout << "categories[" << j << "]:\t";
						for (int l = 0; l < (int)categories[j].size(); l++) {
							cout << categories[j][l] << "\t";
						}
						cout << endl;
						cout << "temp_X_opt[" << j << "]:\t";
						for (int i = 0; i < n[j]; i++) {
							cout << temp_X_opt[j][i] << "\t";
						}
						cout << endl;
					}
					cout << "temp_T:\t" << temp_T << endl;
				}

				T_null.push_back(temp_T);
//					G_null.push_back(temp_G);
				sort(T_null.begin(), T_null.end());
//					sort(G_null.begin(), G_null.end());

				if (do_stop) {
					double p_value_T = (double)(T_null.end() - lower_bound(T_null.begin(), T_null.end(), T)) / (itr + 1);
					if ((itr + 1) * (p_value_T - c * (1 - p_value_T)) > a) {
						break;
					}
				}
			}
			sort(T_null.begin(), T_null.end());
//				sort(G_null.begin(), G_null.end());
			if (do_stop && itr < num_itr) {
				p_value_T = (double)(T_null.end() - lower_bound(T_null.begin(), T_null.end(), T)) / (itr + 1);
//					p_value_G = (double)(G_null.end() - lower_bound(G_null.begin(), G_null.end(), G)) / (itr + 1);
			} else {
				p_value_T = (double)(T_null.end() - lower_bound(T_null.begin(), T_null.end(), T)) / num_itr;
//					p_value_G = (double)(G_null.end() - lower_bound(G_null.begin(), G_null.end(), G)) / num_itr;
			}
			if (gene_name.size() == 1) {
				cout << "itr: " << itr << endl;
				cout << "p_value_T: " << p_value_T << endl;
//					cout << "p_value_G: " << p_value_G << endl;
			}
			if (gene_name.size() != 1) {
//					ofs << " " << name[0] << "\t" << T << "\t" << p_value_T << "\t" << G << "\t" << p_value_G << "\t" << n[0] << "\t";
				ofs << " " << name[0] << "\t" << T << "\t" << p_value_T << "\t" << n[0] << "\t";
/*				for (int i = 0; i < n[0]; i++) {
					ofs << combined_X_opt[i] << ",";
				}
				ofs << "\t" << sum(combined_X_opt) << "\t";*/
				for (int j = 0; j < 2; j++) {
					for (int i = 0; i < n[j]; i++) {
						ofs << X_opt[j][i] << ",";
					}
					ofs << "\t" << sum(X_opt[j]) << "\t";
				}
				ofs << endl;
			}
		}
	}

finished:

	for (int j = 0; j < 2; j++) {
		ifs[j].close();
		ifs[j].clear();
	}
	if (gene_name.size() != 1) {
		ofs.close();
		ofs.clear();
	}
}

void scan_reads_handler(vector<map_result_struct> &results) {
	for (int i = 0; i < (int)results.size(); i++) {
		string read = trim_space(toupper(results[i].probe_seq));
		int read_len = (int)read.length();
		num_total_read++;
		if (read_len > MAX_READ_LEN) {
			printf("read is too long, failed.\n");
			return;
		}
		min_read_len = min(min_read_len, read_len);
		max_read_len = max(max_read_len, read_len);
		bool is_good_read = true;
		for (int i = 0; i < read_len; i++) {
			char ch = read[i];
			if (ch == 'A') {
				local->nt_counts[i][0]++;
			} else if (ch == 'C') {
				local->nt_counts[i][1]++;
				total_gc++;
			} else if (ch == 'G') {
				local->nt_counts[i][2]++;
				total_gc++;
			} else if (ch == 'T') {
				local->nt_counts[i][3]++;
			} else {
				is_good_read = false;
				local->nt_counts[i][4]++;
			}
			local->nt_counts[i][5]++;
			total_nt++;
		}
		if (is_good_read) num_good_read++;
	}
}

void scan_reads() {
	__ASSERT(file_exists(local->other_args[0]), string("error: file ") + local->other_args[0] + " does not exist.\n");
	int file_type;
	__ASSERT(is_read_file(local->other_args[0], file_type), "error: input file format unrecognizable.\n");

	local->nt_counts.resize(MAX_READ_LEN);
	for (int i = 0; i < MAX_READ_LEN; i++) {
		local->nt_counts[i].resize(6);
		for (int j = 0; j < 6; j++) {
			local->nt_counts[i][j] = 0;
		}
	}

	printf ("reading input file %s...\n", local->other_args[0].c_str());
	ifstream ifs(local->other_args[0].c_str());
	string read;
	while(get_read(ifs, read, file_type)) {
		read = trim_space(toupper(read));
		int read_len = (int)read.length();
		num_total_read++;
		if (read_len > MAX_READ_LEN) {
			printf("read is too long, failed.\n");
			return;
		}
		min_read_len = min(min_read_len, read_len);
		max_read_len = max(max_read_len, read_len);
		bool is_good_read = true;
		for (int i = 0; i < read_len; i++) {
			char ch = read[i];
			if (ch == 'A') {
				local->nt_counts[i][0]++;
			} else if (ch == 'C') {
				local->nt_counts[i][1]++;
				total_gc++;
			} else if (ch == 'G') {
				local->nt_counts[i][2]++;
				total_gc++;
			} else if (ch == 'T') {
				local->nt_counts[i][3]++;
			} else {
				is_good_read = false;
				local->nt_counts[i][4]++;
			}
			local->nt_counts[i][5]++;
			total_nt++;
		}
		if (is_good_read) num_good_read++;
	}
	ifs.close();
	ifs.clear();
	printf ("done.\n");

	printf("total %d reads, %d (%d%%) good reads, %d (%d%%) bad reads.\n", num_total_read, num_good_read, round_double(num_good_read*100.0/num_total_read), num_total_read-num_good_read, 100-round_double(num_good_read*100.0/num_total_read));
	printf("minimum read length %d, maximum read length %d.\n", min_read_len, max_read_len);
	if (min_read_len != max_read_len) printf("INCONSISTENT READ LENGTHS!!!\n");
	if (max_read_len == 0 || num_total_read == 0) {
		printf("NO READ, FAILED!!!\n");
		return;
	}
	printf("GC content %d%%.\n", round_double(total_gc*100.0/total_nt));

	printf("nucleotide distribution by position:\n");
	printf("A\t\tC\t\tG\t\tT\t\tN\n");
	for (int i = 0; i < max_read_len; i++) {
		for (int j = 0; j < 5; j++) {
			printf("%d (%d%%)\t", local->nt_counts[i][j], round_double(local->nt_counts[i][j]*100.0/local->nt_counts[i][5]));
		}
		printf("\n");
	}
}

void check_rep() {
	printf("read input files.\n");
    
	bar bar1, bar2;
	
	__ASSERT(bar1.read_from_text_file(local->other_args[0]), string("error reading text file ") + local->other_args[0]);
	__ASSERT(bar2.read_from_text_file(local->other_args[1]), string("error reading text file ") + local->other_args[1]);
    
	if (!t_test) {
	// do modified K-S test

		printf("computing modified Kolmogorov - Smirnov tests.\n");
	    
	//	__ASSERT(bar1.sequences.size() == bar2.sequences.size(), "unequal number of sequences");

		printf("sequence\tK_S_p_value\n");
		int i = 0, k = 0;
		while (i < (int)bar1.sequences.size() && k < (int)bar2.sequences.size()) {
			if (bar1.sequences[i].name < bar2.sequences[k].name) {
				printf("%s\tNA\n", bar1.sequences[i].name.c_str());
				i++;
			} else if (bar1.sequences[i].name > bar2.sequences[k].name) {
				printf("%s\tNA\n", bar2.sequences[k].name.c_str());
				k++;
			} else {
				__ASSERT(bar1.sequences[i].name == bar2.sequences[k].name, "different sequences names");
				double p = -1.0;
				if (bar1.sequences[i].data.size() >= 20 && bar2.sequences[k].data.size() >= 20) {
					vector<double> a, b;
					for (int j = 0; j < (int)bar1.sequences[i].data.size(); j++) {
						__ASSERT(bar1.sequences[i].data[j].size() == 2, "ERROR: wrong data column");
						__ASSERT(bar1.sequences[i].data[j][1].data_float == 1.0, "ERROR: count != 1");
						a.push_back((double)bar1.sequences[i].data[j][0].data_int + rand_double());
					}
					for (int j = 0; j < (int)bar2.sequences[k].data.size(); j++) {
						__ASSERT(bar2.sequences[k].data[j].size() == 2, "ERROR: wrong data column");
						__ASSERT(bar2.sequences[k].data[j][1].data_float == 1.0, "ERROR: count != 1");
						b.push_back((double)bar2.sequences[k].data[j][0].data_int + rand_double());
					}
					p = K_S_Test(a, b);
				}
				if (p < 0) {
					printf("%s\tNA\n", bar1.sequences[i].name.c_str());
				} else {
					printf("%s\t%e\n", bar1.sequences[i].name.c_str(), p);
				}
				i++;
				k++;
			}
		}
	} else if (t_test) {

	// do t-test
		printf("computing t test.\n");
		vector<double> a;
		double q11, q12, q21, q22;
		a.clear();
		int i, k;
		for (i = 0; i < (int)bar1.sequences.size(); i++) {
			for (int j = 0; j < (int)bar1.sequences[i].data.size(); j++) {
				a.push_back((double)bar1.sequences[i].data[j][0].data_int);
			}
		}
		sort(a.begin(), a.end());
		q11 = get_quantile(a, 0.05);
		q12 = get_quantile(a, 0.95);

		a.clear();
		for (i = 0; i < (int)bar2.sequences.size(); i++) {
			for (int j = 0; j < (int)bar2.sequences[i].data.size(); j++) {
				a.push_back((double)bar2.sequences[i].data[j][0].data_int);
			}
		}
		sort(a.begin(), a.end());
		q21 = get_quantile(a, 0.05);
		q22 = get_quantile(a, 0.95);
		a.clear();

		printf("sequence\tt_p_value\n");
		i = 0; k = 0;
		while (i < (int)bar1.sequences.size() && k < (int)bar2.sequences.size()) {
			if (bar1.sequences[i].name < bar2.sequences[k].name) {
				printf("%s\tNA\n", bar1.sequences[i].name.c_str());
				i++;
			} else if (bar1.sequences[i].name > bar2.sequences[k].name) {
				printf("%s\tNA\n", bar2.sequences[k].name.c_str());
				k++;
			} else {
				__ASSERT(bar1.sequences[i].name == bar2.sequences[k].name, "different sequences names");
				double p = -1.0;
				if (bar1.sequences[i].data.size() >= 20 && bar2.sequences[k].data.size() >= 20) {
					int c11 = 0, c12 = 0, c21 = 0, c22 = 0;
					for (int j = 0; j < (int)bar1.sequences[i].data.size(); j++) {
						__ASSERT(bar1.sequences[i].data[j].size() == 2, "ERROR: wrong data column");
						__ASSERT(bar1.sequences[i].data[j][1].data_float == 1.0, "ERROR: count != 1");
						if (bar1.sequences[i].data[j][0].data_int < q11) c11++;
						if (bar1.sequences[i].data[j][0].data_int > q12) c12++;
					}
					for (int j = 0; j < (int)bar2.sequences[k].data.size(); j++) {
						__ASSERT(bar2.sequences[k].data[j].size() == 2, "ERROR: wrong data column");
						__ASSERT(bar2.sequences[k].data[j][1].data_float == 1.0, "ERROR: count != 1");
						if (bar2.sequences[k].data[j][0].data_int < q21) c21++;
						if (bar2.sequences[k].data[j][0].data_int > q22) c22++;
					}
					int c1 = (int)bar1.sequences[i].data.size(), c2 = (int)bar2.sequences[k].data.size();
					double p1 = (double)(c11 + c12) / c1, p2 = (double)(c21 + c22) / c2;
					double sd = sqrt(p1*(1-p1)/c1 + p2*(1-p2)/c2);
					if (sd > 0) {
						double t = (p1 - p2) / sd;
						p = NormalProb(fabs(t)) * 2;
					}
				}
				if (p < 0) {
					printf("%s\tNA\n", bar1.sequences[i].name.c_str());
				} else {
					printf("%s\t%e\n", bar1.sequences[i].name.c_str(), p);
				}
				i++;
				k++;
			}
		}
	}
}

void trim_UTR() {
	genefile mygenefile;
	mygenefile.read_from_text_file(local->other_args[0]);
	mygenefile.trim_UTR();
	if (local->output_prefix == "") local->output_prefix = local->other_args[0];
	mygenefile.write_to_text_file(local->output_prefix + ".trim_UTR");
}

void enumerate_fasta() {
	fasta myfasta;
	myfasta.read_from_file(local->other_args[0]);
	if (local->output_prefix == "") local->output_prefix = local->other_args[0];
	myfasta.enumerate(local->output_prefix + "." + int2str(read_length) + ".enum", read_length);
}

void gen_trans(bool reverse_transcripts = false) { // new version, which reverse transcripts if necessary
	genefile mygenefile;
	__ASSERT(mygenefile.read_from_text_file(local->other_args[0]), string("error reading annotation file: ") + local->other_args[0] + ".\n");
	mygenefile.build_map();
	fasta myfasta, outputfasta;
	outputfasta.clear();
	bool is_fasta = is_fasta_file(local->other_args[1], false);
	ifstream ifs(local->other_args[1].c_str());
	string readline;
	while (getline(ifs, readline)) {
		myfasta.clear();
		if (is_fasta) {
			myfasta.read_from_file(local->other_args[1]);
		}
		else {
			myfasta.read_from_file(readline);
		}
		for (map<string, pair<int, int> >::iterator mi = mygenefile.gene_map.begin(); mi != mygenefile.gene_map.end(); mi++) {
			int i = mi->second.first;
			int j = mi->second.second;
			string sequence = myfasta.get_sequence(mygenefile.sequences[i].name, mygenefile.sequences[i].genes[j].exons);
			if (sequence != "") {
				if (reverse_transcripts && !mygenefile.sequences[i].genes[j].strand) {
					sequence = get_reverse2(sequence);
				}
				outputfasta.tags.push_back(mygenefile.sequences[i].genes[j].geneName + "$$" + mygenefile.sequences[i].genes[j].name);
				outputfasta.sequences.push_back(sequence);
			} 
		}
		if (is_fasta) break;
	}
	ifs.close();
	ifs.clear();
	if (local->output_prefix == "") local->output_prefix = local->other_args[0];
	outputfasta.write_to_file(local->output_prefix + ".fa");
}

void gen_exons_junctions() {
	exon_junction_extractor my_exon_junction_extractor;
	my_exon_junction_extractor.load_refFlat_file(local->other_args[0]);
	if (local->output_prefix == "") local->output_prefix = local->other_args[0];
	my_exon_junction_extractor.write_to_text_file(local->output_prefix + ".subexons.txt");

	fasta myfasta, outputfasta;
	outputfasta.clear();
	ifstream ifs(local->other_args[1].c_str());
	string readline;
	string output_filename = local->output_prefix + ".exons_junctions.fa";
	FILE *file=fopen(output_filename.c_str(), "wt");
	if (!file){
		cout << "error opening file for writing:" << output_filename << endl;
		return;
	}
	while (getline(ifs, readline)) {
		myfasta.clear();
		myfasta.read_from_file(readline);
		for (int i = 0; i < (int)my_exon_junction_extractor.genes.size(); i++) {
			string sequence;
			if (version == 2)
				my_exon_junction_extractor.get_sequence(sequence, i, myfasta, read_length - 5);
			else if (version == 3)
				my_exon_junction_extractor.get_sequence(sequence, i, myfasta, read_length - 5, 5, true, false, false);
			else if (version == 4)
				my_exon_junction_extractor.get_sequence(sequence, i, myfasta, read_length - 5, 5, true, true, false);
			//sequence = myfasta.get_sequence(my_exon_junction_extractor.genes[i].trans[0].chr, my_exon_junction_extractor.genes[i].subexons, 5);
			if (sequence != "") {
				outputfasta.write_to_file(file, my_exon_junction_extractor.genes[i].name, sequence);
			}
		}
	}
	fclose(file);
	ifs.close();
	ifs.clear();
}

void gene_exp_handler(vector<target_struct> &_targets) {
	if (version == 1) {
		set<string> target_genes;
		set<string> target_trans;
		target_genes.clear();
		target_trans.clear();
		for (int i = 0; i < (int)_targets.size(); i++) {
			size_t pos = _targets[i].trans_id.find_first_of("$$");
			string gene_name;
			if (pos == string::npos) {
//				panic("no : in target");
//				return; //ribosomal
				gene_name = _targets[i].trans_id;
			} else {
				gene_name = _targets[i].trans_id.substr(0,pos);
			}
			
			if (genes.find(gene_name) != genes.end()) { 
				target_genes.insert(gene_name);
			} else { //Ribosomal etc
				return;
			}
//			string trans_name = _targets[i].trans_id.substr(0, pos);
//			target_trans.insert(trans_name);
		}

		if (first_phase) {
			for (set<string>::iterator it = target_genes.begin(); it != target_genes.end(); it++) {
				__ASSERT(genes.find(*it) != genes.end(), string("gene not found, gene: " + *it));
				if (target_genes.size() == 1) {
					genes[*it].unique_count += 1;
				}
			}
/*			for (set<string>::iterator it = target_trans.begin(); it != target_trans.end(); it++) {
				__ASSERT(transcripts.find(*it) != transcripts.end(), string("transcript not found, transcript: " + *it));
				if (target_trans.size() == 1) {
					transcripts[*it].unique_count += 1;
				}
			}*/
		} else {
			double total_ucount_density = 0;
			for (set<string>::iterator it = target_genes.begin(); it != target_genes.end(); it++) {
				__ASSERT(genes.find(*it) != genes.end(), string("gene not found, gene: " + *it));
				total_ucount_density += (double)genes[*it].unique_count / genes[*it].length;
			}
			for (set<string>::iterator it = target_genes.begin(); it != target_genes.end(); it++) {
				__ASSERT(genes.find(*it) != genes.end(), string("gene not found, gene: " + *it));
				if (total_ucount_density == 0) {
					genes[*it].total_count += 1.0 / target_genes.size();
				} else {
					genes[*it].total_count += (double)genes[*it].unique_count / genes[*it].length / total_ucount_density;
				}
			}

/*			total_ucount_density = 0;
			for (set<string>::iterator it = target_trans.begin(); it != target_trans.end(); it++) {
				__ASSERT(transcripts.find(*it) != transcripts.end(), "transcript not found");
				total_ucount_density += (double)transcripts[*it].unique_count / transcripts[*it].length;
			}
			for (set<string>::iterator it = target_trans.begin(); it != target_trans.end(); it++) {
				__ASSERT(transcripts.find(*it) != transcripts.end(), "transcript not found");
				if (total_ucount_density == 0) {
					transcripts[*it].total_count += 1.0 / target_trans.size();
				} else {
					transcripts[*it].total_count += (double)transcripts[*it].unique_count / transcripts[*it].length / total_ucount_density;
				}
			}*/
		}
	} else if (version >= 2) {
		vector<pair<string, int> > target_genes;
		target_genes.clear();
		for (int i = 0; i < (int)_targets.size(); i++) {
			if (eje.genes_map.count(_targets[i].trans_id) == 0) {
				continue;
			}
			ej_gene &gene = eje.genes[eje.genes_map[_targets[i].trans_id]];
			if (direction == "b" || (direction == "f" && gene.trans[0].strand != _targets[i].reverse_strand) || (direction == "r" && gene.trans[0].strand == _targets[i].reverse_strand)) {
				target_genes.push_back(pair<string, int>(_targets[i].trans_id, _targets[i].trans_coord));
			}
		}

		if (first_phase) {
			if (target_genes.size() == 1) {
				eje.handler(target_genes[0].first, target_genes[0].second, true);
			} else {
				for (int i = 0; i < (int)target_genes.size(); i++) {
					if (multiple) eje.handler(target_genes[i].first, target_genes[i].second, true);
					else eje.handler(target_genes[i].first, target_genes[i].second, false);
				}
			}
		} else {
		}
	}
}

void aggregate(const string output_exp_filename, const vector<string> input_exp_filenames) {
	if (input_exp_filenames.size() == 0) {
		panic("wrong arguments!");
		return;
	}

	mapped_reads = 0;
	genes.clear();
	vector<int> total_counts;
	vector<map<string, double> > counts;

	for (int i = 0; i < (int)input_exp_filenames.size(); i++) {
		if (!file_exists(input_exp_filenames[i])) {
			cout << "file does not exist: " << input_exp_filenames[i] << endl;
			exit(1);
		}
		ifstream ifs(input_exp_filenames[i].c_str());
		string readline;
		bool first_line = true;
		
		int total_count = 0;

		cout << "reading expression file: " << input_exp_filenames[i] << "...\n";

		int num_genes = 0;
		map<string, double> temp_counts;
		while(getline(ifs, readline)) {
			vector<string> tokens = string_tokenize(readline);
			if (tokens.size() != 9) {
				cout << "bad format.\n";
				exit(1);
			}
			if (first_line) {
				first_line = false;
				if (tokens[0] != "name") {
					cout << "bad format.\n";
					exit(1);
				}
				continue;
			}
			num_genes++;
			string gene_name = tokens[0];
			int ucount = str2int(tokens[1]);
			double count = str2double(tokens[2]);
			int total = str2int(tokens[4]);
			if (total_count == 0) {
				total_count = total;
				mapped_reads += total_count;
				total_counts.push_back(total_count);
			} else if (total != total_count) {
				cout << "wrong total count.\n";
				exit(1);
			}

			int len = str2int(tokens[5]);

			if (genes.find(gene_name) == genes.end()) {
				genes[gene_name].length = len;
				genes[gene_name].total_count = count;
				genes[gene_name].unique_count = ucount;
				vector<double> temp;
				temp.push_back(count);
			} else {
				if (genes[gene_name].length != len) {
					cout << "wrong gene len.\n";
					exit(1);
				}
				genes[gene_name].total_count += count;
				genes[gene_name].unique_count += ucount;
			}
			temp_counts[gene_name] = count;
		}
		counts.push_back(temp_counts);
		ifs.close();
		ifs.clear();
		cout << num_genes << " genes read.\n";
	}

	printf ("writing aggragated files...\n");
	ofstream ofs(output_exp_filename.c_str());
	ofstream ofs1((output_exp_filename+".dev").c_str());
	ofs << "name\tucount\tcount\tratio\ttotal\tlen\texp\tlower\tupper\n";
	ofs1 << "name\tp\tdeviances\n";
	vector<vector<double> > deviances;
	vector<vector<double> > deviances_excluded;
	deviances.resize(total_counts.size());
	deviances_excluded.resize(total_counts.size());
	vector<double> probs;
	vector<vector<double> > normal_fit;
	normal_fit.resize(total_counts.size());
	for (map<string, gene_info_struct>::iterator it = genes.begin(); it != genes.end(); it++) {
		ofs << it->first << "\t" << it->second.unique_count << "\t" << it->second.total_count << "\t" << ((it->second.total_count>0)?(it->second.unique_count/it->second.total_count):0) << "\t" << mapped_reads << "\t" << it->second.length;
		double exp = it->second.total_count / mapped_reads / it->second.length * 1000 * 1000000;
		double p = it->second.total_count / mapped_reads;
		double lower = p - 1.96 * sqrt(p * (1-p) / mapped_reads);
		double upper = p + 1.96 * sqrt(p * (1-p) / mapped_reads);
		ofs << "\t" << exp << "\t" << lower / it->second.length * 1000 * 1000000 << "\t" << upper / it->second.length * 1000 * 1000000 << endl;

		ofs1 << it->first << "\t" << p;
		probs.push_back(p);
		for (int i = 0; i < (int)total_counts.size(); i++) {
			ofs1 << "\t" << (counts[i][it->first] - p * total_counts[i]) / (p * total_counts[i]);
			if (p == 0.0) {
				deviances[i].push_back(0);
			} else {
				deviances[i].push_back((counts[i][it->first] - p * total_counts[i]) / (p * total_counts[i]));
			}
			double p_excluded = (it->second.total_count - counts[i][it->first]) / (mapped_reads - total_counts[i]);
			if (p_excluded == 0.0) {
				deviances_excluded[i].push_back(0);
			} else {
				deviances_excluded[i].push_back((counts[i][it->first] - p_excluded * total_counts[i]) / (p_excluded * total_counts[i]));
			}
			normal_fit[i].push_back((counts[i][it->first] - p * total_counts[i]) / sqrt((1-p) * p * total_counts[i]));
		}
		ofs1 << endl;
	}
	ofs.close();
	ofs.clear();
	ofs1.close();
	ofs1.clear();

	ofstream ofs2((output_exp_filename+".stat").c_str());
	ofs2 << "computation based on all genes:\n";
	ofs2 << "mean of deviances:\tstddev of deviances\tfile\n";
	for (int i = 0; i < (int)deviances.size(); i++) {
		ofs2 << mean(deviances[i]) << "\t" << stddev(deviances[i]) << "\t" << input_exp_filenames[i] << endl;
	}
	ofs2 << endl;

	vector<vector<double> > temp;
	temp.resize(deviances.size());
	for (int i = 0; i < (int)probs.size(); i++) {
		if (probs[i] > 1e-4) {
			for (int j = 0; j < (int)deviances.size(); j++) {
				temp[j].push_back(deviances[j][i]);
			}
		}
	}

	ofs2 << "computation based on genes with abundance > 1e-4, totally " << (int)temp[0].size() << " genes:\n";
	ofs2 << "mean of deviances:\tstddev of deviances\tfile\n";
	for (int i = 0; i < (int)temp.size(); i++) {
		ofs2 << mean(temp[i]) << "\t" << stddev(temp[i]) << "\t" << input_exp_filenames[i] << endl;
	}
	ofs2 << endl;

	temp.clear();
	temp.resize(deviances_excluded.size());
	for (int i = 0; i < (int)probs.size(); i++) {
		if (probs[i] > 1e-4) {
			for (int j = 0; j < (int)deviances_excluded.size(); j++) {
				temp[j].push_back(deviances_excluded[j][i]);
			}
		}
	}

	ofs2 << "computation based on genes with abundance > 1e-4, totally " << (int)temp[0].size() << " genes:\n";
	ofs2 << "mean of excluded deviances:\tstddev of excluded deviances\tfile\n";
	for (int i = 0; i < (int)temp.size(); i++) {
		ofs2 << mean(temp[i]) << "\t" << stddev(temp[i]) << "\t" << input_exp_filenames[i] << endl;
	}
	ofs2 << endl;

	temp.clear();
	temp.resize(normal_fit.size());
	for (int i = 0; i < (int)probs.size(); i++) {
		if (probs[i] > 1e-4) {
			for (int j = 0; j < (int)normal_fit.size(); j++) {
				temp[j].push_back(normal_fit[j][i]);
			}
		}
	}

	ofs2 << "computation based on genes with abundance > 1e-4, totally " << (int)temp[0].size() << " genes:\n";
	ofs2 << "mean of normal fit:\tstddev of normal fit\tfile\n";
	for (int i = 0; i < (int)temp.size(); i++) {
		ofs2 << mean(temp[i]) << "\t" << stddev(temp[i]) << "\t" << input_exp_filenames[i] << endl;
	}
	ofs2 << endl;

	ofs2.close();
	ofs2.clear();
	printf("done.\n");
}

void comp_exp() {
//	if (version == 2) {
	if (version >= 2) {
		printf("computing expressions...\n");
		for (int i = 0; i < (int)eje.genes.size(); i++) {
			eje.compute_expression(i, mapped_reads);
		}
		printf("done.\n");
	}

	printf ("writing outputs...\n");
	if (version == 1) {
		ofstream ofs((local->output_prefix + "." + int2str(num_mismatch) + ".gene.exp").c_str());
		ofs << "name\tucount\tcount\tratio\ttotal\tlen\texp\tlower\tupper\treduced_len\n";
		for (map<string, gene_info_struct>::iterator it = genes.begin(); it != genes.end(); it++) {
			if (it->second.total_count > 0 && it->second.length > 1000000000) {
				printf("WARNING: gene length < %d but has reads. gene: %s\n", read_length + it->second.reduce_len, it->first.c_str());
			}
			ofs << it->first << "\t" << it->second.unique_count << "\t" << it->second.total_count << "\t" << ((it->second.total_count>0)?(it->second.unique_count/it->second.total_count):0) << "\t" << mapped_reads << "\t" << it->second.length;
			double exp = it->second.total_count / mapped_reads / it->second.length * 1000 * 1000000;
			double p = it->second.total_count / mapped_reads;
			double lower = p - 1.96 * sqrt(p * (1-p) / mapped_reads);
			double upper = p + 1.96 * sqrt(p * (1-p) / mapped_reads);
			ofs << "\t" << exp << "\t" << lower / it->second.length * 1000 * 1000000 << "\t" << upper / it->second.length * 1000 * 1000000;
			ofs << "\t" << it->second.reduce_len << endl;
		}
		ofs.close();
		ofs.clear();
/*		ofs.open((output_prefix + "." + int2str(num_mismatch) + ".trans.exp").c_str());
		ofs << "name\tucount\tcount\tratio\ttotal\tlen\texp\tlower\tupper\treduced_len\n";
		for (map<string, gene_info_struct>::iterator it = transcripts.begin(); it != transcripts.end(); it++) {
			ofs << it->first << "\t" << it->second.unique_count << "\t" << it->second.total_count << "\t" << ((it->second.total_count>0)?(it->second.unique_count/it->second.total_count):0) << "\t" << mapped_reads << "\t" << it->second.length;
			double exp = it->second.total_count / mapped_reads / it->second.length * 1000 * 1000000;
			double p = it->second.total_count / mapped_reads;
			double lower = p - 1.96 * sqrt(p * (1-p) / mapped_reads);
			double upper = p + 1.96 * sqrt(p * (1-p) / mapped_reads);
			ofs << "\t" << exp << "\t" << lower / it->second.length * 1000 * 1000000 << "\t" << upper / it->second.length * 1000 * 1000000;
			ofs << "\t" << it->second.reduce_len << endl;
		}
		ofs.close();
		ofs.clear();*/
//	} else if (version == 2) {
	} else if (version >= 2) {
		ofstream ofs((local->output_prefix + "." + int2str(num_mismatch) + "." + direction + ".categories").c_str());
		for (int i = 0; i < (int)eje.genes.size(); i++) {
				ofs << eje.genes[i].name << "\t" << (int)eje.genes[i].trans.size() << "\t" << (int)eje.genes[i].rates.size() << endl;
				for (int k = 0; k < (int)eje.genes[i].trans.size(); k++) {
					ofs << eje.genes[i].trans[k].name << "\t";
				}
				ofs << endl;
				for (int j = 0; j < (int)eje.genes[i].rates.size(); j++) {
					for (int k = 0; k < (int)eje.genes[i].trans.size(); k++) {
						ofs << eje.genes[i].rates[j][k] << "\t";
					}
					ofs << eje.genes[i].counts[j];
					ofs << endl;
				}
		}
		ofs.close();
		ofs.clear();

		ofs.open((local->output_prefix + "." + int2str(num_mismatch) + "." + direction + ".info").c_str());
		for (int i = 0; i < (int)eje.genes.size(); i++) {
			ofs << "gene\t" << eje.genes[i].name;
			for (int j = 0; j < (int)eje.genes[i].subexon_sets.size(); j++) {
				ofs << "\t" << eje.genes[i].len_subexon_sets[j];
			}
			ofs << endl;
			ofs << "u_count\t";
			for (int j = 0; j < (int)eje.genes[i].subexon_sets.size(); j++) {
				ofs << "\t" << eje.genes[i].subexon_set_unique_count[j];
			}
			ofs << endl;
			if (eje.genes[i].subexon_sets.size() > 1) {
				ofs << "u_count_junc\t";
				for (int j = 0; j < (int)eje.genes[i].subexon_sets.size() - 1; j++) {
					for (int k = j + 1; k < (int)eje.genes[i].subexon_sets.size(); k++) {
						ofs << eje.genes[i].subexon_set_junc_unique_count[j][k] << ",";
					}
					ofs << "\t";
				}
				ofs << endl;
			}
			ofs << "exp\t";
			for (int j = 0; j < (int)eje.genes[i].subexon_sets.size(); j++) {
				ofs << "\t" << (double)eje.genes[i].subexon_set_unique_count[j] / mapped_reads / eje.genes[i].len_subexon_sets[j] * 1000 * 1000000;
			}
			ofs << endl;
			for (int k = 0; k < (int)eje.genes[i].trans.size(); k++) {
				ofs << "isoform\t" << eje.genes[i].trans[k].name;
				for (int j = 0; j < (int)eje.genes[i].subexon_sets.size(); j++) {
					ofs << "\t" << eje.genes[i].trans[k].subexon_set_vector[j];
				}
				ofs << endl;
			}
			ofs << "gene_exp\t" << eje.genes[i].gene_expression /*<< "\tgene_var\t" << eje.genes[i].gene_var*/ << endl;
			for (int k = 0; k < (int)eje.genes[i].trans.size(); k++) {
				ofs << "iso_exp\t" << eje.genes[i].tran_expression[k];
				//if (eje.genes[i].trans.size() > 1) ofs << "\tiso_var\t" << eje.genes[i].covar[k][k];
				ofs << endl;
			}
			if (eje.genes[i].trans.size() > 1) {
				/*ofs << "covar:\n";
				for (int k = 0; k < (int)eje.genes[i].trans.size(); k++) {
					for (int l = 0; l < (int)eje.genes[i].trans.size(); l++) {
						ofs << eje.genes[i].covar[k][l] << "\t";
					}
					ofs << endl;
				}*/
	/*			ofs << "obs_covar:\n";
				for (int k = 0; k < (int)eje.genes[i].trans.size(); k++) {
					for (int l = 0; l < (int)eje.genes[i].trans.size(); l++) {
						ofs << eje.genes[i].obs_covar[k][l] << "\t";
					}
					ofs << endl;
				}*/
			}
		}
		ofs.close();
		ofs.clear();

		ofstream ofs1((local->output_prefix + "." + int2str(num_mismatch) + "." + direction + ".exp").c_str());
		for (int i = 0; i < (int)eje.genes.size(); i++) {
			ofs1 << eje.genes[i].name << "\t" << eje.genes[i].gene_expression;
			ofs1 << "\t" << (int)eje.genes[i].trans.size() << "\t";
			for (int k = 0; k < (int)eje.genes[i].trans.size(); k++) {
				if (k > 0) ofs1 << ",";
				ofs1 <<  eje.genes[i].trans[k].name ;
			}
			ofs1 << "\t";
			for (int k = 0; k < (int)eje.genes[i].trans.size(); k++) {
				if (k > 0) ofs1 << ",";
				ofs1 <<  eje.genes[i].tran_expression[k];
			}
			ofs1 << endl;
		}
		ofs1.close();
		ofs1.clear();

		ofstream ofs2((local->output_prefix + "." + int2str(num_mismatch) + "." + direction + ".matlab").c_str());
		ofs2 << mapped_reads << "\t" << read_length << endl;
		for (int i = 0; i < (int)eje.genes.size(); i++) {
			ofs2 << eje.genes[i].name << "\t" << (int)eje.genes[i].subexon_sets.size() << "\t" << (int)eje.genes[i].trans.size() << endl;
			for (int k = 0; k < (int)eje.genes[i].trans.size(); k++) {
				ofs2 << eje.genes[i].trans[k].name << "\t";
			}
			ofs2 << endl;
			for (int j = 0; j < (int)eje.genes[i].subexon_sets.size(); j++) {
				ofs2 << eje.genes[i].len_subexon_sets[j] << "\t";
			}
			ofs2 << endl;
			for (int j = 0; j < (int)eje.genes[i].subexon_sets.size(); j++) {
				ofs2 << eje.genes[i].subexon_set_unique_count[j] << "\t";
			}
			ofs2 << endl;
			if (eje.genes[i].subexon_sets.size() > 1) {
				for (int j = 0; j < (int)eje.genes[i].subexon_sets.size() - 1; j++) {
					for (int k = j + 1; k < (int)eje.genes[i].subexon_sets.size(); k++) {
						ofs2 << eje.genes[i].len_subexon_set_junc[j] + eje.genes[i].len_subexon_set_junc[k] << "\t";
					}
				}
				ofs2 << endl;
			}
			if (eje.genes[i].subexon_sets.size() > 1) {
				for (int j = 0; j < (int)eje.genes[i].subexon_sets.size() - 1; j++) {
					for (int k = j + 1; k < (int)eje.genes[i].subexon_sets.size(); k++) {
						ofs2 << eje.genes[i].subexon_set_junc_unique_count[j][k] << "\t";
					}
				}
				ofs2 << endl;
			}
			for (int k = 0; k < (int)eje.genes[i].trans.size(); k++) {
				for (int j = 0; j < (int)eje.genes[i].subexon_sets.size(); j++) {
					ofs2 << eje.genes[i].trans[k].subexon_set_vector[j] << "\t";
				}
				ofs2 << endl;
			}

		}
		ofs2.close();
		ofs2.clear();
//	} else if (version == 3 || version == 4) {
		if (version == 3 || version == 4) {
			ofstream ofs((local->output_prefix + "." + int2str(num_mismatch) + "." + direction + ".count").c_str());
			ofs << "read_length\t" << read_length << endl;
			ofs <<"mapped_reads\t" << mapped_reads << endl;
			for (int i = 0; i < (int)eje.genes.size(); i++) {
				ofs << "gene\t" << eje.genes[i].name << endl;
				ofs << "#exons\t" << (int)eje.genes[i].subexon_sets.size() << endl;
				ofs << "#isoforms\t" << (int)eje.genes[i].trans.size() << endl;
				ofs << "exon_length\t";
				for (int j = 0; j < (int)eje.genes[i].subexon_sets.size(); j++) {
					ofs << eje.genes[i].len_subexon_sets[j] << "\t";
				}
				ofs << endl;
				ofs << "exon_#reads\t";
				for (int j = 0; j < (int)eje.genes[i].subexon_sets.size(); j++) {
					ofs << eje.genes[i].subexon_set_unique_count[j] << "\t";
				}
				ofs << endl;
				if (version == 4) {
					ofs << "junction_length\t";
					if (eje.genes[i].subexon_sets.size() > 1) {
						for (int j = 0; j < (int)eje.genes[i].subexon_sets.size() - 1; j++) {
							for (int k = j + 1; k < (int)eje.genes[i].subexon_sets.size(); k++) {
								ofs << eje.genes[i].len_subexon_set_junc[j] + eje.genes[i].len_subexon_set_junc[k] << "\t";
							}
						}
						ofs << endl;
					}
					ofs << "junction_#reads\t";
					if (eje.genes[i].subexon_sets.size() > 1) {
						for (int j = 0; j < (int)eje.genes[i].subexon_sets.size() - 1; j++) {
							for (int k = j + 1; k < (int)eje.genes[i].subexon_sets.size(); k++) {
								ofs << eje.genes[i].subexon_set_junc_unique_count[j][k] << "\t";
							}
						}
						ofs << endl;
					}
					ofs << "gene_structure\t";
					for (int j = 0; j < (int)eje.genes[i].subexon_sets.size(); j++) {
						ofs << "exon_" << j + 1 << "\t";
					}
					ofs << endl;
					for (int k = 0; k < (int)eje.genes[i].trans.size(); k++) {
						ofs << eje.genes[i].trans[k].name << "\t";
						for (int j = 0; j < (int)eje.genes[i].subexon_sets.size(); j++) {
							ofs << eje.genes[i].trans[k].subexon_set_vector[j] << "\t";
						}
						ofs << endl;
					}
				}
			}
			ofs.close();
			ofs.clear();
		}
	}

	if (version == 1 && local->other_args.size() > 2) {
		vector<string> gene_exp_files;
		for (int i = 1; i < (int)local->other_args.size(); i++) gene_exp_files.push_back(local->other_args[i] + "." + int2str(num_mismatch) + ".gene.exp");
		aggregate("all.exp", gene_exp_files);
/*		vector<string> trans_exp_files;
		for (int i = 1; i < (int)other_args.size(); i++) trans_exp_files.push_back(other_args[i] + "." + int2str(num_mismatch) + ".trans.exp");
		aggregate("all.exp", trans_exp_files);*/
	}

	printf ("done.\n");
}

void extract_handler(vector<target_struct> &_targets) {
	if (version == 1) {
		for (int i = 0; i < (int)_targets.size(); i++) {
			vector<string> tokens = string_tokenize(_targets[i].trans_id, "/:");
			string _gene_name = tokens[2];
			if (gene_name.empty() || gene_name.count(tolower(_gene_name)) > 0) {
				if (_targets[i].reverse_strand) {
					ofs_F << _gene_name << "\t" << _targets[i].trans_coord << "\t1\n";
				} else {
					ofs_R << _gene_name << "\t" << _targets[i].trans_coord << "\t1\n";
				}
			}
		}
	} else if (version == 0) {
		vector<pair<pair<string, int>, bool> > targets;

		for (int i = 0; i < (int)_targets.size(); i++) {
			targets.push_back(pair<pair<string, int>, bool>(pair<string, int>(_targets[i].trans_id, _targets[i].trans_coord), _targets[i].reverse_strand));
		}

		if (targets.size() == 1) {
			if (targets[0].second) {
				ofs_R << targets[0].first.first << "\t" << targets[0].first.second << "\t1\n";
			} else {
				ofs_F << targets[0].first.first << "\t" << targets[0].first.second << "\t1\n";
			}
		} else if (targets.size() > 1){
			for (int i = 0; i < (int)targets.size(); i++) {
				ofs_M << targets[i].first.first << "\t" << targets[i].first.second << "\t1\n";
			}
		}
	}
}

void convert_handler(vector<target_struct> &_targets) {
	if (version == 0) {
		if (_targets.size() == 1) {
			if (_targets[0].reverse_strand) {
				ofs_R << _targets[0].trans_id << "\t" << _targets[0].trans_coord << "\t1\n";
			} else {
				ofs_F << _targets[0].trans_id << "\t" << _targets[0].trans_coord << "\t1\n";
			}
		} else if (_targets.size() > 1){
			for (int i = 0; i < (int)_targets.size(); i++) {
				ofs_M << _targets[i].trans_id << "\t" << _targets[i].trans_coord << "\t1\n";
			}
		}
	} else if (version == 1) {
		vector<pair<string, pair<int, bool> > > targets;
		for (int i = 0; i < (int)_targets.size(); i++) {
			vector<string> tokens = string_tokenize(_targets[i].trans_id, "$");
			__ASSERT(tokens.size() >= 1, "ERROR: no gene name!");
			string name;
			if (tokens.size() == 1) {
//				panic("no : in target");
//				return; //ribosomal
				name = _targets[i].trans_id;
			} else {
				name = tokens[1];
			}
			string seq_name;
			gene_struct gene;
			local->mygenefile.search_by_name(name, seq_name, gene);
			//__ASSERT(seq_name != "", string("gene not found, gene: " + name));
			if (seq_name != "")
				targets.push_back(pair<string, pair<int, bool> >(seq_name, pair<int, bool>(map_gene_to_chr(_targets[i].trans_coord, gene.exons), _targets[i].reverse_strand)));
			else //ribosomal, etc
				continue;
		}

		sort(targets.begin(), targets.end());
		vector<pair<string, pair<int, bool> > > reduced_targets;
		for (int i = 0; i < (int)targets.size(); i++) {
			if (i < (int)targets.size() - 1 && targets[i] == targets[i+1]) continue;
			reduced_targets.push_back(targets[i]);
		}

		//__ASSERT(reduced_targets.size() > 0, "reduced_targets.size() == 0");
		if (reduced_targets.size() == 0) { //ribosomal, etc
		} else if (reduced_targets.size() == 1) {
			if (reduced_targets[0].second.second) {
				ofs_R << reduced_targets[0].first << "\t" << reduced_targets[0].second.first << "\t1\n";
			} else {
				ofs_F << reduced_targets[0].first << "\t" << reduced_targets[0].second.first << "\t1\n";
			}
		} else {
			for (int i = 0; i < (int)reduced_targets.size(); i++) {
				ofs_M << reduced_targets[i].first << "\t" << reduced_targets[i].second.first << "\t1\n";
			}
		}
	} else if (version >= 2) {
		vector<pair<pair<string, int>, int> > targets;
		for (int i = 0; i < (int)_targets.size(); i++) {
			pair<string, int> chr_pos = eje.map_coord_to_chr(_targets[i].trans_id, _targets[i].trans_coord);
			if (chr_pos.second >= 0) targets.push_back(pair<pair<string, int>, int>(chr_pos, i));
		}

		if (targets.size() == 1) {
			if (_targets[targets[0].second].reverse_strand) {
				ofs_R << targets[0].first.first << "\t" << targets[0].first.second << "\t1\n";
				local->ofs << _targets[targets[0].second].trans_id << "\t" << _targets[targets[0].second].trans_coord << "\t1\n";
			} else {
				ofs_F << targets[0].first.first << "\t" << targets[0].first.second << "\t1\n";
				local->ofs << _targets[targets[0].second].trans_id << "\t" << _targets[targets[0].second].trans_coord << "\t1\n";
			}
		} else if (targets.size() > 1){
			for (int i = 0; i < (int)targets.size(); i++) {
				ofs_M << targets[i].first.first << "\t" << targets[i].first.second << "\t1\n";
			}
		}
	}
}

void exon_usage_handler(vector<target_struct> &targets) {
//	ofs << targets[0].trans_id << "\t" << targets[0].trans_coord << "\t" << (targets[0].reverse_strand?"-":"+") << "\t";
	gene_struct *gene;
	int exon;
	local->mygenefile.search_interval(targets[0].trans_id, targets[0].trans_coord, targets[0].trans_coord, gene, exon);
	if (gene == NULL) {
//		ofs << "NA\tNA\tintergene";
		intergene_reads++;
	} else {
//		ofs << gene->geneName << "\t" << gene->name << "\t";
		if (exon == -1) {
//			ofs << "intron";
			intron_reads++;
		} else if ((targets[0].trans_coord < gene->cdsStart && gene->strand) || (targets[0].trans_coord > gene->cdsEnd && !gene->strand)) {
//			ofs << "5UTR";
			UTR5_reads++;
		} else if ((targets[0].trans_coord < gene->cdsStart && !gene->strand) || (targets[0].trans_coord > gene->cdsEnd && gene->strand)) {
//			ofs << "3UTR";
			UTR3_reads++;
		} else {
//			ofs << "exon\t" << (gene->strand?(exon+1):((int)gene->exons.size()-exon));
			exon_reads++;
		}
	}
//	ofs << endl;
	if (exon != -1) {
		vector<pair<int, pair<int, int> > > exons;
		local->mygenefile.search_interval(targets[0].trans_id, targets[0].trans_coord, targets[0].trans_coord, exons);
		for (int i = 0; i < (int)exons.size(); i++) {
			if (local->exon_counts.count(exons[i]) == 0) {
				local->exon_counts[exons[i]] = 1;
			} else {
				local->exon_counts[exons[i]]++;
			}
		}
	}
}

//compute the log-likelihood for binomial r.v. Y with parameters N and p
double loglikelihood_binom(int Y, int N, double p) {
	return (Y * log(p) + (N-Y) * log(1-p));
}

//perform SPRT for binomial samples with H0:p=theta-dtheta and H1:p=theta+dtheta
//compute test thresholds A and B using Wald's approximation given type-I and type-II error levels alpha and beta
//return 1 if H0 is accepted, 2 if H1 is accepted, 0 if uncertain
int test_sequential_binom_SPRT(int Y, int N, double theta = 0.5, double dtheta = 0.01, double alpha = 0.001, double beta = 0.001) {
	double A = log(beta/(1-alpha));
	double B = log((1-beta)/alpha);
	double theta1 = theta - dtheta;
	double theta2 = theta + dtheta;
  
	double S = loglikelihood_binom(Y, N, theta1) - loglikelihood_binom(Y, N, theta2);
	if (S > B) {
		return 1; // p < theta1
	} else if (S < A) {
		return 2; // p > theta2
	} else {
		return 0; // not sure
	}
}

//perform GLRT for binomial samples with H0:p<theta and H1:p>theta
//compute test threshold A using Wald's approximation if given type-I and type-II error levels alpha
//return 1 if H0 is accepted, 2 if H1 is accepted, 0 if uncertain
int test_sequential_binom_GLRT(int Y, int N, double theta = 0.5, double A = -1.0, double alpha = 0.001) {
	if (A < 0) {
		A = log((1-alpha)/alpha);
	}

	double theta_hat = (Y+0.01)/(N+0.02);
	double S = loglikelihood_binom(Y, N, theta_hat) - loglikelihood_binom(Y, N, theta);
	if (S > A && theta_hat < theta) {
		return 1; // p < theta
	} else if (S > A && theta_hat > theta) {
		return 2; // p > theta
	} else {
		return 0; // not sure
	}
}

void two_sample_diff_gene_exp_handler(vector<target_struct> &targets) {
	string gene_name = targets[0].trans_id;
	if (local->gene_counts.count(gene_name) == 0) {
		local->gene_counts[gene_name] = 1;
	} else {
		local->gene_counts[gene_name]++;
	}
	if (active_samples == 1) return;
	if (diff_gene_rejected.count(gene_name) == 0) {
		int count0 = samples[0]->gene_counts[gene_name];
		int count1 = samples[1]->gene_counts[gene_name];
		int min_count = min(count0, count1);
		int test_result = test_sequential_binom_GLRT(min_count, count0 + count1, 0.4);
		if (test_result > 0) {
			diff_gene_rejected[gene_name] = (test_result == 1);
			cout << gene_name << "\t" << count0 << "\t" << count1 << "\t" << test_result << endl;
		}
	}
}

bool read_exon_file() {
	table mytable;
	if (!mytable.read_from_file(exon_file)) {
		printf("error loading exon file.\n");
		return false;
	}
	if (mytable.data.size() == 0) {
		printf("no data in exon file.\n");
		return false;
	}
	if (mytable.data[0].size() != 3) {
		printf("wrong format in exon file.\n");
		return false;
	}
	for (int i = 0; i < (int)mytable.data.size(); i++) {
		string gene_name = mytable.data[i][0];
		if (!is_int(mytable.data[i][1]) || !is_int(mytable.data[i][2])) {
			printf("wrong format in exon file.\n");
			return false;
		}
		int start = str2int(mytable.data[i][1]);
		int end = str2int(mytable.data[i][2]);
		if (start > end || start < 0) {
			printf("wrong exon number.\n");
			return false;
		}
		string chr;
		gene_struct gene;
		local->mygenefile.search_by_name(gene_name, chr, gene);
		if (chr == "") {
			printf("gene not found.\n");
			return false;
		}
		if ((int)gene.exons.size() <= end) {
			printf("wrong exon number.\n");
			return false;
		}
		if (!gene.strand) {
			int temp = start;
			start = gene.exons.size() - end - 1;
			end = gene.exons.size() - temp - 1;
		}
//		printf("%s\t%d\t%d\n", gene_name.c_str(), gene.exons[start].first, gene.exons[end].second);
		if (local->exon_usage_sets.count(gene_name) == 0) {
			local->exon_usage_sets[gene_name] = vector<pair<set<int>, set<int> > >();
			local->exon_usage_set_lengths[gene_name] = vector<pair<int, int> >();
			local->exon_usage_set_counts[gene_name] = vector<pair<int, int> >();
		}
		set<int> exon_set1;
		set<int> exon_set2;
		int exon_set_length1 = 0;
		int exon_set_length2 = 0;
		for (int i = 0; i < (int)gene.exons.size(); i++) {
			if (start <= i && i <= end) {
				exon_set1.insert(i);
				exon_set_length1 += gene.exons[i].second - gene.exons[i].first;
			} else {
				exon_set2.insert(i);
				exon_set_length2 += gene.exons[i].second - gene.exons[i].first;
		}
		}
		local->exon_usage_sets[gene_name].push_back(pair<set<int>, set<int> >(exon_set1, exon_set2));
		local->exon_usage_set_lengths[gene_name].push_back(pair<int, int>(exon_set_length1, exon_set_length2));
		local->exon_usage_set_counts[gene_name].push_back(pair<int, int>(0, 0));
	}
	return true;
}

//perform GLRT for two binomial samples with H0:p0/p1<theta and H1:p0/p1>theta, with theta < 1.0
//compute test threshold A using Wald's approximation if given type-I and type-II error levels alpha
//return 1 if H0 is accepted, 2 if H1 is accepted, 0 if uncertain
int test_sequential_diff_binom_GLRT(int Y0, int N0, int Y1, int N1, double theta = 1.0, double A = -1.0, double alpha = 0.001) {
	if (A < 0) {
		A = log((1-alpha)/alpha);
	}

	double p0_hat = (Y0+0.01)/(N0+0.02);
	double p1_hat = (Y1+0.01)/(N1+0.02);
	double a = theta * (N0 + N1 + 0.04);
	double b = -(theta * (N0 + Y1 + 0.03) + N1 + Y0 + 0.03);
	double c = Y0 + Y1 + 0.02;
	double p1_hat1 = (-b - sqrt(b*b-4*a*c))/2/a;
	double p1_hat2 = (-b + sqrt(b*b-4*a*c))/2/a;
	__ASSERT(0 < p1_hat1 && p1_hat1 < 1 && 1 < p1_hat2, "internal error in test_sequential_diff_binom_GLRT.\n");
	double p1_hat_new = p1_hat1;
	double p0_hat_new = p1_hat_new * theta;
	double S = loglikelihood_binom(Y0, N0, p0_hat) + loglikelihood_binom(Y1, N1, p1_hat) - loglikelihood_binom(Y0, N0, p0_hat_new) - loglikelihood_binom(Y1, N1, p1_hat_new);
	if (S > A && p0_hat/p1_hat < theta) {
		return 1; // p0/p1 < theta
	} else if (S > A && p0_hat/p1_hat > theta) {
		return 2; // p0/p1 > theta
	} else {
		return 0; // not sure
	}
}

void two_sample_diff_exon_usage_summary(string gene = "", int index = -1) {
	if (!start_output) {
		start_output = true;
		printf("gene\texon_set0\tes0_len\texon_set1\tes1_len\tgene_len\tsample0_count\tes_count00\trpkm00\tes_count01\trpkm01\tes_count0\trpkm0\tratio0\tsample1_count\tes_count10\trpkm10\tes_count11\trpkm11\tes_count1\trpkm1\tratio1\ttest_result\n");
	}
	for (auto mi = local->exon_usage_sets.begin(); mi != local->exon_usage_sets.end(); mi++) {
		if (gene != "" && mi->first != gene) continue;
		string gene_name = mi->first;
		auto &exons = local->exon_usage_sets[gene_name];
		auto &exon_lengths = local->exon_usage_set_lengths[gene_name];
		if (diff_exon_usage_results.count(gene_name) == 0) {
			vector<int> results;
			for (int j = 0; j < (int)exons.size(); j++) {
				results.push_back(0);
			}
			diff_exon_usage_results[gene_name] = results;
		}
		for (int i = 0; i < (int)exons.size(); i++) {
			if (index != -1 && i != index) continue;
			cout << gene_name << "\t";
			for (auto mi = exons[i].first.begin(); mi != exons[i].first.end(); mi++) {
				cout << *mi << ",";
			}
			printf("\t%d\t", exon_lengths[i].first);
			for (auto mi = exons[i].second.begin(); mi != exons[i].second.end(); mi++) {
				cout << *mi << ",";
			}
			printf("\t%d", exon_lengths[i].second);
			printf("\t%d", exon_lengths[i].first + exon_lengths[i].second);
			pair<int, int> counts0 = samples[0]->exon_usage_set_counts[gene_name][i];
			pair<int, int> counts1 = samples[1]->exon_usage_set_counts[gene_name][i];
			int sum0 = counts0.first + counts0.second;
			int sum1 = counts1.first + counts1.second;
			int total_counts0 = samples[0]->total_read_counts;
			int total_counts1 = samples[0]->total_read_counts;
			double ratio0 = (sum0==0)?0.5:((double)counts0.first/sum0);
			double ratio1 = (sum0==0)?0.5:((double)counts1.first/sum1);
			double rpkm00 = (double)counts0.first/total_counts0/exon_lengths[i].first*1000*1000000;
			double rpkm01 = (double)counts0.second/total_counts0/exon_lengths[i].second*1000*1000000;
			double rpkm10 = (double)counts1.first/total_counts1/exon_lengths[i].first*1000*1000000;
			double rpkm11 = (double)counts1.second/total_counts1/exon_lengths[i].second*1000*1000000;
			double rpkm0 = (double)sum0/total_counts0/(exon_lengths[i].first + exon_lengths[i].second)*1000*1000000;
			double rpkm1 = (double)sum1/total_counts1/(exon_lengths[i].first + exon_lengths[i].second)*1000*1000000;
			printf("\t%d\t%d\t%lf\t%d\t%lf\t%d\t%lf\t%lf", total_counts0, counts0.first, rpkm00, counts0.second, rpkm01, sum0, rpkm0, ratio0);
			printf("\t%d\t%d\t%lf\t%d\t%lf\t%d\t%lf\t%lf", total_counts1, counts1.first, rpkm10, counts1.second, rpkm11, sum1, rpkm1, ratio1);
			printf("\t%d\n", diff_exon_usage_results[gene_name][i]);
		}
	}
}

void two_sample_diff_exon_usage_handler(vector<target_struct> &targets) {
	if (active_samples == 1) return;
	string gene_name = targets[0].trans_id;
	size_t pos = gene_name.find_first_of("$$");
	if (pos != string::npos) {
		gene_name = gene_name.substr(pos+2);
	}
	int coord = targets[0].trans_coord;
	string chr;
	gene_struct gene;
	local->mygenefile.search_by_name(gene_name, chr, gene);
	__ASSERT(chr != "", "gene not found.\n");
	int exon = map_gene_to_chr(coord, gene.exons, gene.strand, true);
	__ASSERT(exon >= 0 && exon < (int)gene.exons.size(), "error locating exons");
	if (local->exon_usage_sets.count(gene_name) == 0) return;
	__ASSERT(local->exon_usage_set_lengths.count(gene_name) == 1 && local->exon_usage_set_counts.count(gene_name) == 1, "internal error, inconsiscent exon_usage data.\n");
	auto &exons = local->exon_usage_sets[gene_name];
	auto &exon_lengths = local->exon_usage_set_lengths[gene_name];
	auto &exon_counts = local->exon_usage_set_counts[gene_name];
	for (int i = 0; i < (int)exons.size(); i++) {
		if (exons[i].first.count(exon)) {
			exon_counts[i].first++;
		}
		if (exons[i].second.count(exon)) {
			exon_counts[i].second++;
		}
	}
	for (int i = 0; i < (int)exons.size(); i++) {
		if (exons[i].first.count(exon) || exons[i].second.count(exon)) {
			if (diff_exon_usage_results.count(gene_name) > 0 && diff_exon_usage_results[gene_name][i] > 0) continue;
			pair<int, int> counts0 = samples[0]->exon_usage_set_counts[gene_name][i];
			pair<int, int> counts1 = samples[1]->exon_usage_set_counts[gene_name][i];
			if (counts0.first*counts1.second  > counts0.second*counts1.first) {
				pair<int, int> temp = counts0;
				counts0 = counts1;
				counts1 = temp;
			}
			int test_result = test_sequential_diff_binom_GLRT(counts0.first, counts0.first+counts0.second, counts1.first, counts1.first+counts1.second, 0.8);
			if (test_result > 0) {
				if (diff_exon_usage_results.count(gene_name) == 0) {
					vector<int> results;
					for (int j = 0; j < (int)exons.size(); j++) {
						results.push_back(0);
					}
					diff_exon_usage_results[gene_name] = results;
				}
				diff_exon_usage_results[gene_name][i] = test_result;
			} else {
				int test_result1 = test_sequential_binom_GLRT(counts0.first + counts0.second, local->total_read_counts, (double)(exon_lengths[i].first + exon_lengths[i].second)/1000/1000000);
				int test_result2 = test_sequential_binom_GLRT(counts1.first + counts1.second, local->total_read_counts, (double)(exon_lengths[i].first + exon_lengths[i].second)/1000/1000000);
				if (test_result1 == 1 && test_result2 == 1) {
					if (diff_exon_usage_results.count(gene_name) == 0) {
						vector<int> results;
						for (int j = 0; j < (int)exons.size(); j++) {
							results.push_back(0);
						}
						diff_exon_usage_results[gene_name] = results;
					}
					diff_exon_usage_results[gene_name][i] = 3;
				}
			}
			if (diff_exon_usage_results.count(gene_name) > 0 && diff_exon_usage_results[gene_name][i] > 0) {
				two_sample_diff_exon_usage_summary(gene_name, i);
			}
		}
	}
}

void map_stat_handler(map_result_struct &result, vector<target_struct> &targets) {
	int read_length = (int)result.probe_seq.length();
	real_read_length = max(real_read_length, read_length);
	__ASSERT(real_read_length <= MAX_READ_LEN, "real_read_length > MAX_READ_LEN\n");
	string seq1 = local->myfasta.get_sequence(targets[0].trans_id, targets[0].trans_coord, targets[0].trans_coord + read_length);
	if (seq1 == "") return;

	__ASSERT((int)result.probe_seq.length() >= read_length, "result.probe_seq.length() < read_length");

	for (int i = 0; i < read_length; i++) {
		int ch1 = nt2int(result.probe_seq[i]);
		if (ch1 < 0) ch1 = 4;
		int pos = (targets[0].reverse_strand) ? read_length - 1 - i : i;
		int ch2 = (pos < (int)seq1.length()) ? nt2int(seq1[pos]) : -1;
		if (targets[0].reverse_strand) ch2 = 3 - ch2;
		if (ch2 < 0) ch2 = 4;
		map_stat_counts[i][ch1][ch2]++;
	}
}

void map_stat() {
	cout << "output mapping statistics...\n";
	ofstream stat_file((local->output_prefix + ".stat").c_str());

	stat_file << "Nucleotides Sequencing Transition Count Matrix:\n";
	int i, j, k;
	for (i = 0; i < real_read_length; i++) {
		stat_file << "\nNucleotide At Position " << i+1 << ":\n";
		stat_file << "From\\To\tA\tC\tG\tT\tN\n";
		for (j = 0; j < 5; j++) {
			stat_file << ((j < 4) ? DNA_C[j] : 'N');
			for (k = 0; k < 5; k++) {
				stat_file << "\t" << map_stat_counts[i][j][k];
			}
			stat_file << endl;
		}
	}

	stat_file << "\nNucleotides Sequencing Transition Ratio Matrix:\n";
	for (i = 0; i < real_read_length; i++) {
		stat_file << "\nNucleotide At Position " << i+1 << ":\n";
		stat_file << "From\\To\tA\tC\tG\tT\tN\n";
		for (j = 0; j < 5; j++) {
			stat_file << ((j < 4) ? DNA_C[j] : 'N');
			int total = map_stat_counts[i][j][0] + map_stat_counts[i][j][1] + map_stat_counts[i][j][2] + map_stat_counts[i][j][3] +  + map_stat_counts[i][j][4];
			if (total == 0) total = 1;
			for (k = 0; k < 5; k++) {
				stat_file << "\t" << (double)map_stat_counts[i][j][k] / total;
			}
			stat_file << endl;
		}
	}

	stat_file << "\nPos\tError_Rate\n";
	int total_nt = 0, error_nt = 0;
	for (i = 0; i < real_read_length; i++) {
		int total = 0;
		int error_total = 0;
		for (j = 0; j < 5; j++) {
			for (k = 0; k < 5; k++) {
				total_nt += map_stat_counts[i][j][k];
				total += map_stat_counts[i][j][k];
				if (j != k) { 
					error_nt += map_stat_counts[i][j][k];
					error_total += map_stat_counts[i][j][k];
				}
			}
		}
		stat_file << i+1 << "\t" << (double)100.0*error_total/total << "%" << endl;
	}
	stat_file << "\ntotal\terror\n";
	stat_file << total_nt << "\t" << error_nt << "\t" << endl;

	stat_file.close();
	stat_file.clear();
	cout << "done.\n";
}

void paired_end_handler(vector<vector<target_struct> > &targets, vector<map_result_struct> &results) {
	if (targets[0][0].trans_id != targets[1][0].trans_id) {
		diff_chr_reads++;
		return;
	}
	if (targets[0][0].reverse_strand == targets[1][0].reverse_strand) {
		same_dir_reads++;
		return;
	}
	bool reverse_strand = false;
	if (targets[0][0].reverse_strand) {
		reverse_strand = true;
		target_struct temp = targets[0][0];
		targets[0][0] = targets[1][0];
		targets[1][0] = temp;
	}

	int start[2] = {-1, -1}, end[2] = {-1, -1};

	ej_gene *p_gene = NULL;
	int subexon_set_id[2][2][2] = {{{-1, -1}, {-1, -1}}, {{-1, -1}, {-1, -1}}};
	int chr_pos_start[2] = {-1, -1}, chr_pos_end[2] = {-1, -1};
	vector<int> tran_pos_start[2], tran_pos_end[2];
	vector<int> insert_lens;

	if (version <= 1) {
		for (int k = 0; k < 2; k++) {
			start[k] = targets[k][0].trans_coord;
			end[k] = targets[k][0].trans_coord + read_length - 1;
		}
	} else { // version >= 2
		if (eje.genes_map.count(targets[0][0].trans_id) > 0) {
			for (int k = 0; k < 2; k++) {
				//__ASSERT(eje.locate_coord(targets[k][0].trans_id, targets[k][0].trans_coord - 1, p_gene, subexon_set_id[k][0][0], subexon_set_id[k][0][1], start[k], chr_pos_start[k], tran_pos_start[k]), "internal error: failed locate coordinate");
				if (!eje.locate_coord(targets[k][0].trans_id, targets[k][0].trans_coord - 1, p_gene, subexon_set_id[k][0][0], subexon_set_id[k][0][1], start[k], chr_pos_start[k], tran_pos_start[k])) {
					printf("warning: failed locate coordinate\n");
					return;
				}
				//__ASSERT(eje.locate_coord(targets[k][0].trans_id, targets[k][0].trans_coord + read_length - 4, p_gene, subexon_set_id[k][1][0], subexon_set_id[k][1][1], end[k], chr_pos_end[k], tran_pos_end[k]), "internal error: failed locate coordinate");
				if (!eje.locate_coord(targets[k][0].trans_id, targets[k][0].trans_coord + read_length - 4, p_gene, subexon_set_id[k][1][0], subexon_set_id[k][1][1], end[k], chr_pos_end[k], tran_pos_end[k])) {
					printf("warning: failed locate coordinate\n");
					return;
				}
				if (!(subexon_set_id[k][0][0] == subexon_set_id[k][1][0] && subexon_set_id[k][0][1] == subexon_set_id[k][1][1] && start[k] + read_length - 3 <= end[k] && 0 <= chr_pos_start[k] && chr_pos_start[k] + read_length - 3 <= chr_pos_end[k])) {
					printf("warning: failed locate coordinate\n");
					return;
				}
				if (!(tran_pos_start[k].size() == tran_pos_end[k].size())) {
					printf("warning: failed locate coordinate\n");
					return;
				}
				for (int l = 0; l < (int)tran_pos_start[k].size(); l++) {
					if (!((tran_pos_start[k][l] == -1 && tran_pos_end[k][l] == -1) || (tran_pos_start[k][l] >= 0 && tran_pos_start[k][l] + read_length - 3 == tran_pos_end[k][l]))) {
						printf("warning: failed locate coordinate\n");
						return;
					}
				}
//					ejes[k].handler(targets[k][0].trans_id, targets[k][0].trans_coord - 1, true);
			}
			for (int k = 0; k < (int)tran_pos_start[0].size(); k++) {
				if (!(tran_pos_start[0].size() == tran_pos_start[1].size())) {
					printf("warning: failed locate coordinate\n");
					return;
				}
				if (tran_pos_start[0][k] >= 0 && tran_pos_start[1][k] >= 0) {
					insert_lens.push_back(tran_pos_start[1][k] - tran_pos_start[0][k] + read_length);
				}
			}
		} else { // ribosomal, etc
			for (int k = 0; k < 2; k++) {
				start[k] = targets[k][0].trans_coord;
				end[k] = targets[k][0].trans_coord + read_length - 1;
			}
			insert_lens.push_back(start[1] - start[0] + read_length);
		}
	}

	if (insert_lens.size() > 0) {
		sort(insert_lens.begin(), insert_lens.end());
		int insert_len = insert_lens[0];

		if (insert_len == insert_lens[insert_lens.size() - 1]) {
			if (local->insert_len_hist.count(insert_len) == 0) {
				local->insert_len_hist[insert_len] = 1;
			} else {
				local->insert_len_hist[insert_len]++;
			}
		}

		if (insert_len <= read_length) {
			neg_dist_reads++;
		} else if (insert_len > MAX_INS_LEN) {
			large_dist_reads++;
		} else {
			local->ofs << targets[0][0].trans_id << "\t" << start[0] << "\t" << start[1] - start[0] + read_length << endl;

			if (version <= 1) {
				paired_reads++;
			} else { // version >= 2
				if (eje.genes_map.count(targets[0][0].trans_id) > 0) {
					paired_reads++;
//					p_gene->iem->accumulate_read(subexon_set_id11, subexon_set_id12);
//					p_gene->iem->accumulate_read(subexon_set_id[0][0][0], subexon_set_id[0][0][1], subexon_set_id[1][0][0], subexon_set_id[1][0][1], insert_len);
					p_gene->iem->accumulate_read(tran_pos_start[0], tran_pos_start[1]);

					chr_pos_end[0] += 3;
					chr_pos_end[1] += 3;
					ofs_F << p_gene->trans[0].chr << "\t" << chr_pos_start[0] << "\t" << chr_pos_end[1] << "\t \t0\t" << (reverse_strand?"-":"+") << "\t-2\t0\t";
					ofs_F << "0\t2\t" << chr_pos_end[0] - chr_pos_start[0] << "," << chr_pos_end[1] - chr_pos_start[1] << "\t0," << chr_pos_start[1] - chr_pos_start[0] << endl;
					for (int i = 0; i < 2; i++) {
						string temp = "";
						for (int j = 0; j <= num_mismatch; j++) {
							temp += (results[i].targets[0].num_mismatch==j?"1":"0");
							if (j < num_mismatch ) temp += ":";
						}
						ofs_R << results[i].probe_id+"_"+int2str(i+1) << "\t" << results[i].probe_seq << "\t" << temp << "\t" << results[i].targets[0].trans_id + ":" + int2str(results[i].targets[0].trans_coord) + (results[i].targets[0].reverse_strand?"R":"F") + int2str(results[i].targets[0].num_mismatch) << endl;
					}
				} else {
					other_reads++;
				}
			}
		}
	}
}

void expression_analysis_handler(vector<vector<target_struct> > &_targets, vector<map_result_struct> &results) {
	__ASSERT(_targets.size() == 1 || _targets.size() == 2, "error: inconsistent read.\n");
	vector<paired_target_struct> concordant_targets;
	concordant_targets.clear();
	int *read_counter = &unmapped_reads;
	for (int i = 0; i < (int)_targets[0].size(); i++) {
		paired_target_struct temp;
		size_t pos = _targets[0][i].trans_id.find_first_of("$$");
		temp.gene_name = (pos == string::npos) ? _targets[0][i].trans_id : _targets[0][i].trans_id.substr(0,pos);
		string isoform_name = (pos == string::npos) ? _targets[0][i].trans_id : _targets[0][i].trans_id.substr(pos+2);
		if (local->iem2.gene_map.count(temp.gene_name) == 0) {
			read_counter = &other_reads;
			continue; //Ribosomal, etc
		}
		iso_exp_model_2::_gene_str &gene = local->iem2.genes[local->iem2.gene_map[temp.gene_name]];
		__ASSERT (gene.isoform_map.count(isoform_name) == 1, string("error: unknown isoform ") + isoform_name);
		temp.isoform_index = gene.isoform_map[isoform_name];
		temp.isoform_name = isoform_name;
		iso_exp_model_2::_isoform_str &isoform = gene.isoforms[temp.isoform_index];
		temp.pos1 = _targets[0][i].trans_coord;
		temp.read_len1 = read_length;
		temp.read_id = results[0].probe_id;
		temp.reverse_strand = _targets[0][i].reverse_strand;
		if (results[0].probe_seq != "") temp.read_len1 = (int)results[0].probe_seq.length();
		__ASSERT(in_range(temp.pos1, 0, isoform.length), string("error: mapped location out of range for isoform ") + isoform_name);
		if (_targets.size() == 2) {
			for (int j = 0; j < (int)_targets[1].size(); j++) {
				if (_targets[1][j].trans_id != _targets[0][i].trans_id) {
					read_counter = &diff_chr_reads;
					continue;
				}
				if (_targets[1][j].reverse_strand == _targets[0][i].reverse_strand) {
					read_counter = &same_dir_reads;
					continue;
				}
				temp.pos2 = _targets[1][j].trans_coord;
				temp.read_len2 = read_length;
				if (results[1].probe_seq != "") temp.read_len2 = (int)results[1].probe_seq.length();
				__ASSERT(in_range(temp.pos2, 0, isoform.length), string("error: mapped location out of range for isoform ") + isoform_name);
				if (_targets[0][i].reverse_strand) {
					int temp_pos = temp.pos1;
					temp.pos1 = temp.pos2;
					temp.pos2 = temp_pos;
					int temp_len = temp.read_len1;
					temp.read_len1 = temp.read_len2;
					temp.read_len2 = temp_len;
				}
				if (temp.pos2 < temp.pos1) {
					read_counter = &neg_dist_reads;
					continue;
				}
				if (temp.pos2 - temp.pos1 > MAX_INS_LEN) {
					read_counter = &large_dist_reads;
					continue;
				}
				temp.insert_len = temp.pos2 - temp.pos1 + temp.read_len2;
				concordant_targets.push_back(temp);
			}
		} else {
			temp.insert_len = temp.read_len1;
			temp.pos2 = -1;
			temp.read_len2 = 0;
			concordant_targets.push_back(temp);
		}
	}
	set<string> target_genes;
	target_genes.clear();
	set<int> ins_lens;
	for (int i = 0; i < (int)concordant_targets.size(); i++) {
		target_genes.insert(concordant_targets[i].gene_name);
		ins_lens.insert(concordant_targets[i].insert_len);
	}
	if (first_phase) {
		if (concordant_targets.size() == 0) {
			(*read_counter)++;
			return;
		}
		paired_reads++;
		total_targets += (int)concordant_targets.size();
		__ASSERT(ins_lens.size() > 0, "internal error: no insert length.\n");
		if (ins_lens.size() == 1) {
			int insert_len = *ins_lens.begin();
			if (local->insert_len_hist.count(insert_len) == 0) {
				local->insert_len_hist[insert_len] = 1;
			} else {
				local->insert_len_hist[insert_len]++;
			}
		}
		if (target_genes.size() > 1) return;
		string gene_name =  *(target_genes.begin());
		unique_reads++;
		vector<int> isoform_indices;
		vector<int> start_pos;
		vector<int> insert_lens;
		for (int i = 0; i < (int)concordant_targets.size(); i++) {
			isoform_indices.push_back(concordant_targets[i].isoform_index);
			start_pos.push_back(concordant_targets[i].pos1);
			insert_lens.push_back(concordant_targets[i].insert_len);
		}
		__ASSERT (local->iem2.gene_map.count(gene_name) > 0, string("internal error: gene not found ") + gene_name);
		iso_exp_model_2::_gene_str &gene = local->iem2.genes[local->iem2.gene_map[gene_name]];
		gene.accumulate_read(isoform_indices, start_pos, insert_lens);

		if (annotation_file != "") {
			string chr;
			gene_struct gene_s;
			local->mygenefile.search_by_name(concordant_targets[0].isoform_name, chr, gene_s);
			if (chr != "") {
				int pos1 = map_gene_to_chr(concordant_targets[0].pos1, gene_s.exons, gene_s.strand);
				int pos2 = map_gene_to_chr(concordant_targets[0].pos1 + concordant_targets[0].read_len1 - 1, gene_s.exons, gene_s.strand);
				interval_set read_interval_set(min(pos1, pos2), max(pos1, pos2) + 1);
				if (concordant_targets[0].pos2 >= 0) {
					int pos3 = map_gene_to_chr(concordant_targets[0].pos2, gene_s.exons, gene_s.strand);
					int pos4 = map_gene_to_chr(concordant_targets[0].pos2 + concordant_targets[0].read_len2 - 1, gene_s.exons, gene_s.strand);
					read_interval_set.union_with(interval_set(min(pos3, pos4), max(pos3, pos4) + 1));
				}
				if (!gene_s.strand) {
						concordant_targets[0].reverse_strand = !concordant_targets[0].reverse_strand;
				}
				read_interval_set.intersect_with(interval_set(gene_s.exons));
				read_interval_set.remove_empty_intervals();
				read_interval_set.remove_redundant_inner_points();
				vector<pair<int, int> > int_pairs;
				__ASSERT(read_interval_set.convert_to_int_pairs(int_pairs) && int_pairs.size() > 0, "internal error: convert to int pairs failed.\n");
				local->ofs << chr << "\t" << int_pairs[0].first << "\t" << int_pairs[int_pairs.size()-1].second << "\t" << concordant_targets[0].read_id << "\t0\t" << (concordant_targets[0].reverse_strand?"-":"+");
				local->ofs << "\t" << int_pairs[0].first << "\t" << int_pairs[int_pairs.size()-1].second << "\t0\t" << int_pairs.size() << "\t";
				for (int i = 0; i < (int)int_pairs.size(); i++) {
					local->ofs << int_pairs[i].second - int_pairs[i].first;
					if (i < (int)int_pairs.size() - 1) local->ofs << ",";
					else local->ofs << "\t";
				}
				for (int i = 0; i < (int)int_pairs.size(); i++) {
					local->ofs << int_pairs[i].first - int_pairs[0].first;
					if (i < (int)int_pairs.size() - 1) local->ofs << ",";
					else local->ofs << endl;
				}
			} else ;// ribosomal, etc, gene not found
		}
	} else { // second_phase, finish later
		if (target_genes.size() == 1) return;
	}
}

/*vector<double> q_table;

double filtering_q(int len) {
	if (len <= read_length || len >= (int)q_table.size()) return 0;
	else return q_table[len];
}

double lambda = 100;*/

//for paired_end or expression_analysis
void paired_end() {
	printf ("writing fragment length...\n");
	ofstream ofs((local->output_prefix + ".ins_len").c_str());
	int min_ins_len = 100000000, max_ins_len = -100000000, ave_ins_len = 0, ave_ins_len_num = 0, mode_ins_len = 0, mode_ins_len_num = -1;
	double total_ins_len = 0.0;
	for (map<int, int>::iterator mi = local->insert_len_hist.begin(); mi != local->insert_len_hist.end(); mi++) {
		ofs << mi->first << "\t" << mi->second << endl;
		min_ins_len = min(mi->first, min_ins_len);
		max_ins_len = max(mi->first, max_ins_len);
		total_ins_len += (double)mi->first * mi->second;
		ave_ins_len_num += mi->second;
		if (mi->second > mode_ins_len_num) {
			mode_ins_len_num = mi->second;
			mode_ins_len = mi->first;
		}
	}
	ave_ins_len = round_double(total_ins_len / ((double)ave_ins_len_num + 1e-6));
	ofs.close();
	ofs.clear();
	printf ("done.\n");

	printf("fragment length:\n");
	printf("   minimum: %d\n", min_ins_len);
	printf("   maximum: %d\n", max_ins_len);
	printf("   mean:    %d\n", ave_ins_len);
	printf("   mode:    %d\n", mode_ins_len);
	if (version >= 2 || task == "expression_analysis") {

		__ASSERT(mode_ins_len >= read_length, "internal error: mode_ins_len < read_length");

		printf("computing and writing insert length probabilities...\n");
		local->p_table.clear();
		int sum_reads = 0;
		for (int i = read_length; i <= max_ins_len; i++) {
			if (local->insert_len_hist[i] > mode_ins_len_num * .01) {
				local->p_table[i] = local->insert_len_hist[i];
				sum_reads += local->insert_len_hist[i];
			}
		}
		ofs.open((local->output_prefix + ".ins_len_p").c_str());
		for (map<int, double>::iterator mi = local->p_table.begin(); mi != local->p_table.end(); mi++) {
			mi->second /= sum_reads;
			ofs << mi->first << "\t" << mi->second << endl;
		}
		ofs.close();
		ofs.clear();
		printf ("done.\n");

		if (ins_len_p_file != "") {
			printf("reading insert length probabilities...\n");
			local->p_table.clear();
			ifstream ifs(ins_len_p_file.c_str());
			string readline;
			while (getline(ifs, readline)) {
				vector<string> tokens = string_tokenize(readline);
				__ASSERT(tokens.size() == 2 && is_int(tokens[0]) && is_num(tokens[1]), "ERROR: wrong format.");
				local->p_table[str2int(tokens[0])] = str2double(tokens[1]);
			}
			ifs.close();
			ifs.clear();
			printf ("done.\n");
		}

		if (task == "expression_analysis") {
			printf("computing gene expression values...\n");
			local->iem2.get_rates(local->p_table, read_length);
			local->iem2.output_data(local->output_prefix + ".sampling_rates");
			if (!quick) {
				local->iem2.compute();
				local->iem2.output_exp(local->output_prefix + ".exp.xls");
			}
			printf("done.\n");
		}

		if (task == "paired_end" && (!quick)) {
			ofstream ofs3((local->output_prefix + ".sampling_rates").c_str());
			ofstream ofs4((local->output_prefix + ".exp").c_str());
			int total_supporting_reads = 0;
			for (int i = 0; i < (int)eje.genes.size(); i++) {
				eje.genes[i].iem->get_rates();
				total_supporting_reads += (int)eje.genes[i].iem->sampling_rate.size();
			}
			ofs3 << total_supporting_reads << "\t" << read_length << endl;
			for (int i = 0; i < (int)eje.genes.size(); i++) {
				eje.genes[i].iem->compute(total_supporting_reads);
				ofs3 << eje.genes[i].name << "\t" << (int)eje.genes[i].subexon_sets.size() << "\t" << (int)eje.genes[i].trans.size() << "\t" << (int)eje.genes[i].iem->collapsed_sampling_rate.size()  << "\t" << (int)eje.genes[i].iem->collapsed_possible_sampling_rate.size() << endl;
				for (int k = 0; k < (int)eje.genes[i].trans.size(); k++) {
					ofs3 << eje.genes[i].trans[k].name << "\t";
				}
				ofs3 << endl;
				for (int j = 0; j < (int)eje.genes[i].subexon_sets.size(); j++) {
					ofs3 << eje.genes[i].len_subexon_sets[j] << "\t";
				}
				ofs3 << endl;
				ofs3 << "RPKM\t" << eje.genes[i].gene_expression << endl;
				if (eje.genes[i].iem->collapsed_possible_sampling_rate.size() > 0) {
					ofs3 << "MLE\t";
					for (int j = 0; j < (int)eje.genes[i].tran_expression.size(); j++) {
						ofs3 << eje.genes[i].tran_expression[j] << "\t";
					}
					ofs3 << endl;
					ofs3 << "SUM\t" << sum(eje.genes[i].tran_expression) << endl;;
				}
				for (int j = 0; j < (int)eje.genes[i].iem->collapsed_sampling_rate.size(); j++) {
					for (int k = 0; k < (int)eje.genes[i].trans.size(); k++) {
						ofs3 << eje.genes[i].iem->collapsed_sampling_rate[j][k] << "\t";
					}
					ofs3 << eje.genes[i].iem->collapsed_count[j];
					ofs3 << endl;
				}
				for (int j = 0; j < (int)eje.genes[i].iem->collapsed_possible_sampling_rate.size(); j++) {
					for (int k = 0; k < (int)eje.genes[i].trans.size(); k++) {
						ofs3 << eje.genes[i].iem->collapsed_possible_sampling_rate[j][k] << "\t";
					}
					ofs3 << endl;
				}
				for (int k = 0; k < (int)eje.genes[i].trans.size(); k++) {
					for (int j = 0; j < (int)eje.genes[i].subexon_sets.size(); j++) {
						ofs3 << eje.genes[i].trans[k].subexon_set_vector[j] << "\t";
					}
					ofs3 << endl;
				}

				bool single_isoform = (eje.genes[i].trans.size() == 1);
				bool is_zero = (eje.genes[i].iem->collapsed_possible_sampling_rate.size() == 0);
				ofs4 << eje.genes[i].name << "\t" << (single_isoform?eje.genes[i].gene_expression:(is_zero?0:sum(eje.genes[i].tran_expression))) << "\t" << (int)eje.genes[i].trans.size() << "\t";
				for (int j = 0; j < (int)eje.genes[i].trans.size(); j++) {
					if (j > 0) ofs4 << ",";
					ofs4 << eje.genes[i].trans[j].name;
				}
				ofs4 << "\t";
				for (int j = 0; j < (int)eje.genes[i].trans.size(); j++) {
					if (j > 0) ofs4 << ",";
					ofs4 << (single_isoform?eje.genes[i].gene_expression:(is_zero?0:eje.genes[i].tran_expression[j]));
				}
				ofs4 << endl;
				delete eje.genes[i].iem;
			}
			ofs3.close();
			ofs3.clear();
			ofs4.close();
			ofs4.clear();
		}

/*
		__ASSERT(mode_ins_len > read_length, "internal error: mode_ins_len <= read_length");
		printf("setting lambda as %lf.\n", 1.0/mode_ins_len);
		lambda = 1.0/mode_ins_len;

		printf("computing and writing filtering parameters...\n");
		q_table.clear();
		for (int i = 0; i < min(max_ins_len, mode_ins_len * 5); i++) {
			if (i <= read_length) {
				q_table.push_back(0);
			} else {
				q_table.push_back(insert_len_hist[i] / (lambda * exp(-lambda * i)));
			}
		}
		double max_q = *(std::max_element(q_table.begin(), q_table.end()));
		for (int i = 0; i < (int)q_table.size(); i++) {
			q_table[i] /= max_q;
		}
		ofs.open((output_prefix + ".filter_q").c_str());
		for (int i = 0; i < (int)q_table.size(); i++) {
			ofs << i << "\t" << q_table[i] << endl;
		}
		ofs.close();
		ofs.clear();
		printf ("done.\n");

		printf("computing categories...\n");
		ofstream ofs3((output_prefix + ".category").c_str());
		for (int i = 0; i < (int)eje.genes.size(); i++) {
			eje.genes[i].iem->lambda = lambda;
			ofs3 << eje.genes[i].name << "\t" << (int)eje.genes[i].trans.size() << endl;
			if (eje.genes[i].trans.size() > 1) {
				ofs3 << "category\tcount";
				for (int j = 0; j < (int)eje.genes[i].trans.size(); j++) {
					ofs3  << "\t" << eje.genes[i].trans[j].name;
				}
				ofs3 << endl;
				eje.genes[i].iem->compute_category();
				for (int j = 0; j < (int)eje.genes[i].iem->category.size(); j++) {
					ofs3 << j << "\t" << eje.genes[i].iem->category_read_count[j];
					for (int k = 0; k < (int)eje.genes[i].trans.size(); k++) {
						ofs3 << "\t" << (eje.genes[i].iem->category[j][k]?1:0) << ":" << eje.genes[i].iem->p_cf[j][k];
					}
					ofs3 << endl;
				}
			}
		}
		ofs3.close();
		ofs3.clear();
		printf("done.\n");
*/
	}
}

void handler(vector<map_result_struct> &results) {
	__ASSERT(in_range((int)results.size(), 1, 2), "internal error, results.size() != 1 or 2\n");
	if (results.size() == 2) {
		__ASSERT (task == "paired_end" || task == "expression_analysis", "internal error, results.size() == 2\n");
	}
	if (task == "paired_end" || (task == "expression_analysis" && local->other_args.size() == 3)) {
		__ASSERT (results.size() == 2, "internal error, results.size() != 2\n");
	}

	if (first_phase) total_reads++;
	if (report_progress > 0) {
		if (total_reads % report_progress == 0) {
			printf("%d reads processed.\n", total_reads);
		}
	}

	if (task == "scan_reads") scan_reads_handler(results);

	if (results.size() == 2) { //paired end
		if (results[0].targets.size() == 0 && results[1].targets.size() == 0 ) {
			if (first_phase) unmapped_reads++;
			return;
		} else if (results[0].targets.size() == 0 || results[1].targets.size() == 0 ) {
			if (first_phase) half_mapped_reads++;
			return;
		}
	} else {
		if (results[0].targets.size() == 0) {
			if (first_phase) unmapped_reads++;
			return;
		}
	}

	vector<int> min_mismatch;
	min_mismatch.resize(results.size());
	for (int k = 0; k < (int)results.size(); k++) {
		min_mismatch[k] = 10;
		for (int i = 0; i < (int)results[k].targets.size(); i++) {
			min_mismatch[k] = min(min_mismatch[k], results[k].targets[i].num_mismatch);
		}
		__ASSERT(min_mismatch[k] >= 0 && min_mismatch[k] <= read_length, "min_mismatch out of range");
	}

	for (int k = 0; k < (int)results.size(); k++) {
		if (min_mismatch[k] > num_mismatch) {
			if (first_phase) unmapped_reads++;
			return;
		}
	}

	if (first_phase) mapped_reads++; 

	vector<vector<target_struct> > targets;
	targets.resize(results.size());
	for (int k = 0; k < (int)results.size(); k++) {
		for (int i = 0; i < (int)results[k].targets.size(); i++) {
			if (results[k].targets[i].num_mismatch > min_mismatch[k] || results[k].targets[i].num_mismatch > num_mismatch) continue;
			size_t pos = results[k].targets[i].trans_id.find_first_of("/"); //remove file name for eland output
//			size_t pos1 = results[k].targets[i].trans_id.find_first_of(":");
			size_t pos2 = results[k].targets[i].trans_id.find_first_of(".");
// remove the filename before the first "/" in eland format, seq_name may still contain "/", file_name must contain "." and not contain ":"
//			if (pos != string::npos && pos2 != string::npos && pos2 < pos - 1 && (pos1 == string::npos || pos < pos1 - 1)) {
			if (pos != string::npos && pos2 != string::npos && pos2 < pos - 1) {
				results[k].targets[i].trans_id = results[k].targets[i].trans_id.substr(pos+1);
			} else { // no "/", such as humRibosomal
			}
			targets[k].push_back(results[k].targets[i]);
		}
	}

	for (int k = 0; k < (int)results.size(); k++) {
		__ASSERT(targets[k].size() >= 1, "internal error: no target");
	}

	__ASSERT(targets.size() == 1 || targets.size() == 2, "error: targets.size() != 1 or 2.\n");

	bool unique_flag = true;
	for (int k = 0; k < (int)targets.size(); k++) {
		if (task == "convert_coord") convert_handler(targets[k]);
		if (task == "extract") extract_handler(targets[k]);
		if (task == "comp_exp") gene_exp_handler(targets[k]);
		if (task == "map_stat") map_stat_handler(results[k], targets[k]);
		if (active_samples > 0) { // phased operations for sequential tasks
			local->total_read_counts++;
			unique_lock<mutex> lock(samples_lock);
			if (task == "two_sample_diff_gene_exp") two_sample_diff_gene_exp_handler(targets[0]);
			if (task == "two_sample_diff_exon_usage") two_sample_diff_exon_usage_handler(targets[0]);
			waiting_samples++;
			if (waiting_samples < active_samples) {
				samples_cond_var.wait(lock);
			} else {
				waiting_samples = 0;
				samples_cond_var.notify_all();
			}
		}
		if (targets[k].size() > 1) {
			unique_flag = false;
		}
	}

	if (task == "expression_analysis") {
		expression_analysis_handler(targets, results);
	} else {
		for (int k = 0; k < (int)results.size(); k++) {
			if (first_phase) total_targets += (int) targets[k].size();
		}
		if (!unique_flag) return;
		if (first_phase) unique_reads++; // fix later gene_exp
		if (task == "exon_usage") exon_usage_handler(targets[0]); // do not work together with expression_analysis
		if (task == "paired_end") paired_end_handler(targets, results);
	}
}

void do_task() {
	if (task == "pipeline") {
		if (task == "map_stat" || task == "expression_analysis") {
			__ASSERT(reference_file != "", "no reference file specified.\n");
			local->myfasta.clear();
			local->myfasta.read_from_file(reference_file);
			local->myfasta.build_tag_map();
		}
		if (task == "exon_usage" || task == "two_sample_diff_exon_usage") {
			__ASSERT(annotation_file != "", "no annotation file specified.\n");
			__ASSERT(local->mygenefile.read_from_text_file(annotation_file), "reading annotation file failed.\n");
			if (task == "two_sample_diff_exon_usage") {
				__ASSERT(exon_file != "", "no exon file specified.\n");
				__ASSERT(read_exon_file(), "reading exon file failed.\n");
			}
		}
		local->other_args.insert(local->other_args.begin(), "dummy");
	} else {
		if (version == 1) {
			if (task == "convert_coord") {
				local->mygenefile.read_from_text_file(local->other_args[0]);
			} else {
				if (task == "extract") {
					local->other_args.insert(local->other_args.begin(), "dummy");
				} else {
					local->myfasta.clear();
					local->myfasta.read_from_file(local->other_args[0]);
					local->myfasta.build_tag_map();
				}
			}
		} else if (version == 2) {
			eje.get_ready(local->other_args[0], read_length, read_length - 5);
		} else if (version == 3) {
			eje.get_ready(local->other_args[0], read_length, read_length - 5, 5, true, false, false);
		} else if (version == 4) {
			eje.get_ready(local->other_args[0], read_length, read_length - 5, 5, true, true, false);
		} else { //version == 0
			if (task == "exon_usage") {
				local->mygenefile.read_from_text_file(local->other_args[0]);
			} else {
				local->other_args.insert(local->other_args.begin(), "dummy");
			}
		}
	}

	if (local->output_prefix == "") local->output_prefix = local->other_args[1];

	//for map_stat
	if (task == "map_stat") {
		memset(map_stat_counts, 0, sizeof(map_stat_counts));
	}

	//for exon_usage
	if (task == "exon_usage") {
		local->mygenefile.build_interval_lists();
//		ofs.open((output_prefix + ".usage").c_str());
	}

	//for extract
	if (task == "extract") {
		if (version == 1) {
			ofs_F.open((local->output_prefix + "genes.F.coord").c_str());
			ofs_R.open((local->output_prefix + "genes.R.coord").c_str());
		} else { // version == 0
			ofs_F.open((local->output_prefix + ".F.coord").c_str());
			ofs_R.open((local->output_prefix + ".R.coord").c_str());
			ofs_M.open((local->output_prefix + ".M.coord").c_str());
		}
	}

	//for paired_end
	if (version >= 2 && task == "paired_end") {
		for (int i = 0; i < (int)eje.genes.size(); i++) {
//			eje.genes[i].iem = new iso_exp_model(eje.genes[i], lambda, filtering_q, read_length);
			eje.genes[i].iem = new iso_exp_model(eje.genes[i], read_length, local->p_table);
		}
		local->ofs.open((local->output_prefix + ".gene.map").c_str());
		ofs_F.open((local->output_prefix + ".reads.bed").c_str());
		ofs_R.open((local->output_prefix + ".reads.map").c_str());
	}

	//for convert_coord
	if (task == "convert_coord") {
		ofs_F.open((local->output_prefix + ".F.coord").c_str());
		ofs_R.open((local->output_prefix + ".R.coord").c_str());
		ofs_M.open((local->output_prefix + ".M.coord").c_str());
		if (version >= 2) local->ofs.open((local->output_prefix + ".gene.coord").c_str());
	}

	//for comp_exp
	genes.clear();
	if (task == "comp_exp" && version == 1) {
		for (int i = 0; i < (int)local->myfasta.tags.size(); i++) {
			size_t pos = local->myfasta.tags[i].find_first_of("$$");
			string gene_name = (pos == string::npos) ? local->myfasta.tags[i] : local->myfasta.tags[i].substr(0,pos);
//			__ASSERT(tokens.size() == 2, "ERROR: reference sequence tag is in wrong format. Please use option gen_trans to generate reference sequences.");
//			string gene_name = tokens[1];
//			string gene_name = myfasta.tags[i];

			if (genes.find(gene_name) == genes.end()) {
				genes[gene_name].length = (int)local->myfasta.sequences[i].length();
				genes[gene_name].reduce_len = 0;
				genes[gene_name].total_count = 0;
				genes[gene_name].unique_count = 0;
			} else {
				genes[gene_name].length = max(genes[gene_name].length, (int)local->myfasta.sequences[i].length());
			}
/*			//string trans_name = (tokens.size() == 2) ? tokens[0] : myfasta.tags[i];
			string trans_name = tokens[0];
			if (transcripts.find(trans_name) == transcripts.end()) {
				transcripts[trans_name].length = (int)myfasta.sequences[i].length();
				transcripts[trans_name].reduce_len = 0;
				transcripts[trans_name].total_count = 0;
				transcripts[trans_name].unique_count = 0;
			} else {
				printf("warning: duplicate transcript: %s.\n", trans_name.c_str());
				transcripts[trans_name].length = max(transcripts[trans_name].length, (int)myfasta.sequences[i].length());
			}*/
		}

		if (reduce_len_file != "") {
			ifstream ifs(reduce_len_file.c_str());
			string readline;
			while (getline(ifs, readline)) {
				//vector<string> tokens = string_tokenize(readline, ": \t");
				vector<string> tokens = string_tokenize(readline);
				//__ASSERT(tokens.size() == 3, "tokens.size() != 3");
				__ASSERT(tokens.size() == 2, "tokens.size() != 2, wrong format");
				//int count = str2int(tokens[2]);
				int count = str2int(tokens[1]);
				__ASSERT(count > 0, "count <= 0, wrong format");

				//string gene_name = tokens[1];
				string gene_name = tokens[0];
				__ASSERT(genes.find(gene_name) != genes.end(), string("gene not found, gene: " + gene_name));
				genes[gene_name].reduce_len = max(genes[gene_name].reduce_len, count);

/*				string trans_name = tokens[0];
				__ASSERT(transcripts.find(trans_name) != transcripts.end(), string("transcript not found, transcript: " + trans_name));
				if (transcripts[trans_name].reduce_len == 0) {
					transcripts[trans_name].reduce_len = count;
				} else {
					printf("warning: duplicate transcript: %s.\n", trans_name.c_str());
					transcripts[trans_name].reduce_len = max(transcripts[trans_name].reduce_len, count);
				}*/
			}
			ifs.close();
			ifs.clear();
		}

		for (map<string, gene_info_struct>::iterator it = genes.begin(); it != genes.end(); it++) {
			it->second.total_count = 0;
			it->second.unique_count = 0;
			it->second.length -= it->second.reduce_len;
			it->second.length -= read_length;
			if (it->second.length <= 0) {
				//printf("WARNING: gene length <= %d, gene: %s\n", read_length, it->first.c_str());
				it->second.length = 2000000000;
			}
		}
		
/*		for (map<string, gene_info_struct>::iterator it = transcripts.begin(); it != transcripts.end(); it++) {
			it->second.total_count = 0;
			it->second.unique_count = 0;
			it->second.length -= it->second.reduce_len;
			it->second.length -= read_length;
			if (it->second.length <= 0) {
				//printf("WARNING: transcript length <= %d, transcript: %s\n", read_length, it->first.c_str());
				it->second.length = 1;
			}
		}*/
	}

	// for scan_reads
	if (task == "scan_reads") {
		local->nt_counts.resize(MAX_READ_LEN);
		for (int i = 0; i < MAX_READ_LEN; i++) {
			local->nt_counts[i].resize(6);
			for (int j = 0; j < 6; j++) {
				local->nt_counts[i][j] = 0;
			}
		}
	}

	// for expression_analysis
	if (task == "expression_analysis") {
		local->iem2.build_from_fasta(local->myfasta);
		if (annotation_file != "") { // can not work together with -exon_usage, fix later
			__ASSERT(local->mygenefile.read_from_text_file(annotation_file), "reading annotation file failed.\n");
			local->mygenefile.build_map();
			local->iem2.check_genefile(local->mygenefile);
			local->ofs.open((local->output_prefix + ".reads.bed").c_str());
		}
	}

	map_result_reader myreader;

	printf("reading mapping results...\n");

	first_phase = true;

	if (task == "paired_end" || ((task == "pipeline" || task == "expression_analysis") && local->other_args.size() == 3)) {
		__ASSERT(myreader.read_from_file(handler, local->other_args[1], local->other_args[2]), "internal error: read file failed\n");
	} else {
		__ASSERT(myreader.read_from_file(handler, local->other_args[1]), "internal error: read file failed\n");
	}

	//if ((task == "comp_exp" && version == 1) || task == "expression_analysis") { // implmenent second phase for expression_analysis later
	if (task == "comp_exp") { //for comp_exp
		first_phase = false;
		printf("second phase...\n");
		__ASSERT(myreader.read_from_file(handler, local->other_args[1]), "internal error: read file failed\n");
	}

	if (local->other_args.size() == 3 || working_mode == MR_SAM_PAIRED) { 
		printf ("reads:\n");
		printf( "   total:                %d\n", total_reads);
		printf( "   one end mapped:       %d\n", half_mapped_reads);
		printf( "   different chromosome: %d\n", diff_chr_reads);
		printf( "   same direction:       %d\n", same_dir_reads);
		printf( "   negative distance:    %d\n", neg_dist_reads);
		printf( "   large distance:       %d\n", large_dist_reads);
		printf( "   both ends mapped:     %d\n", mapped_reads);
		printf( "   others (rRNA, etc):   %d\n", other_reads);
		printf( "   not mapped:           %d\n", unmapped_reads);
		printf( "   correctly paired:     %d\n", paired_reads);
		printf( "   unique:               %d\n", unique_reads);
		printf( "   total targets:        %d\n", total_targets);
	} else {
		printf ("reads:\n");
		printf( "   total:            %d\n", total_reads);
		printf( "   not mapped:       %d\n", unmapped_reads);
		printf( "   mapped:           %d\n", mapped_reads);
		printf( "   unique:           %d\n", unique_reads);
		printf( "   total targets:    %d\n", total_targets);
	}

	if (task == "map_stat") {
		map_stat();
	}
	
	if (task == "paired_end") {
		local->ofs.close();
		local->ofs.clear();
		ofs_F.close();
		ofs_F.clear();
		ofs_R.close();
		ofs_R.close();
		paired_end();
	}
	
	if (task == "comp_exp") {
		comp_exp();
	}
	
	if (task == "scan_reads") {
		cout << "output read statistics...\n";
		ofstream scan_read_file((local->output_prefix + ".scan").c_str());

		scan_read_file << "total " << num_total_read << " reads, " << num_good_read << "(" << round_double(num_good_read*100.0/num_total_read) << "%) good reads, ";
		scan_read_file << num_total_read-num_good_read << " (" << 100-round_double(num_good_read*100.0/num_total_read) << "%) bad reads.\n";
		scan_read_file << "minimum read length " << min_read_len << ", maximum read length " << max_read_len << ".\n";
		if (min_read_len != max_read_len) scan_read_file << "INCONSISTENT READ LENGTHS!!!\n";
		if (max_read_len == 0 || num_total_read == 0) {
			scan_read_file << "NO READ, FAILED!!!\n";
		} else {
			scan_read_file << "GC content " << round_double(total_gc*100.0/total_nt) << "%.\n";

			scan_read_file << "nucleotide distribution by position:\n";
			scan_read_file << "A\t\tC\t\tG\t\tT\t\tN\n";
			for (int i = 0; i < max_read_len; i++) {
				for (int j = 0; j < 5; j++) {
					scan_read_file << local->nt_counts[i][j] << " (" << round_double(local->nt_counts[i][j]*100.0/local->nt_counts[i][5]) << "%)\t";
				}
				scan_read_file << "\n";
			}
		}
		scan_read_file.close();
		scan_read_file.clear();
		cout << "done.\n";
	}

	if (task == "expression_analysis") {
		if (annotation_file != "") {
			local->ofs.close();
			local->ofs.clear();
		}
		paired_end();//
	}

	if (task == "two_sample_diff_exon_usage") {
		lock_guard<mutex> lk(samples_lock);
		if (active_samples == 2) { // first sample finished
			two_sample_diff_exon_usage_summary();
		}
	}
	
	if (task == "exon_usage") {
		printf("reads:\n");
		if (unique_reads == 0) unique_reads = 1;
		printf("   intergene: %d (%d%%)\n", intergene_reads, intergene_reads*100/unique_reads);
		printf("   intron:    %d (%d%%)\n", intron_reads, intron_reads*100/unique_reads);
		printf("   5UTR:      %d (%d%%)\n", UTR5_reads, UTR5_reads*100/unique_reads);
		printf("   3UTR:      %d (%d%%)\n", UTR3_reads, UTR3_reads*100/unique_reads);
        printf("   exon:      %d (%d%%)\n", exon_reads, exon_reads*100/unique_reads);
//		ofs.close();
//		ofs.clear();
		if (local->other_args[1].substr(0, 4) == "run:" && local->output_prefix == local->other_args[1]) {
			local->ofs.open(string("run.exp.xls").c_str());
		} else {
			local->ofs.open((local->output_prefix + ".exp.xls").c_str());
		}
		local->ofs << "GeneName\tGeneID\tchr\tstrand\tstart\tend\tgene_length\tgene_count\tgene_rpkm\texons\texon_coords\texon_lengths\texon_counts\texon_rpkms\n";
		for (int i = 0; i < (int)local->mygenefile.sequences.size(); i++) {
			for (int j = 0; j < (int)local->mygenefile.sequences[i].genes.size(); j++) {
				gene_struct& gene = local->mygenefile.sequences[i].genes[j];
				local->ofs << " " << gene.geneName << "\t" << gene.name << "\t" << gene.chrom << "\t" << (gene.strand?"+":"-") << "\t" << gene.txStart << "\t" << gene.txEnd << "\t";
				int gene_length = 0;
				int gene_count = 0;
				string str_exon_coords = "";
				string str_exon_lengths = "";
				string str_exon_counts = "";
				string str_exon_rpkms = "";
				for (int k = 0; k < (int)gene.exons.size(); k++) {
					str_exon_coords += string("[") + int2str(gene.exons[k].first) + "," + int2str(gene.exons[k].second) + "],";
					int exon_length = gene.exons[k].second - gene.exons[k].first;
					str_exon_lengths += int2str(exon_length) + ",";
					gene_length += exon_length;
					int exon_count = local->exon_counts[pair<int, pair<int, int> >(i, pair<int, int>(j, k))];
					str_exon_counts += int2str(exon_count) + ",";
					str_exon_rpkms += double2str((double)exon_count/exon_length/unique_reads*1000.0*1000000.0) + ",";
					gene_count += exon_count;
				}
				local->ofs << gene_length << "\t" << gene_count << "\t" << double2str((double)gene_count/gene_length/unique_reads*1000.0*1000000.0) << "\t";
				local->ofs << gene.exons.size() << "\t" << str_exon_coords << "\t" << str_exon_lengths << "\t" << str_exon_counts << "\t" << str_exon_rpkms << endl;
			}
		}
		local->ofs.close();
		local->ofs.clear();
	}
	
	if (task == "convert_coord") {
		ofs_F.close();
		ofs_F.clear();
		ofs_R.close();
		ofs_R.clear();
		ofs_M.close();
		ofs_M.clear();
		if (version >= 2) {
			local->ofs.close();
			local->ofs.clear();
		}
	}
	
	if (task == "extract") {
		ofs_F.close();
		ofs_F.clear();
		ofs_R.close();
		ofs_R.clear();
		if (version == 0) {
			ofs_M.close();
			ofs_M.clear();
		}
	}
}

void task_thread(int sample, string filename) {
	sample_id = sample;
	local = new local_data();
	{
		lock_guard<mutex> lk(samples_lock);
		samples.push_back(local);
	}
	__ASSERT(local->other_args.size() == 0, "error: other_args.size() != 0\n");
	local->other_args.push_back(filename);
	do_task();
	lock_guard<mutex> lk(samples_lock);
	active_samples--;
	samples_cond_var.notify_all();
}

void do_sequential_tasks() {
	vector<thread*> threads;
	active_samples = (int)local->other_args.size();
	waiting_samples = 0;
	for (int i = 0; i < (int)local->other_args.size(); i++) {
		threads.push_back(new thread(task_thread, i, local->other_args[i]));
	}
	for (int i = 0; i < (int)threads.size(); i++) {
		threads[i]->join();
		delete threads[i];
		threads[i] = NULL;
	}
	threads.clear();
	for (int i = 0; i < (int)samples.size(); i++) {
		delete samples[i];
		samples[i] = NULL;
	}
	samples.clear();
}

#endif //#ifdef ISO_EXP_H
