/*
fasta.h - Header file for fasta file manipulation
written by JIANG Hui, 
Institute for Computational and Mathematical Engineering, Stanford University
March, 2008 -
*/

#ifndef FASTA_H
///Define this macro to prevent from including this header file more than once.
#define FASTA_H

#include "stl.h"
#include "file_operation.h"

class fasta{
public:
	bool verbose;
	vector<string> tags;
	vector<string> sequences;

	map<string, int> tag_map;

	string _id;
	ifstream* _ifs;

	fasta(){verbose = true; clear(); _id = ""; _ifs = NULL;};
	void clear(){tags.clear(); sequences.clear();};
	void build_tag_map() {
		tag_map.clear();
		for (int i = 0; i < (int)tags.size(); i++) {
			__ASSERT(tag_map.count(tags[i])==0, string("error: duplicated transcript ") + tags[i]);
			tag_map[tags[i]] = i;
		}
	}
	string get_sequence(const string tag){
		if (tag_map.size() > 0) {
			if (tag_map.count(tag) == 1) return sequences[tag_map[tag]];
			else return "";
		}
		for (int i = 0; i < (int)tags.size(); i++) {
			if (tag == tags[i]) return sequences[i];
		}
		return "";
	}
	string get_sequence(const string tag, int start, int end){
		if (start < 0 || start > end) return "";
		if (tag_map.size() > 0) {
			if (tag_map.count(tag) == 1) {
				if (end > (int)sequences[tag_map[tag]].length()) return "";
				return sequences[tag_map[tag]].substr(start, end - start);
			}
			else return "";
		}
		for (int i = 0; i < (int)tags.size(); i++) {
			if (tag == tags[i]) {
				if (end > (int)sequences[i].length()) return "";
				return sequences[i].substr(start, end - start);
			}
		}
		return "";
	}
	string get_sequence(const string tag, vector<pair<int, int> > &regions, int separate_N = 0){
		string result = "";
		string seperator = "";
		seperator.insert(seperator.begin(), separate_N, 'N');
		for (int i = 0; i < (int)tags.size(); i++) {
			if (tag == tags[i]) {
				for (int j = 0; j < (int)regions.size(); j++) {
					result += sequences[i].substr(regions[j].first, regions[j].second - regions[j].first);
					if (j < (int)regions.size() - 1 && separate_N > 0) result += seperator;
				}
				break;
			}
		}
		return result;
	}
	bool read_from_file(const string filename, bool check_sequence = true);
	bool write_to_file(FILE *file, const string tag, const string &sequence);
	bool write_to_file(const string filename);
	bool enumerate(const string filename, int read_length = 25);
	bool get_first_sequence(const string &filename, string &id, string &seq);
	bool get_next_sequence(string &id, string &seq);
};

inline bool is_fasta_file(const string file_name, bool allow_fastq = true){
	if (!file_exists(file_name)) return false;
	ifstream ifs(file_name.c_str());
	string readline;
	bool result = true;

	while(getline(ifs, readline)) {
		if (readline == "") continue;
		if (!(readline[0] == '>') && !(allow_fastq && readline[0] == '@')) result = false;
		break;
	}

	ifs.close();
	return result;
}

inline bool fasta::read_from_file(const string filename, bool check_sequence){
	if (!is_fasta_file(filename)) {
		if (verbose) cout << filename << " is not a vaid fasta file.\n";
		return false;
	}

	if (verbose) printf("Loading fasta file %s.\n", filename.c_str());

	string tag = "";
	string sequence = "";
	string readline;
	ifstream ifs(filename.c_str());
	int total_len = 0;
	bool is_sequence = false;
	while(getline(ifs, readline)) {
		if (readline == "") continue;
		if (readline[0] == '>' || readline[0] == '@') {
			if (sequence != "") {
				tags.push_back(tag);
				sequences.push_back(sequence);
				total_len += (int)sequence.length();
			}
			vector<string> tokens = string_tokenize(readline.substr(1));
			if (tokens.size() == 0) tag = ""; else tag = tokens[0];
			sequence = "";
			is_sequence = true;
		} else if (readline[0] == '+') {
			is_sequence = false;
		} else {
			if (is_sequence) {
				if (check_sequence) {
					for (int i = 0; i < (int)readline.length(); i++) {
						if (!isalpha(readline[i])) {
							if (verbose) printf("failed: error letter %c occurred when reading fasta file %s.\n", readline[i], filename.c_str());
							ifs.close();
							return false;
						}
					}
				}
				sequence += readline;
			}
		}
	}

	if (sequence != "") {
		tags.push_back(tag);
		sequences.push_back(sequence);
		total_len += (int)sequence.length();
	}
	ifs.close();
	if (verbose) printf("%d sequences read, total length is %d.\n", (int)sequences.size(), total_len);
	return true;
}

inline bool fasta::write_to_file(FILE *file, const string tag, const string &sequence) {
	fprintf(file, ">%s\n", tag.c_str());
	for (int j = 0; j < (int)sequence.length(); j+=50) {
		fprintf(file, "%s\n", sequence.substr(j,50).c_str());
	}
	return true;
}

inline bool fasta::write_to_file(const string filename) {
	FILE *file=fopen(filename.c_str(), "wt");
	if (!file){
		cout << "error opening file for writing:" << filename << endl;
		return false;
	}
	if (verbose) printf("writing fasta file %s.\n", filename.c_str());
	int total_len = 0;
	for (int i = 0; i < (int)tags.size(); i++) {
		if (!write_to_file(file, tags[i], sequences[i])) return false;
		total_len += (int)sequences[i].length();
	}
	if (verbose) printf("%d sequences wrote, total length is %d.\n", (int)sequences.size(), total_len);
	fclose(file);
	return true;
}

inline bool is_dna(const char letter) {
	char temp = letter;
	if (temp >= 'a' && temp <= 'z') temp += 'A' - 'a';
	if (temp == 'A' || temp == 'C' || temp == 'G' || temp == 'T' || temp == 'N' || temp == '.') return true;
	else return false;
}

inline bool is_read(const string &read) {
	string temp = trim_space(toupper(read));
	if (temp == "") return false;
	for (int i = 0; i < (int)temp.length(); i++) {
		if (!is_dna(temp[i])) return false;
	}
	return true;
}

inline bool is_dna_file(const string file_name){
	if (!file_exists(file_name)) return false;
	ifstream ifs(file_name.c_str());
	string readline;
	if (!getline(ifs, readline)) goto failure;
	if (!is_read(readline)) goto failure;
	ifs.close();
	return true;
failure:
	ifs.close();
	return false;
}

inline bool read_dna_from_file(const string filename, vector<string> &sequences) {
	if (!is_dna_file(filename)) {
		cout << filename << " is not a vaid DNA sequence file.\n";
		return false;
	}

	printf("Loading DNA sequence file %s.\n", filename.c_str());

	string readline;
	ifstream ifs(filename.c_str());
	int total_len = 0;
	while(getline(ifs, readline)) {
		if (readline == "") continue;
		sequences.push_back(readline);
		total_len += (int)readline.length();
	}

	ifs.close();
	printf("%d sequences read, total length is %d.\n", (int)sequences.size(), total_len);
	return true;
}

//file_type = -1: unrecognizable, 0:DNA, 1:fasta, 2:fastq
inline bool is_read_file(const string file_name, int &file_type) {
	file_type = -1;
	if (!file_exists(file_name)) return false;
	if (is_dna_file(file_name)) file_type = 0;
	else if (is_fasta_file(file_name, false)) file_type = 1;
	else if (is_fasta_file(file_name, true)) file_type = 2;
	else panic("error: input file format unrecognizable.\n");
	if (file_type >= 0) return true;
	else return false;
}

//file_type = -1: any type, 0:DNA, 1:fasta, 2:fastq
//return if a read is got
inline bool get_read(ifstream &ifs, string &read, int file_type = -1) {
	string read_line;
	if (!getline(ifs, read_line)) return false;
	if (read_line.length() > 0 && read_line[0] == '>') {
		__ASSERT((file_type == -1 || file_type == 1), string("error: file format error in line ") + read_line);
		__ASSERT(getline(ifs, read_line), string("error: file read error in line ") + read_line);
		__ASSERT(is_read(read_line), string("error: file read error in line ") + read_line);
		read = read_line;
		return true;
	} else if (read_line.length() > 0 && read_line[0] == '@') {
		__ASSERT((file_type == -1 || file_type == 2), string("error: file format error in line ") + read_line);
		__ASSERT(getline(ifs, read_line), string("error: file read error in line ") + read_line);
		__ASSERT(is_read(read_line), string("error: file read error in line ") + read_line);
		read = read_line;
		__ASSERT(getline(ifs, read_line), string("error: file read error in line ") + read_line);
		__ASSERT(read_line.length() > 0 && read_line[0] == '+', string("error: file format error in line ") + read_line);
		__ASSERT(getline(ifs, read_line), string("error: file read error in line ") + read_line);
		__ASSERT(read_line.length() == read.length(), string("error: file format error in line ") + read_line);
		return true;
	} else {
		__ASSERT((file_type == -1 || file_type == 0), string("error: file format error in line ") + read_line);
		__ASSERT(is_read(read_line), string("error: file read error in line ") + read_line);
		read = read_line;
		return true;
	}
}

inline bool fasta::enumerate(const string filename, int read_length) {
	FILE *file=fopen(filename.c_str(), "wt");
	if (!file){
		cout << "error opening file for writing:" << filename << endl;
		return false;
	}
	if (verbose) printf("writing fasta file %s.\n", filename.c_str());
	int total_sequences = 0, total_len = 0;
	for (int i = 0; i < (int)tags.size(); i++) {
		int l = (int)sequences[i].length();
		for (int j = 0; j <= l - read_length; j++) {
			if (!write_to_file(file, tags[i] + ":" + int2str(j+1), sequences[i].substr(j, read_length))) return false;
			total_len += read_length;
			total_sequences++;
		}
	}
	if (verbose) printf("%d sequences wrote, total length is %d.\n", total_sequences, total_len);
	fclose(file);
	return true;
}

inline bool fasta::get_first_sequence(const string &filename, string &id, string &seq) {
	if (!is_fasta_file(filename)) {
		if (verbose) cout << filename << " is not a vaid fasta file.\n";
		return false;
	}

	if (verbose) printf("Loading fasta file %s.\n", filename.c_str());

	_id = "";
	_ifs = new ifstream(filename.c_str());
	return get_next_sequence(id, seq);
}

inline bool fasta::get_next_sequence(string &id, string &seq) {
	string readline;
	id = "";
	seq = "";
	bool is_sequence = true;
	if (_ifs == NULL) return false;
	while(getline(*_ifs, readline)) {
		if (readline == "") continue;
		if (readline[0] == '>' || readline[0] == '@') {
			id = _id;
			_id = readline.substr(1);
			if (id != "") return true;
		} else if (readline[0] == '+') {
			is_sequence = false;
		} else {
			if (is_sequence) {
				seq += readline;
			}
		}
	}

	id = _id;
	_ifs->close();
	delete _ifs;
	_ifs = NULL;
	return true;
}

#endif //FASTA_H
