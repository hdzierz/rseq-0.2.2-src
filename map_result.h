/*
map_result.h - Header file for map result file manipulation
written by JIANG Hui, 
Institute for Computational and Mathematical Engineering, Stanford University
March, 2008 -
*/

#ifndef MAP_RESULT_H
///Define this macro to prevent from including this header file more than once.
#define MAP_RESULT_H

#include "stl.h"
#include "file_operation.h"
#include "probe.h"

class target_struct{
public:
	string trans_id;
	int trans_coord;
	int num_mismatch;
	bool reverse_strand;
};

class map_result_struct{
public:
	string probe_id;
	string probe_seq;
	vector<target_struct> targets;
};

typedef void map_result_handler(vector<map_result_struct> &result);

#define MR_UNKNOWN 0
#define MR_ELAND_MULTI 1
#define MR_SAM 2
#define MR_SAM_PAIRED 3

thread_local int working_mode = MR_UNKNOWN; //MR_UNKNOWN, MR_ELAND_MULTI, MR_SAM, MR_SAM_PAIRED

class map_result_reader{
public:
	map_result_reader(){};
	inline bool read_from_file(map_result_handler handler, const string filename, const string filename2 = "");
	inline bool get_read(istream &ifs, vector<map_result_struct> &results, string &sam_readline, bool &sam_header);
	inline bool get_read_eland_multi(istream &ifs, vector<map_result_struct> &results, bool silent = false);
	inline bool get_read_sam(istream &ifs, vector<map_result_struct> &results, string &sam_readline, bool &sam_header, bool paired = false, bool silent = false);
	inline bool detect_file_format(const string filename);
};

bool map_result_reader::detect_file_format(const string filename) {
	if (working_mode != MR_UNKNOWN) return true;
	if (filename == "stdin" || filename.substr(0, 4) == "run:") {
		working_mode = MR_SAM;
		return true;
	}
	if (!file_exists(filename)) return false;
	ifstream ifs(filename);
	vector<map_result_struct> result;
	string sam_readline = "";
	bool sam_header = true;
	if (get_read_sam(ifs, result, sam_readline, sam_header, true, true)) working_mode = MR_SAM_PAIRED;
	else if (get_read_sam(ifs, result, sam_readline, sam_header, false, true)) working_mode = MR_SAM;
	else if (get_read_eland_multi(ifs, result, true)) working_mode = MR_ELAND_MULTI;
	else {
		ifs.close();
		return false;
	}
	ifs.close();
	return true;
}

bool map_result_reader::get_read(istream &ifs, vector<map_result_struct> &results, string &sam_readline, bool &sam_header){
	results.clear();
	switch (working_mode) {
		case MR_ELAND_MULTI:
			return get_read_eland_multi(ifs, results);
			break;
		case MR_SAM_PAIRED:
			return get_read_sam(ifs, results, sam_readline, sam_header, true);
			break;
		case MR_SAM:
			return get_read_sam(ifs, results, sam_readline, sam_header);
			break;
		default:
			printf("unsupported format.\n");
			return false;
			break;
	}
}

bool map_result_reader::get_read_sam(istream &ifs, vector<map_result_struct> &results, string &sam_readline, bool &sam_header, bool paired, bool silent){
	string probe_id = "";
	results.clear();
	while (true) {
		if (sam_readline == "") {
			if (!getline(ifs, sam_readline)) {
				if (probe_id == "") goto failure;
				else goto success;
			}
		}
		vector<string> tokens = string_tokenize(sam_readline, "\t", false);
		if (tokens.size() < 2) goto failure;
		if (sam_header && (tokens[0] == "@HD" || tokens[0] == "@SQ" || tokens[0] == "@RG" || tokens[0] == "@PG" || tokens[0] == "@CO")) {
			sam_readline = "";
			continue;
		}
		sam_header = false;

		if (tokens.size() < 11) goto failure;
		if (tokens[0] == "") goto failure;
		if (probe_id == "") probe_id = tokens[0];
		if (probe_id != tokens[0]) goto success;
		int index = -1;
		string probe_seq = tokens[9];
		if (!is_int(tokens[1])) goto failure;
		int FLAG = str2int(tokens[1]);
		//if ((paired && (0x1 & FLAG) == 0) || (!paired && (0x1 & FLAG) != 0)) goto failure;
		if (paired && (0x1 & FLAG) == 0) goto failure;
		if ((0x10 & FLAG) != 0) probe_seq = get_reverse2(probe_seq);
		if (!paired && results.size() > 0 && results[0].probe_seq != probe_seq) goto success;
		for (int i = 0; i < (int)results.size(); i++) {
			if (results[i].probe_seq == probe_seq) {
				index = i;
				break;
			}
		}
		if (index == -1) {
			map_result_struct temp;
			temp.probe_id = probe_id;
			temp.probe_seq = probe_seq;
			temp.targets.clear();			
			results.push_back(temp);
			index = (int)results.size() - 1;
		}
		
		if ((0x4 & FLAG) == 0) {
			target_struct temp_target;
			if (tokens[2] == "") goto failure; 
			temp_target.trans_id = tokens[2];
			if (!is_int(tokens[3])) goto failure;
			temp_target.trans_coord = str2int(tokens[3])-1;
			temp_target.num_mismatch = 0;
			temp_target.reverse_strand = ((0x10 & FLAG) != 0);

			for (uint i = 11; i < tokens.size(); i++) {
				vector<string> tokens1 = string_tokenize(tokens[i], ":");
				if (tokens1.size() != 3) goto failure;
				if (tokens1[0] == "NM") {
					if (tokens1[1] != "i" || !is_int(tokens1[2])) goto failure;
					temp_target.num_mismatch = str2int(tokens1[2]);
				} else if (tokens1[0] == "XA") {
					if (tokens1[2] == "") goto failure;
					if (tokens1[1] != "Z" || tokens1[2][tokens1[2].length() - 1] != ';') continue; //not bwa formart
					vector<string> tokens2 = string_tokenize(tokens1[2], ";");
					if (tokens2.size() == 0) goto failure;
					for (uint j = 0; j < tokens2.size(); j++) {
						vector<string> tokens3 = string_tokenize(tokens2[j], ",");
						if (tokens3.size() != 4) goto failure;
						target_struct temp_target1;
						if (tokens3[0] == "") goto failure;
						temp_target1.trans_id = tokens3[0];
						if (tokens3[1].length() < 2) goto failure;
						if (tokens3[1][0] != '+' && tokens3[1][0] != '-') goto failure;
						temp_target1.reverse_strand = (tokens3[1][0] == '-');
						if (!is_int(tokens3[1].substr(1))) goto failure;
						temp_target1.trans_coord = str2int(tokens3[1].substr(1))-1;
						if (!is_int(tokens3[3])) goto failure;
						temp_target1.num_mismatch = str2int(tokens3[3]);
						results[index].targets.push_back(temp_target1);
					}
				}
			}
			results[index].targets.push_back(temp_target);
		}
		sam_readline = "";
	}

success:
	if ((paired && results.size() != 2) || (!paired && results.size() != 1)) goto failure;
	return true;
failure:
	if (!silent && sam_readline != "") printf("bad format in line: %s\n", sam_readline.c_str());
	return false;
}

bool map_result_reader::get_read_eland_multi(istream &ifs, vector<map_result_struct> &results, bool silent){
	results.clear();
	map_result_struct result;
	string readline;
	vector<string> tokens;
	if (!getline(ifs, readline)) goto failure;
	tokens = string_tokenize(readline, " \t");
	if (tokens.size() < 3 || tokens.size() > 4) goto failure;
	result.probe_id = tokens[0];
	if (result.probe_id == "") goto failure;
	if (result.probe_id[0] == '>') result.probe_id = result.probe_id.substr(1);
	result.probe_seq = tokens[1];
	if (tokens.size() == 4) {
		vector<string> tokens1 = string_tokenize(tokens[3], ",");
		if (tokens1.size() == 0) goto failure;
		target_struct temp_target;
		temp_target.trans_id = "";
		temp_target.trans_coord = -1;
		temp_target.num_mismatch = -1;
		temp_target.reverse_strand = false;
		for (int i = 0; i < (int)tokens1.size(); i++) {
			string first = "", second = "";
			size_t pos = tokens1[i].find_last_of(':');
			if (pos == string::npos) {
				second = tokens1[i];
			} else {
				first = tokens1[i].substr(0, pos);
				second = tokens1[i].substr(pos + 1);
			}
			if (first != "") temp_target.trans_id = first;
			size_t len = second.length();
			if(len < 3) goto failure;
			if (second[len - 2] != 'R' && second[len - 2] != 'F') goto failure;
			temp_target.reverse_strand = (second[len - 2] == 'R');
//			if (second[len - 1] != '0' && second[len - 1] != '1'&& second[len - 1] != '2') goto failure; //enable >2bp mismatch
			temp_target.num_mismatch = second[len - 1] - '0';
			if (!is_int(second.substr(0, len - 2))) goto failure;
			temp_target.trans_coord = str2int(second.substr(0, len - 2)) - 1;
			if (temp_target.trans_coord < 0) goto failure;
			if (temp_target.trans_id == "") goto failure;
			result.targets.push_back(temp_target);
		}
	}
	results.push_back(result);
	return true;
failure:
	if (!silent && readline != "") printf("bad format in line: %s\n", readline.c_str());
	return false;
}

class stdiobuf
    : public std::streambuf
{
private:
    FILE* d_file;
    char  d_buffer[8192];
public:
    stdiobuf(FILE* file): d_file(file) {}
    ~stdiobuf() {/* if (this->d_file) fclose(this->d_file); */} //can not close file twice
    int underflow() {
        if (this->gptr() == this->egptr() && this->d_file) {
			size_t size = fread(this->d_buffer, 1, 8192, this->d_file);
            this->setg(this->d_buffer, this->d_buffer, this->d_buffer+size);
        }
        return this->gptr() == this->egptr()
            ? traits_type::eof()
            : traits_type::to_int_type(*this->gptr());
    }
};

bool map_result_reader::read_from_file(map_result_handler handler, const string filename1, const string filename2) {
	if(filename1 == "" || (filename1 != "stdin" && filename1.substr(0, 4) != "run:" && !file_exists(filename1))) {
		printf("error: file %s does not exist!\n", filename1.c_str());
		return false;
	}
	if(filename2 != "" && filename2 != "stdin" && filename2.substr(0, 4) != "run:" && !file_exists(filename2)) {
		printf("error: file %s does not exist!\n", filename2.c_str());
		return false;
	}
	if (filename1 == filename2) {
		printf("error: file names can not be the same!\n");
		return false;
	}
	if (!detect_file_format(filename1)) {
		printf("unsupported file format: %s\n", filename1.c_str());
		return false;
	}

	ifstream ifs1;
	istream* is1 = NULL;
	FILE *fp1 = NULL;
	stdiobuf *buf1 = NULL;
	if (filename1 == "stdin") {
		is1 = &cin;
	} else if (filename1.substr(0, 4) == "run:") {
		fp1 = popen(filename1.substr(4).c_str(), "r");
		if (fp1 == NULL) {
			printf("open file 1 failed.\n");
			return false;
		}
		buf1 = new stdiobuf(fp1);
		is1 = new istream(buf1);
	} else {
		ifs1.open(filename1.c_str());
		is1 = &ifs1;
	}
	ifstream ifs2;
	istream *is2 = NULL;
	FILE *fp2 = NULL;
	stdiobuf *buf2 = NULL;
	if (filename2 == "stdin") {
		is2 = &cin;
	} else if (filename2.substr(0, 4) == "run:") {
		fp2 = popen(filename2.substr(4).c_str(), "r");
		if (fp2 == NULL) {
			printf("open file 2 failed.\n");
			return false;
		}
		buf2 = new stdiobuf(fp2);
		is2 = new istream(buf2);
	} else {
		ifs2.open(filename2.c_str());
		is2 = &ifs2;
	}
	if (filename2 == "") 	
		cout << "reading file " << filename1 << "...\n";
	else 
		cout << "reading files " << filename1 << " and " << filename2 << "...\n";

	string sam_readline = "";
	bool sam_header = true;
	string sam_readline2 = "";
	bool sam_header2 = true;
	while (true) {
		vector<map_result_struct> temp_result;
		if (!get_read(*is1, temp_result, sam_readline, sam_header)) break;
		vector<map_result_struct> results = temp_result;
		if (results.size() == 0) panic("error: no mapped reads in the first mapped file.\n");
		if (filename2 != "") {
			__ASSERT(working_mode != MR_SAM_PAIRED, "error: only one mapped file is needed for SAM paired format.\n");
			if (results.size() > 1) panic("error: the first mapped file seems to be in SAM paired format.\n");
			vector<map_result_struct> temp_result2;
			if (!get_read(*is2, temp_result2, sam_readline2, sam_header2)) {
				printf("inconsistent paired files.\n");
				return false;
			}
			if (temp_result2.size() == 0) panic("error: no mapped reads in the second mapped file.\n");
			else if (temp_result2.size() > 1) panic("error: the second mapped file seems to be in SAM paired format.\n");
			results.push_back(temp_result2[0]);
		}
		__ASSERT(results.size() <= 2, "error: more than two reads in a pair.\n");
		handler(results);
	}

	vector<map_result_struct> temp_result;
	if (filename2 != "" && get_read(*is2, temp_result, sam_readline2, sam_header2)) {
		printf("inconsistent paired files.\n");
		return false;
	}

	if (filename1 == "stdin") {
	} else if (filename1.substr(0, 4) == "run:") {
		pclose(fp1);
		delete buf1;
		delete is1;
	} else {
		ifs1.close();
	}
	is1 = NULL;
	fp1 = NULL;
	buf1 = NULL;
	if (filename2 != "") {
		if (filename2 == "stdin")  {
		} else if (filename2.substr(0, 4) == "run:") {
			pclose(fp2);
			delete buf2;
			delete is2;
		} else {
			ifs2.close();
		}
	}
	is2 = NULL;
	fp2 = NULL;
	buf2 = NULL;
	printf("done.\n");
	return true;
}

#endif //#ifdef MAP_RESULT_H
