/*
bar.h - Header file for bar file manipulation
written by JIANG Hui, 
Institute for Computational and Mathematical Engineering, Stanford University
May, 2007 -
*/

#ifndef BAR_H
///Define this macro to prevent from including this header file more than once.
#define BAR_H

#include "type.h"
#include "endian.h"
#include "string_operation.h"
#include "file_operation.h"

#define DATA_DOUBLE 0
#define DATA_FLOAT 1
#define DATA_INT 2
#define DATA_SHORT 3
#define DATA_CHAR 4
#define DATA_UINT 5
#define DATA_USHORT 6
#define DATA_UCHAR 7

typedef union{
	double data_double;
	float data_float;
	int data_int;
	short data_short;
	char data_char;
	unsigned int data_uint;
	unsigned short data_ushort;
	unsigned char data_uchar;
} bar_column;

class bar_parameter{
public:
	uint len_name;
	string name;
	uint len_value;
	string value;
};

class bar_sequence{
public:
	uint len_seq_name;
	string name;
	uint len_group_name;
	string group_name;
	uint len_version_num;
	string version_num;
	uint num_parameter;
	vector <bar_parameter> parameters;
	uint num_data;
	vector<vector<bar_column> > data;
};

class bar{
public:
	float version;
	bool verbose;
	uint num_seq;
	vector<bar_sequence> sequences;
	uint num_col;
	vector <int> data_type;
	int size_data;
	uint num_parameter;
	vector <bar_parameter> parameters;

	bar();
	bool read_from_text_file(string filename, bool sorted = false);
	bool read_from_file(string file_name);
	bool read_from_file_header(FILE *file);
	bool read_from_file_parameter(FILE *file, bar_parameter &par);
	bool read_from_file_sequence_info(FILE *file, bar_sequence &seq);
	bool read_from_file_data(FILE *file, vector<bar_column> &data);
	bool read_from_file_region(string file_name, string chr, int start, int end, vector<vector<bar_column> > &data);
	bool write_to_file(string file_name);
	bool write_to_text_file(string file_name);
	bool dump_data(FILE *file, string chr, vector<vector<bar_column> > &data);
};

inline bool is_bar_file(string file_name){
	if (!file_exists(file_name)) return false;

	FILE *file=fopen(file_name.c_str(), "rb");
	if (!file){
//		if (verbose) cout << "error opening bar file:" << file_name << endl;
		return false;
	}
	size_t num_read;
	char buf[65536];

//	cout << "reading bar file header.\n";
//read magic num
	if (8 != (num_read = big_endian_fread(buf, 1, 8, file))) {
//		cout << "error reading bar file: no magic\n";
		fclose(file);
		return false;
	}
	buf[8] = 0;
	if (strcmp(buf, "barr\r\n\032\n")) {
//		cout << "error reading bar file: bad magic\n";
		fclose(file);
		return false;
	}
	fclose(file);
	return true;
}

inline bool is_bar_text_file(string file_name){
	if (!file_exists(file_name)) return false;
	ifstream ifs(file_name.c_str());
	string readline;

	int line_num = 0;
	while(getline(ifs, readline)) {
		if (readline == "") continue;
		line_num++;
		vector<string> tokens = string_tokenize(readline, " \t");
		if (tokens.size() >= 2 && is_int(tokens[1])) break;
		if (line_num == 2) {
			ifs.close();
			return false;
		}
	}

	ifs.close();
	return true;
}

inline bool convert_from_text_to_bar(const string text_file_name, const string bar_file_name){
	bar mybar;
	if (!mybar.read_from_text_file(text_file_name, true)) {
		if (!mybar.read_from_text_file(text_file_name)) {
			cout << "error reading text file " << text_file_name << endl;
			return false;
		}
	}
	if (!mybar.write_to_file(bar_file_name)) {
		cout << "error writing bar file " << bar_file_name << endl;
		return false;
	}
	return true;
}

inline bar::bar(){
	verbose = false;
}

inline bool bar::dump_data(FILE *file, string chr, vector<vector<bar_column> > &data){
	int j, k;
	for (j = 0; j < (int)data.size(); j++) {
		fprintf(file, "%s\t", chr.c_str());
		for (k = 0; k < (int)num_col; k++) {
			switch (data_type[k]) {
				case DATA_DOUBLE:
					fprintf(file, "%lf\t", data[j][k].data_double);
					break;
				case DATA_FLOAT:
					fprintf(file, "%lf\t", data[j][k].data_float);
					break;
				case DATA_INT:
					fprintf(file, "%d\t", data[j][k].data_int);
					break;
				case DATA_SHORT:
					fprintf(file, "%d\t", data[j][k].data_short);
					break;
				case DATA_CHAR:
					fprintf(file, "%d\t", data[j][k].data_char);
					break;
				case DATA_UINT:
					fprintf(file, "%d\t", data[j][k].data_uint);
					break;
				case DATA_USHORT:
					fprintf(file, "%d\t", data[j][k].data_ushort);
					break;
				case DATA_UCHAR:
					fprintf(file, "%d\t", data[j][k].data_uchar);
					break;
				default:				
					cout << "UNSUPPORTED DATA TYPE" << endl;
					return false;;
			}	
		}
		fprintf(file, "\n");
	}
	return true;
}

inline bool bar::write_to_text_file(string file_name){
	int i;
	FILE *file=fopen(file_name.c_str(), "wt");
	if (!file){
		cout << "error opening output file:" << file_name << endl;
		return false;
	}
	cout << "dumping to file " << file_name << endl;
	for (i = 0; i < (int)num_seq; i++) {
		if (!dump_data(file, sequences[i].name, sequences[i].data)) {
			cout << "error dumping data" << endl;
			return false;
		}
	}
	fclose(file);
	return true;
}

inline bool bar::read_from_file_header(FILE *file){
	size_t num_read;
	int i;
	char buf[65536];

	if (verbose) cout << "reading bar file header.\n";
//read magic num
	if (8 != (num_read = big_endian_fread(buf, 1, 8, file))) {
		cout << "error reading bar file: no magic\n";
		return false;
	}
	buf[8] = 0;
	if (strcmp(buf, "barr\r\n\032\n")) {
		cout << "error reading bar file: bad magic\n";
		return false;
	}
//read version
	if (1 != (num_read = big_endian_fread(&version, 4, 1, file))) {
		cout << "error reading bar file: no version\n";
		return false;
	}
	if (version != 1.0f && version != 2.0f) {
		cout << "error reading bar file: bad version\n";
		return false;
	}

	if (verbose) cout << "bar file version:" << (int) version << ".0" << endl;
//read number of sequences
	if (1 != (num_read = big_endian_fread(&num_seq, 4, 1, file))) {
		cout << "error reading bar file: # of sequences\n";
		return false;
	}
	if (verbose) cout << "# of sequences:" << num_seq << endl;
	sequences.resize(num_seq);
//read number of columns
	if (1 != (num_read = big_endian_fread(&num_col, 4, 1, file))) {
		cout << "error reading bar file: # of columns\n";
		return false;
	}
	if (verbose) cout << "# of columns:" << num_col << endl;
	data_type.resize(num_col);
//read data types
	if (verbose) cout << "reading data type description\n";
	size_data = 0;
	for (i = 0; i < (int)num_col; i++) {
		if (1 != (num_read = big_endian_fread(&data_type[i], 4, 1, file))) {
			cout << "error reading bar file: data type\n";
			return false;
		}
		switch (data_type[i]) {
			case DATA_DOUBLE:
				if (verbose) cout << "DOUBLE ";
				size_data += 8;
				break;
			case DATA_FLOAT:
				if (verbose) cout << "FLOAT ";
				size_data += 4;
				break;
			case DATA_INT:
				if (verbose) cout << "INT ";
				size_data += 4;
				break;
			case DATA_SHORT:
				if (verbose) cout << "SHORT ";
				size_data += 2;
				break;
			case DATA_CHAR:
				if (verbose) cout << "CHAR ";
				size_data += 1;
				break;
			case DATA_UINT:
				if (verbose) cout << "UINT ";
				size_data += 4;
				break;
			case DATA_USHORT:
				if (verbose) cout << "USHORT ";
				size_data += 2;
				break;
			case DATA_UCHAR:
				if (verbose) cout << "UCHAR ";
				size_data += 1;
				break;
			default:				
				cout << "UNSUPPORTED DATA TYPE" << endl;
				return false;;
		}
	}
	if (verbose) cout << endl;

//read number of parameters
	if (1 != (num_read = big_endian_fread(&num_parameter, 4, 1, file))) {
		cout << "error reading bar file: # of parameters " << endl;
		return false;
	}
	if (verbose) cout << "# of parameters:" << num_parameter << endl;
	parameters.resize(num_parameter);
//read parameters
	for (i = 0; i < (int)num_parameter; i++) {				
		if (!read_from_file_parameter(file, parameters[i])) {
			cout << "error reading bar file: bad parameter" << endl;
			return false;
		}
	}
	return true;
}

inline bool bar::read_from_file_parameter(FILE *file, bar_parameter &par){
	size_t num_read;
	char buf[65536];
//read parameter name
	if (1 != (num_read = big_endian_fread(&par.len_name, 4, 1, file))) {
		cout << "error reading bar file: length of name" << endl;
		return false;
	}
	if (par.len_name != (num_read = big_endian_fread(buf, 1, par.len_name, file))) {
		cout << "error reading bar file: parameter name " << endl;
		return false;
	}
	buf[par.len_name] = 0;
	par.name = buf;
//read parameter value;
	if (1 != (num_read = big_endian_fread(&par.len_value, 4, 1, file))) {
		cout << "error reading bar file: length of value" << endl;
		return false;
	}
	if (par.len_value != (num_read = big_endian_fread(buf, 1, par.len_value, file))) {
		cout << "error reading bar file: parameter value" << endl;
		return false;
	}
	buf[par.len_value] = 0;
	par.value = buf;
	if (verbose) cout << par.name << " = " << par.value << endl;
	return true;
}

inline bool bar::read_from_file_sequence_info(FILE *file, bar_sequence &seq){
	size_t num_read;
	int i;
	char buf[65536];
//read sequence name
	if (1 != (num_read = big_endian_fread(&seq.len_seq_name, 4, 1, file))) {
		cout << "error reading bar file: length of sequence name " << endl;
		return false;
	}
	if (seq.len_seq_name != (num_read = big_endian_fread(buf, 1, seq.len_seq_name, file))) {
		cout << "error reading bar file: sequence name " << endl;
		return false;
	}
	buf[seq.len_seq_name] = 0;
	seq.name = buf;
	if (verbose) cout << "name:" << buf << endl;
	if (version >= 2.0f) {
//read group name
		if (1 != (num_read = big_endian_fread(&seq.len_group_name, 4, 1, file))) {
			cout << "error reading bar file: length of group name for sequence " << endl;
			return false;
		}
		if (seq.len_group_name != (num_read = big_endian_fread(buf, 1, seq.len_group_name, file))) {
			cout << "error reading bar file: sequence group name " << endl;
			return false;
		}
		buf[seq.len_group_name] = 0;
		seq.group_name = buf;
		if (verbose) cout << "group name:" << buf << endl;
	}
//read version number
	if (1 != (num_read = big_endian_fread(&seq.len_version_num, 4, 1, file))) {
		cout << "error reading bar file: length of version number" << endl;
		return false;
	}
	if (seq.len_version_num != (num_read = big_endian_fread(buf, 1, seq.len_version_num, file))) {
		cout << "error reading bar file: sequence version number " << endl;
		return false;
	}
	buf[seq.len_version_num] = 0;
	seq.version_num = buf;
	if (verbose) cout << "version number:" << buf << endl;
	if (version >= 2.0f) {
//read number of parameters
		if (1 != (num_read = big_endian_fread(&seq.num_parameter, 4, 1, file))) {
			cout << "error reading bar file: # of parameters" << endl;
			return false;
		}
		if (verbose) cout << "# of parameters:" << seq.num_parameter << endl;
		seq.parameters.resize(seq.num_parameter);
//read parameters
		for (i = 0; i < (int)seq.num_parameter; i++) {				
			if (!read_from_file_parameter(file, seq.parameters[i])) {
				cout << "error reading bar file: bad parameter" << endl;
				return false;
			}
		}
	}
//read number of data points
	if (1 != (num_read = big_endian_fread(&seq.num_data, 4, 1, file))) {
		cout << "error reading bar file: # of data points" << endl;
		return false;
	}
	if (verbose) cout << "# of data points:" << seq.num_data << endl;
	return true;
}

inline bool bar::read_from_file_data(FILE *file, vector<bar_column> &data){
	size_t num_read;
	int k;
	data.resize(num_col);
	for (k = 0; k < (int) num_col; k++) {
		switch (data_type[k]) {
			case DATA_DOUBLE:
				if (1 != (num_read = big_endian_fread(&data[k].data_double, 8, 1, file))) {
					cout << "error reading bar file: bad data" << endl;
					return false;
				}
				break;
			case DATA_FLOAT:
				if (1 != (num_read = big_endian_fread(&data[k].data_float, 4, 1, file))) {
					cout << "error reading bar file: bad data" << endl;
					return false;
				}
				break;
			case DATA_INT:
				if (1 != (num_read = big_endian_fread(&data[k].data_int, 4, 1, file))) {
					cout << "error reading bar file: bad data" << endl;
					return false;
				}
				break;
			case DATA_SHORT:
				if (1 != (num_read = big_endian_fread(&data[k].data_short, 2, 1, file))) {
					cout << "error reading bar file: bad data" << endl;
					return false;
				}
				break;
			case DATA_CHAR:
				if (1 != (num_read = big_endian_fread(&data[k].data_char, 1, 1, file))) {
					cout << "error reading bar file: bad data" << endl;
					return false;
				}
				break;
			case DATA_UINT:
				if (1 != (num_read = big_endian_fread(&data[k].data_uint, 4, 1, file))) {
					cout << "error reading bar file: bad data" << endl;
					return false;
				}
				break;
			case DATA_USHORT:
				if (1 != (num_read = big_endian_fread(&data[k].data_ushort, 2, 1, file))) {
					cout << "error reading bar file: bad data" << endl;
					return false;
				}
				break;
			case DATA_UCHAR:
				if (1 != (num_read = big_endian_fread(&data[k].data_uchar, 1, 1, file))) {
					cout << "error reading bar file: bad data" << endl;
					return false;
				}
				break;
			default:				
				cout << "UNSUPPORTED DATA TYPE" << endl;
				return false;;
		}
	}
	return true;
}

inline bool bar::read_from_file(string file_name){
	if (!is_bar_file(file_name)) {
		if (verbose) cout << file_name << " is not a vaid bar file.\n";
		return false;
	}

	int i, j;
	FILE *file=fopen(file_name.c_str(), "rb");
	if (!file){
		if (verbose) cout << "error opening bar file:" << file_name << endl;
		return false;
	}
	if (verbose) cout << "reading bar file:" << file_name << endl;

	if (!read_from_file_header(file)) {
		if (verbose) cout << "reading bar file header failed.\n";
		return false;
	}	

//read sequences
	if (verbose) cout << "reading sequences\n";
	for (i = 0; i < (int)num_seq; i++) {
		if (verbose && i==0) cout << "sequence " << i+1 << endl;
		if (!read_from_file_sequence_info(file, sequences[i])) {
			cout << "error reading bar file: bad sequence info" << endl;
			return false;
		}
		sequences[i].data.resize(sequences[i].num_data);
//read data points
		if (verbose) cout << "reading data points" << endl;
		for (j = 0; j < (int) sequences[i].num_data; j++) {
			if (!read_from_file_data(file, sequences[i].data[j])) {
				cout << "error reading bar file: bad data" << endl;
				return false;
			}
		}
	}

	fclose(file);
	return true;
}

inline bool bar::read_from_file_region(string file_name, string chr, int start, int end, vector<vector<bar_column> > &data){
	if (!is_bar_file(file_name)) {
		if (verbose) cout << file_name << " is not a vaid bar file.\n";
		return false;
	}

	FILE *file;
	int i, j;
	file=fopen(file_name.c_str(), "rb");
	if (!file){
		cout << "error opening bar file:" << file_name << endl;
		return false;
	}
	if (verbose) cout << "reading bar file:" << file_name << endl;

	if (!read_from_file_header(file)) {
		if (verbose) cout << "reading bar file header failed.\n";
		return false;
	}	
	if (data_type[0] != DATA_INT) {
		cout << "first field is not coordinate" << endl;
		return false;
	}

//read sequences
	if (verbose) cout << "reading sequences\n";
	for (i = 0; i < (int)num_seq; i++) {
		if (verbose && i==0) cout << "sequence " << i+1 << endl;
		if (!read_from_file_sequence_info(file, sequences[i])) {
			cout << "error reading bar file: bad sequence info" << endl;
			return false;
		}
		if (sequences[i].name != chr || (sequences[i].version_num.length() > 0 && tolower(sequences[i].version_num.substr(0,4)) != "ncbi"))  {
			if (0 != fseek(file, size_data * sequences[i].num_data, SEEK_CUR)){
				cout << "error fseek" << endl;
				return false;
			}
		} else {
			long pos = ftell(file);
			if (pos<0) {
				cout << "error ftell" << endl;
				return false;
			}
			int low = 0, high = sequences[i].num_data - 1, value = start, mid = 0;
			vector<bar_column> temp_data;
			temp_data.resize(num_col);
			while (low <= high) {
				mid = (low + high) / 2;
				if (fseek(file, pos+size_data*mid, SEEK_SET)) {
					cout << "error fseek" << endl;
					return false;
				}
				if (!read_from_file_data(file, temp_data)) {
					cout << "error read data" << endl;
					return false;
				}
				if (temp_data[0].data_int > value)
					high = mid - 1;
				else if (temp_data[0].data_int < value)
					low = mid + 1;
				else
					break;
			}
			j = mid;
			while (temp_data[0].data_int <= end) {
				if (temp_data[0].data_int >= start) data.push_back(temp_data);
				j++;
				if (j == (int)sequences[i].num_data) break;
				if (!read_from_file_data(file, temp_data)) {
					cout << "error read data" << endl;
					return false;
				}
			}
			break; //skip other sequences;
		}
	}

	fclose(file);
	return true;
}

class bar_data_struct{
public:
	string chr;
	int coord;
	vector<float> intensities;
};

inline bool operator < (const bar_data_struct &data1, const bar_data_struct &data2) {
	if (data1.chr < data2.chr) return true;
	if (data1.chr == data2.chr && data1.coord < data2.coord) return true;
	return false;
}

inline bool bar::read_from_text_file(string filename, bool sorted){
	if (!is_bar_text_file(filename)) {
		cout << filename << " is not a vaid bar text file.\n";
		return false;
	}

	printf("Loading text file %s.\n", filename.c_str());
	version = 1.0f;
	num_seq = 0;
	sequences.clear();
	num_parameter = 0;
	parameters.clear();

	vector<bar_data_struct> bar_data;
	ifstream ifs(filename.c_str());
	string readline;
	int i, j;
	int line_num = 0;
	int data_field_num = -1;

	string old_chr = "";
	int old_coord = -1;
	vector<string> tags;

	int total_data = 0;
	while (getline(ifs, readline)){	
		if (readline == "") continue;
		line_num++;
		vector<string> tokens = string_tokenize(readline, " \t");
		if (tokens.size() < 2 || !is_int(tokens[1])) {
			if (line_num == 1) continue;
			else {
				cout << "bad format in bar text file, bad coordinate.\n";
				return false;
			}
		}
		if (data_field_num == -1) {
			data_field_num = (int)tokens.size() - 2;
			num_col = data_field_num+1;
			if (num_col == 1) num_col = 2;
		}
		if (data_field_num != (int)tokens.size() - 2) {
				cout << "bad format in bar text file, bad number of data field.\n";
				return false;
		}
		for (i = 2; i < (int)tokens.size(); i++) {
			if (!is_num(tokens[i])) {
				cout << "bad format in bar text file, bad data.\n";
				return false;
			}
		}

		bar_data_struct data1;
		data1.chr = tokens[0];
		data1.coord = str2int(tokens[1]);
		if (data1.coord < 0) {
			cout << "bad format in bar text file, bad coordinate.\n";
			return false;
		}

		for (i = 2; i < (int)tokens.size(); i++) {
			data1.intensities.push_back((float)str2double(tokens[i]));
		}

		if (sorted) {
			if (data1.chr != old_chr) {	
				if (find(tags.begin(), tags.end(), data1.chr) != tags.end()) {
					cout << "file is not sorted.\n";
					return false;
				} else {
					num_seq++;
					bar_sequence seq;
					seq.name = data1.chr;
					seq.len_seq_name = (int)seq.name.length();
					seq.len_group_name = 0;
					seq.group_name = "";
					seq.len_version_num = 0;
					seq.version_num = "";
					seq.num_parameter = 0;
					seq.parameters.clear();
					seq.num_data = 0;
					sequences.push_back(seq);

					tags.push_back(data1.chr);
					old_chr = data1.chr;
					old_coord = data1.coord;
				}
			}

			if (old_coord > data1.coord) {
				cout << "file is not sorted.\n";
				return false;
			} 
			
			vector<bar_column> column;
			column.resize(num_col);

			column[0].data_int = data1.coord;
			for (int j = 0; j < data_field_num; j++) column[j + 1].data_float = data1.intensities[j];
			if (data_field_num == 0) column[1].data_float = 1.0f;
			sequences[num_seq - 1].data.push_back(column);
			sequences[num_seq - 1].num_data++;
			total_data++;
		} else {
			bar_data.push_back(data1);
			total_data++;
		}
	}

	printf("%d data units loaded.\n", total_data);

	data_type.clear();
	data_type.push_back(DATA_INT);
	for (i = 0; i < (int)num_col - 1; i++) data_type.push_back(DATA_FLOAT);
	size_data = 4 * num_col;

	if (!sorted) {
		sort(bar_data.begin(), bar_data.end());
		string current_chr = "";
		for (i = 0; i < (int)bar_data.size(); i++) {
			if (bar_data[i].chr != current_chr) {
				num_seq++;
				bar_sequence seq;
				seq.name = bar_data[i].chr;
				seq.len_seq_name = (int)seq.name.length();
				seq.len_group_name = 0;
				seq.group_name = "";
				seq.len_version_num = 0;
				seq.version_num = "";
				seq.num_parameter = 0;
				seq.parameters.clear();
				seq.num_data = 0;
				sequences.push_back(seq);
				current_chr = bar_data[i].chr;
			} 
			vector<bar_column> column;
			column.resize(num_col);
			column[0].data_int = bar_data[i].coord;
			for (j = 0; j < data_field_num; j++) column[j + 1].data_float = bar_data[i].intensities[j];
			if (data_field_num == 0) column[1].data_float = 1.0f;
			sequences[num_seq - 1].data.push_back(column);
			sequences[num_seq - 1].num_data++;
		}
	}
	printf("%d sequences loaded.\n", num_seq);
	printf("Successfully loaded input file %s.\n", filename.c_str());
	return true;
}

inline bool bar::write_to_file(string file_name){
	int i, j, k;
	FILE *file=fopen(file_name.c_str(), "wb");
	if (!file){
		cout << "error opening bar file:" << file_name << endl;
		return false;
	}
	if (verbose) cout << "writing bar file:" << file_name << endl;
	big_endian_fwrite("barr\r\n\032\n", 1, 8, file);
	big_endian_fwrite(&version, 4, 1, file);
	big_endian_fwrite(&num_seq, 4, 1, file);
	big_endian_fwrite(&num_col, 4, 1, file);
	for (i = 0; i < (int)num_col; i++) {
		big_endian_fwrite(&data_type[i], 4, 1, file);
	}
	big_endian_fwrite(&num_parameter, 4, 1, file);
	for (i = 0; i < (int)num_parameter; i++) {
		big_endian_fwrite(&parameters[i].len_name, 4, 1, file);
		big_endian_fwrite(parameters[i].name.c_str(), 1, parameters[i].len_name, file);
		big_endian_fwrite(&parameters[i].len_value, 4, 1, file);
		big_endian_fwrite(parameters[i].value.c_str(), 1, parameters[i].len_value, file);
	}
	for (i = 0; i < (int)num_seq; i++) {
		big_endian_fwrite(&sequences[i].len_seq_name, 4, 1, file);
		big_endian_fwrite(sequences[i].name.c_str(), 1, sequences[i].len_seq_name, file);
		if (version >= 2.0f) {
			big_endian_fwrite(&sequences[i].len_group_name, 4, 1, file);
			big_endian_fwrite(sequences[i].group_name.c_str(), 1, sequences[i].len_group_name, file);
		}
		big_endian_fwrite(&sequences[i].len_version_num, 4, 1, file);
		big_endian_fwrite(sequences[i].version_num.c_str(), 1, sequences[i].len_version_num, file);
		if (version >= 2.0f) {
			big_endian_fwrite(&sequences[i].num_parameter, 4, 1, file);
			for (j = 0; j < (int)num_parameter; j++) {
				big_endian_fwrite(&sequences[i].parameters[j].len_name, 4, 1, file);
				big_endian_fwrite(sequences[i].parameters[j].name.c_str(), 1, sequences[i].parameters[j].len_name, file);
				big_endian_fwrite(&sequences[i].parameters[j].len_value, 4, 1, file);
				big_endian_fwrite(sequences[i].parameters[j].value.c_str(), 1, sequences[i].parameters[j].len_value, file);
			}
		}
		big_endian_fwrite(&sequences[i].num_data, 4, 1, file);
		for (j = 0; j < (int)sequences[i].num_data; j++) {
			for (k = 0; k < (int)num_col; k++) {
				switch (data_type[k]) {
					case DATA_DOUBLE:
						big_endian_fwrite(&sequences[i].data[j][k].data_double, 8, 1, file);
						break;
					case DATA_FLOAT:
						big_endian_fwrite(&sequences[i].data[j][k].data_float, 4, 1, file);
						break;
					case DATA_INT:
						big_endian_fwrite(&sequences[i].data[j][k].data_int, 4, 1, file);
						break;
					case DATA_SHORT:
						big_endian_fwrite(&sequences[i].data[j][k].data_short, 2, 1, file);
						break;
					case DATA_CHAR:
						big_endian_fwrite(&sequences[i].data[j][k].data_char, 1, 1, file);
						break;
					case DATA_UINT:
						big_endian_fwrite(&sequences[i].data[j][k].data_uint, 4, 1, file);
						break;
					case DATA_USHORT:
						big_endian_fwrite(&sequences[i].data[j][k].data_ushort, 2, 1, file);
						break;
					case DATA_UCHAR:
						big_endian_fwrite(&sequences[i].data[j][k].data_uchar, 1, 1, file);
						break;
					default:				
						cout << "UNSUPPORTED DATA TYPE" << endl;
						return false;;
				}
			}
		}
	}

	fclose(file);

	return true;
}

#endif //BAR_H
