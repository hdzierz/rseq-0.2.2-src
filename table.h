/*
table.h - Header file for table and table file manipulation
written by Hui Jiang
Department of Biostatistics, University of Michigan
July, 2012 -
*/

#ifndef TABLE_H
///Define this macro to prevent from including this header file more than once.
#define TABLE_H

#include "stl.h"
#include "file_operation.h"
#include "string_operation.h"

class table{
public:
	bool verbose;
	string separator;
	int num_skip_header_lines;
	vector<vector<string> > data;
	map<string, vector<pair<int, int> > > content_map;

	table(){
		separator = "\t";
		num_skip_header_lines = 0;
		verbose = true;
		clear();
	};
	void clear(){data.clear();};
	void build_content_map() {
		content_map.clear();
		for (int i = 0; i < (int)data.size(); i++) {
			for (int j = 0; j < (int)data[i].size(); j++) {
				if (data[i][j] == "") continue;
				if (content_map.count(data[i][j]) == 0) {
					vector<pair<int, int> > temp;
					temp.push_back(pair<int, int>(i, j));
					content_map[data[i][j]] = temp;
				} else {
					content_map[data[i][j]].push_back(pair<int, int>(i, j));
				}
			}
		}
	}
	bool is_table_file(const string file_name);
	bool read_from_file(const string filename);
	bool write_to_file(const string filename);
};

inline bool table::is_table_file(const string file_name){
	if (!file_exists(file_name)) return false;
	ifstream ifs(file_name.c_str());
	string readline;
	bool result = true;

	for (int i = 0; i < num_skip_header_lines; i++) {
		if (!getline(ifs, readline)) {
			result = false;
			break;
		}
	}

	if (!getline(ifs, readline)) result = false;
	vector<string> tokens1 = string_tokenize(readline, separator, false);
	if (!getline(ifs, readline)) result = false;
	vector<string> tokens2 = string_tokenize(readline, separator, false);
	if (tokens1.size() == 0 || tokens1.size() != tokens2.size()) result = false;

	ifs.close();
	return result;
}

inline bool table::read_from_file(const string filename){
	if (!file_exists(filename)) {
		if (verbose) cout << filename << " does not exist.\n";
		return false;
	}
	if (!is_table_file(filename)) {
		if (verbose) cout << filename << " is not a vaid table file.\n";
		return false;
	}

	if (verbose) printf("Loading table file %s.\n", filename.c_str());

	data.clear();
	string readline;
	ifstream ifs(filename.c_str());
	int total_rows = 0;
	int num_cols = 0;
	while(getline(ifs, readline)) {
		total_rows++;
		if (total_rows <= num_skip_header_lines) continue;
		vector<string> tokens = string_tokenize(readline, separator, false);
		if (total_rows == num_skip_header_lines + 1) {
			num_cols = (int)tokens.size();
			if (num_cols == 0) {
				if (verbose) printf("warning: no content at line %d\n", total_rows);
			}
		}
		if (num_cols != (int)tokens.size()) {
			if (verbose) printf("warning: unequal number of columns at line %d : %s\n", total_rows, readline.c_str());
		}
		data.push_back(tokens);
	}
	ifs.close();
	if (verbose) printf("%d rows read, %d columns each row.\n", total_rows - num_skip_header_lines, (int)num_cols);
	return true;
}

inline bool table::write_to_file(const string filename) {
	FILE *file=fopen(filename.c_str(), "wt");
	if (!file){
		cout << "error opening file for writing:" << filename << endl;
		return false;
	}
	if (verbose) printf("writing table file %s.\n", filename.c_str());
	for (int i = 0; i < (int)data.size(); i++) {
		fprintf(file, "%s\n", join_strings(data[i], separator).c_str());
	}
	if (verbose) printf("%d rows wrote.\n", (int)data.size());
	return true;
}

#endif //TABLE_H
