#ifndef CHR_REGION_H
#define CHR_REGION_H
#include "string_operation.h"

inline bool parse_region(string region, string &chr, int &start, int &end){
	region = substitute(region, ",", "");
	region = substitute(region, ":", " ");
	region = substitute(region, "-", " ");

	char buf[1024];

	if (3 != sscanf(region.c_str(), " %s %d %d", buf, &start, &end)) {
		//cout << "error parsing region information" << endl;
		return false;
	}
	chr = buf;
	if (tolower(chr.substr(0,3)) == "chr") chr = "chr" + chr.substr(3);
	if (start <= 0 || start > end) {
		//cout << "bad region" << endl;
		return false;
	}
	return true;
}

inline bool is_region(string region) {
	string chr;
	int start, end;
	return parse_region(region, chr, start, end);
}

inline int chr2num(string chr, string species = "human"){
	if (chr.length() <= 3 || chr.length() > 5 || chr.substr(0, 3) != "chr") {
		//cout << "bad chr information" << endl;
		return -1;
	}
	chr = chr.substr(3);
	species = tolower(species);
	int num_autosome = 0;
	if (species == "mouse") {
		num_autosome = 19;
	} else if (species == "human") {
		num_autosome = 22;
	} else {
		//cout << "unknown species" << endl;
		return -1;
	}
	if (chr == "X") {
		return num_autosome + 1;
	} else if (chr == "Y") {
		return num_autosome + 2;
	} else {
		int temp;
		if (1 != sscanf(chr.c_str(), "%d", &temp)) {
			//cout << "bad chr information" << endl;
			return -1;
		}
		if (temp <= 0 || temp > num_autosome) {
			//cout << "bad chr information" << endl;
			return -1;
		}
		return temp;
	}
}

inline string num2chr(int num, string species = "human") {
	int num_autosome = 0;
	if (species == "mouse") {
		num_autosome = 19;
	} else if (species == "human") {
		num_autosome = 22;
	} else {
		//cout << "unknown species" << endl;
		return "";
	}
	if (num == num_autosome + 2) {
		return "chrY";
	} else if (num == num_autosome + 1) {
		return "chrX";
	} else if (num > 0 && num <= num_autosome) {
		char buf[1024];
		sprintf(buf, "chr%d", num);
		return string(buf);
	} else {
		//cout << "bad chr information" << endl;
		return "";
	}
}

class chr_region{
public:
	string chr;
	int start, end;
	chr_region(){};
	chr_region(string chrom, int startpos, int endpos){
		chr = chrom;
		start = startpos; 
		end = endpos;
	};
	chr_region(string region){
		parse_region(region, chr, start, end);
	};
	string get_region(){
		return chr+":"+int2str(start)+"-"+int2str(end);
	}
	int length(){
		return end - start + 1;
	}
	void correct(){
		if (start <= 0) move(1-start);
		if (end <= start) end = start + 1;
	}
	void move(int dist){
		start += dist; end += dist;
	}
	void move(double scale){
		move(int(scale*length()));
	}
	void resize(int length){
		int mid = (start + end)/2; 
		start = mid - length/2; 
		end = mid + length/2;
	}
	void resize(double scale){
		resize(int(length()*scale));
	}
};

#endif //CHR_REGION_H
