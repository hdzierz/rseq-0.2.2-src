/*
genefile.h - Header file for gene file manipulation
written by JIANG Hui, 
Institute for Computational and Mathematical Engineering, Stanford University
May, 2007 -
*/

#ifndef GENEFILE_H
///Define this macro to prevent from including this header file more than once.
#define GENEFILE_H

#include "type.h"
#include "endian.h"
#include "chr_region.h"
#include "file_operation.h"
#include "math_utils.h"

/*
gene file format: all encoded in little endian

file header:
1)magic number:"genf\r\n\032\n":8 bytes
2)gene file version:float:4 bytes, version == 1.0f for refFlat format, version == 2.0f for BED format
3)length of genome version name:int:4 bytes
4)genome version name
5)number of sequences:int:4 bytes

sequence information:(xN)
1)length of sequence name:int:4 bytes
2)sequence name
3)number of genes:int:4 bytes
4)gene information:(xN)
1)txStart:int:4 bytes
2)txEnd:int:4 bytes
3)offset:int:4 bytes

detailed gene information:(xN)
1)length of geneName:int:4 bytes
2)geneName
3)length of name:int:4 bytes
4)name
5)length of chrom:int:4 bytes
6)chrom
7)strand:char:1 byte
8)txStart:int:4 bytes
9)txEnd:int:4 bytes
10)cdsStart:int:4 bytes
11)cdsEnd:int:4 bytes
12)exonCount:int:4 bytes
13)exon information(xN)
exonStart:int:4 bytes
exonEnd:int:4 bytes
14)score: int 4 bytes, only for version 2.0
15)color: uchar 1x3=3 bytes, only for version 2.0
*/

class gene_struct{
public:
	string  geneName;           //"Name of gene as it appears in Genome Browser."
	string  name;               //"Name of gene, or geneID, suppose to be unique"
	string  chrom;              //"Chromosome name"
	bool strand;             //"+ or - for strand"
	int    txStart;            //"Transcription start position"
	int    txEnd;              //"Transcription end position"
	int    cdsStart;           //"Coding region start"
	int    cdsEnd;             //"Coding region end"
	vector<pair<int, int> > exons; //"Exon positions"
	int		offset;
//only for version 2.0
	int score; //score, integer 0-1000 for bed files
	uchar colorR, colorG, colorB; //color in RGB
	
	gene_struct() {
		geneName = name = chrom = "";
		strand = true;
		txStart = txEnd = cdsStart = cdsEnd = 0;
		exons.clear();
		offset = score = 0;
		colorR = colorG = colorB = 0;
	};
};

inline bool operator < (const gene_struct &gene1, const gene_struct &gene2) {
	if (gene1.chrom < gene2.chrom) return true;
	if (gene1.chrom == gene2.chrom && gene1.txStart < gene2.txStart) return true;
	if (gene1.chrom == gene2.chrom && gene1.txStart == gene2.txStart && gene1.txEnd < gene2.txEnd) return true;
	return false;
}

class gene_seq{
public:
	string name;
	uint num_genes;
	vector<gene_struct> genes;
	interval_list<pair<int, int> > exon_intervals;
	interval_list<int> gene_intervals;

	gene_seq(){
		name = "";
		num_genes = 0;
		genes.clear();
	};
};

class genefile{
public:
	float version;
	bool verbose;
	string version_name;
	uint num_seq;
	vector<gene_seq> sequences;
	genefile() {
		version = 1.0f;
		verbose = false;
		version_name = "";
		num_seq = 0;
		sequences.clear();
	}

	bool read_from_file(string filename);
	bool read_from_file_header(FILE *file);
	bool read_from_file_sequence_info(FILE *file, gene_seq &seq);
	bool read_from_file_sequence_gene_info(FILE *file, gene_seq &seq);
	bool read_from_file_gene(FILE *file, gene_struct &gene);
	bool read_from_file_string(FILE *file, string &str);
	bool read_from_file_gene_name(string file_name, const string gene_name, vector<gene_struct> &genes);
	bool read_from_file_region(string file_name, string chr, int start, int end, vector<gene_struct> &genes);
	bool read_from_text_file(string filename);
	bool read_from_GTF_file(string filename);
	bool read_from_refFlat_sorted(string filename, string chr, int startpos, int endpos, vector<gene_struct> &genes);
	bool write_to_file(string filename);
	bool write_to_file_string(FILE *file, string &str);
	bool write_to_text_file(string filename);
	bool write_to_text_file(FILE *file, vector<gene_struct> &genes);
	void set_offsets(vector<uint> &offsets);
	void get_offsets(vector<uint> &offsets);
	void calc_offsets(vector<uint> &offsets);

	map<string, pair<int, int> > gene_map;
	void build_map();
	void search_by_name(const string &name, string &seq_name, gene_struct &gene);

	void build_interval_lists();
	void search_interval(const string chr, const int start, const int end, vector<pair<int, pair<int, int> > > &exons);
	void search_interval(const string chr, const int start, const int end, gene_struct *&gene, int &exon);

	bool search_by_name(const string gene_name, vector<gene_struct> &genes);
	bool search_by_region(string chr, int start, int end, vector<gene_struct> &genes);

	void trim_UTR();
};

inline bool is_genefile(string file_name){
	if (!file_exists(file_name)) return false;

	FILE *file=fopen(file_name.c_str(), "rb");
	if (!file){
		return false;
	}
	size_t num_read;
	char buf[65536];

//	cout << "reading gene file header.\n";
//read magic num
	if (8 != (num_read = little_endian_fread(buf, 1, 8, file))) {
		cout << "error reading gene file: no magic\n";
		fclose(file);
		return false;
	}
	buf[8] = 0;
	if (strcmp(buf, "genf\r\n\032\n")) {
//		cout << "error reading gene file: bad magic\n";
		fclose(file);
		return false;
	}
	fclose(file);
	return true;
}

//1: refFlat, 2: refFlat without geneName, 3:BED, 4:GTF
inline int is_gene_text_file(string file_name){
	if (!file_exists(file_name)) return 0;
	ifstream ifs(file_name.c_str());
	string readline;

	int result = 0;
	int line_num = 0;
	while(getline(ifs, readline)) {
		if (readline == "") continue;
		vector<string> tokens = string_tokenize(readline, "\t");
		if (tokens.size() == 1 && tokens[0].length() >= 7 && tokens[0].substr(0, 7) == "browser") continue; //bed format
		if (tokens.size() == 1 && tokens[0].length() >= 5 && tokens[0].substr(0, 5) == "track") continue; //bed format
		if (tokens.size() == 1 && tokens[0].length() >= 1 && tokens[0].substr(0, 1) == "#") continue; //gtf format
		line_num++;
		if (tokens.size() == 11 && (tokens[3] == "+" || tokens[3] == "-")) { result = 1; break; } //refFlat format
		if (tokens.size() == 10 && (tokens[2] == "+" || tokens[2] == "-")) { result = 2; break; } //refFlat format without geneName
		if (tokens.size() >= 3 && is_int(tokens[1]) && is_int(tokens[2])) { result = 3; break; } //bed format
		if (tokens.size() == 9 && is_int(tokens[3]) && is_int(tokens[4])) { result = 4; break; } //gtf format
		if (line_num == 2) {
			ifs.close();
			return 0;
		}
	}

	ifs.close();
	return result;
}

inline bool convert_from_text_to_gene(const string text_file_name, const string genefile_name){
	genefile mygenefile;
	if (!mygenefile.read_from_text_file(text_file_name)) {
		cout << "error reading text file " << text_file_name << endl;
		return false;
	}
	if (!mygenefile.write_to_file(genefile_name)) {
		cout << "error writing genefile " << genefile_name << endl;
		return false;
	}
	return true;
}

inline bool genefile::read_from_file(string filename){
	if (!is_genefile(filename)) {
		if (verbose) cout << filename << " is not a vaid genefile.\n";
		return false;
	}

	int i, j;
	FILE *file=fopen(filename.c_str(), "rb");
	if (!file){
		cout << "error opening gene file:" << filename << endl;
		return false;
	}
	if (verbose) cout << "reading gene file:" << filename << endl;

	if (!read_from_file_header(file)) {
		if (verbose) cout << "reading filename file header failed.\n";
		return false;
	}	

	//read sequence informaiton
	if (verbose) cout << "reading sequence information\n";
	sequences.resize(num_seq);
	for (i = 0; i < (int)num_seq; i++) {
		if (verbose) cout << "sequence " << i+1 << endl;
		if (!read_from_file_sequence_info(file, sequences[i])) {
			cout << "error reading gene file: bad sequence info" << endl;
			return false;
		}
		sequences[i].genes.resize(sequences[i].num_genes);
		//skip gene information
		if (0 != fseek(file, sequences[i].num_genes * 12, SEEK_CUR)) {
			cout << "error fseek" << endl;
			return false;
		}
	}
	//read genes
	if (verbose) cout << "reading genes" << endl;
	for (i = 0; i < (int)num_seq; i++) {
		for (j = 0; j < (int)sequences[i].num_genes; j++) {
			if (!read_from_file_gene(file, sequences[i].genes[j])) {
				cout << "error reading gene file: bad gene" << endl;
				return false;
			}
		}
	}

	fclose(file);
	return true;
}

inline bool genefile::read_from_file_header(FILE *file){
	size_t num_read;
	char buf[65536];

	if (verbose) cout << "reading gene file header.\n";
	//read magic num
	if (8 != (num_read = little_endian_fread(buf, 1, 8, file))) {
		cout << "error reading gene file: no magic\n";
		return false;
	}
	buf[8] = 0;
	if (strcmp(buf, "genf\r\n\032\n")) {
		cout << "error reading gene file: bad magic\n";
		return false;
	}
	//read version
	if (1 != (num_read = little_endian_fread(&version, 4, 1, file))) {
		cout << "error reading gene file: no version\n";
		return false;
	}
	if (version != 1.0f && version != 2.0f) {
		cout << "error reading gene file: bad version\n";
		return false;
	}

	if (verbose) cout << "gene file version:" << (int) version << ".0" << endl;
	//read gene version name
	if (!read_from_file_string(file, version_name)) {
		cout << "error reading gene file: no version name" << endl;
		return false;
	}
	if (verbose) cout << "gene file version name:" << version_name << endl;
	//read number of sequences
	if (1 != (num_read = little_endian_fread(&num_seq, 4, 1, file))) {
		cout << "error reading gene file: # of sequences\n";
		return false;
	}
	if (verbose) cout << "# of sequences:" << num_seq << endl;
	return true;
}

inline bool genefile::read_from_file_sequence_info(FILE *file, gene_seq &seq){
	size_t num_read;
	if (!read_from_file_string(file, seq.name)) {
		cout << "error reading gene file: sequence name" << endl;
		return false;
	}
	if (1 != (num_read = little_endian_fread(&seq.num_genes, 4, 1, file))) {
		cout << "error reading gene file: number of genes" << endl;
		return false;
	}
	return true;
}

inline bool genefile::read_from_file_sequence_gene_info(FILE *file, gene_seq &seq){
	size_t num_read;
	seq.genes.resize(seq.num_genes);
	int* buf = new int[seq.num_genes * 12];
	if (seq.num_genes * 3 != (num_read = little_endian_fread(buf, 4, seq.num_genes * 3, file))) {
		cout << "error reading gene file: gene information" << endl;
		return false;
	}
	int i;
	for (i = 0; i < (int)seq.num_genes; i++) {
		seq.genes[i].txStart = buf[i*3];
		seq.genes[i].txEnd = buf[i*3+1];
		seq.genes[i].offset = buf[i*3+2];
	}
	delete[] buf;
	return true;
}

inline bool genefile::read_from_file_gene(FILE *file, gene_struct &gene){
	size_t num_read;
	if (!read_from_file_string(file, gene.geneName)) {
		cout << "error reading gene file: geneName" << endl;
		return false;
	}
	if (!read_from_file_string(file, gene.name)) {
		cout << "error reading gene file: gene name" << endl;
		return false;
	}
	if (!read_from_file_string(file, gene.chrom)) {
		cout << "error reading gene file: gene chromosome name" << endl;
		return false;
	}
	char c;
	if (1 != (num_read = little_endian_fread(&c, 1, 1, file))) {
		cout << "error reading gene file: no strand" << endl;
		return false;
	}
	if (c == '+') gene.strand = true;
	else if (c == '-') gene.strand = false;
	else {
		cout << "error reading gene file: bad strand" << endl;
		return false;
	}
	if (1 != (num_read = little_endian_fread(&gene.txStart, 4, 1, file))) {
		cout << "error reading gene file: txStart" << endl;
		return false;
	}
	if (1 != (num_read = little_endian_fread(&gene.txEnd, 4, 1, file))) {
		cout << "error reading gene file: txEnd" << endl;
		return false;
	}
	if (1 != (num_read = little_endian_fread(&gene.cdsStart, 4, 1, file))) {
		cout << "error reading gene file: cdsStart" << endl;
		return false;
	}
	if (1 != (num_read = little_endian_fread(&gene.cdsEnd, 4, 1, file))) {
		cout << "error reading gene file: cdsEnd" << endl;
		return false;
	}
	int exonCount;
	if (1 != (num_read = little_endian_fread(&exonCount, 4, 1, file))) {
		cout << "error reading gene file: exonCount" << endl;
		return false;
	}
	int i;
	for (i = 0; i < exonCount; i++) {
		int exonStart, exonEnd;
		if (1 != (num_read = little_endian_fread(&exonStart, 4, 1, file))) {
			cout << "error reading gene file: exonStart" << endl;
			return false;
		}
		if (1 != (num_read = little_endian_fread(&exonEnd, 4, 1, file))) {
			cout << "error reading gene file: exonEnd" << endl;
			return false;
		}
		gene.exons.push_back(pair<int, int>(exonStart, exonEnd));
	}
	if (version == 2.0f) {
		if (1 != (num_read = little_endian_fread(&gene.score, 4, 1, file))) {
			cout << "error reading gene file: score";
			return false;
		}
		if (1 != (num_read = little_endian_fread(&gene.colorR, 1, 1, file))) {
			cout << "error reading gene file: colorR" << endl;
			return false;
		}
		if (1 != (num_read = little_endian_fread(&gene.colorG, 1, 1, file))) {
			cout << "error reading gene file: colorG" << endl;
			return false;
		}
		if (1 != (num_read = little_endian_fread(&gene.colorB, 1, 1, file))) {
			cout << "error reading gene file: colorB" << endl;
			return false;
		}
	}
	return true;
}

inline bool genefile::read_from_file_string(FILE *file, string &str){
	size_t num_read;
	int len;
	char buf[65536];
	if (1 != (num_read = little_endian_fread(&len, 4, 1, file))) {
		cout << "error reading gene file: length of string" << endl;
		return false;
	}
	if (len > 1024) {
		cout << "error reading gene file: length too long" << endl;
		return false;
	}
	if (len != (int)(num_read = little_endian_fread(buf, 1, len, file))) {
		cout << "error reading bar file: parameter value" << endl;
		return false;
	}
	buf[len] = 0;
	str = buf;
	return true;
}

inline bool genefile::read_from_file_gene_name(string file_name, const string gene_name, vector<gene_struct> &genes){
	if (!read_from_file(file_name)) return false;
	return search_by_name(gene_name, genes);
}

inline bool genefile::read_from_file_region(string file_name, string chr, int start, int end, vector<gene_struct> &genes){
	if (!is_genefile(file_name)) {
		if (verbose) cout << file_name << " is not a vaid genefile.\n";
		return false;
	}

	int i, j;
	FILE *file=fopen(file_name.c_str(), "rb");
	if (!file){
		cout << "error opening gene file:" << file_name << endl;
		return false;
	}
	if (verbose) cout << "reading gene file:" << file_name << endl;

	if (!read_from_file_header(file)) {
		if (verbose) cout << "reading file header failed.\n";
		return false;
	}	

	//read sequence informaiton
	if (verbose) cout << "reading sequence information\n";
	sequences.resize(num_seq);
	for (i = 0; i < (int)num_seq; i++) {
		if (verbose) cout << "sequence " << i+1 << endl;
		if (!read_from_file_sequence_info(file, sequences[i])) {
			cout << "error reading gene file: bad sequence info" << endl;
			return false;
		}
		//		sequences[i].genes.resize(sequences[i].num_genes);
		if (sequences[i].name == chr) {
			//read gene informations
			if (!read_from_file_sequence_gene_info(file, sequences[i])) {
				cout << "error reading gene information" << endl;
				return false;
			}
			/*			for (j = 0; j < (int)sequences[i].num_genes; j++) {
			if (1 != (num_read = little_endian_fread(&sequences[i].genes[j].txStart, 4, 1, file))) {
			cout << "error reading gene file: txStart" << endl;
			return  false;
			}
			if (1 != (num_read = little_endian_fread(&sequences[i].genes[j].txEnd, 4, 1, file))) {
			cout << "error reading gene file: txEnd" << endl;
			return  false;
			}
			if (1 != (num_read = little_endian_fread(&sequences[i].genes[j].offset, 4, 1, file))) {
			cout << "error reading gene file: offset" << endl;
			return  false;
			}
			}*/
			//read genes
			if (verbose) cout << "reading genes" << endl;
			for (j = 0; j < (int)sequences[i].num_genes; j++) {
				if (sequences[i].genes[j].txStart > end || sequences[i].genes[j].txEnd < start) continue;
				if (0 != fseek(file, sequences[i].genes[j].offset, SEEK_SET)) {
					cout << "error fseek" << endl;
					return false;
				}
				gene_struct gene;
				if (!read_from_file_gene(file, gene)) {
					cout << "error reading gene file: bad gene" << endl;
					return false;
				}
				genes.push_back(gene);
			}
			break;
		} else {
			//skip gene information
			if (0 != fseek(file, sequences[i].num_genes * 12, SEEK_CUR)) {
				cout << "error fseek" << endl;
				return false;
			}
		}
	}

	fclose(file);
	return true;
}

class GTF_record{
public:
	string seqname;
	string source;
	string feature;
	int start;
	int end;
	string score;
	bool forward_strand;
	string phase;
	map<string, string> key_values;
};

inline bool genefile::read_from_GTF_file(string filename) {
	if (4 != is_gene_text_file(filename)) {
		if (verbose) cout << filename << " is not a vaid GTF file.\n";
		return false;
	}

	printf("Loading GTF file %s.\n", filename.c_str());
	ifstream ifs(filename.c_str());
	string readline;
	int line_num = 0;
	vector<GTF_record> records;
	vector<gene_struct> genes;
	map<string, int> gene_map;
	genes.clear();
	version_name = "";
	num_seq = 0;
	version = 1.0f;
	string current_chr = "";
	while (getline(ifs, readline)){
		line_num++;
		size_t pos = readline.find("#");
		if (pos != string::npos) readline = readline.substr(0, pos);
		if (readline == "") continue;
		vector<string> tokens = string_tokenize(readline, "\t");
		if (tokens.size() != 9) goto fail;
		for (int i = 0; i < 9; i++) {
			if (tokens[i] == "") goto fail;
		}
		GTF_record record;
		record.seqname = tokens[0];
		record.source = tokens[1];
		record.feature = tokens[2];
		if (!is_int(tokens[3]) || str2int(tokens[3]) < 1) goto fail;
		record.start = str2int(tokens[3]) - 1;
		if (!is_int(tokens[4]) || str2int(tokens[4]) < str2int(tokens[3])) goto fail;
		record.end = str2int(tokens[4]);
		if (tokens[5] != "." && !is_num(tokens[5])) goto fail;
		record.score = tokens[5];
		if (tokens[6] != "+" && tokens[6] != "-") goto fail;
		record.forward_strand = (tokens[6] == "+");
		if (tokens[7] != "." && tokens[7] != "0" && tokens[7] != "1" && tokens[7] != "2") goto fail;
		record.phase = tokens[7];
		vector<string> tokens1 = string_tokenize(tokens[8], ";");
		if (tokens1.size() < 2) goto fail;
		for (int i = 0; i < (int)tokens1.size(); i++) {
			vector<string> tokens2 = string_tokenize(tokens1[i], " ");
			if (tokens2.size() != 2) goto fail;
			string key = tokens2[0];
			if (key == "") goto fail;
			if (i == 0 && key != "gene_id") goto fail;
			//if (i == 1 && key != "transcript_id") goto fail;
			string value = tokens2[1];
			if (value.length() >= 2 && value[0] == '"' && value[value.length()-1] == '"') {
				value = value.substr(1, value.length() - 2);
			}
			if (value == "") goto fail;
			if (key == "gene_id" || key == "transcript_id" || key == "gene_name") { // to save memory
				record.key_values[key] = value;
			}
		}
		if (record.feature == "transcript" || record.feature == "exon") { // to save memory
			records.push_back(record);
		}
	}
	printf("Successfully loaded GTF file %s.\n", filename.c_str());
	printf("%d records loaded.\n", (int)records.size());

	for (int round = 1; round <= 2; round++) {
		for (int i = 0; i < (int)records.size(); i++) {
			GTF_record &record = records[i];
			if (record.feature == "transcript" && round == 1) {
				gene_struct gene;
				gene.name = record.key_values["transcript_id"];
				gene.geneName = "";
				if (record.key_values.count("gene_id") == 1) {
					gene.geneName = record.key_values["gene_id"];
				}
				if (record.key_values.count("gene_name") == 1) {
					if (gene.geneName == "") {
						gene.geneName = record.key_values["gene_name"];
					} else {
						gene.geneName = gene.geneName + "|" + record.key_values["gene_name"];
					}
				}
				if (gene.geneName == "") gene.geneName = gene.name;
				gene.chrom = record.seqname;
				gene.strand = record.forward_strand;
				gene.txStart = record.start;
				gene.txEnd = record.end;
				gene.cdsStart = record.start;
				gene.cdsEnd = record.end;
				gene.exons.clear();
				genes.push_back(gene);
				if (gene_map.count(gene.name) > 0) {
					printf("error, duplicated transcript %s. \n", gene.name.c_str());
					return false;
				}
				gene_map[gene.name] = genes.size() - 1;
			} else if (record.feature == "exon" && round == 2) {
				string transcript_id = record.key_values["transcript_id"];
				if (gene_map.count(transcript_id) == 0) {
					printf("error, transcript %s not found. \n", transcript_id.c_str());
					return false;
				}
				gene_struct &gene = genes[gene_map[transcript_id]];
				gene.exons.push_back(pair<int, int>(record.start, record.end));
			}
		}
	}
	printf("%d genes loaded.\n", (int)genes.size());
	printf("checking exon coordinates...");
	for (int i = 0; i < (int)genes.size(); i++) {
		sort(genes[i].exons.begin(), genes[i].exons.end());
		for (int j = 0; j < (int)genes[i].exons.size(); j++) {
			if (j == 0 && genes[i].exons[j].first != genes[i].txStart || j > 0 && genes[i].exons[j].first <= genes[i].exons[j-1].second 
					|| genes[i].exons[j].second <= genes[i].exons[j].first || j == (int)genes[i].exons.size() - 1 && genes[i].exons[j].second != genes[i].txEnd) {
				printf("exon coordinate error, transcript %s.\n", genes[i].name.c_str());
				return false;
			}
		}
	}
	printf("done.\n");

	sort(genes.begin(), genes.end());
	for (int i = 0; i < (int)genes.size(); i++) {
		if (genes[i].chrom != current_chr) {
			num_seq++;
			gene_seq seq;
			seq.name = genes[i].chrom;
			seq.num_genes = 0;			
			sequences.push_back(seq);
			current_chr = genes[i].chrom;
		} 
		sequences[num_seq - 1].genes.push_back(genes[i]);
		sequences[num_seq - 1].num_genes++;
	}
	printf("%d sequences loaded.\n", num_seq);
	return true;

fail:
	printf("error reading GTF file %s, line %d.\n", filename.c_str(), line_num);
	return false;
}

inline bool genefile::read_from_text_file(string filename){
	int format = is_gene_text_file(filename);
	if (!format) {
		if (verbose) cout << filename << " is not a vaid gene text file.\n";
		return false;
	}

	if (format == 1 || format == 2) {
		printf("Loading refFlat file %s.\n", filename.c_str());
		version = 1.0f;
	} else if (format == 3) {
		printf("Loading BED file %s.\n", filename.c_str());
		version = 2.0f;
	} else if (format == 4) {
		return read_from_GTF_file(filename);
	} else panic("internal error: wrong format.");

	version_name = "";
	num_seq = 0;
	string current_chr = "";
	ifstream ifs(filename.c_str());
	vector<gene_struct> genes;
	bool hasGeneName = (format == 1);
	int line_num = 0;
	while (!ifs.bad()){
		string readline;
		getline(ifs, readline);
		if (readline == "") break;
		line_num++;

		gene_struct gene;

		if (format == 1 || format == 2) {
			string buf1, buf2, buf3, buf4, buf5, chr;
			//		int int1;
			int int2, int3, int4, int5, int6;
			vector<string> tokens = string_tokenize(readline, "\t");
			if (hasGeneName) {
	//			if (11 != sscanf(readline.c_str(), "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s", buf1, buf2, chr, buf3, &int2, &int3, &int4, &int5, &int6, buf4, buf5)) {
				if (11 != tokens.size()) {
					cout << "error format in line " << line_num << endl;
					return false;
				}
				buf1 = tokens[0];
				buf2 = tokens[1];
				chr = tokens[2];
				buf3 = tokens[3];
				int2 = str2int(tokens[4]);
				int3 = str2int(tokens[5]);
				int4 = str2int(tokens[6]);
				int5 = str2int(tokens[7]);
				int6 = str2int(tokens[8]);
				buf4 = tokens[9];
				buf5 = tokens[10];
			} else {
	//			if (10 != sscanf(readline.c_str(), "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s", buf2, chr, buf3, &int2, &int3, &int4, &int5, &int6, buf4, buf5)) {
				if (10 != tokens.size()) {
					cout << "error format in line " << line_num << endl;
					return false;
				}
				buf2 = tokens[0];
				chr = tokens[1];
				buf3 = tokens[2];
				int2 = str2int(tokens[3]);
				int3 = str2int(tokens[4]);
				int4 = str2int(tokens[5]);
				int5 = str2int(tokens[6]);
				int6 = str2int(tokens[7]);
				buf4 = tokens[8];
				buf5 = tokens[9];
			}

			//		string chr = num2chr(int1);
			/*		if (chr != current_chr) {
			num_seq++;
			gene_seq seq;
			seq.name = chr;
			seq.num_genes = 0;			
			sequences.push_back(seq);
			current_chr = chr;
			}*/
			gene.name = buf2;
			if (hasGeneName) gene.geneName = buf1; else gene.geneName = gene.name;
			gene.chrom = chr;
			gene.strand = (buf3 == "+");
			gene.txStart = int2;
			gene.txEnd = int3;
			gene.cdsStart = int4;
			gene.cdsEnd = int5;
			string str1 = buf4, str2 = buf5;
			for (int i=0; i<int6; i++) {
				int temp1, temp2;
				sscanf(str1.c_str(), "%d", &temp1);
				sscanf(str2.c_str(), "%d", &temp2);
				gene.exons.push_back(pair<int, int>(temp1, temp2));
				size_t index;
				index = str1.find (",");
				if (index != string::npos ) 
					str1 = str1.substr(index + 1);
				index = str2.find (",");
				if (index != string::npos ) 
					str2 = str2.substr(index + 1);
			}
			//sequences[num_seq - 1].genes.push_back(gene);
			//sequences[num_seq - 1].num_genes++;

		} else { //format == 3
			vector<string> tokens = string_tokenize(readline, "\t");
			if (readline.substr(0,5) == "track") continue;
			if (tokens.size() >= 3) {
				gene.chrom = tokens[0];
				gene.txStart = str2int(tokens[1]);
				gene.txEnd = str2int(tokens[2]);
			} else {
				cout << "error format in line " << line_num << endl;
				return false;
			}
			if (tokens.size() >= 4) {
				gene.geneName = gene.name = tokens[3];
			} else gene.name = "";
			if (tokens.size() >= 5) {
				gene.score = str2int(tokens[4]);
			} else gene.score = 0;
			if (tokens.size() >= 6) {
				if (tokens[5] == "+") gene.strand = true;
				else if (tokens[5] == "-") gene.strand = false;
				else {
					cout << "error format in line " << line_num << endl;
					return false;
				}
			} else gene.strand = true;
			if (tokens.size() >= 7) {
				gene.cdsStart = str2int(tokens[6]);
			} else gene.cdsStart = gene.txStart;
			if (tokens.size() >= 8) {
				gene.cdsEnd = str2int(tokens[7]);
			} else gene.cdsEnd = gene.txEnd;
			if (tokens.size() >= 9) {
				vector<string> tokens1 = string_tokenize(tokens[8], ",");
				if (tokens1.size() == 3) {
					gene.colorR = uchar(str2int(tokens1[0]));
					gene.colorG = uchar(str2int(tokens1[1]));
					gene.colorB = uchar(str2int(tokens1[2]));
				} else if (tokens1.size() == 1 && tokens1[0] == "0") {
					gene.colorR = gene.colorG = gene.colorB = 0;
				} else {
					cout << "error format in line " << line_num << endl;
					return false;
				}
			} else gene.cdsEnd = gene.txEnd;
			if (tokens.size() == 9) {
				gene.exons.resize(1);
				gene.exons[0].first = gene.txStart;
				gene.exons[0].second = gene.txEnd;
			} else if (tokens.size() == 12) {
				int exon_count = str2int(tokens[9]);
				gene.exons.resize(exon_count);
				vector<string> tokens1 = string_tokenize(tokens[10], ",");
				vector<string> tokens2 = string_tokenize(tokens[11], ",");
				if (tokens1.size() != tokens2.size() || (int)tokens1.size() < exon_count || (int)tokens1.size() > exon_count + 1) {
					cout << "error format in line " << line_num << endl;
					return false;
				}
				for (int i = 0; i < exon_count; i++) {
					gene.exons[i].first = gene.txStart + str2int(tokens2[i]);
					gene.exons[i].second = gene.exons[i].first + str2int(tokens1[i]);
				}
			} else {
				cout << "error format in line " << line_num << endl;
				return false;
			}
		}
		genes.push_back(gene);
	}
	printf("%d genes loaded.\n", (int)genes.size());
	sort(genes.begin(), genes.end());
	int i;
	for (i = 0; i < (int)genes.size(); i++) {
		if (genes[i].chrom != current_chr) {
			num_seq++;
			gene_seq seq;
			seq.name = genes[i].chrom;
			seq.num_genes = 0;			
			sequences.push_back(seq);
			current_chr = genes[i].chrom;
		} 
		sequences[num_seq - 1].genes.push_back(genes[i]);
		sequences[num_seq - 1].num_genes++;
	}
	printf("%d sequences loaded.\n", num_seq);
	printf("Successfully loaded input file %s.\n", filename.c_str());
	return true;
}

inline bool genefile::read_from_refFlat_sorted(string filename, string chr, int startpos, int endpos, vector<gene_struct> &genes){
	if (1 != is_gene_text_file(filename)) {
		if (verbose) cout << filename << " is not a vaid refFlat_sorted text file.\n";
		return false;
	}

	printf("Loading input file %s.\n", filename.c_str());
	version = 1.0f;
	//fflush(stdout);	
	genes.clear();
	ifstream ifs(filename.c_str());
	while (!ifs.bad()){
		string readline;
		getline(ifs, readline);
		if (readline == "") break;
		char buf1[10240], buf2[10240], buf3[10240], buf4[10240], buf5[10240];
		int int1, int2, int3, int4, int5, int6;
		if (11 != sscanf(readline.c_str(), "%s %s %d %s %d %d %d %d %d %s %s", buf1, buf2, &int1, buf3, &int2, &int3, &int4, &int5, &int6, buf4, buf5)) {
			return false;
		}
		int nchr = chr2num(chr);
		if (nchr < 0) {
			cout << "bad chr information" << endl;
			return false;
		}
		if (int1 < nchr) continue;
		if (int1 > nchr) break;
		if (int3 < startpos) continue;
		if (int2 > endpos) break;
		gene_struct gene;
		gene.geneName = buf1;
		gene.name = buf2;
		gene.chrom = chr;
		gene.strand = (string(buf3) == "+");
		gene.txStart = int2;
		gene.txEnd = int3;
		gene.cdsStart = int4;
		gene.cdsEnd = int5;
		string str1 = buf4, str2 = buf5;
		for (int i=0; i<int6; i++) {
			int temp1, temp2;
			sscanf(str1.c_str(), "%d", &temp1);
			sscanf(str2.c_str(), "%d", &temp2);
			gene.exons.push_back(pair<int, int>(temp1, temp2));
			size_t index;
			index = str1.find (",");
			if (index != string::npos ) 
				str1 = str1.substr(index + 1);
			index = str2.find (",");
			if (index != string::npos ) 
				str2 = str2.substr(index + 1);
		}
		genes.push_back(gene);
	}
	printf("Successfully loaded input file %s.\n", filename.c_str());
	return true;
}

inline bool genefile::write_to_file(string filename){
	vector<uint> offsets;
	calc_offsets(offsets);
	set_offsets(offsets);

	int i, j, k;
	FILE *file=fopen(filename.c_str(), "wb");
	if (!file){
		cout << "error opening gene file:" << filename << endl;
		return false;
	}
	if (verbose) cout << "writing gene file:" << filename << endl;
	little_endian_fwrite("genf\r\n\032\n", 1, 8, file);
	little_endian_fwrite(&version, 4, 1, file);
	if (verbose) cout << "gene file version:" << (int) version << ".0" << endl;
	if (!write_to_file_string(file, version_name)) {
		cout << "error writing gene file: bad verison_name" << endl;
		return false;
	}
	little_endian_fwrite(&num_seq, 4, 1, file);
	if (verbose) cout << "# of sequences:" << num_seq << endl;
	if (verbose) cout << "writing sequence description\n";
	for (i = 0; i < (int)num_seq; i++) {
		if (!write_to_file_string(file, sequences[i].name)) {
			cout << "error writing gene file: bad sequence name" << endl;
			return false;
		}
		if (verbose) cout << "sequence name:" << sequences[i].name << endl;
		little_endian_fwrite(&sequences[i].num_genes, 4, 1, file);
		if (verbose) cout << "# of genes:" << sequences[i].num_genes << endl;
		for (j = 0; j < (int)sequences[i].num_genes; j++) {
			little_endian_fwrite(&sequences[i].genes[j].txStart, 4, 1, file);
			little_endian_fwrite(&sequences[i].genes[j].txEnd, 4, 1, file);
			little_endian_fwrite(&sequences[i].genes[j].offset, 4, 1, file);
		}
	}
	if (verbose) cout << "writing genes\n";
	for (i = 0; i < (int)num_seq; i++) {
		for (j = 0; j < (int)sequences[i].num_genes; j++) {
			if (!write_to_file_string(file, sequences[i].genes[j].geneName)) {
				cout << "error writing gene file: bad geneName" << endl;
				return false;
			}
			if (!write_to_file_string(file, sequences[i].genes[j].name)) {
				cout << "error writing gene file: bad gene name" << endl;
				return false;
			}
			if (!write_to_file_string(file, sequences[i].genes[j].chrom)) {
				cout << "error writing gene file: bad gene chromosome" << endl;
				return false;
			}
			if (sequences[i].genes[j].strand)
				little_endian_fwrite("+", 1, 1, file);
			else
				little_endian_fwrite("-", 1, 1, file);
			little_endian_fwrite(&sequences[i].genes[j].txStart, 4, 1, file);
			little_endian_fwrite(&sequences[i].genes[j].txEnd, 4, 1, file);
			little_endian_fwrite(&sequences[i].genes[j].cdsStart, 4, 1, file);
			little_endian_fwrite(&sequences[i].genes[j].cdsEnd, 4, 1, file);
			size_t temp = sequences[i].genes[j].exons.size();
			little_endian_fwrite(&temp, 4, 1, file);
			for (k = 0; k < (int)sequences[i].genes[j].exons.size(); k++) {
				little_endian_fwrite(&sequences[i].genes[j].exons[k].first, 4, 1, file);
				little_endian_fwrite(&sequences[i].genes[j].exons[k].second, 4, 1, file);
			}
			if (version == 2.0f) {
				little_endian_fwrite(&sequences[i].genes[j].score, 4, 1, file);
				little_endian_fwrite(&sequences[i].genes[j].colorR, 1, 1, file);
				little_endian_fwrite(&sequences[i].genes[j].colorG, 1, 1, file);
				little_endian_fwrite(&sequences[i].genes[j].colorB, 1, 1, file);
			}
		}
	}

	fclose(file);

	return true;
}

inline bool genefile::write_to_file_string(FILE *file, string &str){
	int len = (int)str.length();
	if (len > 1024) {
		cout << "error writeing gene file: length too long" << endl;
		return false;
	}
	little_endian_fwrite(&len, 4, 1, file);
	little_endian_fwrite(str.c_str(), 1, len, file);
	return true;
}

inline bool genefile::write_to_text_file(string filename){
	int i;
	FILE *file=fopen(filename.c_str(), "wt");
	if (!file){
		cout << "error opening gene file:" << filename << endl;
		return false;
	}
	for (i = 0; i < (int)num_seq; i++) {
		if (!write_to_text_file(file, sequences[i].genes)) {
			cout << "error write gene" << endl;
			return false;
		}
	}
	fclose(file);
	return true;
}

inline bool genefile::write_to_text_file(FILE *file, vector<gene_struct> &genes){
	int j, k;
	for (j = 0; j < (int)genes.size(); j++) {
		string c;
		if (genes[j].strand) c = "+"; else c = "-";

		if (version == 1.0f) {
			if (genes[j].geneName != "") fprintf(file, "%s\t", genes[j].geneName.c_str());
			fprintf(file, "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t", genes[j].name.c_str(), genes[j].chrom.c_str(), c.c_str(), genes[j].txStart, genes[j].txEnd, genes[j].cdsStart, genes[j].cdsEnd, (int)genes[j].exons.size());
			for (k = 0; k < (int)genes[j].exons.size(); k++) {
				fprintf(file, "%d,", genes[j].exons[k].first);
			}
			fprintf(file, "\t");
			for (k = 0; k < (int)genes[j].exons.size(); k++) {
				fprintf(file, "%d,", genes[j].exons[k].second);
			}
		} else if (version == 2.0f) {
			fprintf(file, "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%d,%d,%d\t%d\t", genes[j].chrom.c_str(), genes[j].txStart, genes[j].txEnd, genes[j].name.c_str(), genes[j].score, c.c_str(), genes[j].cdsStart, genes[j].cdsEnd, (int)genes[j].colorR, (int)genes[j].colorG, (int)genes[j].colorB, (int)genes[j].exons.size());
			for (k = 0; k < (int)genes[j].exons.size(); k++) {
				fprintf(file, "%d,", genes[j].exons[k].second - genes[j].exons[k].first);
			}
			fprintf(file, "\t");
			for (k = 0; k < (int)genes[j].exons.size(); k++) {
				fprintf(file, "%d,", genes[j].exons[k].first - genes[j].txStart);
			}
		} else {
			panic("internal error: wrong version.");
		}

		fprintf(file, "\n");
	}
	return true;
}

inline void genefile::get_offsets(vector<uint> &offsets){
	int i, j;
	for (i = 0; i < (int)num_seq; i++) {
		for (j = 0; j < (int)sequences[i].num_genes; j++) {
			offsets.push_back(sequences[i].genes[j].offset);
		}
	}
}

inline void genefile::set_offsets(vector<uint> &offsets){
	int i, j, k = 0;
	for (i = 0; i < (int)num_seq; i++) {
		for (j = 0; j < (int)sequences[i].num_genes; j++) {
			sequences[i].genes[j].offset = offsets[k];
			k++;
		}
	}
}

inline void genefile::calc_offsets(vector<uint> &offsets){
	uint curr_offset = 0;
	curr_offset += 8; //magic
	curr_offset += 4; //version
	curr_offset += 4; //length of genome version name
	curr_offset += (uint)version_name.length(); //version_name
	curr_offset += 4; //num_seq
	int i, j;
	for (i = 0; i < (int)num_seq; i++) {
		curr_offset += 4; //len_seq_name
		curr_offset += (uint)sequences[i].name.length(); //name
		curr_offset += 4; // num_gene
		curr_offset += (uint)sequences[i].num_genes * 12; //gene info
	}
	for (i = 0; i < (int)num_seq; i++) {
		for (j = 0; j < (int)sequences[i].num_genes; j++) {
			offsets.push_back(curr_offset);
			curr_offset += 4; //length of geneName
			curr_offset += (uint)sequences[i].genes[j].geneName.length(); //geneName
			curr_offset += 4; //length of name
			curr_offset += (uint)sequences[i].genes[j].name.length(); //name
			curr_offset += 4; //length of chrom
			curr_offset += (uint)sequences[i].genes[j].chrom.length(); //chrom
			curr_offset += 1; //strand
			curr_offset += 4; //txStart
			curr_offset += 4; //txEnd
			curr_offset += 4; //cdsStart
			curr_offset += 4; //cdsEnd
			curr_offset += 4; //exonCount
			curr_offset += (uint)sequences[i].genes[j].exons.size() * 8; //exons
			if (version == 2.0f) {
				curr_offset += 4; //score
				curr_offset += 3; //colorR, colorG, colorB
			}
		}
	}
}

inline void genefile::build_map() {
	gene_map.clear();
	int num_warnings = 0;
	for (int i = 0; i < (int)num_seq; i++) {
		for (int j = 0; j < (int)sequences[i].num_genes; j++) {
			if (sequences[i].genes[j].name == "") {
				num_warnings++;
				if (num_warnings <= 5) printf("warning: transcript id is empty, skipped.\n");
			} else  if (gene_map.count(sequences[i].genes[j].name) != 0) {
				num_warnings++;
//				if (num_warnings <= 5) printf("warning: found duplicated transcript %s, renamed as %s.\n", sequences[i].genes[j].name.c_str(), (sequences[i].genes[j].name + "_" + int2str(i) + "_" + int2str(j)).c_str());
//				sequences[i].genes[j].name = sequences[i].genes[j].name + "_" + int2str(i) + "_" + int2str(j);
				if (num_warnings <= 5) printf("warning: duplicated transcript %s, skipped.\n", sequences[i].genes[j].name.c_str());
			} else if (sequences[i].genes[j].geneName == "") {
				num_warnings++;
				if (num_warnings <= 5) printf("warning: gene name is empty for transcript %s, skipped.\n", sequences[i].genes[j].name.c_str());
			} else {
				gene_map[sequences[i].genes[j].name] = pair<int, int>(i, j);
			}
		}
	}
	if (num_warnings > 5) printf("%d warnings suppressed.\n", num_warnings - 5);
}

inline void genefile::search_by_name(const string &name, string &seq_name, gene_struct &gene) {
	if (gene_map.empty()) build_map();
	if (gene_map.find(name) == gene_map.end()) {
		seq_name = "";
		return;
	}
	pair<int, int> index = gene_map[name];
	seq_name = sequences[index.first].name;
	gene = sequences[index.first].genes[index.second];
}

inline int map_gene_to_chr(int coord, vector<pair<int, int> > &regions, bool forward_strand = true, bool get_exon = false) {
	if (forward_strand) {
		int i = 0;
		while (i < (int)regions.size()) {
			if (coord < regions[i].second - regions[i].first) {
				return get_exon?i:(regions[i].first + coord);
			}
			coord -= (regions[i].second - regions[i].first);
			i++;
		}
		return get_exon?i:(coord + regions[regions.size() - 1].second - 1);
	} else {
		int i = (int)regions.size() - 1;
		while (i >= 0) {
			if (coord < regions[i].second - regions[i].first) return get_exon?i:(regions[i].second - 1 - coord);
			coord -= (regions[i].second - regions[i].first);
			i--;
		}
		return get_exon?i:(regions[0].first - coord);
	}
}

inline bool genefile::search_by_name(const string gene_name, vector<gene_struct> &genes) {
	int i, j;
	for (i = 0; i < (int)sequences.size(); i++) {
		for (j = 0; j < (int)sequences[i].genes.size(); j++) {
			if (tolower(sequences[i].genes[j].geneName).find(tolower(gene_name)) == string::npos && tolower(sequences[i].genes[j].name).find(tolower(gene_name)) == string::npos) continue;
			genes.push_back(sequences[i].genes[j]);
		}
	}
	return true;
}

inline bool genefile::search_by_region(string chr, int start, int end, vector<gene_struct> &genes) {
	int i, j;
	for (i = 0; i < (int)sequences.size(); i++) {
		if (sequences[i].name == chr) {
			for (j = 0; j < (int)sequences[i].num_genes; j++) {
				if (sequences[i].genes[j].txStart > end || sequences[i].genes[j].txEnd < start) continue;
				genes.push_back(sequences[i].genes[j]);
			}
			break;
		}
	}
	return true;
}

inline void genefile::trim_UTR() {
	int i, j;
	for (i = 0; i < (int)sequences.size(); i++) {
		for (j = 0; j < (int)sequences[i].genes.size(); j++) {
			gene_struct gene = sequences[i].genes[j];
			interval_set is(gene.exons);
			is.intersect_with(interval_set(gene.cdsStart, gene.cdsEnd));
			gene.exons.clear();
			if (!is.convert_to_int_pairs(gene.exons)) {
				panic(string("internal error: failed converting to int pair, gene ") + gene.name + "\n");
			}
			gene.txStart = gene.cdsStart;
			gene.txEnd = gene.cdsEnd;
			sequences[i].genes[j] = gene;
		}
	}
}

inline void genefile::build_interval_lists() {
	for (int i = 0; i < (int)sequences.size(); i++) {
		for (int j = 0; j < (int)sequences[i].genes.size(); j++) {
			sequences[i].gene_intervals.add_interval(sequences[i].genes[j].txStart, sequences[i].genes[j].txEnd, j);
			for (int k = 0; k < (int)sequences[i].genes[j].exons.size(); k++) {
				sequences[i].exon_intervals.add_interval(sequences[i].genes[j].exons[k].first, sequences[i].genes[j].exons[k].second, pair<int, int>(j, k));
			}
		}
		sequences[i].gene_intervals.prepare();
		sequences[i].exon_intervals.prepare();
	}
}

inline void genefile::search_interval(const string chr, const int start, const int end, vector<pair<int, pair<int, int> > > &exons) {
	exons.clear();
	for (int i = 0; i < (int)sequences.size(); i++) {
		if (sequences[i].name == chr) {
			vector<pair<int, int> > exon_results;
			sequences[i].exon_intervals.search_interval(start, end, exon_results);
			for (int j = 0; j < (int)exon_results.size(); j++) {
				exons.push_back(pair<int, pair<int, int> >(i, pair<int, int>(exon_results[j].first, exon_results[j].second)));
			}
		}
	}
}

inline void genefile::search_interval(const string chr, const int start, const int end, gene_struct *&gene, int &exon) {
	gene = NULL;
	exon = -1;
	for (int i = 0; i < (int)sequences.size(); i++) {
		if (sequences[i].name == chr) {
			vector<pair<int, int> > exon_results;
			sequences[i].exon_intervals.search_interval(start, end, exon_results);
			if (exon_results.size() > 0) {
				gene = &(sequences[i].genes[exon_results[0].first]);
				exon = exon_results[0].second;
			} else {
				vector<int> gene_results;
				sequences[i].gene_intervals.search_interval(start, end, gene_results);
				if (gene_results.size() > 0) {
					gene = &(sequences[i].genes[gene_results[0]]);
					return;
				}
			}
		}
	}
}


#endif //GENEFILE_H
