#include "iso_exp.h"

int main(int argc, char* argv[]){

	local = new local_data();

	int arg_options = 1;

	if (argc <= 2) {
		goto bad_options;
	} 

	local->other_args.clear();

	if (string(argv[1]).length() < 1 || string(argv[1])[0] == '-') {
		task = "pipeline";		
	} else {
		task = argv[1];
		arg_options = 2;
	}

	for (int i = arg_options; i < argc; i++) {
		string arg = argv[i];
		__ASSERT(arg.length() > 0, "arg.length() == 0");
		if (arg[0] == '-') {
			if (arg == "-r" || arg == "--read_length") {
				if (i + 1 >= argc) {
					cout << "no read_length given\n";
					goto bad_options;
				} else {
					read_length = str2int(argv[i+1]);
					if (read_length > MAX_READ_LEN || read_length <= 0) {
						panic("invalid read length");
					}
					i++;
				}
			} else if (arg == "-n" || arg == "--num_mismatch") {
				if (i + 1 >= argc) {
					cout << "no num_mismatch specified\n";
					goto bad_options;
				} else {
					num_mismatch = str2int(argv[i+1]);
					i++;
				}
			} else if (arg == "-d" || arg == "--direction") {
				if (i + 1 >= argc) {
					cout << "no direction specified\n";
					goto bad_options;
				} else {
					direction = tolower(string(argv[i+1]));
					if (direction != "f" && direction != "r" && direction != "b") {
						cout << "invalid direction specified\n";
						goto bad_options;
					}
					i++;
				}
			} else if (arg == "-f" || arg == "--reduce_len_file") {
				if (i + 1 >= argc) {
					cout << "no reduce_len_file specified\n";
					goto bad_options;
				} else {
					reduce_len_file = argv[i+1];
					i++;
				}
			} else if (arg == "-a" || arg == "--annotation") {
				if (i + 1 >= argc) {
					cout << "no annotation file specified\n";
					goto bad_options;
				} else {
					annotation_file = argv[i+1];
					i++;
				}
			} else if (arg == "-ref" || arg == "--reference") {
				if (i + 1 >= argc) {
					cout << "no reference file specified\n";
					goto bad_options;
				} else {
					reference_file = argv[i+1];
					i++;
				}
			} else if (arg == "-exons") {
				if (i + 1 >= argc) {
					cout << "no exon file specified\n";
					goto bad_options;
				} else {
					exon_file = argv[i+1];
					i++;
				}
			} else if (arg == "-prefix") {
				if (i + 1 >= argc) {
					cout << "no prefix specified\n";
					goto bad_options;
				} else {
					local->output_prefix = argv[i+1];
					i++;
				}
			} else if (arg == "-p" || arg == "--ins_len_p_file") {
				if (i + 1 >= argc) {
					cout << "no ins_len_p_file specified\n";
					goto bad_options;
				} else {
					ins_len_p_file = argv[i+1];
					i++;
				}
			} else if (arg == "-g" || arg == "--gene") {
				if (i + 1 >= argc) {
					cout << "no gene specified\n";
					goto bad_options;
				} else {
					if (file_exists(argv[i+1])) {
						vector<string> names;
						if (!load_file(argv[i+1], names)) {
							cout << "read gene list file failed\n";
							goto bad_options;
						}
						for (int i = 0; i < (int)names.size(); i++) {
							gene_name.insert(tolower(names[i]));
						}
					} else {
						gene_name.insert(tolower(argv[i+1]));
					}
					i++;
				}
			} else if (arg == "-v" || arg == "--version") {
				if (i + 1 >= argc) {
					cout << "no version specified\n";
					goto bad_options;
				} else {
					version = str2int(argv[i+1]);
					i++;
				}
			} else if (arg == "-progress") {
				if (i + 1 >= argc) {
					cout << "no num_reads specified\n";
					goto bad_options;
				} else {
					report_progress = str2int(argv[i+1]);
					i++;
				}
			} else if (arg == "-t" || arg == "--t-test") {
				t_test = true;
			} else if (arg == "-m" || arg == "--multiple") {
				multiple = true;
			} else if (arg == "-q" || arg == "--quick") {
				quick = true;
			} else if (arg == "-ns" || arg == "--no_early_stop") {
				do_stop = false;
			} else if (arg == "-it" || arg == "-iteration") {
				if (i + 1 > argc) {
					cout << "no num_iteration specified\n";
					goto bad_options;
				} else {
					num_itr = str2int(argv[i+1]);
					i++;
				}
			} else if (arg == "-rs" || arg == "-reset-srand") {
				reset_srand = true;
			} else if (arg == "-o" || arg == "--old") {
				old = true;
			} else if (arg == "-em") {
				do_EM = true;
			} else if (arg == "-l" || arg == "--lasso") {
				lambda = str2double(argv[i+1]);
				i++;
			} else if (arg == "-i" || arg == "--isoform") {
				if (i + 1 >= argc) {
					cout << "no num_max_isoform specified\n";
					goto bad_options;
				} else {
					num_max_isoform = str2int(argv[i+1]);
					i++;
				}
			} else if (arg == "--num_chr") {
				num_random_chr = str2int(argv[i+1]);
				i++;
			} else if (arg == "--chr_size") {
				chr_random_size = str2int(argv[i+1]);
				i++;
			} else if (arg == "--num_gene") {
				num_random_gene = str2int(argv[i+1]);
				i++;
			} else if (arg == "--num_iso") {
				num_random_iso = str2int(argv[i+1]);
				i++;
			} else if (arg == "--num_reads") {
				num_random_reads = str2int(argv[i+1]);
				i++;
			} else if (arg == "-samse") {
				working_mode = MR_SAM;
			} else if (arg == "-sampe") {
				working_mode = MR_SAM_PAIRED;
			} else if (arg == "-eland") {
				working_mode = MR_ELAND_MULTI;
			} else if (arg == "-scan_reads") {
				task = "scan_reads";
			} else if (arg == "-map_stat") {
				task = "map_stat";
			} else if (arg == "-exon_usage") {
				task = "exon_usage";
			} else if (arg == "-expression_analysis") {
				task = "expression_analysis";
			} else if (arg == "-sequential") {
				task = "sequential";
			} else if (arg == "-two_sample_diff_gene_exp") {
				task = "two_sample_diff_gene_exp";
			} else if (arg == "-two_sample_diff_exon_usage") {
				task = "two_sample_diff_exon_usage";
			} else {
				cout << "error: unrecognized option: " << arg << endl;
				goto bad_options;
			}
		} else {
			local->other_args.push_back(arg);
		}
	}

	if (task == "pipeline") {
		if (!in_range(version, 0, 1)) {
			panic("version is not between 0 and 1!");
		}		
		if (task == "exon_usage") {
			if (version != 0) {
				panic("version is not 0!");
			}
		}
		if (task == "expression_analysis" || task == "two-sample-diff-gene-exp") {
			if (version != 1) {
				panic("version is not 1!");
			}
		}
	} else {
		if (task == "paired_end" || task == "convert_coord") {
			if (!in_range(version, 0, 4)) {
				panic("version is not between 0 and 4!");
			}
		} else if (task == "comp_exp") {
			if (!in_range(version, 1, 4)) {
				panic("version is not between 1 and 4!");
			}
		} else if (task == "gen_exons_junctions") {
			if (!in_range(version, 2, 4)) {
				panic("version is not between 2 and 4!");
			}
		} else if (task == "extract") {
			if (!in_range(version, 0, 1)) {
				panic("version is not between 0 and 1!");
			}
		} else if (task == "exon_usage") {
			if (version != 0) {
				panic("version is not 0!");
			}
		} else {
			if (!in_range(version, 1)) {
				panic("version is not 1!");
			}
		}
	}

	if (task == "pipeline") {
		if (task == "sequential") {
			 if (task == "two-sample-diff-gene-exp") {
				if (!in_range((int)local->other_args.size(), 2)) {
					panic("number of arguments is not 2!");
				}
			 } else {
				if (local->other_args.size() < 1) {
					panic("number of arguments < 1!");
				}
			 }
		} else {
			if (!in_range((int)local->other_args.size(), 1, 2)) {
				panic("number of arguments is not 1 or 2!");
			}
		}		
	} else {
			if (task == "mapability" || task == "scan_reads" || task == "trim_UTR" || task == "enumerate_fasta" 
				|| task == "extract" || (task == "convert_coord" && version == 0)
				|| task == "random_genome" || task == "random_annotation" || task == "random_expression" || task == "gtf2refFlat" || task == "gtf2bed") {
			if (!in_range((int)local->other_args.size(), 1)) {
				panic("number of arguments is not 1!");
			}
		} else if (task == "paired_end" && version >= 1) {
			if (!in_range((int)local->other_args.size(), 3)) {
				panic("number of arguments is not 3!");
			}
		} else if (task == "expression_analysis") {
			if (!in_range((int)local->other_args.size(), 2, 3)) {
				panic("number of arguments is not 2 or 3!");
			}
		} else {
			if (!in_range((int)local->other_args.size(), 2)) {
				panic("number of arguments is not 2!");
			}
		}
	}

	if (task == "sequential") {
		do_sequential_tasks();
	} else if (task == "pipeline" || task == "extract" || task == "map_stat" || task == "paired_end" || task == "comp_exp" 
			|| task == "convert_coord" || task == "exon_usage" || task == "expression_analysis") {
		do_task();	
	} else if (task == "trim_UTR") trim_UTR();
	else if (task == "enumerate_fasta") enumerate_fasta();
	else if (task == "gen_trans") gen_trans();
	else if (task == "gen_exons_junctions") gen_exons_junctions();
	else if (task == "scan_reads") scan_reads();
	else if (task == "check_rep") check_rep();
	else if (task == "differential") differential();
	else if (task == "mapability") mapability();
	else if (task == "denovo") denovo();
	else if (task == "random_genome") random_genome();
	else if (task == "random_annotation") random_annotation();
	else if (task == "random_expression") random_expression();
	else if (task == "random_reads") random_reads();
	else if (task == "annotate_transcripts") annotate_transcripts();
	else if (task == "generate_transcripts") gen_trans(true);
	else if (task == "gtf2refFlat") gtf_convert(1);
	else if (task == "gtf2bed") gtf_convert(2);
	else {
		cout << "unrecognized task: " << task << endl;
		goto bad_options;
	}

	delete local;
	local = NULL;
	return 0;

bad_options:
	delete local;
	local = NULL;

//#define FULL_VERSION

#ifndef FULL_VERSION //trimmed version
	printf("Usage: %s <TASK> [OPTION] ... [FILE] ...\n", argv[0]);
	printf("rSeq: RNA-Seq Analyzer.\n");
	printf("\nOPTIONs:\n");
	printf("  -a, --annotation refFlat.txt\n");
	printf("\nTASKs, OPTIONs and FILEs:\n");
	printf("  scan_reads input.reads\n");
	printf("  annotate_transcripts refMrna.fa kgXref.txt\n");
	printf("  generate_transcripts refFlat.txt chrfilelist.txt\n");
	printf("  expression_analysis [-a] refMrna.fa.new.fa mapped.reads.output1 [mapped.reads.output2]\n");
	printf("\nReport bugs to <jiangh@stanford.edu>.\n");
	return 1;
#else //full version
	printf("Usage: %s <TASK> [OPTION] ... [FILE] ...\n", argv[0]);
	printf("rSeq: RNA-Seq Analyzer (full version).\n");

	printf("\nOPTIONs:\n");
	printf("  -r, --read_length READ_LENGTH (default is 25)\n");
	printf("  -n, --num_mismatch NUM_MISMATCH (default is 3)\n");
	printf("  -d, --direction [f|r|b] (default is b)\n");
	printf("  -f, --reduce_len_file FILE_NAME\n");
	printf("  -g, --gene GENE_NAME or GENE_LIST_FILENAME\n");
	printf("  -q, --quick (for paired_end or differential)\n");
	printf("  -ns, --no_early_stop (for differential, default is stop)\n");
	printf("  -it, --iteration MUM_ITERATION (for differential, default is 1000)\n");
	printf("  -rs, --reset-srand (for differential, default is false)\n");
	printf("  -samse\n");
	printf("  -sampe\n");
	printf("  -eland\n");
	printf("  -i, --isoform NUM_MAX_ISOFORM (for differential, default is 3, 0 for infinity)\n");
	printf("  -t, --t-test\n");
	printf("  -p, --ins_len_p FILENAME\n");
	printf("  -progress NUM_READS\n");
	printf("  -m, --multiple\n");
	printf("  -o, --old\n");
	printf("  -em\n");
	printf("  -l, --lasso LAMBDA (default is 0, -1 for auto-selection)\n");
	printf("  -a, --annotation refFlat.txt\n");
	printf("  -ref, --reference reference.fa\n");
	printf("  -prefix output_prefix\n");
	printf("  -v, --version VERSION (default is 1)\n");
	printf("     0: work with genome\n");
	printf("     1: work with transcripts\n");
	printf("     2: work with exons and junctions, combine consecutive exons\n");
	printf("     3: work with exons only, do not combine consecutive exons\n");
	printf("     4: work with exons and junctions, do not combine consecutive exons\n");

	printf("\nTASKs, OPTIONs and FILEs:\n");
	printf("  scan_reads input.reads\n");
	printf("  check_rep [-t] coord.1 coord.2\n");
	printf("  trim_UTR  refFlat.txt\n");
	printf("  enumerate_fasta [-r] input.fasta\n");
	printf("  gen_trans refFlat.txt <genome.fasta or chrfilelist.txt>\n");
	printf("  gen_exons_junctions <-v 2, 3 or 4> [-r] refFlat.txt chrfilelist.txt\n");
	printf("  mapability [-r] chrfilelist.txt\n");
	printf("  denovo <region> reads.refFlat.txt\n");

	printf("  map_stat [-r] [-n] refFlat.txt.fa eland.multi.out\n");
	printf("  exon_usage <-v 0> [-n] refFlat.txt eland.multi.out\n");

	printf("  extract <-v 0> [-n] eland.multi.out\n");
	printf("  extract [-g] [-n] eland.multi.out\n");
	printf("  convert_coord <-v 0> [-n] eland.multi.out\n");
	printf("  convert_coord [-n] refFlat.txt eland.multi.out\n");
	printf("  convert_coord <-v 2, 3 or 4> [-r] [-n] subexons.txt eland.multi.out\n");
	printf("  comp_exp [-r] [-n] [-f] refFlat.txt.fa eland.multi.out\n");
	printf("  comp_exp <-v 2, 3 or 4> [-r] [-n] [-d] [-m] subexons.txt eland.multi.out\n");
	printf("  paired_end <-v 0> [-r] [-n] eland.multi.out.1 eland.multi.out.2\n");
	printf("  paired_end [-r] [-n] refFlat.txt.fa eland.multi.out.1 eland.multi.out.2\n");
	printf("  paired_end <-v 2, 3 or 4> [-r] [-n] [-p] [-q] subexons.txt eland.multi.out.1 eland.multi.out.2\n");

	printf("  differential <-o> [-g] [-q] [-i] condition1.sampling_rates condition2.sampling_rates\n");
	printf("  differential [-g] [-q] [-i] condition1.categories condition2.categories\n");

	printf("  random_genome [--num_chr (default = 1)] [--chr_size (default = 1000000)] prefix\n");
	printf("  random_annotation [--num_chr (default = 1)] [--chr_size (default = 1000000)] [--num_gene (default = 1)] [--num_iso (default = 2)] genes.refFlat.txt\n");
	printf("  random_expression genes.refFlat.txt\n");
	printf("  random_reads [--num_reads (default = 10000)] [-r] genes.refFlat.txt.fa genes.refFlat.txt.exp\n");

	printf("  gtf2refFlat genes.gtf\n");
	printf("  gtf2bed genes.gtf\n");
	printf("  annotate_transcripts refMrna.fa kgXref.txt\n");
	printf("  generate_transcripts refFlat.txt <genome.fasta or chrfilelist.txt>\n");
	printf("  expression_analysis [-l] [-a] refMrna.fa.new.fa mapped.reads.output1 [mapped.reads.output2]\n");

	printf("  pipeline [task1, task2, ...] mapped.reads.output1 [mapped.reads.output2]\n");
	printf("    -scan_reads\n"); // -v 0 or -v 1
	printf("    -map_stat <-ref references.fa>\n"); // -v 0 or -v 1
	printf("    -exon_usage <-v 0> <-a refFlat.txt>\n"); // -v 0
	printf("    -expression_analysis [-a] <-ref references.fa>\n"); // -v 1
	printf("	-sequential\n"); // -v 0 or -v 1
	printf("	-two_sample_diff_gene_exp\n"); // -v 1
	printf("	-two_sample_diff_exon_usage <-a refFlat.txt> <-exons exon_list>\n"); // -v 1

	printf("\nReport bugs to <jiangh@stanford.edu>.\n");
	return 1;
#endif
}
