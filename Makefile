# rSeq Makefile

.PHONY: all clean

.DELETE_ON_ERROR:

CXXFLAGS = -O3 -std=c++0x

all: rseq seqmap

rseq: RNASeq.o
	g++ -o $@ $(CXXFLAGS) RNASeq.o

seqmap: match.o
	g++ -o $@ $(CXXFLAGS) match.o

clean:
	rm -f rseq seqmap *.o
		
# generated via 'g++ -MM'
RNASeq.o: RNASeq.cpp iso_exp.h stl.h genefile.h type.h endian.h \
  system_utils.h timer.h chr_region.h string_operation.h math_utils.h \
  file_operation.h bar.h fasta.h map_result.h exon_junction_extractor.h \
  isoform_opt.h probe.h iso_exp_model.h table.h
match.o: match.cpp probe_match.h stl.h probe.h type.h string_operation.h \
  math_utils.h system_utils.h timer.h file_operation.h fasta.h
