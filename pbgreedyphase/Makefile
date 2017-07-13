all: partitionByPhasedSNVs readToSNVList

SEQAN=../seqan/include
BOOSTLIB=../lib
BOOST=../boost/include
BLASR=../blasr/common
VCFLIB=../vcflib
HTSLIB=$(VCFLIB)/tabixpp/htslib
CPPOPTS=  -g
CPP=g++ 
partitionByPhasedSNVs: PartitionByPhasedSNVs.cpp FastaIndex.h
	$(CPP) -g -static $(CPPOPTS) $< \
     -I $(SEQAN) \
     -I $(BLASR) \
     -I $(BOOST) \
     -I $(VCFLIB)/include \
     -I $(HTSLIB) \
     -L $(BOOSTLIB) -l boost_program_options \
     -L $(VCFLIB)/lib -l vcflib  -L $(HTSLIB) -l hts  -lpthread -lz \
     -o $@ 


readToSNVList: ReadToSNVList.cpp
	$(CPP) -static $(CPPOPTS) $^ \
     -I $(SEQAN) \
     -I $(BLASR) \
     -I $(BOOST) \
     -L $(BOOSTLIB) -l boost_program_options \
     -o $@ 

