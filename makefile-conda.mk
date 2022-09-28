INCLUDE=-I$(PREFIX)/include/htslib
LIB=-ldl $(PREFIX)/lib/libhts.a -lpthread -lz -lcurl -llzma -lbz2 -ldeflate -lcrypto -lm 

all: cvlr-meth-of-bam cvlr-cluster 

cvlr-cluster: cvlr-cluster.c clustering.c
	$(GCC) -Wall -pedantic -o $@ $< -lm			

cvlr-meth-of-bam: cvlr-meth-of-bam.c cvlrlib.c tree.c
	$(GCC) -o $@ $(INCLUDE) -L$(PREFIX)/lib -fPIC  $< $(LIB)



