INCLUDE=-I$(PREFIX)/include/htslib
LIB= $(PREFIX)/lib/libhts.a -lpthread -lz -lcurl -llzma -lbz2 -lcrypto -lm -ldeflate

all: cvlr-meth-of-bam cvlr-cluster 

cvlr-cluster: cvlr-cluster.c clustering.c
	gcc -Wall -pedantic -o $@ $< -lm			

cvlr-meth-of-bam: cvlr-meth-of-bam.c cvlrlib.c tree.c
	gcc -o $@ $(INCLUDE) -L$(PREFIX)/lib -fPIC  $< $(LIB)

