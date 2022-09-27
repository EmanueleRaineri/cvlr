INCLUDE=-I$(PREFIX)/include/htslib
LIB=$(PREFIX)/lib/libhts.a -lpthread -lz -lcurl -lcrypto -lm -lbz2 -llzma

all: cvlr-meth-of-bam cvlr-cluster 

cvlr-cluster: cvlr-cluster.c clustering.c
	$(CC) -Wall -pedantic -o $@ $< -lm			

cvlr-meth-of-bam: cvlr-meth-of-bam.c cvlrlib.c tree.c
	$(CC) -o $@ $(INCLUDE) -L$(PREFIX)/lib -fPIC  $< $(LIB)



