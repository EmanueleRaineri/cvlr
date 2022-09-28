INCLUDE=-I$(PREFIX)/include/htslib
LIB=$(PREFIX)/lib/libhts.a -ldl -lbz2 -llzma -lpthread -lz -lcurl -lcrypto -lm 

all: cvlr-meth-of-bam cvlr-cluster 

cvlr-cluster: cvlr-cluster.c clustering.c
	$(GCC) -Wall -pedantic -o $@ $< -lm			

cvlr-meth-of-bam: cvlr-meth-of-bam.c cvlrlib.c tree.c
	$(GCC) -o $@ $(INCLUDE) -L$(PREFIX)/lib -fPIC  $< $(LIB)



