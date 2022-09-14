## odin
##INCLUDE=-I/home/emanuele/src/htslib/htslib
##LIB= /home/emanuele/src/htslib/libhts.a -lpthread -lz -lcurl -llzma -lbz2 -lcrypto -lm
## minosse
INCLUDE=-I/home/emanuele/src/htslib/htslib
LIB= /home/emanuele/src/htslib/libhts.a -lpthread -lz -lcurl -llzma -lbz2 -lcrypto -lm

all: cvlr-meth-of-bam cvlr-cluster cvlr-clean-matrix make-matrix genome-scan\
	cvlr-stats bernoulli print-genomic-positions print-cg-positions\
	check-positions-vs-cpg-list cvlr-randomize-clusters 

cvlr-cluster: cvlr-cluster.c clustering.c
	gcc -Wall -pedantic -o $@ $< -lm			

cvlr-meth-of-bam: cvlr-meth-of-bam.c cvlrlib.c tree.c
	gcc -o $@ $(INCLUDE) -fPIC  $< $(LIB)

cvlr-clean-matrix: cvlr-clean-matrix.c clustering.c tree.c
	gcc -o $@ $(INCLUDE) -fPIC  $< $(LIB)

make-matrix: make-matrix.ml
	ocamlopt str.cmxa make-matrix.ml -o make-matrix

genome-scan: genome-scan.ml
	ocamlopt $< -o $@

binomlib.cmx: binomlib.ml
	ocamlopt -c $< -o $@

cvlr-stats: cvlr-stats.ml binomlib.cmx
	ocamlopt str.cmxa binomlib.cmx $< -o $@

bernoulli: bernoulli.ml
	ocamlopt str.cmxa $< -o $@

print-genomic-positions: print-genomic-positions.c
	gcc -o $@   $< 

print-cg-positions: print-cg-positions.c
	gcc -o $@   $< 

check-positions-vs-cpg-list: check-positions-vs-cpg-list.ml
	ocamlopt.opt str.cmxa $< -o $@

cvlr-randomize-clusters: cvlr-randomize-clusters.ml 
	ocamlopt.opt str.cmxa $< -o $@

