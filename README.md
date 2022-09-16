# Installation

cvlr depends on [htslib](https://github.com/samtools/htslib).
You need to install it and write the relevant path in the Makefile.
Specifically look at the lines starting with INCLUDE and LIB.
If you have installed htslib in the directory
/home/user/src/htslib`
those two lines will have to be:

INCLUDE=-I/home/user/src/htslib/htslib
LIB=/home/user/src/htslib/libhts.a -lpthread -lz -lcurl -llzma -lbz2 -lcrypto -lm

Once you have edited the Makefile, just type `make`.

# Usage

You need to have a BAM/CRAM file containing methylation data
encoded with the Mm/Ml tages as explained in the SAM documentation
(https://samtools.github.io/hts-specs/SAMtags.pdf). This is also the standard output
format (at the moment) of Nanopore's Megalodon.
First you create the matrix from the BAM file with

`cvlr-meth-of-bam NA12878.cram chr20:58839718-5891119 > GNAS-matrix.txt`

then you cluster the reads with

`cvlr-cluster GNAS-matrix.txt 2 1 100 > GNAS-clusters.txt`

(That is you are creating 2 clusters, with seed 1 for the random numbers generator and
a maximum of 100 EM iterations).

You can then look at the clusters (for example for plotting the methylation values) with

`cvlr-stats GNAS-clusters.txt GNAS-matrix.txt > GNAS-stats.txt`

If you have 2 clusters columns 1,4,7 of contain the genomic position,
methylation in cluster 0 and methylation in cluster 1 respectively.



