
#include "sam.h"
#include "hts.h"
#include "bgzf.h"

#include "tree.c"

#define MAXREADS 65536
#define RNAMESIZE 1000
#define MAXCIGARSIZE 65536
#define MAXCG 500000


int binarize(int meth, int lb, int ub ){
  if ((meth >= 0) & (meth <= lb)) return 0;
  if ((meth > lb) & (meth <= ub)) return -1;
  if ((meth > ub) & (meth <= 255)) return 1;
  fprintf(stderr, "error: illegal meth score:%d\n", meth);
  exit(EXIT_FAILURE);
}



uint32_t get_cigar(bam1_t* b, char* scigar, unsigned int* lcigar ){
  /* 
     scigar and lcigar are allocated by the caller. 
     fills scigar with the operations and
     lcigar with the corresponding lengths 
  */
  uint32_t* cigar = bam_get_cigar(b);
  uint32_t nop = b->core.n_cigar;
  for( uint32_t i =0; i< nop; i++){
    scigar[i] = bam_cigar_opchr(cigar[i]);
    lcigar[i] = bam_cigar_oplen(cigar[i]);
  }
  return nop;
}

uint32_t get_mc_tag(bam1_t* b, int64_t* cgmeth ){
  /* returns methylation values */
  uint8_t* mctag = bam_aux_get(b, "MC");
  uint32_t cgmethl;
  if ( NULL == mctag ){
    fprintf(stderr, "MC tag not present\n");
    cgmethl = 0;
  } else {
    cgmethl = bam_auxB_len(mctag);
    for(uint32_t i =0; i < cgmethl ; i++){
      cgmeth[i] = bam_auxB2i(mctag, i);
    }
  }
  return cgmethl;
}

void get_query_seq(uint8_t* seq, int32_t le, char* qseq){
  // Transform an array of uint8_t into an array
  // characters and store it in qseq
  // you have to get a pointer to seq and compute its
  // length outside of this function.
  int32_t i;
  uint8_t tmp;
  for(i = 0; i < le; i++){
    tmp = bam_seqi(seq, i);
    if ( 1 == tmp ) {qseq[i] = 'A'; continue;}
    if ( 2 == tmp ) {qseq[i] = 'C'; continue;}
    if ( 4 == tmp ) { qseq[i] = 'G'; continue;}
    if ( 8 == tmp ) { qseq[i] ='T'; continue;}
    if ( 15 == tmp ) { qseq[i] = 'N'; continue;}
  }
  qseq[i] = '\0';
}

unsigned int get_aligned_pairs(bam1_t* b, int64_t* gpos, int32_t* rpos ){
  /*
    fill gpos and rpos
   */
  char* cigarops = malloc(MAXCIGARSIZE);
  unsigned int* opslen = malloc(MAXCIGARSIZE * sizeof(unsigned int));
  uint32_t nops = get_cigar(b, cigarops, opslen);
  int64_t g = b->core.pos + 1; int32_t r = 0;
  unsigned int j;
  unsigned int nmatches=0;
  for(uint32_t op = 0; op < nops; op++){
    char c = cigarops[op];
    switch(c) {
    case 'M': case '=': case 'X':
      for (j = 0; j < opslen[op]; j++){
	gpos[nmatches  +  j] = g; rpos[nmatches   + j] = r;
	g++; r++;
      }
      nmatches += opslen[op];
      break;
    case 'I':case'S':
      r += opslen[op];
      break;
    case 'D': case 'N':
      g += opslen[op];
      break;
    case 'H': case 'P':case 'B':
      break;
    }
  }
  return nmatches;
}


void create_rname_pos_trees(char* fname,
			    struct pos_tnode **pos_root,
			    struct rname_tnode **rname_root){
  FILE* fp = fopen(fname, "r");
  char *rname = malloc(1000);
  char isrev;
  uint32_t rpos;
  int64_t gpos;
  int64_t meth;
  while(1){
    fscanf(fp, "%s\t%c\t%u\t%"SCNd64"\t%"SCNd64"\n",
	   rname, &isrev, &rpos, &gpos, &meth);
    *pos_root = pos_addtree(*pos_root, gpos);
    *rname_root = rname_addtree(*rname_root, rname);
    if (feof(fp)) break;
  }
  fclose(fp);
}

void print_mm_matrix(int* mm, size_t nr, size_t npos){
  for(size_t i =0; i< nr; i++){
    for(size_t j=0;j<npos; j++){
      printf("%d,", mm[i*npos+j]);
    }
    printf("\n");
  }
}

void make_matrix_from_file_trees(int32_t* mm, size_t nr, size_t npos,
				 char* fname,
				 struct rname_tnode *rname_root,
				 struct pos_tnode *pos_root){
  mm = memset(mm, -1, nr * npos * sizeof(int32_t));
  FILE *fp = fopen(fname, "r");
  char *rname = malloc(1000);
  char isrev;
  uint32_t rpos;
  int64_t gpos;
  int64_t meth;
  int i, j;
  while(1){
    fscanf(fp, "%s\t%c\t%u\t%"SCNd64"\t%"SCNd64"\n",
	   rname, &isrev, &rpos, &gpos, &meth);
    i = rname_index_of(rname_root, rname);
    if ( -1 == i ) {
      fprintf(stderr, "read:%s not found\n", rname);
      exit(1);
    }
    j = pos_index_of(pos_root, gpos);
    if ( -1 == j ) {
      fprintf(stderr, "pos:%"PRId64" not found\n", gpos);
      exit(1);
    }
    mm[i * npos + j ] = (int)(meth > 127)?1:0;
    if (feof(fp)) break; 
  }
  fclose(fp); free(rname);
}

void make_matrix_from_file(char *fname, uint64_t* n, uint64_t* d){
  struct pos_tnode *pos_root = NULL;
  struct rname_tnode *rname_root = NULL;
  create_rname_pos_trees(fname, &pos_root, &rname_root );
  *n = rname_index(rname_root, 0);
  *d = pos_index(pos_root, 0);
  printf("%"PRId64" unique reads\n", *n);
  printf("%"PRId64" unique position(s)\n", *d);
  int32_t* mm = calloc( (*n) * (*d) , sizeof(int32_t));
  make_matrix_from_file_trees(mm, *n, *d, fname, rname_root, pos_root);
}

size_t filter_columns(int* mm, double* mean, size_t n, size_t d,
		    int** cmm, uint64_t* gpos, uint64_t** gpos2){
  size_t i,j, fd =  0;
  for(j =0 ; j < d; j++) {if ( -1 != mean[j] ) fd++;}
  *cmm = malloc( n * fd * sizeof(int));
  *gpos2 = malloc( fd * sizeof(uint64_t));
  size_t j2 = 0;
  for( j = 0; j < d; j++ ){
    if ( -1 == mean[j] ) {
      continue;
    } else { /*copy the column*/
      (*gpos2)[j2] = gpos[j];
      for( i = 0; i < n; i++ ){
    	(*cmm)[i * fd + j2] = mm[i * d + j];
      }
      j2++;
    }
  }
  return fd;
}
