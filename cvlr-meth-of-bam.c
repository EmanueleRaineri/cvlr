/* cvlr-meth-of-bam

Emanuele Raineri

parse Mm Ml tags to extract methylation matrix 

*/

#include <stdio.h>
#include <stdbool.h>
#include "cvlrlib.c"

#define MAXCGCOUNT 262144
#define RNAMESIZE 1000
#define MAXREADLE 1048576
#define DEBUG 0

static char *code(int id) {
    static char code[20];
    if (id > 0) {
        code[0] = id;
        code[1] = 0;
    } else {
        sprintf(code, "(%d)", -id);
    }
    return code;
}

void splash(){
  fprintf(stderr, "cvlr-meth-of-bam source code timestamp: %s\n", __TIMESTAMP__);
  fprintf(stderr, "compiled against htslib version: %d\n", HTS_VERSION );
}

int main (int argc, char* argv[]){
  char *bamfn, *region, *rname;
  const char* regout;
  char nuc1, nuc2;
  char** rnames = malloc(MAXCGCOUNT * sizeof(char*));
  htsFile *fp;
  sam_hdr_t *hdr;
  bam1_t* b;
  hts_idx_t* idx;
  hts_pos_t start, end;
  hts_itr_t* iter;
  hts_base_mod_state *m; 
  hts_base_mod mods[5];
  
  uint isrev;
  uint  nmatches /**< number of positions in the reads which match the genome */; 
  /** number of useful CGs in a given read. 
      (useful == in window and matching properly, see below) */
  uint  cg_in_read; 
  uint64_t minpos, maxpos, span, discarded=0;

  //int32_t l_qseq;
  int32_t *rpos = malloc(MAXREADLE * sizeof(int32_t));
  /**< positions on the read which match the genome */
  int64_t *gpos = malloc(MAXREADLE * sizeof(int64_t));;
  /**< positions on the genome matched by rpos */
  int64_t *allgpos = malloc(MAXCGCOUNT * sizeof(int64_t)); 

  bool cgmatch;
  
  int i, j, n, npos, tid_out = -1, ismeth;
  int* meth = malloc(MAXCGCOUNT * sizeof(int));
  
  size_t cgcount = 0, nreads = 0, *idxread = malloc(MAXCGCOUNT*sizeof(size_t));

  struct pos_tnode* gpostree = NULL;
  splash();
  if (argc != 3){
    fprintf(stderr, "usage %s <cram/sam/bam file> <region>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  for(i=0; i < MAXCGCOUNT; i++){
    rnames[i] = malloc(RNAMESIZE * sizeof(char));
  }
  
  bamfn = argv[1]; region = argv[2];
  fprintf(stderr, "parsing %s\n", bamfn);
  fp = hts_open(bamfn, "r");
  hdr = sam_hdr_read(fp);
  if ( !fp ||  !hdr ) {
    fprintf(stderr, "fatal: can't open %s\n", bamfn);
    goto err;
  }
  b = bam_init1();
  m = hts_base_mod_state_alloc();
  if (!m) exit(EXIT_FAILURE);
  idx = sam_index_load(fp, bamfn);
  if (!idx) {
    fprintf(stderr, "fatal: can't find index\n");
    goto err;
  }
  regout = sam_parse_region(hdr, region, &tid_out, &start, &end, 0);
  if (!regout) {
    fprintf(stderr, "fatal: can't parse region\n");
    goto err;
  }
  fprintf(stderr, "region %s (tid_out:%d start:%"PRId64" end:%"PRId64")\n",
	  region, tid_out, start, end);
  iter = sam_itr_querys(idx, hdr, region);
  if ( !iter ){
    fprintf(stderr, "fatal: empty sam_itr_querys\n");
    goto err;
  }
  
  /** extract reads overlapping with the selected region */
  while ( sam_itr_next(fp, iter, b) >= 0 ){
    isrev = bam_is_rev(b); rname = bam_get_qname(b);
    //l_qseq = b->core.l_qseq;
    if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP)) {
      fprintf(stderr, "discarding %s [flag:%d]\n", rname,b->core.flag );
      discarded++;
      continue;
    }
    //fprintf(stderr, "processing %s\n", rname);
    if ((!bam_aux_get(b, "MM")) && (!bam_aux_get(b, "Mm"))){
      fprintf(stderr, "Skipping read %s (no MM/ML aux tag)\n", rname);
      continue;
    }
    
    if ( bam_parse_basemod(b, m) < 0 ) {
      fprintf(stderr, "Failed to parse MM/ML aux tags in %s\n", rname);
            goto err;
    }
    
    nmatches = get_aligned_pairs(b, gpos, rpos);
    //fprintf(stderr, "number of matches on %s:%u\n", rname, nmatches);
    /** look only at matching positions to find usable CGs*/
    cg_in_read=0;
    for(i = 0; i < nmatches; i++) {
      /* check if matched position is inside [start, end] */
      if ( (gpos[i] < start) || (gpos[i] > end) ) continue;
      nuc1 = seq_nt16_str[bam_seqi(bam_get_seq(b), rpos[i])];
      nuc2 = seq_nt16_str[bam_seqi(bam_get_seq(b), rpos[i+1])];
      /** only consider CpGs on the read which match positions on the genome 
	  ie avoid indels */
      cgmatch = ('C' == nuc1) && ('G' == nuc2);
      cgmatch = cgmatch && (rpos[i+1] == rpos[i]+1) && (gpos[i+1] == gpos[i]+1);
      if ( !cgmatch ) continue;
      cg_in_read++;
      if ( isrev ){
	n = bam_mods_at_qpos(b, rpos[i]+1, m, mods, 5);
      } else {
	n = bam_mods_at_qpos(b, rpos[i], m, mods, 5);
      }
      if ( 0 == n ) ismeth = 0; //no mods listed  
      if ( n > 0 ){
	// check if methylation is one of the listed mods
	ismeth = -1;
	for (j = 0; j < n && j < 5; j++){
	  if ( (109 == mods[j].modified_base)  ){
	    ismeth = mods[j].qual;
	    break;
	  }
	}
	if (-1 == ismeth) {ismeth = 0; n=0;} 
      }
      strcpy( rnames[cgcount], rname );
      allgpos[cgcount] = gpos[i];
      gpostree = pos_addtree( gpostree, gpos[i] );
      idxread[cgcount] = nreads;
      meth[cgcount] = binarize(ismeth, 127,127);
      /* debug info */
      if (DEBUG){
	fprintf(stderr, "%ld\t%s\t%d\t%"PRId64"\t%c\t%d\t", cgcount, rname,
		rpos[i], gpos[i], nuc1, n);
	if (n > 0) {
	  fprintf(stderr, "%c\t%c\t%s\t%d\n",
		  mods[j].canonical_base,
		  "+-"[mods[j].strand],
		  code(mods[j].modified_base), mods[j].qual);
	} else{
	  fprintf(stderr,"\n");
	}
      }
      /* -------------- */
      cgcount++;
    } 
    if (cg_in_read > 0 ) nreads++;
  }
  if (!gpostree){
    fprintf(stderr, "warning: no CpGs position(s) found\n");
    goto out;
  }
  pos_index(gpostree, 0);
  npos = count_pos_tnodes(gpostree);
  minpos = min_pos(gpostree);
  maxpos = max_pos(gpostree);
  span = maxpos - minpos + 1;
  /** print out matrix 
   this does not necessarily print out any -1s 
   ( the presence of -1s depends on the boundaries passed to 
   binarize() ) 
   When reading this file and organizing it into a 2d array (rows are reads,
   columns CpG positions the -1s will be present everytime some position is
   not covered by all reads
  */
  fprintf(stdout, "#@N:%ld\n#@D:%d\n", nreads, npos );
  fprintf(stdout, "#@MINPOS:%"PRIu64"\n#@MAXPOS:%"PRIu64"\n#@SPAN:%ld\n",
	  minpos, maxpos, span);
  fprintf(stdout, "#@REGION:%s\n", region);
  fprintf(stdout, "#@CGCOUNT:%lu\n", cgcount);
  fprintf(stdout, "#@DISCARDED_READS:"PRIu64"\n", discarded);
  for(i=0; i< cgcount; i++){
      fprintf(stdout, "%ld\t%s\t%d\t%"PRId64"\t%d\n",
	      idxread[i], rnames[i],
	      pos_index_of(gpostree, allgpos[i]),
	      allgpos[i], meth[i]);
  }
 out:
  free(rpos);
  free(gpos);
  bam_destroy1(b);
  sam_hdr_destroy(hdr);
  hts_base_mod_state_free(m);
  fprintf(stderr, "done.\n");
  return EXIT_SUCCESS;
 err:
  bam_destroy1(b);
  sam_hdr_destroy(hdr);
  hts_base_mod_state_free(m);
  return EXIT_FAILURE;
}
