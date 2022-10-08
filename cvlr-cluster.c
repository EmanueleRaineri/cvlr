/*  cvlr-cluster.c -- cluster reads based on methylation patterns  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXLINE 1000

#include "clustering.c"

int main(int argc, char* argv[]){
  clock_t begin = clock();
  char* splash = malloc(MAXLINE);
  sprintf(splash, "cvlr-cluster source code timestamp: %s\n", __TIMESTAMP__);
  fprintf(stderr, "%s", splash);
  if (argc < 5){
    fprintf(stderr,
	    "usage:%s <matrix-file> ('-' for stdin) <k> <seed> <maxit>\n",
	    argv[0]);
    fprintf(stderr, "<seed> can be -1 for random initialization\n");
    exit(EXIT_FAILURE);
  }
  int k = atoi(argv[2]);
  int seed = atoi(argv[3]);
  int maxit = atoi(argv[4]);
  srand((unsigned int) seed);
  FILE* infile = stdin;
  if (strcmp(argv[1], "-")){
    infile = fopen(argv[1], "r");
    if (NULL == infile ){
      fprintf(stderr, "could not open %s\n", argv[1]);
      exit(EXIT_FAILURE);
    }
    fprintf(stderr, "reading from %s\n", argv[1]);
  } else {
      fprintf(stderr, "reading from stdin\n");
  }
  fprintf(stderr, "K:%d\n", k);
  fprintf(stderr, "seed:%d\n", seed);
  size_t i,j, n,d;
  read_n_d_from_file(infile, &n, &d);
  fprintf(stderr, "n:%lu d:%lu\n", n, d);
  int32_t* m = malloc(n * d * sizeof(int32_t));
  for(i = 0; i < n * d; i++) m[i] = -1;
  /* store names of rows and columns */
  char** rnames = malloc(n * sizeof(char*));
  for(i = 0; i < n; i++) rnames[i] = NULL;
  uint32_t* gpos = malloc(d * sizeof(uint32_t));
  for(i = 0; i < d; i++) gpos[i] = 0;
  read_matrix_from_file(infile, n,d, rnames, gpos, m);
  fclose(infile);
  // number of parameters
  fprintf(stdout, "#@N:%lu\n#@D:%lu\n#@K:%d\n", n, d, k);
  uint32_t df = k * d + k - 1;
  fprintf(stdout, "#@DF:%d\n", df);
  double* mu = calloc(k * d , sizeof(double));
  double* gamma = malloc(n * k * sizeof(double));
  double* pi = calloc(k , sizeof(double));
  double* pop = malloc(k * sizeof(double));
  //double s = 0.0;
  size_t* cluster = malloc(n * sizeof(size_t));
  for (i = 0; i < n; i++){
    cluster[i] = random_cluster(k);
  }
  size_t* n0 = calloc(k*d, sizeof(size_t));
  size_t* n1 = calloc(k*d, sizeof(size_t));
  fill_n0_n1(m, cluster, n, d, k, n0, n1);
  init_pi(cluster, n, k, pi);
  init_mu(n0,n1,k,d,mu);
  fprintf(stdout, "#@INITPI:");
  for(i = 0; i < (k-1); i++) fprintf(stdout, "%.3lf\t", pi[i]);
  fprintf(stdout, "%.3lf\n", pi[k-1]);
  fprintf(stdout, "#@GPOS:");
  for(i=0; i<d-1; i++){
    fprintf(stdout, "%u\t", gpos[i]);
  }
  fprintf(stdout, "%u\n", gpos[d-1]);
  for(i=0; i<k; i++){
    fprintf(stdout, "#@INITMU%ld:",i);
    for(j=0; j < (d - 1); j++){
      fprintf(stdout, "%.3lf\t",mu[i * d + j]);
    }
    fprintf(stdout, "%.3lf\n",mu[i * d + d - 1]);
  }
  uint32_t it = 0;
  double ll;
  /* EM loop */
  while (it < maxit){
    bernoulli_gamma(pi, m, mu, n, d, k, gamma);
    update_mu(m, gamma,  n, d, k, mu, pop);
    for(i=0;i < k; i++) pi[i] = pop[i] / n;
    ll = llbern(m , mu, pi, n, d, k);
    it++;
    fprintf(stdout,"#@IT:%u\tLL:%.3g\n", it, ll);
  }
  double bic = ll - ( df / 2 ) * log(n);
  double aic = ll - df;
  fprintf(stdout,"#@LL:%.3g\tBIC:%.3g\tAIC:%.3g\n", ll, bic, aic);
  /* ***************************** */
  fprintf(stdout, "#@PI:");
  for(i=0; i< (k-1); i++) fprintf(stdout, "%.3lf\t", pi[i]);
  fprintf(stdout, "%.3lf\n", pi[k-1]);
  for(i=0; i<k; i++){
    fprintf(stdout, "#@MU%ld:",i);
    for(j=0; j<(d-1);j++){
      fprintf(stdout, "%.3lf\t",mu[i*d+j]);
    }
    fprintf(stdout, "%.3lf\n",mu[i*d+d-1]);
  }
  for(i=0; i<n; i++){
    fprintf(stdout, "%s\t", rnames[i]);
    double maxg = -1.0;
    for(j=0; j<k;j++){
      if (gamma[i*k+j] > maxg) {
	maxg = gamma[i*k+j];
	cluster[i] = j;
      }
    }
    fprintf(stdout,"%ld\t", cluster[i]);
    for(j=0; j<k-1;j++){
      fprintf(stdout, "%.3g\t", gamma[i*k+j]);
    }
    fprintf(stdout, "%.3g\n", gamma[i*k+k-1]);
  }
  free(m);
  free(mu);
  free(pi);
  free(rnames);
  free(gpos);
  clock_t end = clock();
  fprintf(stderr, "cvlr-cluster done in %.3gs\n", (double)(end-begin)/CLOCKS_PER_SEC);
  exit(EXIT_SUCCESS);
}
