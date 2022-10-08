
#include <float.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

size_t random_cluster(size_t k){
  // generate a random (uniformly) integer
  // between 0 and k-1
  double r = (double) rand()/RAND_MAX, s=0.0;
  size_t l;
  for(l=0; l < k;l++){
    s += 1.0/k;
    if (r <= s) return l;
  }
  return (k-1);
}

void fill_n0_n1(int* m, size_t* cluster, size_t n, size_t d, size_t k,
		size_t* n0, size_t* n1){
  // given the clustering of a binary matrix
  // compute the number of 0s (1s) per cluster, per position
  // NAs are represented as -1s in the input matrix.
  size_t i, j;
  for( i =0; i<n; i++){
    for(j=0; j<d; j++){
      switch(m[i * d + j]){
      case -1: break;
      case 0: n0[cluster[i]*d+j]+=1; break;
      case 1: n1[cluster[i]*d+j]+=1; break;
      default:
	fprintf(stderr, "%s:%d:invalid matrix entry:%d\n",
		__FILE__,__LINE__, m[i*d+j]);
      }
    }
  }
}


void init_pi(size_t* cluster, size_t n, size_t k, double* pi){
  size_t i;
  for(i=0; i<n; i++){
    pi[cluster[i]] += 1.0/n;
  }
}

void init_mu(size_t* n0, size_t* n1, size_t k, size_t d, double* mu ){
  size_t l,j, depth, idx;
  for (l=0; l< k; l++){
    for(j=0; j<d; j++){
      idx=l*d+j;
      depth= n0[idx]+n1[idx];
      if (depth>0) mu[idx]=n1[idx]/depth; else mu[idx]=0.5;
      if ( 0 == mu[idx] ) mu[idx]=0.1;
      if (1 == mu[idx] ) mu[idx]=0.9;
    }
  }
}


void read_n_d_from_file(FILE* infile,  size_t* n, size_t* d){
  char* line = malloc(MAXLINE);
  line = fgets(line, MAXLINE, infile);
  if ( NULL == line) {
    fprintf(stderr, "fatal: empty input file\n");
    exit(EXIT_FAILURE);
  }
  int scstatus = sscanf(line, "#@N:%lu", n);
  if (scstatus != 1) {
    fprintf(stderr, "malformed matrix file (no #@N):%d\n", scstatus);
    exit(EXIT_FAILURE);
  }
  line = fgets(line, MAXLINE, infile);
  scstatus = sscanf(line, "#@D:%lu", d);
  if (scstatus != 1) {
    fprintf(stderr, "malformed matrix file (no #@D):%d\n", scstatus);
    exit(EXIT_FAILURE);
  }
  free(line);
  return;
}

void read_matrix_from_file(FILE* infile, size_t n, size_t d,
			   char** rnames, uint32_t* gpos,
			   int32_t* m){
  size_t readidx, gposidx ;
  char* line = malloc(MAXLINE);
  char* rname = malloc(MAXLINE);
  uint32_t tmpgpos;
  int status;
  uint32_t cgcount=0;
  while(1) {
    line = fgets(line, MAXLINE, infile);
    if ( NULL == line && feof(infile) ) {
      fprintf(stderr, "done\n");
      break;
    }
    if ('#' == line[0]) continue;
    sscanf( line,
	    "%lu %s %lu %u %d",
	    &readidx, rname, &gposidx, &tmpgpos, &status );
    if ((status <0) || (status > 1)) {
      fprintf(stderr, "invalid input matrix. status=%u\n", status);
      exit(EXIT_FAILURE);
    }
    if ( NULL == rnames[readidx] ) {
      rnames[readidx] = (char*)malloc(MAXLINE);
      strcpy(rnames[readidx], rname);
    }
    gpos[gposidx]=tmpgpos;
    m[readidx * d + gposidx] = status;
    cgcount++;
  }
  free(rname);  
  free(line);
  return;
}

double logsumexp(double *v, size_t k){
  size_t i;
  double s = 0, m = -FLT_MAX;
  for(i = 0; i < k; i++){
    if ( v[i] > m ) m = v[i];
  }
  for(i = 0; i < k; i++){
    s += exp( v[i] - m );
  }
  return m + log(s);
}

double log_bernoulli_term(int32_t x, double y){
  double res;
  if ( x > 1 ) x=1;
  if (x < 0 ) x=0;
  if ( y > 1 ) y = 1;
  if (y < 0 ) y = 0;
  if (( (0 == x) && (0 == y) ) || ( (1 == x) && (1 == y) )){
    res =0.;
  } else {
    res = x * log(y) + (1 - x) * log1p(-y); 
  }
  return res;
}

double xlogy(double x, double y){
  double res;
  if ((0 == x) && (0 == y)){
    res =0.;
  } else {
    res= x * log(y);
  }
  return res;
}

void mu0( int32_t* x, double* mu0,
	  size_t n, size_t d ){
  /* 
     compute mean column by column, skipping NAs.
     (modifies mu0 in place).
   */
  size_t i,j;
  unsigned int *depth = malloc(d*sizeof(unsigned int));
  for(j = 0; j < d; j++){
    mu0[j]=0; depth[j]=0;
    for (i = 0; i < n; i++){
      if (-1 == x[i * d + j]) continue;
      mu0[j] += x[i * d + j];
      depth[j] += 1;
    }
    mu0[j] = mu0[j] / depth[j];
  }
  free(depth);
}

double ll0(int32_t* x, double* mu0, size_t n, size_t d){
  double s = 0.0;
  for(size_t i = 0 ; i < n; i++){
    for(size_t j =  0; j < d; j++){
      if ( ( -1 == x[i * d + j] ) ) continue;
      s += xlogy(x[i * d + j], mu0[j]) +
	xlogy((1 - x[i * d + j]), (1 - mu0[j]));
    }
  }
  return s;
}

double logmvbern(int32_t* x, double* mu, size_t d){
  double logp = 0;
  size_t j;
  for(j = 0; j < d; j++){
    if ( (x[j] < 0)   ) continue;
    logp += log_bernoulli_term(x[j], mu[j]);
  }
  return logp;
}

double logphi(double pi, int32_t* x, double* mu, size_t d){
  return log(pi) + logmvbern(x, mu, d);
}

/* loglikelihood log(p(x|theta)) */
double llbern(int32_t* x, double* mu, double* pi,
	      size_t n, size_t d, size_t k){
  size_t i,l;
  double ll=0.;
  double* v = malloc(k*sizeof(double));
  for (i = 0; i < n; i++){
    for (l = 0; l < k; l++){
      v[l] = logphi(pi[l], x + i * d, mu + l * d, d); 
    }
    ll += logsumexp(v, k);
  }
  free(v);
  return ll;
}

double log2(double x){
  return log10(x) / 0.3010299957;
}

void depth_mean_entropy(int32_t *x,
			unsigned int* depth, double* mean, double* entropy,
			unsigned int* n0,
			size_t n, size_t d){
  size_t i,j;
  double p;
  unsigned int n1;
  for(j = 0; j < d; j++){
    n1 = 0; n0[j]=0; depth[j]=0;
    for(i = 0; i < n; i++){
      if ( 0 == x[i * d + j] ) {n0[j]++; depth[j]++;}
      if ( 1 == x[i * d + j] ) {n1++; depth[j]++;} 
    }
    if (  ( 0 == n0[j] ) ||  ( 0 == n1 ) )  {
      mean[j] = -1;
      entropy[j] = -1;
    } else {
      p = (double) n1 / (n0[j] + n1);
      entropy[j] = -p * log2(p) - (1 - p) * log2(1 - p);
      mean[j] = p;
    }
  }
}

void bernoulli_gamma(double* pi, int32_t* x, double* mu, 
		     size_t n, size_t d, size_t k, double* g){
  /* use mu and pi to compute gamma (double* g) */
  size_t i,l;
  for( i = 0; i < n; i++ ){
    for( l = 0; l < k; l++ ){
      g[i * k + l] = log(pi[l]) + logmvbern(x + i * d, mu + l * d, d);
    }
    double denom = logsumexp(g + i * k, k);
    for( l = 0; l < k; l++ ){
      g[i * k + l] = exp(g[i * k + l] - denom);
    }
  }
}

void update_mu(int32_t* x,  double*g, 
	       size_t n, size_t d, size_t k,
	       double* mu, double* pop ){
  /* use g to compute  mu and pop */
  size_t i,j,l;
  double xij;
  // sum columns of gamma
  for(l = 0; l < k; l++){
    pop[l] = 0.0;
    for(i = 0; i < n; i++){
      pop[l] += g[i * k + l];
    }
  }
  double tmp;
  for(l = 0; l < k; l++){
    for(j = 0; j < d; j++){
      tmp = 0;
      for(i = 0; i < n; i++){
	xij = x[i * d + j];
	if ( -1 == xij ) xij = mu[l * d + j];
	tmp += g[i * k + l] * xij ;
      }
      mu[l*d+j]=tmp;
    }
    for(j = 0; j < d; j++) mu[l*d+j] = mu[l*d+j]/pop[l];
  }
}

void bin_meth_matrix( int32_t* m, size_t n, size_t d, int thres){
  /*
    binarize the methylation matrix (Nanopore gives a number from 0 to 255).
    * arguments
    m matrix (changed in place)
    n number of rows
    d number of columns
    thres threshold for binarization
   */
  size_t i,j;
  int me;
  for (i = 0; i < n; i++){
    for(j = 0; j < d; j++){
      me  = m[i * d + j];
      if (-1 == me ) continue;
      if ( me > thres) {
	m[i * d + j] = 1;
      } else   m[i * d + j] = 0;
    }
  }
}


double dist_to_centroid(int32_t* x, double* mu, size_t d, size_t i, size_t l){
  /* 
     computes the mean squared distance of x_i from centroid mu_l.
     Only using observed positions
  */
  size_t j;
  double dist=0;
  int el; unsigned int nobs=0;
  for(j=0; j < d; j++){
    el = x[i * d + j];
    if (-1 == el) continue;
    dist += ((double)el - mu[l * d + j]) * ((double)el - mu[l * d + j]);
    nobs++;
  }
  return (dist / nobs);
}

void kmeans(int32_t* x, double* mu, size_t n, size_t d, size_t k,
	    unsigned int maxit, unsigned int* cl){
  size_t l;
  size_t i, j;
  unsigned int it=0;
  unsigned int* pop = calloc(k, sizeof(unsigned int));
  double dtc, mindtc, ssdold=FLT_MAX, ssd;
  while (it < maxit){
    ssd = 0;
    // init pop
    for(l = 0; l < k; l++) pop[l]=0;
    //
    for(i = 0; i < n; i++){
      mindtc = FLT_MAX;
      for(l = 0; l < k ; l++){
	dtc = dist_to_centroid(x, mu, d, i, l); 
	if ( dtc < mindtc){
	  mindtc = dtc; cl[i] = l; 
	}
      }
      ssd += mindtc;
      fprintf(stdout, "i=%lu  ssd=%lf mindtc=%lf\n",i, ssd, mindtc);
      pop[cl[i]]++;
    }
    for (l = 0; l < k; l++){
      for (j = 0 ; j < d; j++){
	mu[l*d+j] = 0;
	for(i=0; i<n;i++){
	  if (cl[i] != l) continue;
	  if ( -1 == x[i * d + j] ) continue;
	  mu[l * d + j] += x[i * d + j];
	}
      }
    }
    for (l = 0; l < k; l++){
      for (j = 0; j < d; j++){
	mu[l * d + j] = mu[l * d + j] / pop[l];
      }
    }
    double dssd = fabs(ssd - ssdold)/ssdold; 
    fprintf(stdout, "kmeans::iteration %u ssd:%.3g dssd:%.3g\n",
	    it, ssd, dssd);
    if ( dssd < 1e-4) break;
    ssdold = ssd;
    it++;
  }
  free(pop);
}

double ln_beta_function(int alpha, int beta){
    return lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta);
}

double h(int a1, int b1, int a2, int b2){
    double tmp;
    tmp = ln_beta_function((a1 + a2), (b1 + b2)) -
        ln_beta_function(a1, b1) -
        ln_beta_function(a2, b2);
    return exp(tmp);    
}

double g (int a1, int b1, int a2, int b2){
    double res;
    if ( ( a1 == a2 ) && ( b1 == b2 ) ) return 0.5;
    if (a1 > a2){ 
      res = g( ( a1 - 1), b1, a2, b2 ) + 
	h ((a1 - 1), b1, a2, b2) / ( a1 - 1 );  
      return res;
    }
    if (a2 > a1){
      res =  g (a1, b1, (a2 - 1), b2 ) - 
            h (a1, b1, (a2 - 1), b2 ) /  ( a2 - 1 ) ;
      return res;
    }
    if  ( b1 > b2 ) {
      res =  g (a1 , (b1 - 1) ,a2, b2 ) - 
	h (a1 , (b1 - 1) ,a2, b2) /  ( b1 - 1 ) ;
      return res;
    }
    if  ( b2 > b1 ) {
      res =  g( a1, b1, a2, ( b2 - 1 ) ) + 
	h (a1, b1, a2, ( b2 - 1 ) ) /   ( b2 - 1 );
      return res;
    }
    return -1;
}

void r_methyldiff(int* nc1, int* c1,
		 int* nc2, int* c2, double* res){
  *res = g(*nc1 + 1 , *c1 + 1 , *nc2 + 1 , *c2 + 1 );
}

double methyldiff(int nc1, int c1, int nc2, int c2){
  double res = g(nc1 + 1 , c1 + 1 , nc2 + 1 , c2 + 1 );
  return res;
}
