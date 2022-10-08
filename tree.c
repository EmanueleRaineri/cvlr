
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

struct pos_tnode{
  uint64_t pos;
  unsigned int count;
  struct pos_tnode *left;
  struct pos_tnode *right;
  uint64_t index;
};

void pos_treeprint(struct pos_tnode *p){
  if (p != NULL) {
    pos_treeprint(p->left);
    printf("@%lu %u %lu\n", p->pos, p->count, p->index );
    pos_treeprint(p->right);
  }
}


// sort genomic position by walking the tree in preorder
uint64_t pos_index(struct pos_tnode *p, uint64_t idx0){
  if ( p != NULL ){
    idx0 = pos_index(p->left, idx0);
    p->index = idx0 ;
    idx0 = pos_index(p->right, idx0+1);
  }
  return idx0;
}

// returns the index of a gpos
int pos_index_of(struct pos_tnode *p, uint64_t pos){
  int res;
  if ( p != NULL ) {
    if ( pos == p->pos ) res = p-> index;
    if ( pos < p->pos ) res = pos_index_of(p->left, pos);
    if (pos > p->pos ) res = pos_index_of(p->right, pos);
  } else {
    res = -1;
  }
  return res;
}

struct pos_tnode *talloc(void){
  return (struct pos_tnode *) malloc(sizeof(struct pos_tnode));
}

// here p is the root
struct pos_tnode *pos_addtree(struct pos_tnode *p, unsigned int pos){
  if (NULL == p){
    //fprintf(stderr, "allocating\n");
    p = talloc();
    p->pos = pos;
    p->count = 1;
    p->left = p->right = NULL;
  }
  else if ( pos == p->pos ) p->count++; 
  else if ( pos < p->pos) p->left = pos_addtree( p->left, pos );
  else p->right = pos_addtree( p->right, pos );
  return p;
}

void unroll_pos_tree(struct pos_tnode *p, uint64_t* sorted_pos){
  if (NULL == p) return;
  sorted_pos[p->index] = p->pos;
  unroll_pos_tree(p->left, sorted_pos);
  unroll_pos_tree(p->right, sorted_pos);
}


int count_pos_tnodes(struct pos_tnode *p){
  if (NULL == p){return 0;}
  else {
    return 1 + count_pos_tnodes(p->left) + count_pos_tnodes(p->right);
  }
}

uint64_t min_pos(struct pos_tnode *p){
  if (NULL == p->left) return p->pos;
  else return min_pos(p->left);
}

uint64_t max_pos(struct pos_tnode *p){
  if (NULL == p->right) return p->pos;
  else return max_pos(p->right);
}


/* rname */

struct rname_tnode{
  char *rname;
  unsigned int count;
  struct rname_tnode *left;
  struct rname_tnode *right;
  uint64_t index;
};

void rname_treeprint(struct rname_tnode *p){
  if (p != NULL) {
    rname_treeprint(p->left);
    printf("@%s %u %lu\n", p->rname, p->count, p->index );
    rname_treeprint(p->right);
  }
}

uint64_t rname_index(struct rname_tnode *p, uint64_t idx0){
  if ( p != NULL ){
    idx0 = rname_index(p->left, idx0);
    p->index = idx0 ;
    idx0 = rname_index(p->right, idx0+1);
  }
  return idx0;
}

int rname_index_of(struct rname_tnode *p, char* rname){
  int res;
  if ( p != NULL ) {
    int cmp = strcmp(rname, p->rname);
    if ( 0 == cmp) res= p-> index;
    if ( cmp < 0 ) res = rname_index_of(p->left, rname);
    if (cmp >0 ) res = rname_index_of(p->right, rname);
  } else {
    res = -1;
  }
  return res;
}

struct rname_tnode *rname_talloc(void){
  return (struct rname_tnode *) malloc(sizeof(struct rname_tnode));
}

struct rname_tnode *rname_addtree(struct rname_tnode *p, char* rname){
  if (NULL == p){
    p = rname_talloc();
    p->rname = strdup(rname);
    p->count = 1;
    p->left = p->right = NULL;
  }
  else if ( 0 == strcmp( rname, p-> rname) ) p->count++; 
  else if ( strcmp(rname, p-> rname) <0 )  p->left = rname_addtree(p->left, rname);
  else p->right = rname_addtree(p->right, rname);
  return p;
}

void unroll_rname_tree(struct rname_tnode *p, char** sorted_rname){
  if (NULL == p) return;
  sorted_rname[p->index] = p->rname;
  unroll_rname_tree(p->left, sorted_rname);
  unroll_rname_tree(p->right, sorted_rname);
}


int count_rname_tnodes(struct rname_tnode *p){
  if (NULL == p){return 0;}
  else {
    return 1 + count_rname_tnodes(p->left) + count_rname_tnodes(p->right);
  }
}
