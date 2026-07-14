// Copyright (c) 1999, 2000, and 2001 The Institute for Genomic Research


#include <list>
#include <malloc.h>
#include <string>
//#if defined(__cplusplus) 
//extern "C" {
//#endif

#define MAX_COEFFICIENT 	1.0
#define MAX_NO_OF_ATTRIBUTES	100
#define MAX_DT_DEPTH 		50 

#define LINESIZE 		80

#define TRUE 			1
#define FALSE 			0

  using std::list;

typedef struct point
 {
   float *dimension;
   int category;
   double prob[2]; /* New for prob. classification, Added by Xin Chen   */
                   /* probability of point being classified as category */
   
   double val; /*Value obtained by substituting this point in the 
                equation of the hyperplane under consideration.
                This field is maintained to avoid redundant
                computation. */
 }POINT;

struct endpoint
 {
  float x,y;
 };

typedef struct edge
 {
  struct endpoint from,to;
 }EDGE;

struct tree_node
 {
  float *coefficients;
  int *left_count, *right_count;
  int left_total, right_total;  /* New for prob. classification */
                                /* Added by Xin Chen            */
  struct tree_node *parent,*left,*right;
  int left_cat,right_cat;
  char label[MAX_DT_DEPTH];
  float alpha; /* used only in error_complexity pruning. */
  int no_of_points;
  EDGE edge; /* used only in the display module. */
   double rght_prob;
   double lft_prob;
   list<float*> lft_lst;
   list<float> lft_val;
   list<float*> rght_lst;
   list<float> rght_val;
   unsigned lft_cnt,rght_cnt;
   signed lft_label,rght_label;
 };

int error(const std::string&);
void free_ivector(int *,int,int);
void free_vector(float *,int,int);
void free_dvector(double*,int,float);
float myrandom(float,float);
float *c_vector(int,int);
double *dvector(int,int);
int *ivector(int,int);
float average(float*,int);
float sdev(float*,int);
void classify(POINT*,int,struct tree_node**,int);
struct tree_node* read_tree(const char*,int,int&);
void modify (POINT*,int no_of_dimensions, double val, 
	     struct tree_node **roots, int no_of_trees);
int assignLabel(const float *dimension, int no_of_dimensions, 
		struct tree_node *root);
void classify2(POINT *point, int no_of_dimensions, 
	       struct tree_node **roots, int no_of_trees) ;
void write_tree(struct tree_node*,const std::string&,int);
void init_prob(struct tree_node **trees,int);
//#if defined(__cplusplus)
//}
//#endif
