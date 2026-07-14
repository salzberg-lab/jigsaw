//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "moc1.h"
#include <assert.h>
#include <stdio.h>

#define PLEVELH 4
#define PLEVEL 2
#if _RUNTIME_DEBUG >= PLEVEL
extern bool MYDEBUG;
#endif

void init_node(struct tree_node *node) {
  if(!node) return;
  node->lft_prob = 0.0;
  node->rght_prob = 0.0;
  node->lft_lst.clear();
  node->rght_lst.clear();
  node->lft_cnt = 0;
  node->rght_cnt = 0;
  node->lft_label = node->rght_label = -1;
  init_node(node->left);
  init_node(node->right);
}

void init_prob(struct tree_node **trees, int no_of_trees) {
  for (int t=1; t<=no_of_trees; t++) {
    struct tree_node* cur_node = trees[t-1];
    init_node(cur_node);
  }
}

void modify (POINT* point, int no_of_dimensions, double val, 
	       struct tree_node **roots, int no_of_trees) {
  int j, t;
  struct tree_node *cur_node;
  double sum;
  double probtmp = 0;

  /* initialize Exon and Intron prob */
  point->prob[0]=0;
  point->prob[1]=0;
  /* get length and skip if <= 2 */
  for (t=1; t<=no_of_trees; t++) {
    cur_node = roots[t-1];
    while (cur_node != NULL) {
      sum = cur_node->coefficients[no_of_dimensions];
      for (j=1;j<=no_of_dimensions;j++) {
	sum += cur_node->coefficients[j-1] * point->dimension[j-1];
      }
      if (sum < 0) {
	if (cur_node->left != NULL)
	  cur_node = cur_node->left;
	else {
	  /* New for prob. classification, added by Xin Chen */
	  probtmp =
	    (double)cur_node->left_count[(cur_node->left_cat)-1]/
	    cur_node->left_total;
	  if(cur_node->left_cat == 1) {
	    point->prob[0] += probtmp;
	    point->prob[1] += 1-probtmp;
	  } else {
	    point->prob[0] += 1-probtmp;
	    point->prob[1] += probtmp;
	  }
#if _RUNTIME_DEBUG >= PLEVEL
	  if(MYDEBUG) {
   printf("left cat= %d  left_count[%d]= %d no_points= %d\n",
      cur_node->left_cat,
      cur_node->left_cat,
      cur_node->left_count[cur_node->left_cat-1],
      cur_node->left_total);
   for(int j = 0; j < no_of_dimensions; ++j) {
     printf("%.1f ",point->dimension[j]);
   }
   printf("\n");
	  }
#endif
	  cur_node->lft_prob += val;
	  ++cur_node->lft_cnt;
	  cur_node->lft_lst.push_back( point->dimension );
	  cur_node->lft_val.push_back( val );
	  break;
	}
      } else {
	if (cur_node->right != NULL)
	  cur_node = cur_node->right;
	else {
	  /* New for prob. classification, added by Xin Chen */
	  probtmp =
	    (double)cur_node->right_count[(cur_node->right_cat)-1]/
	    cur_node->right_total;
	  if(cur_node->right_cat == 1) {
	    point->prob[0] += probtmp;
	    point->prob[1] += 1-probtmp;
	  } else {
	    point->prob[0] += 1-probtmp;
	    point->prob[1] += probtmp;
	  }
#if _RUNTIME_DEBUG >= PLEVEL
	  if(MYDEBUG) {
   printf("right cat= %d  right_count[%d]= %d no_points= %d\n",
        cur_node->right_cat,
        cur_node->right_cat,
        cur_node->right_count[cur_node->right_cat-1],
        cur_node->right_total);
   for(int j = 0; j < no_of_dimensions; ++j) {
     printf("%.1f ",point->dimension[j]);
   }
   printf("\n");
	  }
#endif
	  cur_node->rght_prob += val;
	  ++cur_node->rght_cnt;
	  cur_node->rght_lst.push_back( point->dimension );
	  cur_node->rght_val.push_back( val );
	  break;
	}
      }
    }
  }

  point->prob[0] /= (double) no_of_trees;
  point->prob[1] /= (double) no_of_trees;
  if(point->prob[0] >= point->prob[1]){
    point->category = 1;
  }
  else{
    point->category = 2;
  }
}

void classify2(POINT *point, int no_of_dimensions, 
	       struct tree_node **roots, int no_of_trees) {
  int j, t;
  struct tree_node *cur_node;
  double sum;

  /* initialize Exon and Intron prob */
  point->prob[0]=0;
  point->prob[1]=0;
  /* get length and skip if <= 2 */
  for (t=1; t<=no_of_trees; t++) {
    cur_node = roots[t-1];
    while (cur_node != NULL) {
      sum = cur_node->coefficients[no_of_dimensions];
      for (j=1;j<=no_of_dimensions;j++) {
	sum += cur_node->coefficients[j-1] * point->dimension[j-1];
#if _RUNTIME_DEBUG >= PLEVELH
	if(MYDEBUG) {
	  printf("%d: %f=%f*%f\n",j-1,sum,cur_node->coefficients[j-1],point->dimension[j-1]);
	}
#endif
      }
#if _RUNTIME_DEBUG >= PLEVELH
      if(MYDEBUG) {
	printf("leftorright? %f\n",sum);
      }
#endif
      if (sum < 0) {
	if (cur_node->left != NULL)
	  cur_node = cur_node->left;
	else {
	  double val = 0;
	  unsigned size = cur_node->lft_cnt;
	  if( size > 0 )
	    val = cur_node->lft_prob / size;
	  point->prob[0] += val;
   /****/
#if _RUNTIME_DEBUG >= PLEVEL
	  if(MYDEBUG) {
   printf("left cat= %d  left_count[%d]= %d no_points= %d\n",
      cur_node->left_cat,
      cur_node->left_cat,
      cur_node->left_count[cur_node->left_cat],
      cur_node->left_total);
	  printf("avg val: %.2f, tot: %d\n",val,size);
   for(int j = 0; j < no_of_dimensions; ++j) {
     printf("%.1f ",point->dimension[j]);
   }
   printf("\n");
	  write_tree(cur_node,"tt.txt",no_of_dimensions);
	  }
#endif
	  break;
	}
      }
      else {
	if (cur_node->right != NULL)
	  cur_node = cur_node->right;
	else {
	  double val = 0;
	  unsigned size = cur_node->rght_cnt;
	  if( size > 0 )
	    val = cur_node->rght_prob / size;
	  point->prob[0] += val;
 
   /****/
#if _RUNTIME_DEBUG >= PLEVEL
	  if(MYDEBUG) {
   printf("right cat= %d  right_count[%d]= %d no_points= %d\n",
        cur_node->right_cat,
        cur_node->right_cat,
        cur_node->right_count[cur_node->right_cat],
        cur_node->right_total);

	  printf("avg val: %.2f, tot: %d\n",val,size);
   for(int j = 0; j < no_of_dimensions; ++j) {
     printf("%.1f ",point->dimension[j]);
   }
   printf("\n");
	  write_tree(cur_node,"tt.txt", no_of_dimensions);
	  }
#endif
	  break;
	}
      }
    }
  }

#if _RUNTIME_DEBUG >= PLEVEL
  if(MYDEBUG) {
    int stop =0;
  }
#endif
  point->prob[0] /= (double) no_of_trees;
  point->prob[1] /= (double) no_of_trees;
  if(point->prob[0] >= point->prob[1]){
    point->category = 1;
  }
  else{
    point->category = 2;
  }
}

int assignLabel(const float *dimension, int no_of_dimensions, 
		struct tree_node *root) { 
  int j, res_label = -1;
  struct tree_node *cur_node = root;
  double sum;

    while (cur_node != NULL) {
      sum = cur_node->coefficients[no_of_dimensions];
      for (j=1;j<=no_of_dimensions;j++) {
	sum += cur_node->coefficients[j-1] * dimension[j-1];
#if _RUNTIME_DEBUG >= PLEVEL
	if(MYDEBUG) {
	  printf("%d: %f=%f*%f\n",j-1,sum,cur_node->coefficients[j-1],dimension[j-1]);
	}
#endif
      }
#if _RUNTIME_DEBUG >= PLEVEL
      if(MYDEBUG) {
	printf("leftorright? %f\n",sum);
      }
#endif
      if (sum < 0) {
	if (cur_node->left != NULL)
	  cur_node = cur_node->left;
	else {
	  res_label = cur_node->lft_label;
#if _RUNTIME_DEBUG >= PLEVEL
	  if(MYDEBUG) {
   printf("left cat= %d  left_count[%d]= %d no_points= %d\n",
      cur_node->left_cat,
      cur_node->left_cat,
      cur_node->left_count[cur_node->left_cat],
      cur_node->left_total);
   /* printf("avg val: %.2f, tot: %d\n",val,size);*/
   for(int j = 0; j < no_of_dimensions; ++j) {
     printf("%.1f ",dimension[j]);
   }
   printf("\n");
	  write_tree(cur_node,"tt.txt",no_of_dimensions);
	  }
#endif
	  break;
	}
      }
      else {
	if (cur_node->right != NULL)
	  cur_node = cur_node->right;
	else {
	  res_label = cur_node->rght_label;
#if _RUNTIME_DEBUG >= PLEVEL
	  if(MYDEBUG) {
   printf("right cat= %d  right_count[%d]= %d no_points= %d\n",
        cur_node->right_cat,
        cur_node->right_cat,
        cur_node->right_count[cur_node->right_cat],
        cur_node->right_total);

   /* printf("avg val: %.2f, tot: %d\n",val,size);*/
   for(int j = 0; j < no_of_dimensions; ++j) {
     printf("%.1f ",dimension[j]);
   }
   printf("\n");
	  write_tree(cur_node,"tt.txt", no_of_dimensions);
	  }
#endif
	  break;
	}
      }
    }
    if( res_label == -1 ) {
      write_tree(root,"tt.txt",no_of_dimensions);
      for(int j = 0; j < no_of_dimensions; ++j) {
	printf("%.1f ",dimension[j]);
      }
      printf("\n");
    }
    assert(res_label != -1);
    return res_label;
}
