//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 1999, 2000, and 2001 The Institute for Genomic Research
#if _RUNTIME_DEBUG >= 4
extern bool MYDEBUG;
#endif

/****************************************************************/
/****************************************************************/
/* File Name :	tree_util_prob.c				*/
/* Contains modules : 	read_tree				*/
/*			read_subtree				*/
/*			read_hp					*/
/*			read_header				*/
/*			isleftchild				*/
/*			isrightchild				*/
/* Uses modules in :	oc1.h					*/
/*			util.c					*/ 
/* Is used by modules in :	substringscores.c		*/
/* Remarks       : 	These routines are mainly used to read	*/
/*			a decision tree from a file, and to     */
/*                      write a tree to a file.			*/
/****************************************************************/		
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "moc1.h"

extern int no_of_categories;

struct tree_node *extra_node;
char train_data[LINESIZE];


/* prototypes */
void read_subtree(struct tree_node *, FILE *, int , int& );
struct tree_node *read_hp(FILE *, int);
int read_header(FILE *, int);
int isleftchild(struct tree_node *, struct tree_node *);
int isrightchild(struct tree_node *x, struct tree_node *y);

/************************************************************************/
/* Module name : vector							*/
/* Functionality :	Allocates a 1-D float array, whose indices	*/
/*			range from "nl" through "nh", and returns a	*/
/*			pointer to this array.				*/
/* Parameters :	nl,nh : lowest and highest indices.			*/
/* Calls modules :	error						*/
/************************************************************************/
float *fvector(int nl, int nh) {
  float *v;
  
  v=(float *)malloc((unsigned)(nh-nl+1)*sizeof(float));
  if (v==NULL) error("Vector : Memory allocation failure.");
  return (v-nl);
}
     
/************************************************************************/
/* Module name :	read_tree					*/ 
/* Functionality : High level routine for reading in a decision tree	*/
/* Parameters :	decision_tree : Name of the file in which the tree is	*/
/*		stored.							*/
/* Returns :	pointer to the root node of the tree.			*/
/* Calls modules :	read_subtree					*/
/*			read_header					*/
/*			read_hp						*/
/*			error (util.c)					*/
/* Is called by modules :	main (mktree.c)				*/
/*				main (gen_data.c)			*/
/* Remarks : 	It is assumed that the file "decision_tree" is		*/
/* 		written in a format similar to the output of the	*/
/*		write_tree module. A sample decision tree is given in   */
/*              the file sample.dt.				        */ 
/************************************************************************/
struct tree_node *read_tree(const char *decision_tree, int num_dim, int& labels)
{
  FILE *dtree;
  struct tree_node *root;

  if ((dtree = fopen(decision_tree,"r")) == NULL) {
	return NULL;
  }
#if _RUNTIME_DEBUG >= 4
  if(MYDEBUG) {
    int stop = 0;
  }
#endif
  if ( !(read_header(dtree,num_dim))) {
	 printf("dt: %s\n",decision_tree);
    error("Decision tree invalid/absent.");
  }
  
  if ((root = read_hp(dtree,num_dim)) == NULL) {
	 printf("dt: %s\n",decision_tree);
    error("Decision tree invalid/absent.");
  }
  root->parent = NULL;
  extra_node = NULL;
  labels = 0;
  read_subtree(root,dtree,num_dim,labels);
  fclose(dtree);
  return(root);
}


/************************************************************************/
/* Module name :	read_subtree					*/ 
/* Functionality :	recursively reads in the hyperplane, left 	*/
/*			subtree and the right subtree at a node of 	*/
/*			the decision tree. 				*/
/* Parameters :	root : node, the subtree at which is to be read.	*/
/*		dtree: file pointer where the tree is available.	*/
/* Returns :	nothing.						*/
/* Calls modules :	read_subtree					*/
/*			read_hp						*/
/*			isleftchild					*/
/*			isrightchild					*/
/* Is called by modules :	read_tree				*/
/*				read_subtree				*/
/* Important Variables Used :	extra_node 				*/
/*	Hyperplanes are read from the file "dtree" in the order "parent,*/
/*	left child, right child". In case a node does not have either a	*/
/* 	left child or a right child or both, this routine reads one 	*/
/*	hyperplane before it is needed. Such hyperplanes, that are read	*/
/*	before they are needed, are stored in extra_node.		*/
/************************************************************************/
void read_subtree(struct tree_node *root, FILE *dtree, int num_dim, int& label)
{
  struct tree_node *cur_node;
  
  if (extra_node != NULL) {
      cur_node = extra_node;
      extra_node = NULL;
  } else {
    cur_node = read_hp(dtree, num_dim);
  }
  
  if (cur_node == NULL) {
    root->lft_label = label;
    ++label;
    root->rght_label = label;
    ++label;
    return;
  }
  if (isleftchild(cur_node,root)) {
      cur_node->parent = root;
      root->left = cur_node;
      
      read_subtree(cur_node,dtree, num_dim,label);
      if (extra_node != NULL) {
	  cur_node = extra_node;
	  extra_node = NULL;
      } else {
	cur_node = read_hp(dtree, num_dim);
      }
      if (cur_node == NULL) {
	root->rght_label = label;
	++label;
	return;
      }
  } else {
    root->lft_label = label;
    ++label;
  }

  if (isrightchild(cur_node,root)) {
      cur_node->parent = root;
      root->right = cur_node;
      read_subtree(cur_node,dtree, num_dim, label);
  } else {
    root->rght_label = label;
    ++label;
    extra_node = cur_node;
  }
}

/************************************************************************/
/* Module name : read_hp						*/
/* Functionality :	Reads a hyperplane (one node of the decision	*/
/*			tree).						*/
/* Parameters :	dtree : file pointer to the decision tree file.		*/
/* Returns : pointer to the decision tree node read.			*/
/* Calls modules :	vector (util.c)					*/
/*			error (util.c)					*/
/* Is called by modules :	read_tree				*/
/*				read_subtree				*/
/* Remarks :	Rather strict adherance to format.			*/
/*		Please carefully follow the format in sample.dt, if	*/
/*		your decision tree files are not produced by "mktree".	*/
/************************************************************************/
struct tree_node *read_hp(FILE *dtree, int no_of_dimensions)
{
  struct tree_node *cur_node;
  float temp;
  char c;
  int i;
  extern float* c_vector(int,int);

  cur_node = new tree_node;
  cur_node->coefficients = new float[no_of_dimensions+1];
  cur_node->left_count = ivector(0,no_of_categories);
  cur_node->right_count = ivector(0,no_of_categories);
  cur_node->lft_prob = 0.0;
  cur_node->rght_prob = 0.0;
  cur_node->lft_lst.clear();
  cur_node->rght_lst.clear();
  cur_node->lft_cnt = cur_node->rght_cnt = 0;
  cur_node->lft_label = cur_node->rght_label = -1;
  
  for (i=1;i<=no_of_dimensions+1;i++) cur_node->coefficients[i-1] = 0;
  
  cur_node->left = cur_node->right = NULL;
  
  while (isspace(c = getc(dtree)));
  ungetc(c,dtree); 
  
#if _RUNTIME_DEBUG >= 4
  if( MYDEBUG ) {
    int stop = 0;
  }
#endif
  if (fscanf(dtree,"%[^' '] Hyperplane: Left = [", cur_node->label) != 1) {
    return(NULL);
  }
  /*
  printf("dt: %s\n",cur_node->label);
  */

  for (i=1;i<no_of_categories;i++) {
    if (fscanf(dtree,"%d,",&cur_node->left_count[i-1]) != 1)
      return(NULL); 
    /*
    printf("left %d: %d\n",i-1,cur_node->left_count[i-1]);
    */
  }
  if (fscanf(dtree,"%d], Right = [",
	     &cur_node->left_count[no_of_categories-1]) != 1) {
    return(NULL); 
  }
  /*
  printf("nc- %d: %d\n",no_of_categories-1,cur_node->left_count[no_of_categories-1]);
  */

  for (i=1;i<no_of_categories;i++) {
    if (fscanf(dtree,"%d,",&cur_node->right_count[i-1]) != 1)
      return(NULL); 
    /*
    printf("right %d: %d\n",i-1,cur_node->right_count[i-1]);
    */
  }
  if (fscanf(dtree,"%d]\n", &cur_node->right_count[no_of_categories-1]) != 1)
    return(NULL); 

  /*
  printf("rc- %d: %d\n",no_of_categories-1,cur_node->right_count[no_of_categories-1]);
  */

  if (!strcmp(cur_node->label,"Root")) strcpy(cur_node->label,"");
  
  while (TRUE)
    {
      if ((fscanf(dtree,"%f %c",&temp,&c)) != 2)
	error("Invalid/Absent hyperplane equation.");
      if (c == 'x')
	{ 
	  if ((fscanf(dtree,"[%d] +",&i)) != 1) 
	    error("Read-Hp: Invalid hyperplane equation.");
	  if (i <= 0 || i > no_of_dimensions+1) 
	    error("Read_Hp: Invalid coefficient index in decision tree.");
	  cur_node->coefficients[i-1] = temp;
	}
      else if (c == '=')
	{
	  fscanf(dtree," 0\n\n");
	  cur_node->coefficients[no_of_dimensions] = temp;
	  break;
	}
    }

  cur_node->no_of_points = cur_node->left_total = cur_node->right_total = 0;
  cur_node->left_cat = cur_node->right_cat = 1;
  for (i=1;i<=no_of_categories;i++)
    {
      cur_node->left_total += cur_node->left_count[i-1];
      cur_node->right_total += cur_node->right_count[i-1];
      if(cur_node->left_count[i-1] > cur_node->left_count[(cur_node->left_cat)-1])
	cur_node->left_cat = i;
      if(cur_node->right_count[i-1] > cur_node->right_count[(cur_node->right_cat)-1])
	cur_node->right_cat = i;
    }
  cur_node->no_of_points =
    cur_node->left_total + cur_node->right_total;
  
  return(cur_node);
}

/************************************************************************/
/* Module name : isleftchild						*/
/* Functionality : 	Checks if node x is a left child of node y.	*/
/*			i.e., checks if the label of node x is the same	*/
/*			as label of y, concatenated with "l".		*/
/* Parameters : x,y : pointers to two decision tree nodes.		*/
/* Returns :	1 : if x is the left child of y				*/
/*		0 : otherwise						*/
/* Is called by modules :	read_subtree				*/
/************************************************************************/
int isleftchild(struct tree_node *x, struct tree_node *y)
{
  char temp[MAX_DT_DEPTH];
  
  strcpy(temp,y->label);
  if (!strcmp(strcat(temp,"l"),x->label)) return(1);
  else return(0);
}

/************************************************************************/
/* Module name : isrightchild						*/
/* Functionality : 	Checks if node x is a right child of node y.	*/
/*			i.e., checks if the label of node x is the same	*/
/*			as label of y, concatenated with "l".		*/
/* Parameters : x,y : pointers to two decision tree nodes.		*/
/* Returns :	1 : if x is the right child of y			*/
/*		0 : otherwise						*/
/* Is called by modules :	read_subtree				*/
/************************************************************************/
int isrightchild(struct tree_node *x, struct tree_node *y)
{
  char temp[MAX_DT_DEPTH];
  
  strcpy(temp,y->label);
  if (!strcmp(strcat(temp,"r"),x->label)) return(1);
  else return(0);
}

/************************************************************************/
/* Module name : read_header						*/
/* Functionality :	Reads the header information in a decision tree	*/
/*			file.						*/
/* Parameters :	dtree : file pointer to the decision tree file.		*/
/* Returns : 	1 : if the header is successfully read.			*/
/*		0 : otherwise.						*/
/* Calls modules : none.						*/
/* Is called by modules :	read_tree				*/
/* Remarks :	Rather strict adherance to format.			*/
/*		Please carefully follow the format in sample.dt, if	*/
/*		your decision tree files are not produced by "mktree".	*/
/************************************************************************/
int read_header(FILE *dtree, int no_of_dimensions)
{
  if ((fscanf(dtree,"Training set: %[^,], ",train_data)) != 1) return(0);
  if ((fscanf(dtree,"Dimensions: %d, Categories: %d\n",
	      &no_of_dimensions,&no_of_categories)) != 2) return(0);
  return(1);
}

/************************************************************************/
/************************************************************************/
