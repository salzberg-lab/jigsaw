//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "moc1.h"
#include <stdio.h>
#include <cstring>

extern int no_of_categories;
void write_header(FILE*,int);
void write_subtree(struct tree_node*,FILE*,int);
void write_hp(struct tree_node *cur_node, FILE* dtree,int);

/************************************************************************/
/* Module name : write_tree                  */
/* Functionality :   High level routine to write a decision tree to  */
/*       a file.                 */
/* Parameters :   root : pointer to the structure containing the root of   */
/*           decision tree.               */
/*    dt_file : Name of the file into which the dt is to be */
/*           stored.                */
/* Returns :   Nothing.                */
/* Calls modules :   error (util.c)             */
/*       write_header               */
/*       write_subtree              */
/* Is called by modules :  main (mktree.c)            */
/*          cross_validate (mktree.c)     */
/* Remarks :   The input/output formats are rather strict at this */
/*    stage for OC1. The read_tree routines                   */
/*    produce run time errors while reading decision trees  */
/*    even if there is the slightest deviation from the  */
/*    format outputted by the "write_*" routines.             */
/************************************************************************/
void write_tree(struct tree_node* root,const std::string& dt_file, int num_dim)
{
  FILE *dtree;

  if ((dtree = fopen(dt_file.c_str(),"w")) == NULL)
    error("Write_Tree: Decision Tree file can not be opened.");

  write_header(dtree,num_dim);
  write_subtree(root,dtree,num_dim);
  fclose(dtree);
}

/************************************************************************/
/* Module name :  write_subtree              */
/* Functionality :   Initiates writing a hyperplane, and recursively */
/*       writes the subtrees on the left and right of the*/
/*       hyperplane.             */
/* Parameters :   cur_node : Pointer to the DT node under consideration.   */
/*    dtree : File pointer to the output file.     */
/* Calls modules :   write_subtree              */
/*       write_hp             */
/* is called by modules :  write_subtree           */
/*          write_tree           */
/************************************************************************/
void write_subtree(struct tree_node* cur_node,FILE* dtree, int num_dim)
{
  if (cur_node == NULL) return;

  write_hp(cur_node,dtree,num_dim);
  write_subtree(cur_node->left,dtree,num_dim);
  write_subtree(cur_node->right,dtree,num_dim);
}


/************************************************************************/
/* Module name :  write_hp             */
/* Functionality : Writes one hyperplane.          */
/* Parameters :   cur_node : Pointer to the DT node under consideration.   */
/*    dtree : File pointer to the output file.     */
/* Is called by modules :  write_subtree           */
/* Important Variables used : cur_node->label : Label is the empty   */
/*          string for the root node. If string x  */
/*          is the label of a node, xl is the label   */
/*          of its left child, and xr of its right */
/*          child. Label plays an important role,  */
/*          because it informs the tree reading */
/*          routines (in classify_util.c) about the   */
/*          structure of the DT.       */
/************************************************************************/
void write_hp(struct tree_node *cur_node, FILE* dtree, int no_of_dimensions) {
  int i;

  if (dtree == NULL) return;

  if (strcmp(cur_node->label,"\0") == 0)
    fprintf(dtree, "Root Hyperplane: ");
  else
    fprintf(dtree, "%s Hyperplane: ",cur_node->label);

  fprintf(dtree,"Left = [");
  for (i=1;i<=no_of_categories;i++)
    if (i == no_of_categories) fprintf(dtree,"%d], ",cur_node->left_count[i-1]);
    else fprintf(dtree,"%d,",cur_node->left_count[i-1]);
  fprintf(dtree,"Right = [");
  for (i=1;i<=no_of_categories;i++)
    if (i == no_of_categories) fprintf(dtree,"%d]\n",cur_node->right_count[i-1]);
    else fprintf(dtree,"%d,",cur_node->right_count[i-1]);
  fprintf(dtree,"%d %d",cur_node->lft_label,cur_node->rght_label);

/*  for (i=1;i<=no_of_dimensions+1;i++)
    if (cur_node->coefficients[i])
      {
   if (i <= no_of_dimensions)
     fprintf(dtree,"%f x[%d] + ",cur_node->coefficients[i],i);
   else
     fprintf(dtree,"%f = 0\n\n",cur_node->coefficients[i]);
      }
*/
  {
    fprintf(dtree,"Left: %f\n",cur_node->lft_prob);
    list<float>::const_iterator iter2 = cur_node->lft_val.begin();
    list<float*>::const_iterator iter = cur_node->lft_lst.begin();
    list<float*>::const_iterator stop = cur_node->lft_lst.end();
    while( iter != stop ) {
      const float* val = *iter;
      for (i=0;i<no_of_dimensions;i++)
	fprintf(dtree,"%.2f ",val[i]);
      fprintf(dtree,"%.2f\n",*iter2);
      ++iter;
      ++iter2;
    }
  }
  {
    fprintf(dtree,"Right: %f\n",cur_node->rght_prob);
    list<float>::const_iterator iter2 = cur_node->rght_val.begin();
    list<float*>::const_iterator iter = cur_node->rght_lst.begin();
    list<float*>::const_iterator stop = cur_node->rght_lst.end();
    while( iter != stop ) {
      const float* val = *iter;
      for (i=0;i<no_of_dimensions;i++)
	fprintf(dtree,"%.2f ",val[i]);
      fprintf(dtree,"%.2f\n",*iter2);
      ++iter;
      ++iter2;
    }
  }

  for (i=1;i<=no_of_dimensions;i++)
    if (cur_node->coefficients[i-1])
     fprintf(dtree,"%f x[%d] + ",cur_node->coefficients[i-1],i);
  fprintf(dtree,"%f = 0\n\n",cur_node->coefficients[no_of_dimensions]);
}

/************************************************************************/
/* Module name :  write_header               */
/* Functionality :   Writes the decision tree header.    */
/* Parameters :   dtree : file pointer to the output file.     */
/* Is called by modules :  write_tree           */
/************************************************************************/
void write_header(FILE* dtree, int no_of_dimensions)
{
  extern int /*no_of_dimensions,*/no_of_categories;
  extern char train_data[LINESIZE];

  if (dtree == NULL) return;

  fprintf(dtree,"Training set: %s, ",train_data);
  fprintf(dtree,"Dimensions: %d, Categories: %d\n",
     no_of_dimensions,no_of_categories);
  fprintf(dtree,"\n\n");
}

