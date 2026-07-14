//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef DTREE_H
#define DTREE_H

#include <string>
#include "FeatVec.h"
#include <vector>
class ostream;
class DataDescrip;

class DTree {
 public:
  /**
   * tree file from oc1
   */
  DTree(const string&,const string&, const DataDescrip&,int);
  double eval(const vector<FeatVec*>&) const;
  int getClusterId(const vector<FeatVec*>&, const vector<double>&) const;
  int getNumClusters() const { return _numClus; }
  struct tree_node** getRoot() { return _root;}
 private:
  unsigned getNumDt() const;

 private:
  
  struct tree_node** _root;
  unsigned _numDt;
  int _numClus;
};

#endif // DTREE_H
