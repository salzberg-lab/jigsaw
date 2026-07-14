//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef ANNOTATIONMATRIX_H
#define ANNOTATIONMATRIX_H

#include <string>
#include <list>
#include "AnnotationScore.h"
//#include "OutputModel.h"

class SeqLoc;
class AnnotationItem;
class DnaStr;
class UnitScore;
class Compare;
class CustomScore;

using std::list;

/**
 * <b>Purpose:</b><br>
 * A glorified linked list of AnnotationItem objects<br> 
 * class acts as a wrapper container class.
 */
class AnnotationMatrix {
 public:
  typedef list<AnnotationItem*>::iterator iterator;
  typedef list<AnnotationItem*>::const_iterator const_iterator;

  AnnotationMatrix();
  ~AnnotationMatrix();

  void del();
  /**
   * finds the best starting point to begin the traceback in the matrix.<br> 
   */
  //trace_back_idx initTraceBack(const CustomScore&, const string&,int,bool,bool = false) const; 
/**
 * <b>Precond:</b><BR>
 * AnnotationItems are inserted in
 * ascending order based on the 5 prime end<br>
 *
 **/
  void insert(AnnotationItem*);

  iterator begin();
  iterator end();
  
  unsigned size() const { return _data.size(); }

  const_iterator begin() const;
  const_iterator end() const;

  bool empty() const { return _data.empty(); }

 protected:
  list<AnnotationItem*> _data;

 private:
  const Compare* _cmp; 
  /**
   * @link aggregation
   * @supplierCardinality 0..* 
   */
  /*# AnnotationItem lnkAnnotationItem; */
};

#endif // ANNOTATIONMATRIX_H
