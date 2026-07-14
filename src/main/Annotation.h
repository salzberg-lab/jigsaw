//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef ANNOTATION_H
#define ANNOTATION_H

#include <vector>
#include <list>
#include "SeqLoc.h"
#include "GffFeature.h"
#include "AnnotationItem.h"
#include "AnnotationMatrix.h"
#include "AnnotationScore.h"
#include "DnaStr.h"

class RegionPrediction;
class CustomScore;

/**
 * <br><b>Purpose:</b> <br>
 * Holds all possible gene annotations
 * given a set of input annotations
 * This structure will be used by the dynamic prog.<br>
 * algorithm to record the optimal annotation
 * with respect to the scoring system.<br>
 * This class is mainly a wrapper around the AnnotationMatrix object
 * and is responsible for inserting into the matrix and extracting the 
 * optimal parse from the matrix
 * <br>
 * <b>Collaborators:</b><br>
 * AnnotationMatrix, AnnotationItem, AnnotationScore, Combine<br>
**/

class Annotation {
public:
  typedef std::pair<AnnotationScore,int> aipair_t;

#if _RUNTIME_DEBUG >= 1
  void scoreExon(int,int) const;
#endif
/**
 * initialize variables
 *<br>
 * @param theSequence: reference to the raw dna sequence 
 * @param forceWholeGeneModel: when set to true, while constructing the parse
 * require all coding sequences to begin with a start and end with a stop
 * when false, only check that no inframe stop codons exist
 */
  Annotation(const string& theSequence, bool forceWholeGeneModel, int cutOff) 
  : _theSeq(theSequence), _posSeq(theSequence,dsu::ePos), _negSeq(theSequence,dsu::eNeg),
    _forceWholeGeneModel(forceWholeGeneModel), _cutOff(cutOff) {} 
  ~Annotation() { }

  //void del() { _matrix.del(); }

  const MatrixItemContents
  getBestStart(const CustomScore&) const;
  
  /**
   * <b>Pre Cond:</b><br>
   * RegionPrediction's are inserted 'left to right'.<br>
   * Both sequences are annotated at the same time , so left to right
   * means 5' to 3' on positive side and 3' to 5' on negative side
   * @param const CustomScore& is used to check splice connections and add additional weight
   * <br>
   * <b>Post Cond:</b><br>
   * New piece of evidence is inserted into the annotation matrix<br> 
   */
  void insert(const std::vector<RegionPrediction*>&, const CustomScore&);

  aipair_t 
    scoreHelp(int end5, int end3, dsu::Strand_t strnd, int mid, const CustomScore& cs) const;
  
  /**
   *
   */
  bool forceWholeGeneModel() const { return _forceWholeGeneModel; }

  /**
   *
   */
  bool empty() const { return _matrix.empty(); }

  /**
   *
   */
  const AnnotationMatrix& getMatrix() const { return _matrix; }

  /**
   * <b>Precond:</b><br>
   * Completed scoring of the annotation matrix<br>
   * <br>
   * <b>Postcond:</b><br>
   * Traces back through the matrix connecting the exons<br>
   * and ouputing complete gene models<br> 
   */
  int printOptimalParse(const std::string&,std::ostream&,int,const CustomScore&,const string&) const;

  /**
   *
   */
  void print(std::ostream&,const CustomScore&) const;

  /**
   * <b>Predcond:</b><br>
   * presumes scoring matrix is complete. identical to printOptimalParse
   * but returns results in list instead of directly to ostream
   */
  std::list<GffFeature*>*
    getGeneList(const std::string&, int&, const CustomScore&) const;

/*
   * score a gene list
   */
  void
    scoreGeneList(std::list<GffFeature*>& flst, const CustomScore& cs) const;

protected:
 AnnotationMatrix _matrix;
 const string& _theSeq;
 const DnaStr _posSeq, _negSeq;

private:
 AnnotationScore score(AnnotationItem*, const AnnotationScore&) const;

  /**
   * <b>Precond:</b><br>
   * Fully Scored annotation matrix<br>
   * <br>
   * <b>Post 
   */
 bool _forceWholeGeneModel;

 const int _cutOff;
};

#endif // ANNOTATION_H
