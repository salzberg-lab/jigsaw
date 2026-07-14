//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef ANNOTATIONITEM_H
#define ANNOTATIONITEM_H

#include "AnnotationScore.h"
#include "RegionPrediction.h"
#include "MatrixItemContents.h"
#include "AnnotationMatrix.h"
#include "localtypedef.h"
#include <vector>

class DnaStr;
class CustomScore;

/**
 * <b>Purpose:</b><br>
 * The annotation matrix is made up of a linked list of this object
 * One instance of this class records all the necessary information<br>
 * for one unique genomic region.  This is a manager class that delegates<br>
 * a lot of responsibility to other objects, for scoring and gene prediction evidence <br>
 * management.  The actual code for recording the optimal parse is done here, nested in 
 * dynamicProgrammingUpdateScore<br>
 * this class contains a RegionPrediction object that refers to the associated evidence
 * and an AnnotationScoreMatrix which records the possible 6 frame scores and their 
 * predecessor pointers
 */
class AnnotationItem {
 public:
  typedef pair<const AnnotationItem*,AnnotationScore> apair_t; 
  typedef pair<int,apair_t> bpair_t;
  typedef pair<dsu::Strand_t,int> strnd_frm_t;
  struct tback {
    tback() : _item(NULL) { }
    tback(const AnnotationItem* i, strnd_frm_t strf, dsu::side_t sd, int type) : _item(i),
	 _strf(strf.first,strf.second), _side(sd), _type(type) { }
    int get_type() const { return _type; }
    const AnnotationItem* _item;
    strnd_frm_t _strf;
    dsu::side_t _side;
    int _type;
  };
public:
  AnnotationItem(const std::vector<RegionPrediction*>&,const DnaStr&, const DnaStr&);
  ~AnnotationItem() { }

  // Accessors
  pair<AnnotationScore,MatrixItemContents>
    endScore(const CustomScore&,const AnnotationScore&,bool) const;
  int getEnd5AdjForStop(int model, dsu::Strand_t strnd) const;
  int getEnd3AdjForStop(int model, dsu::Strand_t strnd) const;
  dsu::Strand_t getStrnd() const { return _strnd; }
  int getEnd5() const { return _seq.getEnd5(); }
  int getEnd3() const { return _seq.getEnd3(); }
  int length() const { return _seq.length(); }
  const RegionPrediction& getRegionPrediction() const { return _seq; }
  const DnaStr& getSeq(dsu::Strand_t strnd) const 
    { assert(strnd != dsu::eEither); return strnd==dsu::ePos? _posSeq : _negSeq; }

  tback
  getPred(int,dsu::Strand_t,dsu::side_t,int) const;

  void
  print(std::ostream&,const CustomScore&) const;

  // Mutators
  /**
   * <b>Precond:</b><br>
   * Takes as input a valid iterator reference to a downstream<br>
   * previously scored region<br>
   * <b>Postcond:</b><br>
   * Score the connection of "this" region with the input region<br>
   * update the best score for "this" region.<br>
   * This can be a fairly involved process, becuase it checks for legal<br>
   * exon types, and scores only legal orfs for each of the possible 3 <br>
   * reading frames.<br> 
   */
  void dpUpdate(const AnnotationItem&, std::vector<scr_pr>&,bool, const CustomScore&);
  void dpUpdate(std::vector<scr_pr>&, bool, bool, const CustomScore&);
  AnnotationScore computeNonCodingScore() const;

  void
  update(const MatrixIdx&,const MatrixIdx&, const int, const AnnotationScore&, const AnnotationItem& ); 
  void
  update(const MatrixIdx&,const MatrixIdx&, const int, const AnnotationScore&); 

  static bool 
  checkOrf(const AnnotationItem*, const MatrixIdx&,int,int,const DnaStr&);

  pair<int,dsu::Strand_t>
    getPrevExonType(const MatrixIdx&) const;

  const MatrixItemContents& getBest(const MatrixIdx&) const;
  
  const RegionPrediction* getLeft() const { return _left; }
  const RegionPrediction* getRght() const { return _rght; }

private:
  MatrixItemContents& getBest(const MatrixIdx&);

  void
    dpUpdate(const AnnotationItem& litem, dsu::Strand_t lstrnd,
	     std::vector<scr_pr>& inv_scores,
	     bool isContig, dsu::Strand_t rstrnd) ;

  void
  updateScore(const AnnotationItem&, bool, std::vector<scr_pr>& inv_scores, int frame,dsu::Strand_t);

  void setStrnd(dsu::Strand_t strnd) { _strnd = strnd; }
  
private:
  RegionPrediction& _seq;
  const DnaStr& _posSeq, _negSeq;
  dsu::Strand_t _strnd;
  std::vector<MatrixItemContents> _best;
  const RegionPrediction* _left;
  const RegionPrediction* _rght;
};

#endif //ANNOTATIONITEM_H
