//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef REGIONPREDICTION_H
#define REGIONPREDICTION_H

#include "DataSetAbstractProduct.h"
#include "AnnotationScore.h"
#include "DataSetAbstractFactory.h"
#include "InputData.h"
#include "EvScores.h"
#if __GNUC__ > 4 || ( __GNUC__ == 4 && __GNUC_MINOR__ >= 3 )
 #include <utility>
 using namespace std;
#else
 #include <pair.h>
#endif

class PlaceHolder;
class SpliceSite;
class DnaStr;
class CustomScore;

/**
 * <b>Responsibilities:</b><br>
 * Manage all available evidence for a given sequence region<br>
 * <br>
 * <b>Collaborators:</b><br>
 * DataSetAbstractProduct, SpliceSite, AnnotationScore<br>
**/

class RegionPrediction : public DataSetAbstractProduct {
public:
#if 0
  /**
   * @stereotype factory 
   */
  class GenePredDataSetConcreteFactory : public DataSetAbstractFactory {
  public:

      virtual DataSetAbstractProduct * createDataSetAbstractProduct(const string&);

  private:

      /**
       * @link
       * @shapeType PatternLink
       * @pattern Abstract Factory
       * @supplierRole Product 
       */
      /*# GenePredDataSet _genePredDataSet; */
  };
#endif

  RegionPrediction(const SeqLoc&, const string&);
  ~RegionPrediction();

  // accessors
  static bool
    is_igr(const RegionPrediction&, const RegionPrediction&,dsu::Strand_t,int);
  static bool
    is_intron(const RegionPrediction&, const RegionPrediction&,dsu::Strand_t);

  const string& getSeq() const { return _seq; }
  bool noData()const;
  bool noStrndData(dsu::Strand_t) const;
  virtual void print(std::ostream&, const CustomScore&) const;
  const string& getTypeStr() const;
  bool isValid(EvScores::Bnd_t, dsu::Strand_t) const; 

  /**
   * returns true if region falls on a valid acceptor boundry
   * @param - strnd = dna strnd of interest
   */
  bool isAccpBndry(dsu::Strand_t strnd) const { return _accpBndry[strnd]; }
  /**
   * returns true if region falls on a valid donor boundry
   * @param - strnd = dna strnd of interest
   */
  bool isDonrBndry(dsu::Strand_t strnd) const { return _donrBndry[strnd]; }

  bool isValidAcceptorSignals(dsu::Strand_t strnd) const { return isAccpBndry(strnd) || hasStart(strnd); }
  bool isValidDonorSignals(dsu::Strand_t strnd) const { return isDonrBndry(strnd) || hasStop(strnd); }

  /**
   * <b>Precondition:</b><br>
   * takes as input a valid reference to the genmoic sequence.<br>
   * <br>
   * <b>Postcondition:</b><br>
   * returns a score for each exon model type.  it is possible that <br>
   * none of the exon model type's will be scored, if the available <br>
   * information doesn't "warrant" providing a score, this is basically<br>
   * saying not enough information to indicate an exon of anytype. 
   * CustomScore does all the work now, essentially now it is plugged into this
   * method which use to do all the work. Now you can test different scoring methods
   * through CustomScore implementations without changing the framework code
   */
  void score(const CustomScore&,const string&);

  /** 
   * Needs to track only the 5' 3' region for<br>
   * which all of the available evidence have in common<br>
   * e.g. only the overlapping region<br>
   **/
  void virtual insert(const PlaceHolder*);
  void storeRelaventData(InputData::const_iterator&, const InputData::const_iterator&);
   
  void setScores(dsu::Strand_t, const EvScores&);
  void setScoresOpp(dsu::Strand_t, const EvScores&);
  /**
   * falls on a valid start boundary
   */
  bool hasStart(dsu::Strand_t strnd) const { 
    return _hasStart[strnd]; 
  }

  /**
   * falls on a valid stop boundary
   */
  bool hasStop(dsu::Strand_t strnd) const { 
    return _hasStop[strnd]; 
  }

  /**
   * check to see if interval is right up against the end of the sequence
   */
  bool isEndOfSeq() const {
    //note "getEnd" starts at index 1
    return getEnd3()==(signed)_seq.length();
  }

  void scoreHelp(const string&, dsu::Strand_t);
  inline const EvScores& getScores(dsu::Strand_t strnd) const; 
  inline const EvScores& getScoresOpp(dsu::Strand_t strnd) const; 
  inline const AnnotationScore& getScore(unsigned,dsu::Strand_t strnd) const; 
  inline const AnnotationScore& getScoreOpp(unsigned,dsu::Strand_t strnd) const; 

private:
  AnnotationScore getProbStore(EvScores::Bnd_t,dsu::Strand_t) const;
  void setProbHelp(dsu::Strand_t);
  inline EvScores& getScores(dsu::Strand_t strnd); 
  inline EvScores& getScoresOpp(dsu::Strand_t strnd); 
  void doubleCheck();
  void setDonorEvidenceBoundary(const string&,dsu::Strand_t);
  void setAcceptorEvidenceBoundary(const string&,dsu::Strand_t);
  void setDonorEvidenceBoundaryFromString(const DnaStr&,dsu::Strand_t);
  void setAcceptorEvidenceBoundaryFromString(const DnaStr&,dsu::Strand_t);
  void checkProteinStops(const string& str, dsu::Strand_t strnd);
  inline int getIdx(dsu::Strand_t) const;

  void setProbVal(unsigned,const AnnotationScore&);
private:
  vector<EvScores> _score;
  vector<EvScores> _scoreOpp;

  // caching these results
  bool _donrBndry[2], _accpBndry[2];
  bool _hasStart[2], _hasStop[2];

  int _negCnt, _posCnt;
  const string& _seq;
};

#include "RegionPrediction.hpp" // inlines

#endif // REGIONPREDICTION_H
