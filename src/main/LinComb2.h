//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef LINCOMB2_H
#define LINCOMB2_H

#include "CustomScore.h"
#include "DataDescrip.h"
#include "dtree.h"
class RegionPrediction;

class LinComb2 : public CustomScore {
 public:
  LinComb2() : CustomScore(EvScores::eNumBnd) { }
  void createTree(const string&); 
  static void readDiscreteVals(const string&);
  static void newInstance() {
    assert(!LinComb2::_instance);
    LinComb2::_instance = new LinComb2();
  }
  static LinComb2* getInstance() {
    assert(LinComb2::_instance);
    return LinComb2::_instance; 
  }
  static void deleteInstance() {
    assert(LinComb2::_instance);
    delete LinComb2::_instance;
  }
  CustomScore::LabelScore score(const RegionPrediction&, const string&) const;
  void score(vector<RegionPrediction*>&, const string&) const;
  double getProb(const AnnotationItem&, EvScores::Bnd_t, dsu::Strand_t) const;

 protected:
  EvScores scoreHelp(const vector<RegionPrediction*>&, const string&, dsu::Strand_t) const;
  EvScores scoreHelpOpp(const vector<RegionPrediction*>&, const string&, dsu::Strand_t) const;

  AnnotationScore getProbHelp(const RegionPrediction&, EvScores::Bnd_t, dsu::Strand_t, int) const;
  AnnotationScore getRawProb(const RegionPrediction&, EvScores::Bnd_t, dsu::Strand_t) const;

 private:
  static LinComb2 * _instance;
};

#endif //LINCOMB2_H
