//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef LINCOMB_H
#define LINCOMB_H

#include "CustomScore.h"
#include "DataDescrip.h"
#include "dtree.h"
class RegionPrediction;

class LinComb : public CustomScore {
 public:
  LinComb() : CustomScore(EvScores::eNumBnd) { }
  void createTree(const string&); 
  static void readDiscreteVals(const string&);
  static void newInstance() {
    assert(!LinComb::_instance);
    LinComb::_instance = new LinComb();
  }
  static LinComb* getInstance() {
    assert(LinComb::_instance);
    return LinComb::_instance; 
  }
  static void deleteInstance() {
    assert(LinComb::_instance);
    delete LinComb::_instance;
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
  static LinComb * _instance;
};

#endif //LINCOMB_H
