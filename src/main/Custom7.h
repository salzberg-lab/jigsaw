//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef CUSTOM7_H
#define CUSTOM7_H

#include "CustomScore.h"
#include "DataDescrip.h"
#include "dtree.h"
class RegionPrediction;

class Custom7 : public CustomScore {
 public:
  Custom7() : CustomScore(EvScores::eNumBnd), _dTree(EvScores::eNumBnd) { }
  void createTree(const string&); 
  static void readDiscreteVals(const string&);
  static void newInstance() {
    assert(!Custom7::_instance);
    Custom7::_instance = new Custom7();
  }
  static Custom7* getInstance() {
    assert(Custom7::_instance);
    return Custom7::_instance; 
  }
  static void deleteInstance() {
    assert(Custom7::_instance);
    delete Custom7::_instance;
  }
  CustomScore::LabelScore score(const RegionPrediction&, const string&) const;
  void score(vector<RegionPrediction*>&, const string&) const;
  //double getProbInterval(const RegionPrediction*, const RegionPrediction*, dsu::Strand_t) const;
  double getProb(const AnnotationItem&, EvScores::Bnd_t, dsu::Strand_t) const;

 protected:
  ~Custom7() { 
    for(unsigned i = 0; i < _dTree.size(); ++i) {
      delete _dTree[i];
    }
  }
  EvScores scoreHelp(const vector<RegionPrediction*>&, const string&, dsu::Strand_t) const;
  //vector<FeatVec*> calcVec(const RegionPrediction*, dsu::Strand_t,bool) const;

 private:
  vector<DTree*> _dTree;

 private:
  static Custom7 * _instance;
};

#endif //CUSTOM7_H
