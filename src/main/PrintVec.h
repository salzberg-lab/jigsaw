//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef PRINTVEC_H
#define PRINTVEC_H

#include "CustomScore.h"
#include "Options.h"
#include <fstream>
#include <stdexcept>

using std::vector;

class RegionPrediction;
class DataSetAbstractProduct;

class PrintVec : public CustomScore {
 public:
  static void newInstance(const Options& opts, 
			  const DataSetAbstractProduct& dsap) {
    assert(!_instance);
    _instance = new PrintVec(opts,dsap);
  }
  static PrintVec* getInstance() {
    assert(_instance);
    return _instance; 
  }
  static void deleteInstance() {
    assert(_instance);
    delete _instance;
  }
  CustomScore::LabelScore score(const RegionPrediction&, const string&) const;
  void score(vector<RegionPrediction*>&, const string&) const;
  bool parseOnly() const { return true;} 
  int getEnd5() const { return _end5; }
  int getEnd3() const { return _end3; }
  
  double getProb(const AnnotationItem&, EvScores::Bnd_t, dsu::Strand_t) const { return -1; }
  
 public:
  vector<std::ostream*> _ostrmVec;

 protected:
  virtual ~PrintVec() throw() { 
    for(unsigned i = 0; i < _ostrmVec.size(); ++i) {
      delete _ostrmVec[i];
    }
  }
  PrintVec(const Options&, const DataSetAbstractProduct&);

 private:
  void scoreHelp(const RegionPrediction&, const string&, dsu::Strand_t) const;

 private:
  static PrintVec * _instance;

 private:
  PrintVec(const PrintVec&);
  const DataSetAbstractProduct& _geneList;
  int _end5, _end3;
  
};

#endif // PRINTVEC_H
