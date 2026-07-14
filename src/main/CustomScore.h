//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef CUSTOMSCORE_H
#define CUSTOMSCORE_H

#include "AnnotationScore.h"
#include "DataDescrip.h"
#include "dsu.h"
#include "EvScores.h"
#include "FeatVec.h"
#if __GNUC__ > 4 || ( __GNUC__ == 4 && __GNUC_MINOR__ >= 3 )
 #include <utility>
 using namespace std;
#else
 #include <pair.h>
#endif

class RegionPrediction;

/**
 * Interface, defines the behavior of the scoring function
 * for the client application. Each version of the combiner must derive an 
 * instance of this class
 * This is the labeling function that takes an evidence set E (RegionPrediction)
 * returns a 2-tuple with a weight for coding, and a weight for noncoding
 *
 * The user is expected to access the feature maps which map evidence into
 * their feature vector representation, so that the approriate weights can
 * be applied
 */
class CustomScore {
 public:
  typedef vector<EvScores> LabelScore;
  friend class InputData; // initializes the feature maps
  virtual ~CustomScore() { 
    for(unsigned i = 0; i < _ddVec.size(); ++i) {
      delete _ddVec[i];
    }
  }

  vector<FeatVec*> calcVec(const RegionPrediction*, dsu::Strand_t) const;
  virtual double getProb(const AnnotationItem&, EvScores::Bnd_t, dsu::Strand_t) const = 0;

  const vector<DataDescrip*>& getFeatMaps() const { return _ddVec; }
  vector<DataDescrip*>& getFeatMaps() { return _ddVec; }
  virtual LabelScore score(const RegionPrediction&, const string&) const = 0;
  virtual void score(vector<RegionPrediction*>&, const string&) const = 0;
  /**
   * an unfortunate hack at the moment, allows the client application to turn off the 
   * dynamic programming scoring, allows the application to quickly create the evidence sets
   * with out the scoring matrix (might want to move to Options)
   */
  virtual bool parseOnly() const { return false; }

 protected:

  CustomScore(unsigned size = (unsigned)EvScores::eNumBnd) : _ddVec(size) { 
    for(unsigned i = 0; i < _ddVec.size(); ++i) {
      _ddVec[i] = new DataDescrip;
    }
  } 
 protected:
  vector<DataDescrip*> _ddVec;

};

#endif //CUSTOMSCORE_H
