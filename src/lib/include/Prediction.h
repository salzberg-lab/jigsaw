//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef PREDICTION_H
#define PREDICTION_H

#include "SeqLoc.h"
#include "AnnotationScore.h"
#include <dsu.h>
#include <list>

class RegionPrediction;

/**
 * <b>Responsibilities:</b><br>
 * sequence location.  GeneModels are a collection of <i>GenePrediction</i> objects. <br>
 * Sequence Match models are a collection of <i>MatchPrediction<i> objects, both of which<br>
 * are just <i>DataSetAbstractProduct</i>, whose behavior are outlined by this object.<br>
 * In other words, one knowledge source might be made up of a collection of predictions about related
 * sequence locations or it might be just one sequence location
 * Presumeably, there will be many types of predictions:
 * CDS, start or stop codons, 5' UTRs, 3' UTRS, introns, exons for starters.
 * <br>
 * <b>Collaborators:</b><br>
 * SeqLoc, AnnotationScore<br>
 * @stereotype interface 
 */
class Prediction : public SeqLoc {
public:
  typedef std::list<Prediction*> DataStorageType;

  /**
   * @param SeqLoc - sequence location this prediction refers to
   * @param AnnotationScore - the evaluation score of this prediction
   * @param int - a counter of ith element in the prediction model (only usefull
   * for models made up of multiple components this should be phased out) 
   * @param string - string identifier, defining the "source" of this prediction
   * this value is gnerally computed automatically by the InputData class, so that the scoring
   * functions can uniquely seperate different prediction sources...(e.g. nap versus gap)
   */
  Prediction(const SeqLoc& sl, AnnotationScore d, int modelId, const string& str) : 
    SeqLoc(sl), _modelId(modelId), _source(str), _score(d) { }

  /**
   * copy constructor
   */
  Prediction(const Prediction& pred) 
    : SeqLoc(pred), _modelId(pred._modelId), _source(pred._source), _score(pred._score) { }

  /**
   * no dynamic memory management, does contain one 'stl string' object.
   */
  virtual ~Prediction() { }

  virtual void del() = 0;
  /**
   * return score associated with this object
   */
  AnnotationScore getScore() const;

  /**
   * set score valuie
   */
  void setScore(const AnnotationScore& s) { _score = s; }

  /**
   * get integer idx of this coding sequence prediction
   */
  int getModelId() const { return _modelId; }

  /**
   * get string identifier of the source of this prediction
   */
  const string& getDataSource() const { return _source; }

  /**
   * all prediction type objects must be able convert themselfs
   * into a list of pointers that contain each coding sequence seperately
   */
  virtual DataStorageType* makePlaceHolders() const  = 0;

  /**
   * print to ostream, information about this prediction
   */
  virtual void print(std::ostream& os, const string& str) const; 

private:
  int _modelId;
  string _source;
  AnnotationScore _score;
};

#endif // PREDICTION_H
