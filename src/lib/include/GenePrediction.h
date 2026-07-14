//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef GENEPREDICTION_H
#define GENEPREDICTION_H

#include "Prediction.h"
class GeneModelType;

/**
 * GenePrediction is derived from Prediction and is 
 * designed to record whole exon genes. (should have 
 * called it ExonPrediction oh well)
 * A GenePrediction is a prediction that is made up
 * of Initial,Internal,Single and Terminal exons.
 * 
 */

class GenePrediction : public Prediction {
public:
  /**
   * same construction as Prediction, but also takes as input an exon type identifier 
   * currently there are 4 types, initial, terminal, internal, single.
   * becuase there currently is not enough information on this, I've kept the
   * GeneModelType private under the assumption that it will not be used for now
   * I want to keep it in the framework, however, in case it can become usefull
   */
  GenePrediction(const SeqLoc& sl, AnnotationScore d, int modelId, const string& source, const GeneModelType& mt, int frame = -1) :
    Prediction(sl,d,modelId,source), _modelType(mt), _frame(frame) { }

  GenePrediction(const GenePrediction& gl, const GeneModelType& mt) : Prediction(gl), _modelType(mt), _frame(gl._frame) {}


  /**
   * no dynamic memory management
   */
  virtual ~GenePrediction() { }

  virtual void del() { }

  /**
   * standard display of this object
   */
  void print(std::ostream&, const std::string&) const;

  /**
   * type identifier for this object
   * @return the exon model type
   */
  const string& getTypeStr() const;

  /**
   * get the unique integer index associated with this
   * exon model.
   */
  const int getScoreIdx() const;

  /**
   * @return list of exons
   */
  DataStorageType* makePlaceHolders() const ;

  /**
   * @return coding frame
   */
  int getFrame() const { return _frame; }

  static const int NUM_ANNOTATION_SCORED_TYPES = 4;
  static const int NOT_GENE_ID = 4;
private:
  const GeneModelType& _modelType;
  int _frame;
};

#endif // GENEPREDICTION_H
