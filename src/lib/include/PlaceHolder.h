//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef PLACEHOLDER_H
#define PLACEHOLDER_H

#include "Prediction.h"

class ostream;

/**
 * <b>Responsibilities:</b><br>
 * represents an individual exon and points to the exons parent genemodel<br>
 * The parent gene models is what gets stored in the big scoring matrix.<br>
 * This class is used as a matter of conveniance for<br>
 * breaking up the genomic sequence into regions made up of unique datasets<br>
 * of annotation data<br>
 * <br>
 * <b>Collaborators:</b><br>
 * Prediction 
 */
class PlaceHolder : public Prediction {
public:
  PlaceHolder(const Prediction& pred, const Prediction* big) : Prediction(pred), _totalSequence(big) { }

  virtual ~PlaceHolder() { }

  virtual void del() { }
  DataStorageType* makePlaceHolders() const { assert(0); return NULL; }

  /**
   * returns a constant string identifier that identifies this
   * as a "PlaceHolder" object
   */
  const string& getTypeStr() const { return _ph; }

  /**
   * return a pointer to the entire gene model data
   */
  const Prediction* getWholeShabang() const { return _totalSequence; }

  void print(std::ostream&,const std::string& ) const;

private:
  const static string _ph;
  const Prediction* _totalSequence;
};


#endif
