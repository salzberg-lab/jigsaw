//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef SPLICESITE_H
#define SPLICESITE_H

#include "Prediction.h"

class RegionPrediction;

/**
 * <b>Responsibilities:</b><br>
 * Represents a splice site which is a specialization of <i>Prediction</i> object.<br>
 * <br>
 * <b>Collabors:</b><br>
 * Prediction 
 */
class SpliceSite : public Prediction {
public:

  /**
   * calls basic Prediction constructor and also
   * stores the splice type, either donor, or acceptor
   */
  SpliceSite(const SeqLoc& sl, AnnotationScore sc,  int modelId, 
	     const string& source, const string& spliceType) :
		Prediction(sl,sc,modelId,source), _spliceType(spliceType) { }
  ~SpliceSite() { }
  virtual void del() { }

  /**
   * print out info on splice prediction
   */
  void print(std::ostream& os, const string&) const;

  /**
   * returns constant string identifier of this object
   * which is called "SpliceSite"
   */
  const string& getTypeStr() const;

  /**
   * return string identifier for splice site
   * type, donor or acceptor.
   */
  const string& getSpliceType() const;

  /**
   * true if donor splice site false otherwise
   */
  bool isDonor() const;

  /**
   * true if acceptor splice site false otherwise
   */
  bool isAcceptor() const;
  
  /**
   * each prediction must implement this method
   */
  DataStorageType* makePlaceHolders() const ;
private:
  string _spliceType;


};

#endif // SPLICESITE_H
