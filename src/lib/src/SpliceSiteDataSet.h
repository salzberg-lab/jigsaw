//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef SPLICESITEDATASET_H
#define SPLICESITEDATASET_H

#include "DataSetAbstractProduct.h"

/**
 * <b>Responsibilities:</b><br>
 * 
 * Contains all the Gene Prediction data read from a Splice Site predictor data source.<BR>
 * Specialized version of DataSetAbstractProduct which specifies the basic format that all input<BR>
 * data must adhere to.<BR>
 * 
 * <b>Collaborators:</b><br>
 * SpliceSiteDataSetConcreteFactory, DataSetAbstractProduct<br>
 */
class SpliceSiteDataSet : public DataSetAbstractProduct {
public:
  SpliceSiteDataSet(const SeqLoc&, int, const string&);
  ~SpliceSiteDataSet() { }
  void print(std::ostream& os, const string&) const;
  /**
   * @param const string& --> source identifier... (e.g.) geneplicer
   * @param const string& --> confidence type (not used, can be null)
   * @param int --> 5' coordinate 
   * @param int --> 3' coordinate
   * @param const string& --> splice type (acceptor or donor)
   */
   void insert(const string&,int,int,const string&,double);

   const string& getTypeStr() const;

   // do not add introns
   void add_introns() { }

private:
  int _geneId;
  string _source;
};

#endif //SPLICESITEDATASET_H
