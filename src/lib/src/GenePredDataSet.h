//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef GENEPREDDATASET_H
#define GENEPREDDATASET_H

#include "DataSetAbstractProduct.h"
#include <map>

/**
 * <b>Responsibilities:</b><br>
 * 
 * Contains all the Gene Prediction data read from a single source. e.g. GLIMMER, GENSCAN etc.<BR>
 * Specialized version of DataSetAbstractProduct which specifies the basic format that all input<BR>
 * data must adhere to.<BR>
 * 
 * <b>Collaborators:</b><br>
 * GenePredDataSetConcreteFactory, DataSetAbstractProduct<br> 
 */
class GenePredDataSet : public DataSetAbstractProduct {
public:
  friend class Gff3Factory;
  typedef std::map<int, GenePredDataSet*> TranscriptHash;
  /**
   ** @param const SeqLoc& --> sequence location of the prediction
   ** @param int --> the gene model id
   ** @param const string& --> string identifying the source (for readability)
   **/
  GenePredDataSet(const SeqLoc&, int, const string&, TranscriptHash* = NULL);

  virtual ~GenePredDataSet();


  /**
   * derived from DataSetAbstractProduct
   * removes temprary hash tables
   */
  void del_temp();

  /** 
   ** standard printing
   **/
  void print(std::ostream& os, const string&) const;

  /**
   ** each SeqLoc derived object can reimplement this 
   ** identify itself when print information to the screen
   ** for readability only, this actually can be phased
   ** out with the use of virtual print method
   **/
  const string& getTypeStr() const;

  /**
   ** One predicted gene model is built up, one exon at a time
   ** by repeatedly calling insert, providing a user defined numeric id of a model
   ** gene model to associate the exon with.
   ** @param const string& --> specifies the identifying source of the prediction
   ** like..glimmer, genscan etc.
   ** @param int --> numeric id associating this exon with a specific gene model
   ** @param const SeqLoc& --> the sequence location of the exon
   ** @param const string& --> a string identifier of the exon type: <br>
   ** internal, terminal, single, initial (this is perhaps a tiresome requirement
   ** that can be relaxed in the future
   ** @param double --> a confidence socre associated with this exon
   **/
  void insert(const string&, int, const SeqLoc&, const string&, double);

 private:
  void add_introns_help();
  void label_exons(const std::string&);

 private:
  int getModelId() const;


  TranscriptHash* _thPtr;

};


#endif //GENEPREDDATASET_H
