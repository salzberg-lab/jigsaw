//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef HOMOLOGYDATASET_H
#define HOMOLOGYDATASET_H

#include "DataSetAbstractProduct.h"
#include <map>
class MatchPrediction;

/**
 * <b>Responsibilities:</b><br>
 * 
 * Contains all the Gene Prediction data from a sequence similarity data source.<BR>
 * Specialized version of DataSetAbstractProduct which specifies the basic format that all input<BR>
 * data must adhere to.<BR>
 * I originally was goinging to have a separate representation for distinc types of homology data
 * like "HOMOLOGY", but this isn't needed right now so the Homology namesake got used accross the board
 * plan to change this name to HomologyDataSet
 * 
 * <b>Collaborators:</b><br>
 * HomologyDataSetConcreteFactory, DataSetAbstractProduct<br>
 */
class HomologyDataSet : public DataSetAbstractProduct {
 public:
  //struct my_hash {
    //size_t operator()(const std::string& __s) const {
      //return __gnu_cxx::__stl_hash_string(__s.c_str()); }
  //};
  typedef std::map<std::string, DataSetAbstractProduct*> TranscriptHash;

public:
  HomologyDataSet(const SeqLoc&, const string&, TranscriptHash* = NULL);
  ~HomologyDataSet();

  void del_temp();
  void print(std::ostream& os, const string&) const;

  void insert(MatchPrediction*,const string&);
  const string& getTypeStr() const;

private:
  void add_introns_help();

 private:
  TranscriptHash* _thPtr;

};

#endif //HOMOLOGYDATASET_H
