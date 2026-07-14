//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef INPUTDATA_H
#define INPUTDATA_H

#include <string>
#include "Prediction.h"
#include "DataDescrip.h"
#include <vector>
#include <list>

using std::vector;
using std::list;

class CustomScore;
class DataSetAbstractProduct;

/**
 * <b>Responsibilities:</b><br>
 * Stores all data input information, in sorted 5' -> 3' order.<br>
 * <br>
 * <b>Collaborators:</b><br>
 * DataSetAbstractFactory, GenePredDataSetFactory, NapDataSetConcreteFactory, SpliceSiteDataSetConcreteFactory,<br>
 * Prediction, Combine<br> 
 */
class InputData {
public:
  typedef vector<Prediction*>::const_iterator const_iterator;

  /**
   * <b>Precond:</b><br>
   * Takes as input a file that specifies a list of flat files containing data.<br>
   * The format of the input file is:
   * flatfile name (fullpath), data type, weights
   * there is a factory for each data type to read and store the data.  A reference to the entire <br>
   * <b>Postcond:</b><br>
   * All data stored in sorted 5'->3' order.  Can now be used by the <br>
   * the Combine object to go through and parse and score<br>
   * sub regions of the genmoic sequence<br> 
   * CustomScore object stores the approriate feature vector descriptions, as determined by the users input file
   * see CustomScore for more details
   *
   * InputData calls factories based on the following keyword string ids:
   * homology -- NapDataSetConcreteFactory
   * geneprediction -- GenePredDataSetConcreteFactory
   * spliceprediction -- SpliceSiteDataSetConcreteFactory
   * curation -- CurationConcreteFactory
   */
  InputData(); 
  ~InputData();

  const_iterator begin() const { return _data.begin(); }
  const_iterator end() const { return _data.end(); }
  int getEnd5() const { assert(_data[0]); return _data[0]->getEnd5(); }
  const string& getTypeStr() const;
  friend std::ostream& operator<<(std::ostream&, const InputData&);    
  
  void setTrainExamples(DataSetAbstractProduct* dp) { _trainExamples = dp; }
  DataSetAbstractProduct* getTrainExamples() { return _trainExamples; }

  void readEvidenceFile(const string&, const string&, vector<DataDescrip*>&);
  void readEvidenceFileLinComb(const string&, const string&, vector<DataDescrip*>&);

 private:
  vector<Prediction*> _data;
  list<DataSetAbstractProduct*> _stopLeak;
  DataSetAbstractProduct* _trainExamples;
};

#endif // INPUTDATA_H
