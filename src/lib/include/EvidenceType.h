//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef EVIDENCETYPE_H
#define EVIDENCETYPE_H

#include "AnnotationScore.h"
#include "dsu.h"
#include <string>

using std::string;

class Prediction;
/** 
 * Interface specifies the requiremnts
 * of how user defined evidence types must behave
 * as a feature in a vector.
 * Note that the instanciated derived types of
 * this class, must correspond to an element in the Prediction 
 * class heirarchy.  In other words types of Prediction
 * objects map to types of evaluation feature objects
 * (there is not an explicit  relationship established in
 * the code)
 */
class EvidenceType {
 public:
  EvidenceType(const string& str, int idx, double val, double wght) : _val(val), _wght(wght),
    _str(str), _idx(idx) { }
  virtual ~EvidenceType() { }

  /**
   * takes a Prediction about a given region of the data 
   * where that Prediction correspons to a particular
   * evidence type, and updates it's values in the feature vector
   * @param end5 = 5' end of region looking at
   * @param end3 = 3' end of region looking at
   */
  virtual void update(const Prediction* pred, int end5, int end3, dsu::Strand_t) =  0; 
  virtual void finish() { }
  
  //virtual void 
  //update(const Prediction*, int, const Prediction*, int, dsu::Strand_t) { }

  virtual unsigned getNum() const = 0;

  /**
   * returns the i th feature value
   * some features might be encoded with dependent values
   * such as homology data that stores best match and total sum of matches
   */
  virtual double getVal(int) const { return _val; }

  /**
   * create a duplicate of this evidence type
   * used to create heterogeneous feature vectors
   * without having to actually know what type of evidence
   * the feature vector is made up of
   */
  virtual EvidenceType* dup() const = 0;
  
  /**
   * set vals
   * @param - const EvidenceType& specifies values
   */
  virtual void setVals(const EvidenceType& et) {
    _val = et._val;
  }

  virtual double eval() const { return _val*_wght; }
  /**
   * during the training mode, method defines how this element
   * of the feature vector should look (in print)
   */
  virtual bool train_print(std::ostream& os) const = 0;

  /**
   * get the string that identifies this type of evidence
   */ 
  const char* getStr() const { return _str.c_str(); }

  /**
   * get a unique integer index that identifies this type of evidence
   */
  const int getIdx() const { return _idx; }
  virtual bool read(std::istream&) = 0;

  double& getValRef() { return _val; }
  virtual double getVal() const { return _val; }
  double getWght() const { return _wght; }

  virtual double getSumVal() const { return getVal(); }
  virtual void normalize(double sum) { getValRef() /= sum; }
  
 protected:
  double _val, _wght;

 private:
  string _str;
  int _idx;
};

#endif
