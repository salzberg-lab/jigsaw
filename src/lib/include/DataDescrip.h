//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef DATADESCRIP_H
#define DATADESCRIP_H

#include "FeatVec.h"
#if __GNUC__<= 2
#include <hash_map>
#else
//#include <map>
#include <ext/hash_map>
#endif
#include <iostream>


/**
 * Defines the mapping function that 
 * takes user supplied evidence sets and
 * maps them into feature vectors
 * for use in identifying and evalutating patterns 
 * by the combiner
 */

class DataDescrip {
 public: 
  DataDescrip() : _size(0), _totWght(0) { }

  /**
   * same as FeatVec::getNumDim()
   */
  int getNumDim() const;

  /**
   * Add an evidence type 
   */
  void add(EvidenceType*);

  /**
   * combined, returns the integer id that maps to hash table
   * idx that records the information about what type of 
   * evidence this is
   * @param const Prediction* = input evidence that is being
   */
  int getIdx(const Prediction*) const;
  int getIdx(const char*) const;

  /**
   * return the number of discrete value types there are
   * used right now for binary encoding
   */
  unsigned getDiscSize() const;

  /**
   * initializes a feature vector to the values
   * that reflect the predefined vector topology
   * as defined by each EvidenceType
   * e.g. homolgy = "max % sim, Sum % sim"
   * @param feature vector
   */
  void initVec(FeatVec&) const;

  FeatVec* read(std::istream&) const;
  FeatVec* read_wght(std::istream&) const;

  /**
   * query how many types of unique evidence sets there are
   */
  int size() const { return _size; }

  void calcTotWght();
  double getTotWght() const { return _totWght; }

  friend std::ostream& operator<<(std::ostream&, const DataDescrip&);
 private:
  struct myeqstr {
    bool operator()(const char* s1, const char* s2) const {
      return strcmp(s1,s2)==0;
    }
  };
#if __GNUC__<=2
  typedef hash_map<const char*, EvidenceType*, hash<const char*>, myeqstr> evHashMap;
#elif defined(__INTEL_COMPILER) && __INTEL_COMPILER <= 810
  typedef std::hash_map<const char*, EvidenceType*, std::hash<const char*>, myeqstr> evHashMap;
#else
  typedef __gnu_cxx::hash_map<const char*, EvidenceType*, __gnu_cxx::hash<const char*>, myeqstr> evHashMap;
  //typedef std::map<const char*, EvidenceType*> evHashMap;
#endif
  evHashMap _evMap;
  int _size;
  double _totWght;
};

#endif // DATADESCRIP_H
