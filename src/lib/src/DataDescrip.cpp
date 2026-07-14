//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "DataDescrip.h"
#include "Prediction.h"
#include <fstream>
#include <string.h>
#include <assert.h>

std::ostream& operator<<(std::ostream& os, const DataDescrip& dd) {
  DataDescrip::evHashMap::const_iterator iter = dd._evMap.begin();
  const DataDescrip::evHashMap::const_iterator stop = dd._evMap.end();
  while( iter != stop ) {
    const EvidenceType* ev = (*iter).second;
    assert(ev); 
    ++iter;
  }
  return os;
}

void DataDescrip::add(EvidenceType* ev) {
  assert(ev);
  const char* str_id = strdup(ev->getStr());
  _evMap[str_id] = ev;
  ++_size;
}

void DataDescrip::initVec(FeatVec& vec) const {
  evHashMap::const_iterator iter = _evMap.begin();
  evHashMap::const_iterator stop = _evMap.end();
  while( iter != stop ) {
    const EvidenceType* ev = (*iter).second;
    assert(ev);
    vec[ev->getIdx()] = ev->dup();
    ++iter;
  }
}

void DataDescrip::calcTotWght() {
  if( _totWght == 0 ) {
    evHashMap::const_iterator iter = _evMap.begin();
    evHashMap::const_iterator stop = _evMap.end();
    double sum = 0;
    while( iter != stop ) {
      const EvidenceType* ev = (*iter).second;
      assert(ev);
      sum += ev->getWght();
      ++iter;
    }
    _totWght = sum;
  }
}

int DataDescrip::getNumDim() const {
  evHashMap::const_iterator iter = _evMap.begin();
  evHashMap::const_iterator stop = _evMap.end();
  int cnt = 0;
  while( iter != stop ) {
    const EvidenceType* ev = (*iter).second;
    assert(ev);
    cnt += ev->getNum();
    ++iter;
  }
  return cnt;
}

FeatVec* DataDescrip::read(std::istream& ifs) const {
  FeatVec* retVec = new FeatVec(size());
  initVec(*retVec);
  FeatVec& vec = *retVec;
  evHashMap::const_iterator iter = _evMap.begin();
  evHashMap::const_iterator stop = _evMap.end();
  int cnt = 0;
  while( iter != stop ) {
    EvidenceType& et = *(vec[cnt]);
    assert(et.getIdx()==cnt);
    if( !et.read(ifs) ) {
      delete retVec;
      return NULL;
    }
    ++iter;
    ++cnt;
  }
  return retVec;
}

int DataDescrip::getIdx(const Prediction* pred) const {
  const char* id = pred->getDataSource().c_str();
  return getIdx(id);
}

int DataDescrip::getIdx(const char* id) const {
  evHashMap::const_iterator iter = _evMap.find(id);
  if(iter == _evMap.end()) {
    return -1;
  }
  const EvidenceType* ev = (*iter).second; 
  return ev->getIdx();
}

