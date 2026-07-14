//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "DataSetAbstractProduct.h"
#include "PlaceHolder.h"
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <list>

using std::list;
using std::endl;

DataSetAbstractProduct::DataSetAbstractProduct(const SeqLoc& sl, int modelId, const string& source,
					       AnnotationScore score) :
  Prediction(sl,score,modelId,source) {
}

DataSetAbstractProduct::~DataSetAbstractProduct() {
  lnkSeqLoc.clear();
}

void 
DataSetAbstractProduct::del() {
  DataStorageType::iterator iter = begin();
  DataStorageType::iterator stop = end();
  while(iter != stop) {
    Prediction* pred = *iter;
    delete pred;
    ++iter;
  }
}

DataSetAbstractProduct::iterator 
DataSetAbstractProduct::end() {
  return lnkSeqLoc.end();
}

DataSetAbstractProduct::const_iterator 
DataSetAbstractProduct::end() const {
  return lnkSeqLoc.end();
}

DataSetAbstractProduct::iterator 
DataSetAbstractProduct::begin() {
  return lnkSeqLoc.begin();
}

DataSetAbstractProduct::const_iterator 
DataSetAbstractProduct::begin() const {
  return lnkSeqLoc.begin();
}

unsigned 
DataSetAbstractProduct::size() const {
  return lnkSeqLoc.size();
}

const Prediction* 
DataSetAbstractProduct::getLast() const {
  return lnkSeqLoc.back();
}

// important precondition is that end5 end3 refers to only one possible prediction
void
DataSetAbstractProduct::getRelatedPrediction(int end5, int end3, list<const Prediction*>& lst) const {
  DataStorageType::const_iterator iter = begin(), stop = end();
  while( iter != stop ) {
    const SeqLoc& sl = **iter;
    if( sl.inRange(end5) || sl.inRange(end3) ) {
      //return *iter;
      lst.push_back(*iter);
    }
    ++iter;
  }
}

void
DataSetAbstractProduct::getRelatedPrediction5(int end5, list<const Prediction*>& lst) const {
  DataStorageType::const_iterator iter = begin(), stop = end();
  while( iter != stop ) {
    const SeqLoc& sl = **iter;
    if( sl.getEnd5() == end5) {
      lst.push_back(*iter);
    }
    ++iter;
  }
}

void
DataSetAbstractProduct::getRelatedPrediction3(int end3, list<const Prediction*>& lst) const {
  DataStorageType::const_iterator iter = begin(), stop = end();
  while( iter != stop ) {
    const SeqLoc& sl = **iter;
    if( sl.getEnd3() == end3) {
      lst.push_back(*iter);
    }
    ++iter;
  }
}

// RIGHT NOW THIS IS SET UP TO FIND INTRON PREDICTIONS!!!!! NOTE THAT THE 5',3' ORDER
// IS REVERSED FOR EXON PREDICTIONS
pair<const Prediction*,const Prediction*>
DataSetAbstractProduct::getRelatedPredictionInterval(int lidx, int ridx) const {
  const DataStorageType::const_iterator stop = end();
  DataStorageType::const_iterator iter1 = begin();
  const Prediction* left = NULL;
  const Prediction* rght = NULL;
  for(; iter1 != stop; ++iter1 ) {
    const Prediction* sl = *iter1;
    if( sl->getEnd3()==lidx ) {
      left=sl;
      DataStorageType::const_iterator iter2 = iter1;
      ++iter2;
      if( iter2 != stop ) {
	const Prediction* sl2 = *iter2;
	assert(sl2);
	if( sl2->getEnd5()==ridx ) {
	  rght=sl;
	  break;
	}
      }
    }
  }
  return pair<const Prediction*,const Prediction*>(left,rght);
}

// 
// looks to see if end5,end3 falls inbtween a gene prediction to indicate that it is
// with in an intron
pair<const Prediction*,const Prediction*>
DataSetAbstractProduct::getIntronBoundaries(int end5, int end3) const {
  const DataStorageType::const_iterator stop = end();
  DataStorageType::const_iterator iter1 = begin();
  const Prediction* left = NULL;
  const Prediction* rght = NULL;
  for(; iter1 != stop; ++iter1 ) {
    const Prediction* lPred = *iter1;
    DataStorageType::const_iterator iter2 = iter1;
    ++iter2;
    if( iter2 != stop ) {
      const Prediction* rPred = *iter2;
      if( end5 > lPred->getEnd3() && end5 < rPred->getEnd5() &&
	  end3 > lPred->getEnd3() && end3 < rPred->getEnd5() ) {
	left = lPred;
	rght = rPred;
	break;
      }
    }
  }
  return pair<const Prediction*,const Prediction*>(left,rght);
}

void
DataSetAbstractProduct::add_introns_help() {
  std::cerr<<"Warning this should only be called by derived data set forms"<<endl;
}

void
DataSetAbstractProduct::add_introns() {
  const DataStorageType::const_iterator stop = end();
  DataStorageType::const_iterator iter1 = begin();
  for(; iter1 != stop; ++iter1 ) {
    DataSetAbstractProduct* model = dynamic_cast<DataSetAbstractProduct*>(*iter1);
    if(model) {
      model->add_introns_help();
    }
  }
}

struct Compare : public std::binary_function<const SeqLoc*, const SeqLoc*,bool> {
  bool operator()(const SeqLoc* a, const SeqLoc* b) const {
    return *a < *b;
  }
};

void
DataSetAbstractProduct::sort() {
  //std::vector<Prediction*> data(size());
  const DataStorageType::const_iterator stop = end();
  DataStorageType::const_iterator iter1 = begin();
  for(unsigned cnt = 0; iter1 != stop; ++iter1, ++cnt ) {
    Prediction* pred = *iter1;
    assert(pred);
    DataSetAbstractProduct* model = dynamic_cast<DataSetAbstractProduct*>(pred);
    if( model ) {
      model->lnkSeqLoc.sort(Compare());
    }
  }
  lnkSeqLoc.sort(Compare());
#if 0
  sort(data.begin(),data.end(),Compare());
  reset();
  for(unsigned cnt = 0; cnt < data.size(); ++cnt ) {
    push_back(data[cnt]);
  }
#endif
}

// important precondition is that end5 end3 refers to only one possible prediction
const Prediction* 
DataSetAbstractProduct::getRelatedPrediction(int endp) const {
  DataStorageType::const_iterator iter = begin(), stop = end();
  while( iter != stop ) {
    const SeqLoc& sl = **iter;
    if( sl.inRange(endp) ) return *iter;
    ++iter;
  }
  return NULL;
}

void 
DataSetAbstractProduct::insert(Prediction* seqLoc) {
  if( getEnd5() > seqLoc->getEnd5() || !isEnd5Set())
    setEnd5(seqLoc->getEnd5());
  if( getEnd3() < seqLoc->getEnd3() || !isEnd3Set())
    setEnd3(seqLoc->getEnd3());
  
  lnkSeqLoc.push_back(seqLoc);
}

void 
DataSetAbstractProduct::sort_insert(Prediction* seqLoc) {
  if( getEnd5() > seqLoc->getEnd5() || !isEnd5Set())
    setEnd5(seqLoc->getEnd5());
  if( getEnd3() < seqLoc->getEnd3() || !isEnd3Set())
    setEnd3(seqLoc->getEnd3());
  
  DataStorageType::iterator iter = begin(), stop = end();
  while( iter != stop ) {
    const SeqLoc* curr = *iter;
    if( seqLoc->getEnd5() < curr->getEnd5() ||
	(seqLoc->getEnd5() == curr->getEnd5() &&
	 seqLoc->getEnd3() < curr->getEnd3() )) {
      break;
    }
    ++iter;
  }
  lnkSeqLoc.insert(iter,seqLoc);
}

Prediction::DataStorageType* 
DataSetAbstractProduct::makePlaceHolders() const {
  DataStorageType* retVal = new DataStorageType;
  DataStorageType::const_iterator iter = begin(), stop = end();
  while( iter != stop ) {
    Prediction* ph = new PlaceHolder(**iter,this); 
    retVal->push_back( ph );
    ++iter;
  }
  return retVal;
}

int DataSetAbstractProduct::getClosestEnd5(int end5val) const {
  DataStorageType::const_iterator iter = begin(), stop = end();
  int best = -1, closestEnd5 = -1;
  while( iter != stop ) {
    const SeqLoc& sl = **iter;
    int curr_val = sl.getEnd5() - end5val;
    if(sl.getEnd5() >= end5val && (curr_val < best || best == -1)) {
      closestEnd5 = sl.getEnd5();
      best = curr_val;
    }
    ++iter;
  }
  return closestEnd5;
}

const Prediction* DataSetAbstractProduct::getClosestEnd5Exon(int end5val) const {
  DataStorageType::const_iterator iter = begin(), stop = end();
  int best = -1; 
  const Prediction* bestPred = NULL;
  while( iter != stop ) {
    const Prediction* pred = *iter;
    const SeqLoc& sl = *pred;
    int curr_val = end5val - sl.getEnd5();
    if(sl.getEnd5() <= end5val && (curr_val < best || best == -1)) {
      best = curr_val;
      bestPred = pred;
    }
    ++iter;
  }
  if(bestPred) {
    const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(bestPred);	
    if(dp) {
      return dp->getClosestEnd5Exon(end5val);
    }
  }
  return bestPred;
}

int DataSetAbstractProduct::getListPos(const Prediction* pred) const {
  DataStorageType::const_iterator iter = begin(), stop = end();
  int cnt = 0;
  while( iter != stop ) {
    const Prediction* p = *iter;
    if( p == pred ) return cnt;
    ++cnt;
    ++iter;
  }
  return -1;
}


const Prediction* DataSetAbstractProduct::getClosestEnd3Exon(int end3val, int aMinEnd5) const {
  DataStorageType::const_iterator iter = begin(), stop = end();
  int best = -1; 
  const Prediction* bestPred = NULL;
  while( iter != stop ) {
    const Prediction* pred = *iter;
    const SeqLoc& sl = *pred;
    if( sl.getEnd3() < aMinEnd5 ) { ++iter; continue; }
    int curr_val = abs(sl.getEnd3() - end3val);
    if( (curr_val < best || best == -1 ) ) {
      best = curr_val;
      bestPred = pred;
    }
    ++iter;
  }
  if(bestPred) {
    const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(bestPred);	
    if(dp) {
      return dp->getClosestEnd3Exon(end3val,aMinEnd5);
    }
  }
  return bestPred;
}

void DataSetAbstractProduct::print(std::ostream& os, const string& str) const {
  SeqLoc::print(os,str);
  os<<endl;
  DataStorageType::const_iterator iter = begin(), stop = end();
  while( iter != stop ) {
    const SeqLoc& sl = **iter;
    sl.print(os,str+"\t");
    os<<endl;
    ++iter;
  }
}

std::ostream& operator<<(std::ostream& os, const DataSetAbstractProduct& dsap) {
  dsap.print(os,"");
  return os;
}
