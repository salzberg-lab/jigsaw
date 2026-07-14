//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "HomologyDataSet.h"
#include "MatchPrediction.h"
#include "ParamFile.h"
#include "Intron.h"
#include <list>
#include <stdexcept>
#include <sstream>
#include <assert.h>

using std::cout;
using std::endl;
using std::list;

HomologyDataSet::HomologyDataSet(const SeqLoc& sl, const string& accessionId, TranscriptHash* th) 
   : DataSetAbstractProduct(sl, -1, accessionId, AnnotationScore::defaultScore()), _thPtr(th) {
  if(!_thPtr) {
    _thPtr = new TranscriptHash;
  }
}

HomologyDataSet::~HomologyDataSet() {
}

void
HomologyDataSet::del_temp() {
  delete _thPtr;
}

void HomologyDataSet::print(std::ostream& os, const string& str) const {
  Prediction::print(os,str);
  os<<endl;
  DataStorageType::const_iterator iter = begin();
  while( iter != end() ) {
    const SeqLoc& sl = **iter;
    sl.print(os,str+"\t");
    os<<endl;
    ++iter;
  }
}

static const string geneStr = "HOMOLOGY";
const string& HomologyDataSet::getTypeStr() const {
  return geneStr;
}

void HomologyDataSet::insert(MatchPrediction* mp, const string& type) {
  /** WARNING THIS SHOULD BE FILTERED ON INPUT NOT HARDCODED!!!**/
  if(mp->length() < 5) {
    delete mp;
    return;
  }
  if( mp->getPcntSim() < 0 || mp->getPcntId() < 0 ) {
    std::ostringstream ostrm;
    ostrm<<"HomologyDataSet::insert error bad % identiy/similarity value: "<<mp->getPcntSim()<<" "<<mp->getPcntId()<<endl;
    throw std::runtime_error(ostrm.str());
  }
  TranscriptHash* thPtr = _thPtr;
  const ParamFile& pf = *(ParamFile::getInstance());
  //Same as GenePredDataSet, can not throw away data unforunately right now
  // without screwing up the intron model
#if 0
  if( (options._subEnd5 != -1 && options._subEnd3 != -1) &&
      (mp->getEnd3() < options._subEnd5 || mp->getEnd5() > options._subEnd3)) {
    delete mp;
    return;
  }
#endif
  const int INTRON_CUTOFF = pf.getAlignmentConnectionCutoff(); 
  DataStorageType::iterator iter = begin();
  const DataStorageType::iterator stop = end();
  bool foundMatch=false;
#warning "right now all homology data types are assumed to provide intron info"
  if( thPtr ) {
    DataSetAbstractProduct* match = (*thPtr)[mp->getDataSource()];
    if( match ) {
      match->insert(mp);
      foundMatch = true;
    }
    
  } else {
    for(;iter != stop; ++iter) {
      Prediction* sl = *iter;
      // existing gene
      DataSetAbstractProduct* currMatch = dynamic_cast<DataSetAbstractProduct*>(sl);
      assert(currMatch);
      DataStorageType::iterator miter = currMatch->begin();
      const DataStorageType::iterator mstop = currMatch->end();
      for(; miter != mstop; ++miter ) {
	const Prediction* pred = *miter;
	const MatchPrediction* smatch = dynamic_cast<const MatchPrediction*>(pred);
	assert(smatch);
	if( mp->getDataSource() != smatch->getDataSource() ) {
	  break;
	}
	if( mp->getStrnd() == smatch->getStrnd() && (abs(mp->getEnd3()-smatch->getEnd5()) < INTRON_CUTOFF ) ) {
	  currMatch->insert( mp );
	  foundMatch = true;
	  break;
	}
      }
      if(foundMatch) {
	break;
      }
    }
  }
  if(!foundMatch) {
    DataSetAbstractProduct* ds = new HomologyDataSet(*mp,type,thPtr);
    ds->insert(mp);
    DataSetAbstractProduct::insert(ds);
    if( thPtr ) {
      (*thPtr)[mp->getDataSource()] = ds;
    }
  }
}

void
HomologyDataSet::add_introns_help() {
  //cout<<"Here I am: "<<*this<<endl;
  const DataStorageType::const_iterator stop = end();
  DataStorageType::const_iterator iter1 = begin();
  const dsu::Strand_t strnd = getStrnd();
  assert(strnd != dsu::eEither);
  list<Prediction*> tlst;
  while(iter1 != stop) {
    MatchPrediction* model1 = dynamic_cast<MatchPrediction*>(*iter1);
    DataStorageType::const_iterator iter2 = iter1;
    ++iter2;
    while(iter2 != stop) {
      MatchPrediction* model2 = dynamic_cast<MatchPrediction*>(*iter2);
      if( !isOverlap(*model1,*model2) ) {
	assert(model2->getScoreIdx() != Intron::getInstance().getScoreIdx());
	const int rghtSide = model1->getEnd3()+1;
	const int leftSide = model2->getEnd5()-1;
	SeqLoc nl(strnd,rghtSide,leftSide);
	const double percid = (model1->getPcntId()+model2->getPcntId())/2.0;
	const string accession = model1->getDataSource();
	Prediction* thePrediction = new MatchPrediction(nl,percid,-1,accession,percid,percid,percid,SeqLoc(dsu::eEither,-1,-1),Intron::getInstance());
	tlst.push_back(thePrediction);
	iter2 = stop;
      } else {
	iter1 = iter2;
	break;
      }
    }
    if(iter1 != iter2) {
      ++iter1;
    }
  }
  list<Prediction*>::iterator iter = tlst.begin();
  const list<Prediction*>::iterator tstop = tlst.end();
  for(; iter != tstop; ++iter) {
    DataSetAbstractProduct::sort_insert(*iter);
  }
}
