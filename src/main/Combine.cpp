//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "Combine.h"
#include "genemodels.h"
#include "RegionPrediction.h"
#include "dsu.h"
#include "CustomScore.h"
#include "Options.h"
#include <stdexcept>

using std::exception;
using std::cerr;
using std::endl;
using std::cout;

Combine::Combine(const InputData & id, const string& theSequence, bool forceWholeGeneModel, const CustomScore& cstmScr, int cutOff) : 
  _inData(id), _theSequence(theSequence), _annotation(NULL), _cstmScr(&cstmScr) {
  _annotation = new Annotation(theSequence,forceWholeGeneModel,cutOff);
}

void
Combine::run()  {
  try {
    InputData::const_iterator iter = _inData.begin();
    //int end5 = _inData.getEnd5();
    const Options& opt = *(Options::getInstance());
    const int end5 = opt._subEnd5 == -1 ? _inData.getEnd5() : opt._subEnd5;
    getAllRelaventData(end5, iter, true);
  } catch(const exception& e) {
    cerr<<"Combine::Combine exception: "<<e.what()<<endl;
    exit(1);
  }
}

Combine::~Combine() {
  delete _annotation;
  list<RegionPrediction*>::iterator iter = _lst.begin(); 
  const list<RegionPrediction*>::iterator stop = _lst.end(); 
  for(; iter != stop; ++iter) {
    RegionPrediction* delme = *iter;
    delete delme;
  }
  _lst.clear(); //automatic
}

bool Combine::empty() const {
  return _annotation ?  _annotation->empty() : false;
}

#if _RUNTIME_DEBUG >= 5
bool COMBINE_DEBUG_ON = true;
#endif

static int getNext(const InputData::const_iterator& data_iter, int bp_idx, const InputData::const_iterator& stop) {
  InputData::const_iterator iter = data_iter;
  int curr_best = -1, closest = -1;
  while( iter != stop ) {
    Prediction * currPred = *iter;
    int baseline = currPred->getClosestEnd5(bp_idx);
    if(baseline < 0) {
      ++iter;
      continue;
    }
    int diff = baseline-bp_idx;
    if( diff < curr_best && diff > 0 || curr_best == -1 ) {
      curr_best = baseline-bp_idx;
      closest = baseline;
    } else if(baseline-bp_idx >= curr_best)
      break;
    ++iter;
  }
  return closest;
}

bool TDEBUG;

void Combine::getAllRelaventData(int bp_idx, InputData::const_iterator& data_iter, bool) {
  const Options& opt = *(Options::getInstance());
  int subEnd3 = opt._subEnd3;
  int subEnd5 = opt._subEnd5;
  if(opt._subEnd5 == -1 || opt._subEnd3 == -1) {
    subEnd3 = _theSequence.size();
    subEnd5 = 0;
  }
  while(bp_idx != -1 || data_iter != _inData.end()) {
    SeqLoc sl(dsu::eEither,bp_idx,bp_idx);
    RegionPrediction* relavent_data = new RegionPrediction(sl,_theSequence);
    //cout<<"hell: "<<relavent_data->getEnd3()<<endl;
    // for each seperate data source (e.g. glimmer, dps, dds etc.)
    int nextBestPossible = -1;
    const InputData::const_iterator& stop = _inData.end();
    if( bp_idx == 41506841 ) {
      TDEBUG=true;
    } else {
      TDEBUG=false;
    }
    relavent_data->storeRelaventData(data_iter,stop);
    // if nothing was added from this data source, find a possible bp idx to start
    if(relavent_data->empty()) {
      delete relavent_data;
      if( data_iter != _inData.end() ) {
	nextBestPossible = getNext(data_iter,bp_idx,stop);
      }
    } else {
#if _RUNTIME_DEBUG >= 5
      if(relavent_data->getEnd5() == 41506841 && relavent_data->getEnd3() == 41508442) {
	cout<<"insertion: "<<*relavent_data<<endl;
	//exit(1);
      }
#endif
      nextBestPossible = relavent_data->getEnd3() + 1;
      if( relavent_data->noData() ) {
	delete relavent_data;
      } else if( !relavent_data->noData() ) {
	if( relavent_data->getEnd3() > subEnd3 && subEnd3 != -1 ) {
	  relavent_data->setEnd3(subEnd3);
	}
	if(!_cstmScr->parseOnly() && ((relavent_data->getEnd5() >= subEnd5 && relavent_data->getEnd5() <= subEnd3) ||
					       (relavent_data->getEnd3() >= subEnd5 && relavent_data->getEnd3() <= subEnd3))) {
	  //relavent_data->scoreHelp(_theSequence,dsu::ePos);
	  //relavent_data->scoreHelp(_theSequence,dsu::eNeg);
	  //cout<<"hmm: "<<relavent_data->getEnd5()<<" "<<relavent_data->getEnd3()<<" "<<relavent_data->size()<<endl;
	  _lst.push_back(relavent_data);
	} else if(((relavent_data->getEnd5() >= subEnd5 && relavent_data->getEnd5() <= subEnd3) ||
		   (relavent_data->getEnd3() >= subEnd5 && relavent_data->getEnd3() <= subEnd3))) {
	  relavent_data->scoreHelp(_theSequence,dsu::ePos);
	  relavent_data->scoreHelp(_theSequence,dsu::eNeg);
	  vector<RegionPrediction*> vec(3);
	  vec[1]=relavent_data;
	  _cstmScr->score(vec,_theSequence);
	  _lst.push_back(relavent_data);
	  //cout<<"hmm: "<<*relavent_data<<endl;
	} else {
	  delete relavent_data;
	}
      }
    }
    bp_idx = nextBestPossible;
    if( bp_idx > subEnd3) {
      bp_idx = -1;
      data_iter = _inData.end();
    }
  }
  const bool doScore = !_cstmScr->parseOnly();
  if( doScore ) {
    score(_lst);
  }
#if 0
  const Annotation& anot = getAnnotation();
  anot.print(cout,*_cstmScr);
  //cout<<"HMMM: "<<getAnnotationMatrix()<<endl;
  exit(1);
#endif
}

void Combine::score(list<RegionPrediction*>& rlst) {
  list<RegionPrediction*>::iterator iter = rlst.begin(), last = rlst.end(), next = rlst.end();
  list<RegionPrediction*>::iterator stop = rlst.end();
  while(iter != stop) {
    vector<RegionPrediction*> vec(3);
    if( last == stop ) {
      vec[0] = NULL;
      if( !_cstmScr->parseOnly() ) {
	(*iter)->scoreHelp(_theSequence,dsu::ePos);
	(*iter)->scoreHelp(_theSequence,dsu::eNeg);
      }
    } else {
      vec[0] = *last;
    }
    vec[1] = *iter;
    next = iter;
    ++next;
    if( next == stop ) {
      vec[2] = NULL;
    } else {
      if( !_cstmScr->parseOnly() ) {
	(*next)->scoreHelp(_theSequence,dsu::ePos);
	(*next)->scoreHelp(_theSequence,dsu::eNeg);
      }
      vec[2] = *next;
    }
    _cstmScr->score(vec,_theSequence);
    if( !_cstmScr->parseOnly() ) {
      _annotation->insert(vec,*_cstmScr);
    }
    last = iter;
    ++iter;
  }
}
