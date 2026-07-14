//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "GenePredDataSet.h"
#include "GenePrediction.h"
#include "Options.h"
#include "Terminal.h"
#include "Initial.h"
#include "Single.h"
#include "Internal.h"
#include "Intron.h"
#include <sstream>
#include <stdexcept>
#include <assert.h>
#include <list>

using std::cout;
using std::endl;
using std::runtime_error;
using std::list;

GenePredDataSet::GenePredDataSet(const SeqLoc& sl, int geneId, const string& str, TranscriptHash* th) : 
  DataSetAbstractProduct(sl, geneId, str,AnnotationScore::defaultScore()), _thPtr(th) {
  if(!_thPtr) {
    _thPtr = new TranscriptHash;
  }
}

GenePredDataSet::~GenePredDataSet() {
  //DataSetAbstractProduct::~DataSetAbstractProduct();
}


void
GenePredDataSet::del_temp() {
  delete _thPtr;
}

void GenePredDataSet::print(std::ostream& os, const string& str) const {
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

void GenePredDataSet::insert(const string& source, int geneId, const SeqLoc& theSl, const string& type, double score) {
  if( theSl.getEnd5() < 1 || theSl.getEnd3() < 1 ) {
    std::ostringstream ostrm;
    ostrm<<"Coordinates are expected to be 1 or greater: "<<theSl;
    throw std::runtime_error(ostrm.str());
  }
  const Options& options = *(Options::getInstance());
  if(options.isLinComb())
    score = 1.0; // prediction scores ignored.
  GenePredDataSet* currLocGene = (*_thPtr)[geneId];
  bool newGene = false;
  if(!currLocGene) {
    currLocGene = new GenePredDataSet(theSl,geneId,source,_thPtr);
    (*_thPtr)[geneId] = currLocGene;
    newGene = true;
  }

  Prediction *thePrediction = NULL;
  const GeneModelType* exon_type = NULL;
  if(type == "Internal" || type == "internal" || type == "internal-exon") {
    exon_type = &(Internal::getInstance());
  } else if(type == "Terminal" || type == "terminal" || type == "terminal-exon")
    exon_type = &(Terminal::getInstance());
  else if(type == "Initial" || type == "initial" || type == "initial-exon")
    exon_type = &(Initial::getInstance());
  else if(type == "Single" || type == "single" || type == "single-exon")
    exon_type = &(Single::getInstance());
  else {
    throw runtime_error("unrecognized exon type: "+type);
  }
  thePrediction = new GenePrediction(theSl,score,-1,"",*exon_type);
  assert( thePrediction );
  currLocGene->DataSetAbstractProduct::insert( thePrediction );
  if(newGene)
    DataSetAbstractProduct::insert(currLocGene);
}

void
GenePredDataSet::label_exons(const string& seq) {
  sort();
  const DataStorageType::iterator stop = end();
  DataStorageType::iterator iter1 = begin();
  for(; iter1 != stop; ++iter1) {
    DataSetAbstractProduct* gPtr = dynamic_cast<DataSetAbstractProduct*>(*iter1);
    const unsigned num_exons = gPtr->size();
    const DataStorageType::iterator initial = gPtr->begin();
    DataStorageType::iterator exon_iter = gPtr->begin();
    const DataStorageType::iterator exon_stop = gPtr->end();
    unsigned ex_cnt = 0; 
    list<GenePrediction*> tlst;
    while(exon_iter != exon_stop ) {
      GenePrediction* exon_val = dynamic_cast<GenePrediction*>(*exon_iter);
      GenePrediction* np = NULL;
      if(num_exons == 1) {
	np  = new GenePrediction(*exon_val,Single::getInstance());
      } else if( exon_iter == initial ) {
	np  = exon_val->getStrnd() == dsu::ePos ? new GenePrediction(*exon_val,Initial::getInstance()) :new GenePrediction(*exon_val,Terminal::getInstance());
	
      } else if( ex_cnt == num_exons-1 ) {
	np  = exon_val->getStrnd() == dsu::ePos ? new GenePrediction(*exon_val,Terminal::getInstance()) :new GenePrediction(*exon_val,Initial::getInstance());
      } else {
	np  = new GenePrediction(*exon_val,Internal::getInstance());
      }
      assert(np);
      delete *exon_iter;
      DataStorageType::iterator next_iter = exon_iter;
      ++next_iter;
      gPtr->erase(exon_iter);
      tlst.push_back(np);
      //gPtr->insert(np);
      exon_iter = next_iter;
      ex_cnt++;
    }
    list<GenePrediction*>::iterator ti = tlst.begin();
    const list<GenePrediction*>::iterator is = tlst.end();
    for(; ti != is; ++ti) gPtr->insert(*ti);
    
  }
}
void
GenePredDataSet::add_introns_help() {
  const DataStorageType::const_iterator stop = end();
  DataStorageType::const_iterator iter1 = begin();
  const dsu::Strand_t strnd = getStrnd();
  assert(strnd != dsu::eEither);
  list<Prediction*> tlst;
  while(iter1 != stop) {
    GenePrediction* model1 = dynamic_cast<GenePrediction*>(*iter1);
    DataStorageType::const_iterator iter2 = iter1;
    ++iter2;
    while(iter2 != stop) {
      GenePrediction* model2 = dynamic_cast<GenePrediction*>(*iter2);
      if( !isOverlap(*model1,*model2) ) {
	assert(model2->getScoreIdx() != Intron::getInstance().getScoreIdx());
	const int rghtSide = model1->getEnd3()+1;
	const int leftSide = model2->getEnd5()-1;
	SeqLoc nl(strnd,rghtSide,leftSide);
	const GeneModelType& exon_type = Intron::getInstance();
	const double score = (model1->getScore().getVal()+model2->getScore().getVal())/2.0;
	Prediction* thePrediction = new GenePrediction(nl,score,-1,"",exon_type);
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

static const string geneStr = "One_Gene_Prediction_Set";
const string& GenePredDataSet::getTypeStr() const {
  return geneStr;
}
