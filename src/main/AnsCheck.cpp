//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "AnsCheck.h"

#if 0
bool AnsCheck::checkCuration(const RegionPrediction& rp, dsu::Strand_t strnd) {
  DataSetAbstractProduct::DataStorageType::const_iterator iter = rp.begin(), stop = rp.end();
  while(iter != stop) {
    const Prediction* pred = *iter;
    if(pred->getStrnd() != strnd) {
      ++iter;
      continue;
    }
    const char* id = pred->getDataSource().c_str();
    if(strstr(id,"curation")) return true;
    ++iter;
  }
  return false;
}
#endif

static void checkGeneModel(const DataSetAbstractProduct& glist, const RegionPrediction& rp, 
			   dsu::Strand_t strnd, vector<bool>& bndryVec) {
  DataSetAbstractProduct::DataStorageType::const_iterator iter = glist.begin(), stop = glist.end();
  bool isFirst = true;
  int cnt = 0;
  while( iter != stop) {
    const Prediction& pred = **iter;
    if( glist.size() == 1 ) {
      if(pred.getEnd5() == rp.getEnd5() && strnd == pred.getStrnd()) {
	if(strnd == dsu::ePos)
	  bndryVec[EvScores::eStart] = true;
	else
	  bndryVec[EvScores::eStop] = true;
      }
      if(pred.getEnd3() == rp.getEnd3() && strnd == pred.getStrnd()) {
	if(strnd == dsu::ePos)
	  bndryVec[EvScores::eStop] = true;
	else
	  bndryVec[EvScores::eStart] = true;
      }
      return;
    }
    if(isFirst) {
      bool hit = false;
      if(glist.size() > 1 && pred.getEnd3() == rp.getEnd3() && strnd == pred.getStrnd()) {
	hit = true;
	if(strnd == dsu::ePos)
	  bndryVec[EvScores::eDon] = true;
	else
	  bndryVec[EvScores::eAcc] = true;
      }
      if(pred.getEnd5() == rp.getEnd5() && strnd == pred.getStrnd()) {
	hit = true;
	if(strnd == dsu::ePos)
	  bndryVec[EvScores::eStart] = true;
	else
	  bndryVec[EvScores::eStop] = true;
      }
      if( hit ) return;
      isFirst = false;
    } else if(cnt == (signed)(glist.size()-1)) {
      bool hit = false;
      if(pred.getEnd3() == rp.getEnd3() && strnd == pred.getStrnd()) {
	hit = true;
	if(strnd == dsu::ePos)
	  bndryVec[EvScores::eStop] = true;
	else
	  bndryVec[EvScores::eStart] = true;
      }
      if(pred.getEnd5() == rp.getEnd5() && strnd == pred.getStrnd()) {
	hit = true;
	if(strnd == dsu::ePos)
	  bndryVec[EvScores::eAcc] = true;
	else
	  bndryVec[EvScores::eDon] = true;
      }
      if( hit ) return;
    } else {
      bool hit = false;
      if(glist.size() > 1 && pred.getEnd3() == rp.getEnd3() && strnd == pred.getStrnd()) {
	hit = true;
	if(strnd == dsu::ePos)
	  bndryVec[EvScores::eDon] = true;
	else
	  bndryVec[EvScores::eAcc] = true;
      }
      if(pred.getEnd5() == rp.getEnd5() && strnd == pred.getStrnd()) {
	hit = true;
	if(strnd == dsu::ePos)
	  bndryVec[EvScores::eAcc] = true;
	else
	  bndryVec[EvScores::eDon] = true; 
      }
      if( hit ) return;
    }
    ++iter;
    ++cnt;
  }
}

void AnsCheck::onBorder(const DataSetAbstractProduct& glist, const RegionPrediction& rp, dsu::Strand_t strnd, vector<bool>& vecB) {
  DataSetAbstractProduct::DataStorageType::const_iterator iter = glist.begin(), stop = glist.end();
  while(iter != stop) {
    const Prediction& pred = **iter;
    const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(&pred);
    assert(dp);
    checkGeneModel(*dp, rp, strnd, vecB);
    for(unsigned i = 0; i < vecB.size(); ++i) {
      if( vecB[i] ) return;
    }
    ++iter;
  }
}
