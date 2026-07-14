#include "TypeManage.h"
#include "DataSetAbstractProduct.h"
#include <list>

using std::list;

pair<const Prediction*,const Prediction*>
TypeManage::checkForIntronPrediction(const Prediction* pred5, int end5, const Prediction* pred3, int end3) {
  assert(end5 < end3);
  // since we're looking for an intron the 5' coord is the 3' of the upstream exon
  // and the 3' coord is the 5' coord of the downstream exon
  pair<const Prediction*,const Prediction*> rval(NULL,NULL);
  if(pred5 == pred3) {
    const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(pred5);
    assert(dp);
    pair<const Prediction*,const Prediction*> mtch = dp->getRelatedPredictionInterval(end5,end3);
    const Prediction* mtch5 = mtch.first;
    const Prediction* mtch3 = mtch.second;
    if( mtch5 && mtch3 ) {
      if( mtch5->getDataSource() == mtch3->getDataSource() ) {
	// what alignment score we use needs to be determined by the database type!!!
	// not sure if this is happening, so default to % identity.
	rval.first = mtch5;
	rval.second = mtch3;
      }
    }
  }
  return rval;
}

/** 
 * this version checks to see if the coordinates end5, end3 fall with in a intronic
 * region predicted by pred
 */
pair<const Prediction*,const Prediction*>
TypeManage::checkForIntronPrediction(const Prediction* pred, int end5, int end3) {
  assert(end5 < end3);
  // since we're looking for an intron the 5' coord is the 3' of the upstream exon
  // and the 3' coord is the 5' coord of the downstream exon
  pair<const Prediction*,const Prediction*> rval(NULL,NULL);
  const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(pred);
  assert(dp);
  //pair<const Prediction*,const Prediction*> mtch = dp->getRelatedPredictionInterval(end5,end3);
  pair<const Prediction*,const Prediction*> mtch = dp->getIntronBoundaries(end5,end3);
  const Prediction* mtch5 = mtch.first;
  const Prediction* mtch3 = mtch.second;
  if( mtch5 && mtch3 ) {
    assert( mtch5->getDataSource() == mtch3->getDataSource() );
    rval.first = mtch5;
    rval.second = mtch3;
  }
  return rval;
}


const MatchPrediction* 
TypeManage::checkForMatchPrediction(const Prediction* pred, int end5, int end3) {
  const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(pred);
  if(!dp) return NULL;
  list<const Prediction*> lst;
  dp->getRelatedPrediction(end5,end3,lst);
  assert(!lst.empty());
  list<const Prediction*>::const_iterator iter = lst.begin();
  list<const Prediction*>::const_iterator stop = lst.end();
  const MatchPrediction* gp = NULL;
  while(iter != stop) {
    const MatchPrediction* curr =  dynamic_cast<const MatchPrediction*>(*iter);
    if( !gp ) {
      gp = curr;
    } else if( curr->getPcntSim() > gp->getPcntSim() ) {
      gp = curr;
    }
    ++iter;
  }
  return gp;
}

const MatchPrediction* 
TypeManage::checkForMatchPrediction5(const Prediction* pred, int end5) {
  const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(pred);
  if(!dp) return NULL;
  list<const Prediction*> lst;
  dp->getRelatedPrediction5(end5,lst);
  list<const Prediction*>::const_iterator iter = lst.begin();
  list<const Prediction*>::const_iterator stop = lst.end();
  const MatchPrediction* gp = NULL;
  while(iter != stop) {
    const MatchPrediction* curr =  dynamic_cast<const MatchPrediction*>(*iter);
    if( !gp ) {
      gp = curr;
    } else if( curr->getPcntSim() > gp->getPcntSim() ) {
      gp = curr;
    }
    ++iter;
  }
  return gp;
}

const MatchPrediction* 
TypeManage::checkForMatchPrediction3(const Prediction* pred, int end3) {
  const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(pred);
  if(!dp) return NULL;
  list<const Prediction*> lst;
  dp->getRelatedPrediction3(end3,lst);
  list<const Prediction*>::const_iterator iter = lst.begin();
  list<const Prediction*>::const_iterator stop = lst.end();
  const MatchPrediction* gp = NULL;
  while(iter != stop) {
    const MatchPrediction* curr =  dynamic_cast<const MatchPrediction*>(*iter);
    if( !gp ) {
      gp = curr;
    } else if( curr->getPcntSim() > gp->getPcntSim() ) {
      gp = curr;
    }
    ++iter;
  }
  return gp;
}

const GenePrediction* 
TypeManage::checkForGenePrediction(const Prediction* pred, int end5, int end3) {
  const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(pred);
  if(!dp) return NULL;
  list<const Prediction*> lst;
  dp->getRelatedPrediction(end5,end3,lst);
  list<const Prediction*>::const_iterator iter = lst.begin();
  list<const Prediction*>::const_iterator stop = lst.end();
  const GenePrediction* gp = NULL;
  while(iter != stop) {
    const GenePrediction* curr =  dynamic_cast<const GenePrediction*>(*iter);
    if( !gp ) {
      gp = curr;
    } else if( curr->getScore() > gp->getScore() ) {
      gp = curr;
    }
    ++iter;
  }
  return gp;
}

const GenePrediction* 
TypeManage::checkForGenePrediction5(const Prediction* pred, int end5) {
  const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(pred);
  if(!dp) return NULL;
  list<const Prediction*> lst;
  dp->getRelatedPrediction5(end5,lst);
  list<const Prediction*>::const_iterator iter = lst.begin();
  list<const Prediction*>::const_iterator stop = lst.end();
  const GenePrediction* gp = NULL;
  while(iter != stop) {
    const GenePrediction* curr =  dynamic_cast<const GenePrediction*>(*iter);
    if( !gp ) {
      gp = curr;
    } else if( curr->getScore() > gp->getScore() ) {
      gp = curr;
    }
    ++iter;
  }
  return gp;
}

const GenePrediction* 
TypeManage::checkForGenePrediction3(const Prediction* pred, int end3) {
  const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(pred);
  if(!dp) return NULL;
  list<const Prediction*> lst;
  dp->getRelatedPrediction3(end3,lst);
  list<const Prediction*>::const_iterator iter = lst.begin();
  list<const Prediction*>::const_iterator stop = lst.end();
  const GenePrediction* gp = NULL;
  while(iter != stop) {
    const GenePrediction* curr =  dynamic_cast<const GenePrediction*>(*iter);
    if( !gp ) {
      gp = curr;
    } else if( curr->getScore() > gp->getScore() ) {
      gp = curr;
    }
    ++iter;
  }
  return gp;
}
