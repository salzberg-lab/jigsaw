#include "CustomScore.h"
#include "RegionPrediction.h"
#include <iostream>

#if _RUNTIME_DEBUG >= 2
extern bool MYDEBUG;
#endif

using std::cout;
using std::endl;

vector<FeatVec*> CustomScore::calcVec(const RegionPrediction* rpPtr, dsu::Strand_t strnd) const {
  const vector<DataDescrip*>& ddVec = getFeatMaps();
  vector<FeatVec*> feats(ddVec.size());
  for(unsigned j = 0; j < ddVec.size(); ++j) {
    const DataDescrip& fm = *(ddVec[j]);
    feats[j] = new FeatVec(fm.size());
    fm.initVec(*(feats[j]));
  }
  if( !rpPtr ) {
    return feats;
  }
  const RegionPrediction& rp = *rpPtr;
  Prediction::DataStorageType::const_iterator iter = rp.begin(), stop = rp.end();
  while(iter != stop) {
    const Prediction* pred = *iter;
    //if(pred->getStrnd() != strnd) {
      //++iter;
      //continue;
    //}
    for(unsigned k = 0; k < ddVec.size(); ++k) {
#if _RUNTIME_DEBUG >= 3
      if( MYDEBUG  && (EvScores::Bnd_t)k == EvScores::eDon) {
	int stop = 0;
	cout<<*pred<<endl;
      }
#endif
      const DataDescrip& fm = *(ddVec[k]);
      const int idx = fm.getIdx(pred);
      if(idx != -1) {
	EvidenceType* type = (*feats[k])[idx];
	type->update(pred,rp.getEnd5(),rp.getEnd3(),strnd);
      }
    }
    ++iter;
  }
  // give objects a chance to do recalculations here
  for(unsigned i = 0; i < feats.size(); ++i) {
    feats[i]->finish();
  }
  return feats;
}
