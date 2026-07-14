//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "PrintVec.h"
#include "DataDescrip.h"
#include "RegionPrediction.h"
#include "Internal.h"
#include "Intron.h"
#include "GenePrediction.h"
#include <algorithm>
#include <vector>
#include <math.h>
#include <cassert>

using std::cerr;
using std::endl;
using std::cout;

PrintVec* PrintVec::_instance = NULL;
#if _RUNTIME_DEBUG >= 2
extern bool MYDEBUG;
#endif
static string idx2str(int idx);
static bool 
checkIntronInterval(const DataSetAbstractProduct&, int, int, dsu::Strand_t); 

PrintVec::PrintVec(const Options& options, const DataSetAbstractProduct& dsap) : 
  _ostrmVec(_ddVec.size()), _geneList(dsap) {
  const string& id = options._vecFile;
  _end5 = options._subEnd5;
  _end3 = options._subEnd3;
  for(unsigned i = 0; i < _ostrmVec.size(); ++i) {
    string fname = id + idx2str(i); 
    _ostrmVec[i] = new std::ofstream(fname.c_str());
    if( !*_ostrmVec[i] ) {
      cerr<<"unable to write to: "<<fname<<endl;
      throw "failed to open file";
    }
  }
}

typedef std::pair<EvScores::Bnd_t,EvScores::Bnd_t> tpair_t;
static tpair_t 
getEnd5Bndry( const GenePrediction& pred) {
  EvScores::Bnd_t leftP = EvScores::eStart;
  EvScores::Bnd_t leftN = EvScores::eStop;
  if( pred.getScoreIdx() == Internal::getInstance().getScoreIdx() ) {
    leftP = EvScores::eAcc;
    leftN = EvScores::eDon;
  }
  // first parameter ("first") is for positive, 2nd is for negative
  return tpair_t(leftP,leftN);
}


static tpair_t 
getEnd3Bndry( const GenePrediction& pred) {
  EvScores::Bnd_t leftP = EvScores::eStop;
  EvScores::Bnd_t leftN = EvScores::eStart;
  if( pred.getScoreIdx() == Internal::getInstance().getScoreIdx() ) {
    leftP = EvScores::eDon;
    leftN = EvScores::eAcc;
  }
  // first parameter ("first") is for positive, 2nd is for negative
  return tpair_t(leftP,leftN);
}

static void checkGeneModel( const DataSetAbstractProduct& glist, 
			    const RegionPrediction& rp, dsu::Strand_t strnd, vector<bool>& bndryVec) {
  DataSetAbstractProduct::DataStorageType::const_iterator iter = glist.begin(), stop = glist.end();
  bool isFirst = true;
  int cnt = 0;
  const int cINTRON = Intron::getInstance().getScoreIdx();
  while( iter != stop) {
    const GenePrediction& pred = dynamic_cast<const GenePrediction&>(**iter);
    if( pred.getScoreIdx() == cINTRON ) {
      ++cnt;
      ++iter;
      continue;
    }
    if( rp.getEnd5() >= pred.getEnd5() && rp.getEnd3() <= pred.getEnd3() && strnd == pred.getStrnd())
      bndryVec[EvScores::eCoding] = true;
    tpair_t pn5 = getEnd5Bndry(pred);
    tpair_t pn3 = getEnd3Bndry(pred);
    if( glist.size() == 1 ) {
      if(pred.getEnd5() == rp.getEnd5() && strnd == pred.getStrnd()) {
	if(strnd == dsu::ePos) {
	  bndryVec[pn5.first] = true;
	} else {
	  bndryVec[pn5.second] = true;
	}
      }
      if(pred.getEnd3() == rp.getEnd3() && strnd == pred.getStrnd()) {
	if(strnd == dsu::ePos)
	  bndryVec[pn3.first] = true;
	else
	  bndryVec[pn3.second] = true;
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
	  bndryVec[pn5.first] = true;
	else
	  bndryVec[pn5.second] = true;
      }
      if( hit ) return;
      isFirst = false;
    } else if(cnt == (signed)(glist.size()-1)) {
      bool hit = false;
      if(pred.getEnd3() == rp.getEnd3() && strnd == pred.getStrnd()) {
	hit = true;
	if(strnd == dsu::ePos)
	  bndryVec[pn3.first] = true;
	else
	  bndryVec[pn3.second] = true;
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

static void onBorder( const DataSetAbstractProduct& glist, const RegionPrediction& rp, dsu::Strand_t strnd, vector<bool>& vecB) {
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


static bool 
isInValidRangeOfTrainGene( const DataSetAbstractProduct& trainGenes, const RegionPrediction& rp, const int cEnd5, const int cEnd3) {
  DataSetAbstractProduct::DataStorageType::const_iterator trainGeneIter = trainGenes.begin(), stop = trainGenes.end();
  while(trainGeneIter != stop) {
    const Prediction& pred = **trainGeneIter;
    const int gEnd5 = pred.getEnd5() - cEnd5;
    const int gEnd3 = pred.getEnd3() + cEnd3;
    if( (rp.getEnd5() >= gEnd5 && rp.getEnd5() <= gEnd3) ||
	(rp.getEnd3() >= gEnd5 && rp.getEnd3() <= gEnd3) ) {
      return true;
    }
    ++trainGeneIter;
  }
  return false;
}

#if 1
//static void printVec(const vector<EvidenceType*>& vec, std::ostream& os) {
static void printVec(const FeatVec& vec, std::ostream& os) {
  for(unsigned i = 0; i < vec.size(); ++i) {
    vec[i]->train_print(os);
    os<<" ";
  }
}
#endif

void PrintVec::score(vector<RegionPrediction*>& vrp, const string& str) const {
  RegionPrediction& rp = *vrp[1];
  score(rp,str);
}

CustomScore::LabelScore 
PrintVec::score(const RegionPrediction& rp, const string& str) const {
  const Options& opt = *(Options::getInstance());
  const int cOffSet = opt._trainOffSet;
  if( isInValidRangeOfTrainGene(_geneList,rp,cOffSet,cOffSet))  {
    scoreHelp(rp,str,dsu::ePos);
    scoreHelp(rp,str,dsu::eNeg);
  }
  return CustomScore::LabelScore(); // ignored during printout
}

void PrintVec::scoreHelp(const RegionPrediction& rp, const string& str, dsu::Strand_t strnd) const {
#if _RUNTIME_DEBUG >= 2
  bool localDebug = false;
  if(rp.getEnd5() == 70248 && strnd == dsu::ePos ) {
    MYDEBUG = localDebug = true;
    cout<<rp<<endl;
  } else 
    MYDEBUG = false;
#endif
  if( getEnd5() != -1 && getEnd3() != -1) {
    if( rp.getEnd5() < getEnd5() || rp.getEnd3() > getEnd3()) 
      return; 
  }
  vector<FeatVec*> vec = calcVec(&rp,strnd);
  assert(vec.size()==EvScores::eNumBnd);
  vector<bool> ans(EvScores::eNumBnd,false);
  onBorder( _geneList, rp, strnd, ans);
  ans[EvScores::eIntron] = checkIntronInterval(_geneList,rp.getEnd5(),rp.getEnd3(),strnd);

  PrintVec *shutup = const_cast<PrintVec*>(this);
  assert(_ostrmVec.size()==vec.size());
  for(unsigned i = 0; i < vec.size(); ++i) {
    std::ostream& ostrm = *(shutup->_ostrmVec[i]);
    ostrm<<rp.getEnd5()<<" "<<rp.getEnd3()<<" "<<dsu::strnd2str(strnd)<<" "; 
    FeatVec* feat = vec[i];
    printVec(*feat,ostrm);
    //printVec(*vec[i],ostrm);
    ostrm<<(ans[i]?1:0)<<endl;
    delete feat;
  }
}

static string idx2str(int idx) {
  string str;
  if( idx == EvScores::eCoding )
      str = ".code";
  else if( idx == EvScores::eAcc )
      str = ".acc";
  else if( idx == EvScores::eDon )
      str = ".don";
  else if( idx == EvScores::eStart )
      str = ".start";
  else if( idx == EvScores::eStop )
      str = ".stop";
  else if( idx == EvScores::eIntron )
      str = ".intron";
  else {
    throw "idx2str failed\n";
  }

  return str;
}

static bool 
checkIntronHelp( const DataSetAbstractProduct& glist, int end5, int end3, dsu::Strand_t strnd) {
  DataSetAbstractProduct::DataStorageType::const_iterator iter1 = glist.begin(), stop = glist.end();
  bool found = false;
  const int cINTRON = Intron::getInstance().getScoreIdx();
  for( ; iter1 != stop; ++iter1) {
    const GenePrediction& pred1 = dynamic_cast<const GenePrediction&>(**iter1);
    if( pred1.getScoreIdx() == cINTRON ) {
      if( strnd == pred1.getStrnd() && 
	  end5 >= pred1.getEnd5() && end5 <= pred1.getEnd3() &&
	  end3 >= pred1.getEnd5() && end3 <= pred1.getEnd3() ) {
	found = true;
	break;
      }
    }
#if 0
    DataSetAbstractProduct::DataStorageType::const_iterator iter2 = iter1;
    ++iter2;
    if(iter2 != stop) {
      const GenePrediction& pred2 = dynamic_cast<const GenePrediction&>(**iter2);
      //if( end5 == pred1.getEnd3() && end3 == pred2.getEnd5() ) {
      if( strnd == pred1.getStrnd() && 
	  end5 > pred1.getEnd3() && end5 < pred2.getEnd5() &&
	  end3 > pred1.getEnd3() && end3 < pred2.getEnd5() ) {
	found = true;
	break;
      }
    }
#endif
  }
  return found;
}

bool 
checkIntronInterval(const DataSetAbstractProduct& glist, int end5, int end3, dsu::Strand_t strnd ) {
  DataSetAbstractProduct::DataStorageType::const_iterator iter = glist.begin(), stop = glist.end();
  while(iter != stop) {
    const Prediction& pred = **iter;
    const DataSetAbstractProduct* dp = dynamic_cast<const DataSetAbstractProduct*>(&pred);
    assert(dp);
    bool val = checkIntronHelp(*dp, end5, end3, strnd);
    if(val)
      return val;
    ++iter;
  }
  return false;
}
