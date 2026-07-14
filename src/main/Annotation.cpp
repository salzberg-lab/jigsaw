//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "Annotation.h"
#include "localtypedef.h"
#include "IntronInterval.h"
#include "Compare.h"
#include "AnnotationItem.h"
#include "ScoreInterval.h"
#include "TraceBack.h"
#include "genemodels.h"
#include "CustomScore.h"
#include "Options.h"
#include <iostream>
#include <assert.h>
#include <stdexcept>
#include <sstream>

using std::cerr;
using std::cout;
using std::endl;

static bool 
isCon(int end5, int end3) {
  return (end5+1)==end3;
}

#if _RUNTIME_DEBUG >= 2
bool DEBUGON = false;
#endif


void Annotation::insert(const vector<RegionPrediction*>& rpVec, const CustomScore& cs) {
  // for each possiple continuation of the parse
  // find the maximal scoring parse for each reading frame
  // for each current piece of stored evidence that falls within 
  // the range of the new piece of evidence, re-evaluate, those
  // stored evidence scores and add new evidence.
  // Check all available predecessor parses
  // e.g. the dynamic algo. parse looks like this:
  // based on notation in MORGAN paper: salzberg et al. (1998)
  // for D[t,n] 
  // t = dynamic type of seqLocPtr (e.g. Terminal, Initial, etc.)
  // this is implicitly stored so you don't really have to worry about
  // this here.  
  // n = end3 (seqLocPtr->getEnd3())
  // for each available end3 "i" 1<= i < n
  // D[n] = max { D[i] + new_score }
  const bool allowPartialGene = !forceWholeGeneModel();
  assert(rpVec.size()==3);
  AnnotationItem* newItem = new AnnotationItem(rpVec, _posSeq, _negSeq);
  //cout<<"Pos: "<<newItem->getEnd5()<<endl;
  AnnotationMatrix::iterator iter = _matrix.end();
  --iter;
#if _RUNTIME_DEBUG >= 2
  bool debug=false;
  if(newItem->getEnd3()==47194 ) {
    debug = true;
  } else {
    debug = false;
  }
#endif
  // This assumes going right to left,
  // starting with the closest predecesor parse first
  AnnotationScore noCodeSum = 0;
  int num_intervals = 0;
  vector<scr_pr> inv_score(2);
  for(unsigned i = 0; i < 2; ++i) {
    inv_score[i].first=inv_score[i].second.first=inv_score[i].second.second = 0;
  }
  bool isContig = true;
  const AnnotationItem* last = newItem;
  while(iter != _matrix.end()) {
    const AnnotationItem& item = **iter;
#if _RUNTIME_DEBUG >= 2
    if((item.getEnd5() == 46992 || item.getEnd5()==47099)
       && debug) {
      int stop =0;
      DEBUGON=true;
    } else
      DEBUGON=false;
#endif
    if(isContig) {
      isContig=isCon(item.getEnd3(),last->getEnd5());
    }
    newItem->dpUpdate(item,inv_score,isContig,cs);
    const AnnotationScore t1 = cs.getProb(item,EvScores::eCoding,dsu::ePos);
    if(t1.isValid()) {
      inv_score[dsu::ePos].first += t1;
    }
    const AnnotationScore t2 = cs.getProb(item,EvScores::eCoding,dsu::eNeg);
    if(t2.isValid()) {
      inv_score[dsu::eNeg].first += t2;
    }
    const AnnotationScore nc = cs.getProb(item,EvScores::eNonCoding,dsu::eEither);
    //cout<<"what the fuck: "<<item.getEnd5()<<" "<<item.getEnd3()<<" "<<t2<<" "<<nc<<endl;
    assert(nc.isValid());
    inv_score[dsu::ePos].second.first += nc;
    inv_score[dsu::eNeg].second.first += nc;
    const AnnotationScore isp = cs.getProb(item,EvScores::eIntron,dsu::ePos);
    if(isp.isValid()) {
      inv_score[dsu::ePos].second.second += isp;
    }
    const AnnotationScore isn = cs.getProb(item,EvScores::eIntron,dsu::eNeg);
    if(isn.isValid()) {
      inv_score[dsu::eNeg].second.second += isn;
    }
    last = &item;
    --iter;
    ++num_intervals;
  }
#if _RUNTIME_DEBUG >= 2
  if( debug ) {
    int stop=0;
  }
#endif
  newItem->dpUpdate(inv_score,isContig,allowPartialGene,cs);
  _matrix.insert(newItem);
}

static GffFeature*
print_gff(const string& seq_id, std::ostream& os, const string& strid, const TraceBack::const_iterator& iter, 
	  const TraceBack::const_iterator& last, unsigned& model, bool saveGff) {
  const GenePrediction& curr_gp = **iter;
  const GenePrediction& last_gp = **last;
  // check if new gene model
  //cout<<"before: "<<model<<endl;
  if( last != iter ) {
    // cases for new gene models
    if( curr_gp.getStrnd() != last_gp.getStrnd() ) {
      ++model;
    } else if(curr_gp.getScoreIdx() == Single::getInstance().getScoreIdx()) {
      ++model;
    } else if(curr_gp.getStrnd() == dsu::ePos && 
	      curr_gp.getScoreIdx() == Initial::getInstance().getScoreIdx()) {
      ++model;
    } else if(curr_gp.getStrnd() == dsu::eNeg &&
	      curr_gp.getScoreIdx() == Terminal::getInstance().getScoreIdx()) {
      ++model;
    } 
  }
  os<<strid<<"\t";
  dsu::Strand_t strnd = curr_gp.getStrnd();
  const string& strndStr = dsu::strnd2str(strnd);
  string typeStr = genemodels::idx2str(curr_gp.getScoreIdx());
  typeStr += "-exon";
  int end5 = curr_gp.getEnd5();
  int end3 = curr_gp.getEnd3();
  int frame = curr_gp.getFrame();

  /** internally phase references left side, must be reoriented to 5' end
   ** which is on the right side on the negative strand
   */
  if( dsu::eNeg == strnd ) {
    const int exonLen = (end3-end5)+1;
    const int diff = (exonLen-frame) % 3;
    frame = diff;
  }

  std::ostringstream ostrm;
  GffFeature* res = NULL;
  if(saveGff) {
    ostrm<<"dum\t"<<"dum"<<"\t"<<typeStr<<"\t"<<end5<<"\t"<<end3<<"\t1\t"<<strndStr<<"\t"<<frame<<"\t"<<seq_id<<"."<<model;
    res = new GffFeature;
    res->fillFromBuffer(ostrm.str());
  } else {
    os<<typeStr<<"\t"<<end5<<"\t"<<end3<<"\t1\t"<<strndStr<<"\t"<<frame<<"\t"<<seq_id<<"."<<model<<endl;
  }
  return res;
}

const MatrixItemContents
Annotation::getBestStart(const CustomScore& cs) const {
  const bool allowPartialGene = !forceWholeGeneModel();
  const Options& opt = *Options::getInstance();
  const Compare* cmp = LtCompare::getInstance();
  if( opt.isLinComb() ) 
    cmp = GtCompare::getInstance();
  MatrixItemContents btb;
  AnnotationScore bscore = AnnotationScore::badVal();
  AnnotationMatrix::const_iterator iter = _matrix.end();
  --iter;
  AnnotationMatrix::const_iterator stop = _matrix.end();
  AnnotationScore t1=0;
  while(iter != stop) {
    const AnnotationItem* item = *iter;
#if _RUNTIME_DEBUG >= 2
    if(item->getEnd5() == 26400 ) {
      int stop = 0;
    }
#endif
    pair<AnnotationScore,MatrixItemContents> val = item->endScore(cs,t1,allowPartialGene);
    if( val.first.isValid() ) {
      AnnotationScore cscore = val.first;
#if _RUNTIME_DEBUG >= 2
      cout<<"getBestStart: "<<cscore<<" "<<item->getEnd5()<<" "<<item->getEnd3()<<" "<<bscore<<endl;
#endif
      if( cmp->compare(cscore,bscore) ) {
	bscore = cscore;
 	btb = val.second;
      }
    }
    --iter;
    /** 
     * currently ignores scoring partial ending introns
     * since we don't have a good representation for this right now
     * we use the generic non-coding score in all cases.
     **/
    t1 += cs.getProb(*item,EvScores::eNonCoding,dsu::eEither);
  }
  if( cmp->compare(t1,bscore) ) {
    btb = MatrixItemContents();
  }
  return btb;
}

int Annotation::printOptimalParse(const std::string& seq_id, std::ostream& os, int startModelIdx, const CustomScore& cs, const string& strid) const {
  const Options& options = *(Options::getInstance());
  if(options.isVerbose()) {
    cout<<"Number of Evidence Vectors: "<<_matrix.size()<<endl;
  }
  unsigned model = startModelIdx;
  if(startModelIdx == 0 ) {
    os<<"# Evidence List File: "<<options.m_evidenceListFile<<endl;
    os<<"# Prediction Model Directory: "<<options._dtFile<<endl;
    model++;
  }
  try 
    {
      const MatrixItemContents bstart = getBestStart(cs);
      TraceBack tb;
      tb.run(bstart);
      TraceBack::const_iterator iter = tb.begin();
      TraceBack::const_iterator last = tb.begin();
      const TraceBack::const_iterator stop = tb.end();
      const bool saveGff=false;
      while(iter != stop) {
	print_gff(seq_id,os,strid,iter,last, model, saveGff);
	last = iter;
	++iter;
      }
      return tb.begin() != stop ? model+1 : model;
    } 
  catch(const std::runtime_error& err) 
    {
      cerr<<err.what()<<endl;
      return -1;
    }
}

list<GffFeature*>*
Annotation::getGeneList(const string& seq_id, int& startModelIdx, const CustomScore& cs) const {
  unsigned model = startModelIdx;
  if(startModelIdx == 0 ) {
    model++;
  }
  list<GffFeature*>* lst = new list<GffFeature*>;
  try 
    {
      const MatrixItemContents bstart = getBestStart(cs);
      TraceBack tb;
      tb.run(bstart);
      TraceBack::const_iterator iter = tb.begin();
      TraceBack::const_iterator last = tb.begin();
      const TraceBack::const_iterator stop = tb.end();
      const bool saveGff=true;
      while(iter != stop) {
	GffFeature* gff = print_gff(seq_id,cout/*ignored*/,"",iter,last, model,saveGff);
	last = iter;
	++iter;
	lst->push_back(gff);
      }
      model = tb.begin() != stop ? model+1 : model;
    } 
  catch(const std::runtime_error& err) 
    {
      cerr<<err.what()<<endl;
      
    }
  startModelIdx = model;
  scoreGeneList(*lst,cs);
  return lst;
}

void Annotation::print(std::ostream& os, const CustomScore& cs) const {
  os<<"Annotation Matrix"<<endl;
  AnnotationMatrix::const_iterator iter = _matrix.begin();
  while(iter != _matrix.end()) {
    const AnnotationItem& matrix_cell = **iter;
    matrix_cell.print(os,cs);
    ++iter;
  }
}

void 
Annotation::scoreGeneList(list<GffFeature*>& flst, const CustomScore& cs) const {
  const int cIntronId = Intron::getInstance().getScoreIdx();
  DataSetAbstractProduct* geneLst = GffFeature::makeDataSet(flst);
  geneLst->add_introns();
  DataSetAbstractProduct::iterator iter = geneLst->begin();
  const DataSetAbstractProduct::iterator stop = geneLst->end();
  list<double> exon_scores;
  for(; iter != stop; ++iter) {
    DataSetAbstractProduct* prod = dynamic_cast<DataSetAbstractProduct*>(*iter);
    assert(prod);
    int num_occ = 0;
    AnnotationScore score_sum = 0;
    DataSetAbstractProduct::const_iterator giter = prod->begin();
    const DataSetAbstractProduct::const_iterator gstop = prod->end();
    for(; giter != gstop; ++giter) {
      const GenePrediction* pred = dynamic_cast<const GenePrediction*>(*giter);
      const int end5 = pred->getEnd5();
      const int end3 = pred->getEnd3();
      const dsu::Strand_t strnd = pred->getStrnd();
      const int mid = pred->getScoreIdx();
      const aipair_t pscore = scoreHelp(end5,end3,strnd,mid,cs);
      score_sum += pscore.first;
      num_occ += pscore.second;
      if( mid != cIntronId ) {
	const AnnotationScore exon_score = pscore.first / pscore.second;
	exon_scores.push_back(exon_score.getVal());
      }
    }
    if( num_occ <= 0 ) {
      cout<<"Warning in scoreGeneList, score set to 0"<<endl;
    }
    const AnnotationScore fval = num_occ > 0 ? score_sum/(AnnotationScore)num_occ : 0;
    prod->setScore(fval);
  }
  list<double>::const_iterator di = exon_scores.begin();
  const list<double>::const_iterator ds = exon_scores.end();
  assert(exon_scores.size()==flst.size());
  list<GffFeature*>::iterator gff_iter = flst.begin();
  const list<GffFeature*>::iterator gff_stop = flst.end();
  iter = geneLst->begin();
  for(; iter != stop; ++iter) {
    DataSetAbstractProduct* prod = dynamic_cast<DataSetAbstractProduct*>(*iter);
    assert(prod);
    AnnotationScore score_sum = 0;
    DataSetAbstractProduct::const_iterator giter = prod->begin();
    const DataSetAbstractProduct::const_iterator gstop = prod->end();
    for(; giter != gstop; ++giter ) {
      const GenePrediction* pred = dynamic_cast<const GenePrediction*>(*giter);
      if( pred->getScoreIdx() != cIntronId ) {
	assert(di != ds);
	assert(gff_iter != gff_stop);
	GffFeature* gff = *gff_iter;
	std::ostrstream os;
	os<<*di<<std::ends;
	gff->setCol5(os.str());

	std::ostrstream os1;
	const double score = prod->getScore().getVal();
	os1<<"; gene_score="<<score<<std::ends;
	gff->addComment(os1.str());
	
	++gff_iter;
	++di;
      }
    }
  }
}

Annotation::aipair_t
Annotation::scoreHelp(int end5, int end3, dsu::Strand_t strnd, int mid, const CustomScore& cs) const {
  const Options& opt = *Options::getInstance();
  const bool isLinComb = opt.isLinComb();
  const int cIgrId = Igr::getInstance().getScoreIdx();
  const int cIntronId = Intron::getInstance().getScoreIdx();
  const int cTerminalId = Terminal::getInstance().getScoreIdx();
  const int cSingleId = Single::getInstance().getScoreIdx();
  if( mid == cIntronId ) {
    end5-=1;
    end3+=1;
  }
  // This assumes going right to left,
  // starting with the closest predecesor parse first
  const AnnotationItem* lmatch = NULL, * rmatch = NULL, *sitem = NULL;
  AnnotationMatrix::const_iterator siter = _matrix.end();
  AnnotationMatrix::const_iterator riter = _matrix.end();
  const AnnotationMatrix::const_iterator stop = _matrix.end();
  AnnotationMatrix::const_iterator iter1 = _matrix.begin();
  for(; iter1 != stop; ++iter1) {
    const AnnotationItem* litem = *iter1;
    if( mid == cIntronId && litem->getEnd5()==end3 ) {
      rmatch = litem;
      riter = iter1;
      break;
    } else if( mid != cIntronId && litem->getEnd3()==end3 ) {
      rmatch = litem;
      riter = iter1;
      break;
    }
  }
  if( !rmatch ) {
    AnnotationMatrix::const_iterator iter1 = _matrix.begin();
    for(; iter1 != stop; ++iter1) {
      const AnnotationItem* litem = *iter1;
      // check that the alignments don't end 1 codon up from the stop codon
      if( sitem && ((dsu::ePos == strnd && mid == cTerminalId) || mid == cSingleId)) {
	if( sitem->getEnd3() == (end3-3) ) {
	  rmatch = sitem;
	  riter = siter;
	  break;
	}
      }
      sitem = litem;
      siter = iter1;
    }
  }
  if( !rmatch ) {
    cout<<"Warning unable to score this model, no 3' match "<<end5<<" "<<end3<<endl;
    return aipair_t(0,1);
  }
  assert(rmatch);
  iter1 = riter;
  AnnotationScore igrScore = 0;
  AnnotationScore sum = 0;
  int count = 0;
  sitem = NULL;
  for(; iter1 != _matrix.end(); --iter1) {
    const AnnotationItem* litem = *iter1;
    ++count;
    if( mid == cIntronId && litem->getEnd3()==end5 ) {
      lmatch = litem;
      break;
    } else if( mid != cIntronId && litem->getEnd5()==end5 ) {
      lmatch = litem;
      break;
    }
    const AnnotationScore cigr = cs.getProb(*litem,EvScores::eNonCoding,dsu::eEither);
    igrScore += cigr;
    if( mid == cIntronId ) {
      const AnnotationScore isp = cs.getProb(*litem,EvScores::eIntron,strnd);
      sum += isp;
    } else {
      assert( mid != cIgrId );
      const AnnotationScore isp = cs.getProb(*litem,EvScores::eCoding,strnd);
      sum += isp;
    }
    sitem = litem;
  }
  if( !lmatch ) {
    iter1 = riter;
    igrScore = 0;
    sum = 0;
    count = 0;
    sitem = NULL;
    for(; iter1 != _matrix.end(); --iter1) {
      const AnnotationItem* litem = *iter1;
      ++count;
      if( sitem && ((dsu::eNeg == strnd && mid == cTerminalId) || mid == cSingleId)) {
	if( sitem->getEnd5() == (end5+3) ) {
	  lmatch = sitem;
	  break;
	}
      }
      const AnnotationScore cigr = cs.getProb(*litem,EvScores::eNonCoding,dsu::eEither);
      igrScore += cigr;
      if( mid == cIntronId ) {
	const AnnotationScore isp = cs.getProb(*litem,EvScores::eIntron,strnd);
	sum += isp;
      } else {
	assert( mid != cIgrId );
	const AnnotationScore isp = cs.getProb(*litem,EvScores::eCoding,strnd);
	sum += isp;
      }
      sitem = litem;
    }
  }
  if( !lmatch ) {
    cout<<"Warning unable to score this model, no 5' match "<<end5<<" "<<end3<<endl;
    return aipair_t(0,1);
  }
  const AnnotationScore is=ScoreInterval::scoreOne(*lmatch,*rmatch,mid,strnd,sum,cs);
  AnnotationScore  ls(igrScore-is);
  if( isLinComb ) ls = AnnotationScore(is-igrScore); // means - max score always best score
  return aipair_t(ls,count);
}
