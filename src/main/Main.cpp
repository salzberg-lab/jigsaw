//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "Main.h"
#include "ParamFile.h"
#include "Options.h"
#include "InputData.h"
#include "Combine.h"
#include "dsu.h" // my dna string utilities (starting over)
#include "CustomScore.h"
#include <stdexcept>
#include <string>
#include <fstream>
#include <time.h>
#include <stdlib.h>

using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;
using std::invalid_argument;
using std::ostream;

static void
fixGffFormat(list<GffFeature*>& glst);
static bool 
isSameGeneModelInPerfectOverlap(const GffFeature& left, const GffFeature& rght); 
static void
printGeneList(list<GffFeature*>& glst, ostream& os, const string& seqId, const string& progId); 
static void
mergeGeneList(const string&,list<GffFeature*>&, list<GffFeature*>&, Options&, CustomScore*, 
	      InputData&, const string&, int&); 

void
Main::run(Options& options, CustomScore* scoringSystem, InputData& inputData, 
	  const string& theSequence, const string& sequenceHeader, const string& program_id) {
  srand( time(NULL));

  if(options.quit) return;
  try {
    ofstream* predOutputStream = NULL;
    ofstream* matrixStream = NULL;
    const string& matrixFileName = options._matrixFile;
    if(matrixFileName.size() > 0) {
      matrixStream = new ofstream(matrixFileName.c_str());
    }
    if(!options.doTraining()) {
      predOutputStream = new ofstream(options._outputGeneModelFile.c_str()); 
      if(! (*predOutputStream) ) { 
	cerr<<"combiner output genemodels file, unable to open: ["<<options._outputGeneModelFile.c_str()<<"]"<<endl; 
	throw invalid_argument("could not continue");
      }
    }
    const bool forceWholeGeneModel = options._forceModels;
    const int cutOff = options._cutOff;
    int modelIdx=0;
    DataSetAbstractProduct* examples = inputData.getTrainExamples();
    if( examples ) {
      const int offset = options._trainOffSet;
      DataSetAbstractProduct::DataStorageType::const_iterator iter = examples->begin();
      const DataSetAbstractProduct::DataStorageType::const_iterator stop = examples->end();
      while(iter != stop) {
	const Prediction* pred1 = *iter;
	++iter;
	const int gEnd5 = pred1->getEnd5() - offset-1;
	int gEnd3 = pred1->getEnd3() + offset+1;
	while( iter != stop ) {
	  const Prediction* pred2 = *iter;
	  if( (pred2->getEnd5()-offset-1) < (pred1->getEnd3()+offset+1) ) {
	    gEnd3 = pred2->getEnd3() + offset+1;
	    ++iter;
	  } else {
	    break;
	  }
	  pred1 = pred2;
	}
	options._subEnd5 = gEnd5;
	options._subEnd3 = gEnd3;
	cout<<"extract examples around region: "<<gEnd5<<" "<<gEnd3<<endl;
	Combine comb(inputData,theSequence,forceWholeGeneModel,*scoringSystem, cutOff);
	comb.run();
	const Annotation* aPtr= NULL;
	if( !options.doTraining() ) {
	  Annotation& anot = comb.getAnnotation();
	  aPtr = &anot;
	  const string seq_id = string(&sequenceHeader.c_str()[1]);
	  const string strid = seq_id + string("\t") + program_id;
	  const unsigned id = anot.printOptimalParse(seq_id,*predOutputStream,modelIdx,*scoringSystem,strid); 
	  modelIdx = id;
	}
	if(aPtr && matrixStream ) {
	  (*matrixStream)<<"Annotation Matrix: "<<endl;
	  aPtr->print(*matrixStream,*scoringSystem);
	  (*matrixStream)<<endl;
	}
	else if(matrixStream) {
	  const list<RegionPrediction*>& rpall = comb.getOutput();
	  list<RegionPrediction*>::const_iterator iter = rpall.begin();
	  const list<RegionPrediction*>::const_iterator stop = rpall.end();
	  for(; iter != stop; ++iter) {
	    (*iter)->print(*matrixStream,*scoringSystem);
	    (*matrixStream)<<endl;
	  }
	}
      }
    } else {
      bool runOverlaps = false;
      if( !predOutputStream ) {
	cerr<<"Prediction output file not specified"<<endl;
	cerr<<"If you are in training, you may have forgotton to specify the 'curation' line in the evidene list file"<<endl;
	throw invalid_argument("could not continue");
      }
      if(! (*predOutputStream) ) { 
	      cerr<<"jigsaw genemodel output file, unable to open: ["<<options._outputGeneModelFile.c_str()<<"]"<<endl; 
	      throw invalid_argument("could not continue");
      }
      int end5 = options._subEnd5;
      int end3 = options._subEnd3;
      unsigned slen = theSequence.length();
      unsigned defInterval=slen;
      const ParamFile& pf = (*ParamFile::getInstance());
      const signed cMaxLen=pf.getMaxSeqLen(); //should be user parameter
      if( (options._subEnd5 != -1 && options._subEnd3 != -1) || cMaxLen == -1 ) {
	defInterval = (options._subEnd3-options._subEnd5) + 1;
      } else {
	end5 = 1;
	end3 = slen-1;
	if( (signed)slen > cMaxLen ) {
	  runOverlaps = true;
	  const int overlapLen = pf.getOverlapLen();
	  defInterval = (unsigned)(cMaxLen - overlapLen); //(unsigned)cMaxLen/10;
	  assert(defInterval>0);
	}
      }
      list<GffFeature*> glst;
      int modelIdx = 0;
      for(signed iter = end5; iter <= end3; iter += defInterval ) {
	if( runOverlaps ) {
	  options._subEnd5 = iter;
	  options._subEnd3 = (iter+(unsigned)cMaxLen>=slen?slen:iter+(unsigned)cMaxLen);
	}
	cout.flush();
	cout<<"proc subsec: "<<options._subEnd5<<" "<<options._subEnd3<<std::ends;
	cout<<endl;
	cout.flush();
	Combine comb(inputData,theSequence,forceWholeGeneModel,*scoringSystem, cutOff);
	comb.run();

	Annotation& anot = comb.getAnnotation();
	const string seq_id = string(&sequenceHeader.c_str()[1]);
	if( !runOverlaps ) {
	  const string strid = seq_id + string("\t") + program_id;
	  list<GffFeature*>* lst = anot.getGeneList(seq_id,modelIdx,*scoringSystem);
	  //anot.printOptimalParse(seq_id,*predOutputStream,modelIdx,*scoringSystem,strid);
	  printGeneList(*lst,*predOutputStream,seq_id,program_id);
	} else {
	  list<GffFeature*>* lst = anot.getGeneList(seq_id,modelIdx,*scoringSystem);
	  mergeGeneList(seq_id,glst,*lst,options,scoringSystem,inputData,theSequence,modelIdx);
	  delete lst;
	}
	if(matrixStream ) {
	  (*matrixStream)<<"Annotation Matrix: "<<endl;
	  anot.print(*matrixStream,*scoringSystem);
	  (*matrixStream)<<endl;
	}
      }
      if( runOverlaps ) {
	const string sid = string(&sequenceHeader.c_str()[1]);
	fixGffFormat(glst);
	printGeneList(glst,*predOutputStream,sid,program_id);
      }
    }
    if( predOutputStream )
      delete predOutputStream;
    if( matrixStream )
      delete matrixStream;
  } catch(const invalid_argument& ia) {
    cerr<<ia.what()<<endl;
  }
}

void
fixGffFormat(list<GffFeature*>& glst) {
  list<GffFeature*>::iterator iter = glst.begin();
  const list<GffFeature*>::iterator stop = glst.end();
  GffFeature* prev = NULL;
  bool fixState = false;
  for(; iter != stop; iter++) {
    GffFeature* gff =*iter;
    if( prev ) {
      //this initiates the gene label fixing state
      const bool isPartial = (gff->getType()=="internal-exon" || (gff->getType()=="terminal-exon" && gff->getStrnd()=='+') || (gff->getType()=="initial-exon" && gff->getStrnd()=='-'));
      if( prev->getCol8() != gff->getCol8() && isPartial) {
	fixState=true;
      }
    }
    if(fixState) {
      gff->setCol(8,prev->getCol8());
    }

    if( gff->getType() != "internal-exon" ) {
      fixState=false;
    }
    prev = gff;
  }
  
}

void
printGeneList(list<GffFeature*>& glst, ostream& os, const string& seqId, const string& progId) {
  const Options& options = *(Options::getInstance());
  os<<"# Evidence List File: "<<options.m_evidenceListFile<<endl;
  os<<"# Prediction Model Directory: "<<options._dtFile<<endl;
  list<GffFeature*>::iterator iter = glst.begin();
  list<GffFeature*>::iterator stop = glst.end();
  for(; iter != stop; iter++) {
    GffFeature* gff =*iter;
    gff->setCol(1,seqId);
    gff->setCol(2,progId);
    os<<*gff<<endl;
  }
}

void
mergeGeneList(const string& seq_id, list<GffFeature*>& mlst, list<GffFeature*>& nlst, Options& options, CustomScore* scoringSystem, 
	      InputData& inputData, const string& theSequence, int& modelIdx) {
  if( nlst.empty() ) {
    return;
  }
  if( mlst.empty() ) {
    list<GffFeature*>::iterator niter = nlst.begin();
    const list<GffFeature*>::iterator nend = nlst.end();
    for(; niter != nend; niter++) {
      GffFeature* ngff = *niter;
      mlst.push_back(ngff);
    }
  }  else {
    // get overlap iterators
    const list<GffFeature*>::iterator nend = nlst.end();
    list<GffFeature*>::iterator add = nlst.end();
    bool perfectOverlap=false;
    const GffFeature* mlast=mlst.back();
    const GffFeature* nfrst=nlst.front();
    if( mlast->getEnd3() < nfrst->getEnd5() ) {
      // you have to check if these divide two complete genes
      // or properly connect a single gene
      perfectOverlap=false;
    } else {
      list<GffFeature*>::iterator miter = mlst.end();
      --miter;
      while( miter != mlst.end() ) {
	const GffFeature* mgff = *miter;
	list<GffFeature*>::iterator niter = nlst.begin();
	bool breakFor=false;
	for(; niter != nend; niter++) {
	  const GffFeature* ngff = *niter;
	  if( ngff->getEnd5() > mgff->getEnd3() ) {
	    breakFor = true;
	    break;
	  }
	  if( //mgff->getEnd5() == ngff->getEnd5() &&
	      mgff->getEnd3() == ngff->getEnd3() &&
	      (mgff->getFrame() == ngff->getFrame() || 
	       mgff->getStrnd() == '+' && mgff->getType()=="terminal-exon")
	       &&
	      mgff->getStrnd() == ngff->getStrnd()
	      ) {
	    perfectOverlap = true;
	    add = niter;
	    ++add;
	    breakFor = true;
	    break;
	  } 
	  if( mgff->getEnd5() == ngff->getEnd5() &&
	      mgff->getFrame() == ngff->getFrame()
	       &&
	      mgff->getStrnd() == ngff->getStrnd()
	      ) {
	    perfectOverlap = true;
	    add = niter;
		 delete *miter;
	    mlst.erase(miter);
	    breakFor = true;
	    break;
	  } 
	  if( //mgff->getEnd5() == ngff->getEnd5() &&
	     mgff->getEnd3() > ngff->getEnd5() && mgff->getEnd3() < ngff->getEnd3() &&
	      mgff->getFrame() == ngff->getFrame() &&
	      mgff->getStrnd() == ngff->getStrnd()
	      ) {
	    GffFeature* gff = new GffFeature(*mgff);
	    bool myassert = (mgff->getEnd5() < ngff->getEnd5() && mgff->getEnd3() < ngff->getEnd3());
	    if(!myassert) {
	      // not handle crap out an use other method
	      breakFor = true;
	      delete gff;
	      break;
	    }
	    perfectOverlap = true;
	    add = niter;
	    ++add;
	    gff->setEnd3(mgff->getEnd3());
	    gff->setEnd5(ngff->getEnd5());
	    delete *niter;
	    delete *miter;
	    mlst.erase(miter);
	    mlst.push_back(gff);
	    breakFor = true;
	    break;
	  } 
	}
	if( !breakFor ) {
	  list<GffFeature*>::iterator del = miter;
	  delete *del;
	  --miter;
	  mlst.erase(del);
	}
	if( breakFor ) {
	  break;
	}
      }
    }
    if(perfectOverlap) {
      const GffFeature* prev = mlst.back();
      const GffFeature* curr = *add;
      const string& pId = prev->getCol8();
      const string& remId = curr->getCol8();
      bool isSame=false;
      if( isSameGeneModelInPerfectOverlap(*prev,*curr) ) {
	isSame=true;
      }
      for(; add != nend; ++add) {
	GffFeature* curr = *add;
	if( curr->getCol8()==remId && isSame ) {
	  curr->setCol8(pId);
	}
	mlst.push_back(curr);
      }
    } else {
      //cout<<"going into ugly zone\n";
      // it get's uglier
      list<GffFeature*>::iterator miter = mlst.end();
      const list<GffFeature*>::iterator mend = mlst.end();
      --miter;
      int geneCnt=0, nStart = -1;
      bool isTermExon = false;
      while(miter != mend) {
	const GffFeature* mgff = *miter;
	if( mgff->getType()=="single-exon" || 
	    (mgff->getType()=="initial-exon" && mgff->getStrnd()=='+') ||
	    (mgff->getType()=="terminal-exon" && mgff->getStrnd()=='-') ) {
	  ++geneCnt;
	  isTermExon = true;
	} else {
	  isTermExon = false;
	}
	nStart = mgff->getEnd5();
	list<GffFeature*>::iterator del = miter;
	--miter;
	mlst.erase(del);
	if( isTermExon && geneCnt>=2 && nStart < nlst.front()->getEnd5() ) {
	  break;
	}
      }
      list<GffFeature*>::iterator niter = nlst.begin();
      const list<GffFeature*>::iterator nend = nlst.end();
      int geneCntN=0, nStop = -1;
      isTermExon = false;
      while(niter != nend) {
	const GffFeature* ngff = *niter;
	if( ngff->getType()=="single-exon" || 
	    (ngff->getType()=="terminal-exon" && ngff->getStrnd()=='+') ||
	    (ngff->getType()=="initial-exon" && ngff->getStrnd()=='-') ) {
	  ++geneCntN;
	  isTermExon = true;
	} else {
	  isTermExon = false;
	}
	nStop = ngff->getEnd3();
	list<GffFeature*>::iterator del = niter;
	++niter;
	nlst.erase(del);
	if( isTermExon && geneCntN>=2 && nStop > mlst.back()->getEnd3() ) {
	  break;
	}
      }
      assert(nStart != -1 && nStop != -1);
      options._subEnd5 = nStart;
      options._subEnd3 = nStop;
      const bool forceWholeGeneModel = true;
      Combine comb(inputData,theSequence,forceWholeGeneModel,*scoringSystem, -1);
      comb.run();
      const Annotation& anot = comb.getAnnotation();
      list<GffFeature*>* lst = anot.getGeneList(seq_id,modelIdx,*scoringSystem);
      list<GffFeature*>::iterator riter = lst->begin();
      const list<GffFeature*>::iterator rend = lst->end();
      for(; riter != rend; ++riter ) {
	mlst.push_back(*riter);
      }
      niter = nlst.begin();
      for(; niter != nlst.end(); ++niter ) {
	mlst.push_back(*niter);
      }
      delete lst;
    }
  }
  nlst.clear();
}


bool 
isSameGeneModelInPerfectOverlap(const GffFeature& left, const GffFeature& rght) {
  const string& leftType = left.getType();
  const char leftStrnd = left.getStrnd();
  //const string& rghtType = rght.getType();
  const char rghtStrnd = rght.getStrnd();

  bool res = false;
  if( leftStrnd == rghtStrnd ) {
    if( (leftType == "internal-exon") || 
	(leftStrnd == '+' && leftType == "initial-exon") ||
	(leftStrnd == '-' && leftType == "terminal-exon") ) {
      res = true;
    }
  }
  return res;
}
