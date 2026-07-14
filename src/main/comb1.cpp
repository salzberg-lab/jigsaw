//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "Main.h"
#include "ParamFile.h"
#include "PrintVec.h"
#include "ExonDistr.h"
#include "Custom7.h"
#include "LinComb.h"
#include "LinComb2.h"
#include "FormatRegistry.h"
#include "InputData.h"
#include <fstream>
#include <stdexcept>
#include <sstream>
#if __GNUC__ > 4 || ( __GNUC__ == 4 && __GNUC_MINOR__ >= 3 )
 #include <utility>
 using namespace std;
#else
 #include <pair.h>
#endif

using std::cerr;
using std::endl;
using std::cout;
using std::exception;

static pair<string,string> getTrainData(const string&); 
static void usage(const char*);

int main(int argc, char** argv) 
{
  string cvs_tag = "3.2.9";
  string id = string(argv[0]) + string("-") + cvs_tag;
  string cid = id + string(": predict gene models using multiple sources of evidence.");
  try 
    {
      Options::createInstance(argc,argv,cid);
      Options& options = *(Options::getInstance());
      if( (options.m_fastaFile.size()==0 || 
	   options.m_evidenceListFile.size()==0)  && 
	  options._vecFile.size()==0 ) {
	usage(argv[0]);
	return 0;
      }
      string pfile = "param.txt";
      if( options._dtFile != "" ) {
	pfile = options._dtFile + string("/param.txt");
      }
     
      ParamFile::createInstance(pfile);
      ExonDistr::createInstance(options._exonDistrFile);
      if(options.quit) return 0;

      const char* procfile = options.m_evidenceListFile.c_str();
      string theSequence, sequenceHeader;

      /** maybe replace with bioinformatics c++ toolkit **/
      dsu::readFastaFile(theSequence,options.m_fastaFile);
      dsu::readFastaHeader(sequenceHeader,options.m_fastaFile);
      /** bioinformatics c++ toolkit 
      bioinfo::fasta_seq fseq(options.m_fastaFile);
      const string& theSequence = fseq.get_seq();
      const string& sequenceHeader = fseq.get_header();
      **/
      CustomScore* cs = NULL;
      DataSetAbstractProduct* geneList = NULL;
      InputData* inputData = new InputData;
      std::pair<string,string> tdata = getTrainData(options.m_evidenceListFile);
      const string format = tdata.first;
      if(format != "" ) {
	const string tfile = tdata.second;
	const string predtype = "curation";
	DataSetAbstractFactory* factory = FormatRegistry::getFactory(format,predtype);
	if(!factory) {
	  cerr<<"format: "<<format<<" for "<<tfile<<" "<<predtype<<" not recognized, can't continue."<<endl;
	  return 1;
	}
	// the "curation" string id
	// value is picked up by the "PrintVec" object
	// to identify which are the true genes for training
	geneList = 
	  factory->createDataSetAbstractProduct(tfile,predtype,predtype,theSequence);
	geneList->add_introns();
	geneList->del_temp();
	inputData->setTrainExamples(geneList);
      }
      if(options.doTraining()) {
	PrintVec::newInstance(options, *geneList);
	cs = PrintVec::getInstance();
	inputData->readEvidenceFile(procfile,theSequence,cs->getFeatMaps());
      } else if( options.isLinComb2() ) {
        LinComb2::newInstance();
        cs = LinComb2::getInstance();
        inputData->readEvidenceFileLinComb(procfile,theSequence,cs->getFeatMaps());
	id += "-lin";
      } else if( options.isLinComb() ) {
        LinComb::newInstance();
        cs = LinComb::getInstance();
        inputData->readEvidenceFileLinComb(procfile,theSequence,cs->getFeatMaps());
	id += "-lin";
      } else {
	Custom7::newInstance();
	cs = Custom7::getInstance(); 
	inputData->readEvidenceFile(procfile,theSequence,cs->getFeatMaps());
	Custom7* c7 = dynamic_cast<Custom7*>(cs);
	if( options._dtFile.size() <= 0  ){
	  cerr<<"Please specify the training directory\n";
	  return 0;
	}
	c7->createTree(options._dtFile);
      }
      Main main;
      main.run(options, cs, *inputData, theSequence,sequenceHeader, id);
      delete cs;
      delete inputData;			
      if(geneList)
	delete geneList;
    } 
  catch(const exception& err) 
    {
      cerr<<"Unhandled exception, unable to continue."<<endl;
      cerr<<err.what()<<endl;
    }
  return 0;
}

static pair<string,string> getTrainData(const string& fname) {
  std::ifstream ifs(fname.c_str());
  if(!ifs) {
    cerr<<"unable to open training data file: ["<<fname<<"]"<<endl;
    exit(1);
  }
  string name,type;
  const int BUFF_SIZE = 15000;
  char buff[BUFF_SIZE];
  while( ifs.getline(buff,BUFF_SIZE) ) {
    if( *buff == '#' ) {
      continue;
    }
    string tmp(buff);
    if(tmp.find("curation") != string::npos) {
      std::istringstream istrm(buff);
      string name,type,format;
      istrm>>name>>format>>type;
      return pair<string,string>(format,name);
    }
  }
  return pair<string,string>("","");
  //no longer a precondition
  //cerr<<"problems finding the training file...didn't find error: ["<<fname<<"]\n";
  //exit(1);
}

void usage(const char* prog) {
  cout<<"----------------------------------------------------------------"<<endl;
  cout<<"Combiner: predict gene models from multiple sources of evidence."<<endl;
  cout<<"----------------------------------------------------------------"<<endl<<endl;
  cout<<endl<<"usage: "<<prog;
  cout<<" -f \"filename\" -e \"filename\" -m \"filename\" -d \"directory\" [-lhxvS] [-t filename] [-b end5 ] [ -n end3]\n\n";
  cout<<endl;
  cout<<"\t-f\tcombiner will construct gene models for the genomic sequence specified in this file\n";
  cout<<endl;
  cout<<"\t-e\tfile lists the location of each evidence file and its associated type.\n";
  cout<<"\t\tformat is: filename (locatable path), file format (default, or btab),\n";
  cout<<"\t\tevidence type (homology, homgene, geneprediction or spliceprediction), prediction type (acc,don,coding,start,start_overlap,stop,intron)\n\n";
  cout<<endl;
  cout<<"\t\texample using statistical option:\n\n";
  cout<<"\t\t/data/estAlignments.ph default homology acc don coding intron\n";
  cout<<"\t\t/data/glimmerM.gp default geneprediction acc don coding start stop intron\n";
  cout<<"\t\t/data/genesplicer.sp default spliceprediction acc don\n\n";
  cout<<"\t\texample using linear (-l) option:\n\n";
  cout<<"\t\t/data/estAlignments.ph default homology acc 1.0 don 1.0 coding 1.0 intron 1.0\n";
  cout<<"\t\t/data/glimmerM.gp default geneprediction acc 1.0 don 1.0 coding 1.0 start 1.0 stop 1.0 intron 1.0\n";
  cout<<"\t\t/data/genesplicer.sp default spliceprediction acc 1.0 don 1.0\n\n";
  cout<<"\t-d\tdirectory where the decision trees are located\n";
  cout<<endl;
  cout<<"\t-m\tname of output file to place predicted genemodel list in\n";
  cout<<endl;
  cout<<"\t-x\tname of output file to place matrix info in. Can be useful for analyzing how predictions are generated.\n";
  cout<<endl;
  cout<<"\t-t\tprints feature vectors to a sepcified file\n";
  cout<<endl;
  cout<<"\t-q\tspecify the number of bases to use surrounding an input gene (e.g. look at +-20,0000 bases surrounding a gene)\n";
  cout<<endl;
  cout<<"\t-p\tuse subset of a specified sequence\n";
  cout<<endl;
  cout<<"\t-S\tlist recognized file formats\n";
  cout<<endl;
  cout<<"\t-l\tuse linear combiner\n";
  cout<<endl;
  cout<<"\t-h\tprints out this usage message\n";
  cout<<endl;
  cout<<"\t-v\tversion information\n";
  cout<<endl;
  cout<<"\n\tA typical execution of the combiner from the command line:\n";
  cout<<endl;
  cout<<"\t"<<prog<<" -f sequence.fasta -e evidence.data -m outputmodel.gff -d training_directory \n\n";
}
