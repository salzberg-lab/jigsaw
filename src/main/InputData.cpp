//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "InputData.h"
#include "Options.h"
#include "CustomScore.h"
#include "DataDescrip.h"
#include "EvidenceTypeRegistry.h"
#include "FormatRegistry.h"
#include "DnaStr.h"
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <algorithm>

using std::list;
using std::cout;
using std::cerr;
using std::endl;

using namespace dsu;

InputData::~InputData() {
  for(unsigned i = 0; i < _data.size(); ++i) {
    _data[i]->del();
    delete _data[i];
  }
  list<DataSetAbstractProduct*>::iterator iter = _stopLeak.begin();
  list<DataSetAbstractProduct*>::iterator stop = _stopLeak.end();
  while(iter != stop) {
    DataSetAbstractProduct* dp = *iter;
    dp->del();
    delete dp;
    ++iter;
  }
}

static string createUniqId(const string& filename, int id) {
  string::size_type pos = filename.rfind("/");
  if(pos == string::npos) {
    pos = filename.rfind("\\");
  }
  std::ostringstream ostrm;
  ostrm<<id<<'\0';
  string retval = filename.substr( pos==string::npos?0:pos+1,filename.size()) + ostrm.str();
  return retval;
}

static string getStr(const string& str, string::size_type& pos) {
  if(pos >= str.size()) {
    cerr<<"error parsing string: ["<<str<<"]"<<endl;
    cerr<<"Unable to continue, did you forgot to tell me what type of evidence this is?"<<endl<<endl<<endl;
    exit(1);
  } 
  string::size_type stop = str.find(" ",pos);
  if(stop==string::npos)
    stop = str.size();
  string ret = str.substr(pos,stop-pos);
  pos = stop+1;
  return ret;
}

static string
parseEvBuff(const string& str, string::size_type& pos) {
  if( pos >= str.size() ) {
    return "";
  }
  assert(pos < str.size());
  string::size_type stop = str.find(" ",pos);
  if(stop==string::npos)
    stop = str.size();
  string sstr = str.substr(pos,stop-pos);
  pos = stop+1;
  return sstr;
}

static std::pair<string,double>
parseEvBuffLin(const string& str, string::size_type& pos) {
  if( pos >= str.size() ) {
    return std::pair<string,double>("",-1.0);
  }
  assert(pos < str.size());
  string::size_type stop = str.find(" ",pos);
  if(stop==string::npos)
    stop = str.size();
  const string sstr = str.substr(pos,stop-pos);
  pos = stop+1;

  stop = str.find(" ",pos);
  if(stop==string::npos)
    stop = str.size();
  const string sstr2 = str.substr(pos,stop-pos);
  std::istrstream is(sstr2.c_str());
  double val = -1.0;
  is>>val;
  pos = stop+1;
  
  return std::pair<string,double>(sstr,val);
}


struct Compare : public std::binary_function<const SeqLoc*, const SeqLoc*,bool> {
  bool operator()(const SeqLoc* a, const SeqLoc* b) const {
    return *a < *b;
  }
};

static void local_insert(DataSetAbstractProduct* data, DataSetAbstractProduct& product) {
  Prediction::DataStorageType::iterator iter = data->begin();
  while( iter != data->end() ) {
    Prediction::DataStorageType* my_list = (*iter)->makePlaceHolders();
    Prediction::DataStorageType::iterator inner_iter = my_list->begin();
    while( inner_iter != my_list->end() ) {
      product.insert( *inner_iter );
      ++inner_iter;
    }
    delete my_list;
    ++iter;
  }
}

static bool
update_data(vector<DataDescrip*>& fms, const vector<EvidenceType*>& evd, const string& type) {
  bool someFound = false;
  for(unsigned i = 0; i < fms.size(); ++i) {
    if(!evd[i] && type != "curation") {
      //someFound = true;
    } else {
      someFound = true;
      DataDescrip& dd = *(fms[i]);
      if( evd[i] ) {
	dd.add(evd[i]);
      }
    }
  }
  return someFound;
}


static unsigned
ev2idx(const string& evidenceType) {
  int idx=0;
  if( ignoreCaseEqual(evidenceType,"coding") ) {
    idx = (unsigned)EvScores::eCoding;
  } else if( ignoreCaseEqual(evidenceType,"acc") ) {
    idx = (unsigned)EvScores::eAcc;
  } else if( ignoreCaseEqual(evidenceType,"don") ) {
    idx = (unsigned)EvScores::eDon;
  } else if( ignoreCaseEqual(evidenceType,"start") || ignoreCaseEqual(evidenceType,"overlap_start") ) {
    idx = (unsigned)EvScores::eStart;
  } else if( ignoreCaseEqual(evidenceType,"stop") ) {
    idx = (unsigned)EvScores::eStop;
  } else if( ignoreCaseEqual(evidenceType,"intron") ) {
    idx = (unsigned)EvScores::eIntron;
  } else {
    cerr<<"homology Evidence type: ["<<evidenceType<<"] not supported, use: coding,acc,don,start,stop and intron "<<endl;
  }
  return idx;
}

InputData::InputData() : _trainExamples(NULL) { }

void
InputData::readEvidenceFile(const string& filename, const std::string& theSeq, vector<DataDescrip*>& fms) {
  DataSetAbstractProduct product(SeqLoc(dsu::eEither,-1,-1),-1,filename,AnnotationScore::badVal()); 
  // file contains a list of filenames
  // that hold the actual annotation data
  std::ifstream ifs(filename.c_str());
  if(!ifs) {
     throw std::invalid_argument(filename + " : not a valid filename");
  }
  // known data formats
  DataSetAbstractFactory* factory = NULL;
  const int BUFF_SIZE = 15000;
  char buff[BUFF_SIZE];
  // loop through each file
  int numDataSources = 0;
  bool someDataFound =false;
  bool useIntron = false;
  while(ifs.getline(buff,BUFF_SIZE)) {
    // let "#" denote a comment, ignore these lines
    if( *buff == '#' ) {
      continue;
    }
    factory = NULL;
    string::size_type pos = 0;
    string parseFileName = getStr(buff,pos);
    string format = getStr(buff,pos);
    string type = getStr(buff,pos);
    string usid = createUniqId(parseFileName,numDataSources); 
    vector<EvidenceType*> evd(fms.size());
    for(unsigned i = 0; i < evd.size(); ++i) {
      evd[i] = NULL;
    }
    for(;;) {
      const string evidenceType = parseEvBuff(buff,pos);
      if( ignoreCaseEqual(evidenceType,"intron") ) {
	useIntron = true;
      } 
      if( evidenceType.length() == 0 )
	break;
      unsigned idx = ev2idx(evidenceType);
      const int uniqId = fms[idx]->size();
      evd[idx] = EvidenceTypeRegistry::createNewEvidenceType(evidenceType,type,usid,uniqId,theSeq);
    }
    if( !ignoreCaseEqual(type,"curation") ) {
      factory = FormatRegistry::getFactory(format,type);
    } else {
      factory = NULL;
    }
    if(factory) {
      try {
	const bool someFound1 = update_data(fms,evd,type);
	if( !someFound1 ) {
	  cerr<<"warning: ["<<type<<"] not described, wont be used in a feature vector"<<endl;
	}
	DataSetAbstractProduct* data = factory->createDataSetAbstractProduct(parseFileName,usid,type,theSeq);
	data->sort();
	data->add_introns();
	data->del_temp();
	local_insert( data , product);
	_stopLeak.push_back(data);
	someDataFound = true;
      } catch(std::exception& e) {
	cerr<<e.what()<<endl;
	cerr<<"Parsing: ["<<parseFileName<<"] failed."<<endl;
	cerr<<"Attempting to continue..."<<endl;
      }
      ++numDataSources;
    }
  }
  if( !someDataFound || product.empty() ) {
    cerr<<"fatal error, no input data read, check the input file: ["<<filename<<"]\n";
    cerr<<"possibly none of the files were readable\n";
    exit(1);
  }
  DataSetAbstractProduct* trainExamples = getTrainExamples();
  if( trainExamples ) {
    local_insert( trainExamples, product );
  }
  _data.resize( product.size() );
  DataSetAbstractProduct::iterator iter = product.begin();
  const DataSetAbstractProduct::iterator end = product.end();
  int cnt = 0;
  while( iter != end ) {
    _data[cnt] = *iter;
    ++iter;
    ++cnt;
  }
  vector<Prediction*>::iterator start = _data.begin();
  const vector<Prediction*>::iterator stop = _data.end();
  std::sort<vector<Prediction*>::iterator,Compare>( start, stop, Compare() );
  if( ! useIntron ) {
    Options& opt = *(Options::getInstance());
    opt.setIntronFalse();
  }
}

std::ostream& operator<<(std::ostream& os, const InputData& id) {
  os<<"implment print";
  return os;
}

const string singleStr = "All Input Data";
const string& InputData::getTypeStr() const {
  return singleStr;
}


void
InputData::readEvidenceFileLinComb(const string& filename, const string& theSeq, 
			    vector<DataDescrip*>& fms) {
  DataSetAbstractProduct product(SeqLoc(dsu::eEither,-1,-1),-1,filename,AnnotationScore::badVal()); 
  // file contains a list of filenames
  // that hold the actual annotation data

  
  std::ifstream ifs(filename.c_str());
  if(!ifs) {
     throw std::invalid_argument(filename + " : not a valid filename");
  }
  // known data formats
  DataSetAbstractFactory* factory = NULL;
  const int BUFF_SIZE = 15000;
  char buff[BUFF_SIZE];
  // loop through each file
  int numDataSources = 0;
  bool someDataFound =false;
  bool useIntron = false;
  while(ifs.getline(buff,BUFF_SIZE)) {
    // let "#" denote a comment, ignore these lines
    if( *buff == '#' ) {
      continue;
    }
    factory = NULL;
    string::size_type pos = 0;
    string parseFileName = getStr(buff,pos);
    string format = getStr(buff,pos);
    string type = getStr(buff,pos);
    string usid = createUniqId(parseFileName,numDataSources); 
    vector<EvidenceType*> evd(fms.size());
    for(unsigned i = 0; i < evd.size(); ++i) {
      evd[i] = NULL;
    }
    for(;;) {
      const pair<string,double> evidenceType = parseEvBuffLin(buff,pos);
      if( dsu::ignoreCaseEqual(evidenceType.first,"intron") ) {
	useIntron = true;
      } 
      if( evidenceType.first.length() == 0 )
	break;
      unsigned idx = ev2idx(evidenceType.first);
      const int uniqId = fms[idx]->size();
      evd[idx] = EvidenceTypeRegistry::createNewEvidenceType(evidenceType.first,type,usid,uniqId,theSeq,evidenceType.second);
    }
    if( !dsu::ignoreCaseEqual(type,"curation") ) {
      factory = FormatRegistry::getFactory(format,type);
    } else {
      factory = NULL;
    }
    if(factory) {
      try {
	const bool someFound1 = update_data(fms,evd,type);
	if( !someFound1 ) {
	  cerr<<"warning: ["<<type<<"] not described, wont be used in a feature vector"<<endl;
	}
	DataSetAbstractProduct* data = factory->createDataSetAbstractProduct(parseFileName,usid,type,theSeq);
	data->sort();
	data->add_introns();
	data->del_temp();
	//cout<<"temp debug: "<<parseFileName<<" "<<usid<<" "<<format<<endl;
	//cout<<*data<<endl;
	local_insert( data , product);
	_stopLeak.push_back(data);
	someDataFound = true;
      } catch(std::exception& e) {
	cerr<<e.what()<<endl;
	cerr<<"Parsing: ["<<parseFileName<<"] failed."<<endl;
	cerr<<"Attempting to continue..."<<endl;
      }
      ++numDataSources;
    }
  }
  if( !someDataFound || product.empty() ) {
    cerr<<"fatal error, no input data read, check the input file: ["<<filename<<"]\n";
    cerr<<"possibly none of the files were readable\n";
    exit(1);
  }
  DataSetAbstractProduct* trainExamples = getTrainExamples();
  if( trainExamples ) {
    local_insert( trainExamples, product );
  }
  _data.resize( product.size() );
  DataSetAbstractProduct::iterator iter = product.begin();
  const DataSetAbstractProduct::iterator end = product.end();
  int cnt = 0;
  while( iter != end ) {
    _data[cnt] = *iter;
    ++iter;
    ++cnt;
  }
  vector<Prediction*>::iterator start = _data.begin();
  const vector<Prediction*>::iterator stop = _data.end();
  std::sort<vector<Prediction*>::iterator,Compare>( start, stop, Compare() );
  if( ! useIntron ) {
    Options& opt = *(Options::getInstance());
    opt.setIntronFalse();
  }
  for(unsigned i = 0; i < fms.size(); ++i) {
    fms[i]->calcTotWght();
  }
}
