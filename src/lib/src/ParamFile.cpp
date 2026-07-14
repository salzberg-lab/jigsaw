#include "ParamFile.h"
#include <fstream>
#include <cstring>
#include <stdexcept>
#include <string>
#include <iostream>
#include <sstream>

using std::string;
using std::endl;
using std::cout;
using std::ifstream;
using std::ofstream;
using std::ostringstream;
using std::istringstream;
using std::runtime_error;


ParamFile* ParamFile::_instance = NULL;

ParamFile::ParamFile(const string& file) {
  load(file);
}

void
ParamFile::createDefault(const string& file) {
  ofstream ofs(file.c_str());
  if( !ofs ) {
    ostringstream emsg;
    emsg<<"Unable to create default parameter file ["<<file<<"]"<<endl;
    throw runtime_error(emsg.str());
  }
  ofs<<"# Intron Length Penalty"<<endl;
  ofs<<"-1 10"<<endl;
  ofs<<"# Intergenic Length Penalty"<<endl;
  ofs<<"0 0"<<endl;
  ofs<<"# Internal Exon Length Penalty"<<endl;
  ofs<<"0 0"<<endl;
  ofs<<"# Minimum Intron Length"<<endl;
  ofs<<"10"<<endl;
  ofs<<"# Minimum Intergenic Length"<<endl;
  ofs<<"20"<<endl;
  ofs<<"# Alignment Connection Cutoff"<<endl;
  ofs<<"-1"<<endl;
  ofs<<"# Donor Consensus Sequence"<<endl;
  ofs<<"gt gc"<<endl;
  ofs<<"# Acceptor Consensus Sequence"<<endl;
  ofs<<"ag"<<endl;
  ofs<<"# Maximum Sequence Length / Overlap Length"<<endl;
  ofs<<"2000000 20000"<<endl;
}

static int justincase=0;
void
ParamFile::load(const string& file) {
  ifstream ifs(file.c_str());
  if( !ifs ) {
    cout<<"Unable to open parameter file["<<file<<"]"<<endl;
    cout<<"Creating default file: ["<<file<<"]"<<endl;
    createDefault(file);
    if( justincase==0 ) {
      ++justincase;
      load(file);
    } else {
      ostringstream emsg;
      emsg<<"Tried creating default parameter file, but still can not read: "<<file<<" sorry, exiting."<<endl;
      throw runtime_error(emsg.str());
    }
  }
  const int buff_size = 4096;
  char buff[buff_size];
  char buff2[buff_size];
  while( ifs.getline(buff,buff_size) ) {
    ifs.getline(buff2,buff_size);
    istringstream istrm(buff2);
    if( strstr(buff,"Intron Length Penalty") ) {
      istrm>>_intronLenPenalty>>_intronPenalty;
    } else if( strstr(buff,"Intergenic Length Penalty") ) {
      istrm>>_igrLenPenalty>>_igrPenalty;
    } else if( strstr(buff,"Internal Exon Length Penalty") ) {
      istrm>>_internalExonLenPenalty>>_internalExonPenalty;
    } else if( strstr(buff,"Minimum Intron Length") ) {
      istrm>>_minIntronLen;
    } else if( strstr(buff,"Minimum Intergenic Length") ) {
      istrm>>_minIgrLen;
    } else if( strstr(buff,"Alignment Connection Cutoff") ) {
      istrm>>_minAlignmentConnectionCutoff;
	 // older versions mispelled consensus
    } else if( strstr(buff,"Donor Consensus Sequence") || strstr(buff,"Donor Concensus Sequence") ) {
      string seq;
      while(istrm>>seq) { 
	_donors.push_back(seq);
      }
    } else if( strstr(buff,"Acceptor Consensus Sequence") || strstr(buff,"Acceptor Concensus Sequence") ) {
      string seq;
      while(istrm>>seq) { 
	_acceptors.push_back(seq);
      }
    } else if( strstr(buff,"Maximum Sequence Length / Overlap Length") ) {
      istrm>>_maxSeqLen>>_overlapLen;
    } else {
      string err = string("Unable to parse string: ") + buff;
      throw runtime_error(err.c_str());
    }
  }
}
