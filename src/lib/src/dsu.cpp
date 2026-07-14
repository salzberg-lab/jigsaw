//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "dsu.h"
#include <stdexcept>
#include <fstream>
#include <sstream>

using std::runtime_error;
using std::endl;

void dsu::readFastaFile(string& sequenceResult, const string& fastaFile) {
	const int BUFF_SIZE = 512000;
   char stra[BUFF_SIZE];
   std::ifstream ifs(fastaFile.c_str());
   if(!ifs) {
		std::ostringstream ostrm;
		ostrm<<"Unable to open fasta file (for reading): ["<<fastaFile<<"]"<<endl;
		throw std::runtime_error(ostrm.str());
   }
	sequenceResult = "";
   // load dna strings
	bool status = ifs.getline(stra,BUFF_SIZE);
	if( !status )
		throw runtime_error("unable to read fasta file");
	if( stra[0] != '>' ) 
		throw runtime_error("fasta does not start with header");

   while( ifs.getline(stra,BUFF_SIZE) ) {
      sequenceResult.append(stra);
   }
	
}

void dsu::readFastaHeader(string& headerResult, const string& fastaFile) {
   std::ifstream ifs(fastaFile.c_str());
   if(!ifs) {
      std::ostringstream ostrm; 
      ostrm<<"Unable to open fasta file (header read): ["<<fastaFile<<"]"<<endl;
      throw std::runtime_error(ostrm.str());
   }
   headerResult = "";
   // load dna strings
	if( ! (ifs>>headerResult) ) {
      throw std::runtime_error("unable to read fasta file, no fasta header");
	} else if( headerResult[0] != '>') {
      throw std::runtime_error("fasta does not start with header");
	}
}


dsu::Strand_t dsu::str2strnd(const string& parseStrand) {
     dsu::Strand_t theStrand = dsu::eEither;
     if(parseStrand == "-")
       theStrand = dsu::eNeg;
     else if(parseStrand == "+")
       theStrand = dsu::ePos;
     return theStrand;
}

dsu::Strand_t dsu::coord2strnd(int end5, int end3) {
  return end5 > end3 ? dsu::eNeg : dsu::ePos;
}

static const string posStr = "+";
static const string negStr = "-";
static const string eitherStr = "either";

const string& dsu::strnd2str(dsu::Strand_t strnd) {
	switch(strnd) {
		case ePos:
			return posStr;
		case eNeg:
			return negStr;
		default:
			return eitherStr;
	}
}

bool dsu::ignoreCaseEqual(const string& a, const string& b) {
  for(unsigned i = 0; i < a.size(); ++i) {
    if(tolower(a[i]) != tolower(b[i])) return false;
  }
  return (true && a.length()==b.length());
}
