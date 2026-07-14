#include "GffFeature.h"
#include "DataSetAbstractProduct.h"
#include "GenePredDataSet.h"
#include <sstream>
#include <string>
#include <iostream>
#include <cassert>

using std::list;
using std::ostream;
using std::cout;
using std::endl;
using std::string;
typedef string::size_type size_type;

GffFeature::GffFeature(const GffFeature& gff) : _c1(gff._c1), _c2(gff._c2), _c3(gff._c3), _end5(gff._end5), _end3(gff._end3),
  _c5(gff._c5), _strnd(gff._strnd), _frame(gff._frame), _c8(gff._c8) {
}

void
GffFeature::setCol(int col, const string& val) {
  switch(col) {
  case 1: _c1 = val; break;
  case 2: _c2 = val; break;
  case 8: _c8 = val; break;
  default:
    std::cerr<<"not used yet "<<col<<" "<<val<<endl;
  }
}

void
GffFeature::printLst(const list<GffFeature*>& lst) {
  list<GffFeature*>::const_iterator iter = lst.begin();
  list<GffFeature*>::const_iterator stop = lst.end();
  for(; iter != stop; ++iter) {
    const GffFeature& gff = **iter;
    cout<<gff<<endl;
  }
}

void
GffFeature::fillFromBuffer(const string& all) {
  unsigned idx = 0;
  size_type last = 0, cursor = 0;
  //while(cursor != string::npos ) {
  bool stop = false;
  while(!stop) {
    cursor = all.find('\t',last);
    if( cursor == string::npos ) {
      cursor = all.size();
      stop = true;
    }
    string col = all.substr(last,cursor-last);
    //cout<<"col: ["<<col<<"]"<<" idx: "<<idx<<" last: "<<last<<" curr: "<<cursor<<endl;
    if( idx == 0 ) {
      this->_c1 = col;
    } else if( idx == 1 ) {
      this->_c2 = col;
    } else if( idx == 2 ) {
      this->_c3 = col;
    } else if( idx == 3 ) {
      std::istringstream istrm(col.c_str());
      istrm>>this->_end5;
    } else if( idx == 4 ) {
      std::istringstream istrm(col.c_str());
      istrm>>this->_end3;
    } else if( idx == 5 ) {
      this->_c5 = col;
    } else if( idx == 6 ) {
      this->_strnd = col[0];
    } else if( idx == 7 ) {
      std::istringstream istrm(col.c_str());
      istrm>>this->_frame;
    } else if( idx == 8 ) {
      this->_c8 = col;
    }
    last = cursor+1;
    if( idx > 8 || stop) {
      break;
    }
    ++idx;
  }
  assert(idx==8);
}

ostream&
operator<<(ostream& os, const GffFeature& gff) {
  os<<gff._c1<<"\t"<<gff._c2<<"\t"<<gff._c3<<"\t";
  os<<gff.getEnd5()<<"\t"<<gff.getEnd3()<<"\t"<<gff._c5<<"\t";
  os<<gff._strnd<<"\t"<<gff._frame<<"\t"<<gff._c8;
  return os;
}

DataSetAbstractProduct* 
GffFeature::makeDataSet(const list<GffFeature*>& lst) {
  SeqLoc sl(dsu::eEither,-1,-1);
  // create a new list of genes
  GenePredDataSet* newGenePredDataSet = new GenePredDataSet(sl,-1,"none");
  int geneId = 0;
  string geneIdStr,lastGeneIdStr;
  string gffid = "none1";
  list<GffFeature*>::const_iterator iter = lst.begin();
  const list<GffFeature*>::const_iterator stop = lst.end();
  for(; iter != stop; ++iter) {
    const GffFeature* gff = *iter;
    const string exon_type = gff->getType();
    const int end5 = gff->getEnd5();
    const int end3 = gff->getEnd3();
    const double score = 1.0;
    const string geneIdStr = gff->getCol8();
    if( geneIdStr != lastGeneIdStr ) {
      ++geneId;
    }
    const char strnd = gff->getStrnd();
    const dsu::Strand_t parseStrnd = strnd == '-' ? dsu::eNeg : dsu::ePos;
    try {
      const SeqLoc theSl(parseStrnd,end5,end3);
      newGenePredDataSet->insert(gffid,geneId,theSl,exon_type,score);
    } catch(std::exception& e) {
      std::cerr<<e.what()<<std::endl;
      std::cerr<<"attempting to continue: ["<<*gff<<"]"<<std::endl;
    }
    lastGeneIdStr = geneIdStr;
  }
  return newGenePredDataSet;
}

