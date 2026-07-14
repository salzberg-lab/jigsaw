//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "SpliceSiteDataSet.h"
#include "SpliceSite.h"
#include "Options.h"
#include <stdexcept>
#include <assert.h>

using std::endl;

static bool mystrstri(const string& a, const char* knownStr);

SpliceSiteDataSet::SpliceSiteDataSet(const SeqLoc& sl, int geneId, const string& str) 
 : DataSetAbstractProduct(sl,geneId,str,AnnotationScore::defaultScore()) {
}

void SpliceSiteDataSet::print(std::ostream& os, const string& str) const {
  os<<str<<"ID: "<<_geneId<<" "<<getTypeStr()<<" "<<_source<<" ";
  SeqLoc::print(os);
  os<<endl;
  DataSetAbstractProduct::DataStorageType::const_iterator iter = begin();
  while( iter != end() ) {
    const Prediction& sl = **iter;
    sl.print(os,str+"\t");
    os<<endl;
    ++iter;
  }
}

void SpliceSiteDataSet::insert(const string& source, int end5, int end3, 
			       const string& splice_type, double score) {
  const Options& options = *(Options::getInstance());
  if( (options._subEnd5 != -1 && options._subEnd3 != -1) &&
      (end3 < options._subEnd5 || end5 > options._subEnd3)) {
    return;
  }
  dsu::Strand_t strnd = end5 > end3 ? dsu::eNeg : dsu::ePos;
  if(mystrstri(splice_type,"acc") && strnd == dsu::ePos ) {
    end5 = end3 = (end3+1);
  } else if(mystrstri(splice_type,"don") && strnd == dsu::eNeg) {
    end5 = end3 = (end5+1);
  } else if(mystrstri(splice_type,"don") && strnd == dsu::ePos) {
    end5 = end3 = (end5-1);
  } else if(mystrstri(splice_type,"acc") && strnd == dsu::eNeg) {
    end5 = end3 = (end3-1);
  } 
  SeqLoc sl(strnd,end5,end3);
  SpliceSite* sp = new SpliceSite(sl,score,-1,source,splice_type);
  DataSetAbstractProduct::insert(sp);
}

static const string geneStr = "One_Predicted_Splice_Site";
const string& SpliceSiteDataSet::getTypeStr() const {
  return geneStr;
}

bool mystrstri(const string& a, const char* knownStr) {
  string tmpl(a.length(),'\0');
  for(unsigned i = 0; i < a.size(); ++i) {
    tmpl[i] = tolower(a[i]);
  }
  return strstr(tmpl.c_str(),knownStr)?true:false;
}
