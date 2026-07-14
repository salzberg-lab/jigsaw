//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef DSU_H
#define DSU_H

#include <string>
#include <iostream>
#include <cstring>

using std::string;

namespace dsu {

  struct myeqstr {
    bool operator()(const char* s1, const char* s2) const {
      return strcmp(s1,s2)==0;
    }
  };

 void readFastaFile(string&, const string&);
 void readFastaHeader(string&, const string&);
 inline char comp(char);

 enum Strand_t { ePos = 0, eNeg, eEither };
 enum side_t { eLeft=0, eRght, eSide };
 enum gbnd_t { eEndOfGene=0, eMdlOfGene, eGbnd };

 Strand_t str2strnd(const string&); 
 Strand_t coord2strnd(int,int);

 const string& strnd2str(Strand_t);

 bool ignoreCaseEqual(const string&, const string&);

}

#include "inline/dsu.hpp"

#endif // DSU_H
