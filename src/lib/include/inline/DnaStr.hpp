//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
inline bool DnaStr::isStop(int bp) const {
    if(length() < static_cast<unsigned>(bp+2) || bp < 1)
        return false;
    char a, b, c;
    if(getStrnd() == dsu::ePos) {
      a = (*this)[bp];
      b = (*this)[bp+1];
      c = (*this)[bp+2];
    } else {
      c = (*this)[bp];
      b = (*this)[bp+1];
      a = (*this)[bp+2];
    }
    bool firstIsU = isU(a);
    bool secIsA = isA(b);
    return  
      (firstIsU && secIsA && isR(c)) ||
      (firstIsU && isR(b) && isA(c));
  }

inline bool DnaStr::isC(char c) {
  switch(c) {
  case 'C':
  case 'c':
    return true;
  default:
    return false;
  }
}

inline bool DnaStr::isA(char c) {
  switch(c) {
  case 'A':
  case 'a':
    return true;
  default:
    return false;
  }
}

inline bool DnaStr::isR(char c) {
  char lc = tolower(c);
  return lc=='r'||lc=='a'||lc=='g'?true:false;
}

inline bool DnaStr::isU(char c) {
  return c=='T'||c=='t'||c=='u'||c=='U'?true:false;
}

inline bool DnaStr::isG(char c) {
  switch(c) {
  case 'g':
  case 'G':
    return true;
  default:
    return false;
  }
}

inline bool DnaStr::isStart(int bp) const {
  if(length() < static_cast<unsigned>(bp+2) || bp < 1)
    return false;
  if(getStrnd()==dsu::ePos) {
    return isA((*this)[bp]) && isU((*this)[bp+1]) && isG((*this)[bp+2]);
  } else {
    return isA((*this)[bp+2]) && isU((*this)[bp+1]) && isG((*this)[bp]);
  }
}

inline bool DnaStr::isStart(const char codon[3]) const {
  if(getStrnd() == dsu::ePos)
    return isA(codon[0]) && isU(codon[1]) && isG(codon[2]);
  else
    return isA(codon[2]) && isU(codon[1]) && isG(codon[0]);
}

inline bool DnaStr::isStop(const char codon[3]) const {
  const char* tmpPtr = codon;
  char a,b,c;
  if(getStrnd()==dsu::ePos) {
    a = *tmpPtr; b = *(++tmpPtr); c = *(++tmpPtr);
  } else {
    c = *tmpPtr; b = *(++tmpPtr); a = *(++tmpPtr);
  }
  bool firstIsU = isU(a);
  bool secIsA = isA(b);

  return 
      (firstIsU && secIsA && isR(c)) ||
      (firstIsU && isR(b) && isA(c));
}
