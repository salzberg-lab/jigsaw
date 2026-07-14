#include "parse_help.h"
#include <string>

using std::string;

static inline bool isWhiteSpace(char c) {
  return (c == ' ' || c == '\t');
}

string get_elem(const string& str, unsigned ith_elem) {
  const char* buff = str.c_str();
  unsigned cnt = 0, offset = 0;
  unsigned start=0,stop=0;
  bool prevIsSpace = true;
  unsigned elem = 0;
  const char* ptr = NULL;
  while(isWhiteSpace(*buff)) {
    ++buff;
    ++offset;
  }
  for(ptr = buff, elem = 0; 
      *ptr && *ptr != '\n' && cnt < ith_elem; 
      ++ptr, ++elem) {
    ++stop;
    //cout<<"ptr: "<<*ptr<<" cnt: "<<cnt<<" start: "<<start<<" stop: "<<stop<<endl;
    if(!isWhiteSpace(*ptr) && prevIsSpace) {
      start = elem;
      prevIsSpace = false;
    }
    else if(!isWhiteSpace(*ptr) && !prevIsSpace) {
      stop = elem;
      prevIsSpace = false;
    }
    if( isWhiteSpace(*ptr) && !prevIsSpace ) {
      stop = elem;
      ++cnt;
      prevIsSpace = true;
    } else if( isWhiteSpace(*ptr) ) {
      prevIsSpace = true;
    }
  }
  // hack
  if( ith_elem == 10 ) {
    ++stop;
  }
  start += offset;
  stop += offset;
  string val = str.substr(start,stop-start);
  //cout<<"i: "<<ith_elem<<" start: "<<start<<" stop: "<<stop<<" val: ["<<val<<"]"<<endl;
  //exit(1);
  return val;
}
