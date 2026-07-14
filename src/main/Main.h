//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#ifndef _MAIN_H
#define _MAIN_H

#include "CustomScore.h"

class InputData;
class Options;

/**
 * Interface to the "main" function
 * runs combiner
 */
class Main {
 public:
  /**
   * @param Options - command line options 
   * @param CustomScore* - user defined scoring function
   * @param InputData& - all the evidence data stored here
   * @param const string& - genomic sequence
   * @param const string& - fasta header
   * @param const string& - program string id output in gff file
   */
  Main() { }
  void
  run(Options&, CustomScore*, InputData&, const string&,const string&, const string&);

 private:
  static void mgenes(const InputData&,const string&, std::ostream&, bool, CustomScore*,const Options&, const string&);
};

#endif // _MAIN_H
