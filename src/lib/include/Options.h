#ifndef _OPTIONS_H
#define _OPTIONS_H
//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include <stdexcept>
#include <strstream>

using std::string;

/**
 * command line interface
 */
struct Options {
  /** allow individual programs (defined by the main)
   ** access private variables to customize program
   ** behavior
   **/
  friend int main(int,char**);

  /**
   * @param argc - number command line arguments
   * @param argv - command line arguments
   * parameters gotton from c's main function.
   */
  bool isExonDistr() const { return _exonDistrFile.size()>0; }
  bool isLinComb() const { return _lnCmb; }
  bool isLinComb2() const { return _lnCmb2; }
  bool _checkSize;
  int _cutOff;
  bool quit, _forceModels;
  string m_fastaFile, m_evidenceListFile, _outputGeneModelFile, _matrixFile, _vecFile, _dtFile, _exonDistrFile;
  bool _train;
  int _subEnd5, _subEnd3, _trainOffSet;
  string _hmmFile;

  // for backward compatibility with previousl combiners that do not use intron model
  void setIntronFalse() { _useIntron = false; }
  
  bool useIntron() const { return _useIntron; }
  bool isVerbose() const { return _verbose; }
  bool doTraining() const { return _train; }
  static void createInstance(int argc, char** argv, const string& id) {
    _instance = new Options(argc,argv,id);
  }
  static Options* getInstance() {
    return _instance; 
  }

protected:
  Options(int argc, char** argv, const string&); 

private:
  void parse(int argc, char** argv, const string&); 
  static void versionStatement(std::ostream&, const string&); 

private:
  bool _lnCmb,_lnCmb2;
  bool _verbose;
  bool _useIntron;
  static Options* _instance;
};

#endif //_OPTIONS_H
