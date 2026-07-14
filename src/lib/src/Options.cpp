//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "Options.h"
#include "FormatRegistry.h"
#include <iostream>
#include <unistd.h>
#include <sstream>

using std::cerr;
using std::cout;
using std::endl;
using std::runtime_error;

Options* Options::_instance = NULL;

Options::Options(int argc, char** argv, const string& id) : 
		_checkSize(false),_cutOff(-1), quit(false), _forceModels(false), _train(false), 
		_subEnd5(-1), _subEnd3(-1), _trainOffSet(250), _lnCmb(false), _lnCmb2(false), _verbose(false), 
		_useIntron(true) { 
  parse(argc,argv,id); 
}

void Options::parse(int argc, char** argv, const string& id) {
	int opt;
	while( (opt=getopt(argc,argv,"x:pSLlVvphe:f:t:m:c:d:H:q:b:n:")) != -1) {
		switch(opt) {
		case 'p':
                  _checkSize=true;
                  break;
		case 'S':
		  FormatRegistry::printFormats(cout);
		  exit(0);
                case 'H':
                  _hmmFile = optarg;
                  cout<<"Using hmm file: "<<_hmmFile<<endl;
                  break;
		case 'L':
		  _lnCmb = _lnCmb2 = true;
		  cout<<"Using the linear combination of weights\n";
		  break;
		case 'l':
		  _lnCmb = true;
		  cout<<"Using the linear combination of weights\n";
		  break;
		case 'd':
		  _dtFile = optarg;
		  cout<<"Using decision tree in: "<<_dtFile<<endl;
		  break;
		case 'V':
		  _verbose = true;
		  break;
		case 'v':
		  versionStatement(cout,id);
		  exit(0);
		case 't': 
		  _train = true;
		  _vecFile = optarg;
		  cout<<"Print Feature Vectors to ["<<_vecFile<<"]"<<endl;
		  break;
		case 'f':
		  assert(optarg);
		  m_fastaFile = optarg;
		  cout<<"Reading genomic sequence from fasta file: "<<optarg<<endl;
		  break;
		case 'e':
		  assert(optarg);
		  m_evidenceListFile = optarg;
		  cout<<"Combining evidence specified in file: "<<optarg<<endl;
		  break;
		case 'm':
		  _outputGeneModelFile = optarg;
		  cout<<"Output gene models to file: "<<_outputGeneModelFile<<endl;
		  break;
		case 'x':
		  _matrixFile = optarg;
		  cout<<"Printing scored evidence matrix to file: "<<_matrixFile<<endl;
		  break;
		case 'c':
		  {
		    std::istringstream istrm(optarg);
		    istrm>>_cutOff;
		    cout<<"Using a connection cutoff set to: "<<_cutOff<<endl;
		  }
		  break;
#if 0
		case 'p':
		  {
		    std::istringstream istrm(optarg);
		    if(!(istrm>>_subEnd5>>_subEnd3)) {
		      throw runtime_error("-p option expects two numbers, got: "+string(optarg));
		    }
		    if( _subEnd5 > _subEnd3 ) {
		      throw runtime_error("-p option failed\n5' coordinate must be less than 3' coordinate: "+string(optarg));
		    }
		    cout<<"use subset of genomic sequence: "<<_subEnd5<<" "<<_subEnd3<<endl;
		    break;
		  }
#endif
		case 'b':
		  {
		    std::istringstream istrm(optarg);
		    if(!(istrm>>_subEnd5)) {
		      throw runtime_error("-b option expects one number, got: "+string(optarg));
		    }
		    cout<<"begin genomic sequence pos : "<<_subEnd5<<endl;
		    break;
		  }

		case 'n':
		  {
		    std::istringstream istrm(optarg);
		    if(!(istrm>>_subEnd3)) {
		      throw runtime_error("-n option one number, got: "+string(optarg));
		    }
		    cout<<"end  genomic sequence pos : "<<_subEnd3<<endl;
		    break;
		  }
		case 'q':
		  {
		    std::istringstream istrm(optarg);
		    if(!(istrm>>_trainOffSet)) {
		      throw runtime_error("-q option expects a single number, got: "+string(optarg));
		    }
		    cout<<"Using a window around input region of +-: "<<_trainOffSet<<" bases"<<endl;
		    break;
		  }
		case 'h':
		  quit = true;
		  break;
		default:
		  cerr<<"Unhandeled case: ["<<(char)opt<<"]\n";
		  quit = true;
		  throw runtime_error("unable to parse the command line arguments sorry");
		  break;
		}
	}
	for(; optind < argc; ++optind) {
	  cout<<"Unhandeled argument: "<<argv[optind]<<endl;
	}
  }

void Options::versionStatement(std::ostream& os, const string& id) {
  os<<id<<endl;
}
