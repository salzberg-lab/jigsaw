#include <stdexcept>
#include "FormatRegistry.h"
#include "GenePredDataSetConcreteFactory.h"
#include "SpliceSiteDataSetConcreteFactory.h"
#include "HomologyDataSetConcreteFactory.h"
#include "Gff3Factory.h"
#include "GtfFactory.h"
#include "GffFactory.h"
#include "BtabFactoryIntron.h"
#include "BtabFactory.h"
#include "PhatFactory.h"
#include "GlimmerMFactory.h"
#include "FgeneshFactory.h"
#include "GenemarkHmmFactory.h"
#include "GenscanFactory.h"
#include "SnapFactory.h"
#include "dsu.h"

using std::cerr;
using std::endl;
using namespace dsu; // ignoreCaseEqual

/**
 ** consider making these singletons
 ** and then add a destructor to the FormatRegistry
 ** that cleans up nicely on program exit
 **
 ** these objects take up so little memory that
 ** for now I just create them all (even if not
 ** used) becuase the memory waste is minimal
 **/
static GtfFactory gtfFactory;
static SnapFactory snapFactory;
static GenscanFactory genscanFactory;
static GenemarkHmmFactory genemarkHmmFactory;
static FgeneshFactory fgeneshFactory;
static GlimmerMFactory glimmermFactory;
static PhatFactory phatFactory;
static Gff3Factory gff3Factory;
static GffFactory gffFactory;
static BtabFactoryIntron btabFactoryIntron;
static BtabFactory btabFactory;
static GenePredDataSetConcreteFactory gpFactory;
static SpliceSiteDataSetConcreteFactory spliceFactory;
static HomologyDataSetConcreteFactory napFactory;

DataSetAbstractFactory* FormatRegistry::getFactory(const string& format, const string& predtype) {
  if(ignoreCaseEqual(format, "default")) {
    if(ignoreCaseEqual(predtype, "spliceprediction")) {
      return &spliceFactory;
    } else if(ignoreCaseEqual(predtype, "geneprediction")) {
      return &gpFactory;
    } else if(ignoreCaseEqual(predtype, "homology")) {
      return &napFactory;
    } else if(ignoreCaseEqual(predtype,"curation")) {
      return &gpFactory;
    } else {
      return NULL;
    }
  } else if(ignoreCaseEqual(format, "gff3")) {
    return &gff3Factory;
  } else if(ignoreCaseEqual(format, "gtf")) {
    return &gtfFactory;
  } else if(ignoreCaseEqual(format, "gff")) {
    return &gffFactory;
  } else if(ignoreCaseEqual(format, "glimmerm")) {
    return &glimmermFactory;
  } else if(ignoreCaseEqual(format, "btab_i")) {
    return &btabFactoryIntron;
  } else if(ignoreCaseEqual(format, "btab")) {
    return &btabFactory;
  } else if(ignoreCaseEqual(format, "phat")) {
    return &phatFactory;
  } else if(ignoreCaseEqual(format, "fgenesh")) {
    return &fgeneshFactory;
  } else if(ignoreCaseEqual(format, "genemarkhmm")) {
    return &genemarkHmmFactory;
  } else if(ignoreCaseEqual(format, "genscan")) {
    return &genscanFactory;
  } else if(ignoreCaseEqual(format, "snap")) {
    return &snapFactory;
  } else {
    cerr<<"Warning, don't recognize the format type: ["<<format<<"]"<<endl;
    return NULL;
  }
  return NULL;
}

void FormatRegistry::printFormats(std::ostream& os) {
  const char* genePredFormats[] = {
    "gtf", "gff3", "gff","glimmerm","phat","fgenesh","genemarkhmm","genscan","snap","default",NULL
  };
  const char* homologyFormats[] = {
    "gtf", "gff3", "btab", "btab_i", "gff","default",NULL
  };
  const char* spliceFormats[] = {
    "gff","default",NULL
  };
  os<<endl<<"\t\tRecognized file formats"<<endl<<endl;
  os<<"\"geneprediction\" formats:";
  for(unsigned i = 0; genePredFormats[i]; ++i) {
    os<<" "<<genePredFormats[i];
  }
  os<<endl;
  os<<"\"homology\" formats:";
  for(unsigned i = 0; homologyFormats[i]; ++i) {
    os<<" "<<homologyFormats[i];
  }
  os<<endl;
  os<<"\"spliceprediction\" formats:";
  for(unsigned i = 0; spliceFormats[i]; ++i) {
    os<<" "<<spliceFormats[i];
  }
  os<<endl<<endl<<endl;
}
