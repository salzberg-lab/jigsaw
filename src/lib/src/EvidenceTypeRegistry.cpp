#include "EvidenceTypeRegistry.h"
#include "SpliceSiteEvidence.h"
#include "HomologyEvidence.h"
#include "HomGeneEvidence.h"
#include "GenePredEvidence.h"
#include <iostream>

using std::endl;
using std::cerr;
using std::cout;
using std::string;
using namespace dsu; // what do with this crap? //ignoreCaseEqual


EvidenceType* 
EvidenceTypeRegistry::createNewEvidenceType(const string& evidenceType, const string& predictionType, 
					    const string& strId, int intId, const string& theSeq, const double wght1) {
  
  if(evidenceType.length() == 0 || predictionType.length() == 0 ) {
    cout<<"Warning: incomplete input: "<<endl;
    return NULL;
  }
  EvidenceType* result = NULL;
  if( ignoreCaseEqual(predictionType,"spliceprediction") ) {
    if( ignoreCaseEqual(evidenceType,"acc") ) {
      result = new SpliceSiteAcc(strId,intId,wght1);
    } else if( ignoreCaseEqual(evidenceType,"don") ) {
      result = new SpliceSiteDon(strId,intId,wght1);
    } else {
      cerr<<"spliceprediction Evidence type: ["<<evidenceType<<"] not supported, only acc and don"<<endl;
    }
  } else if( ignoreCaseEqual(predictionType,"geneprediction") ) {
    //double wght1 = 1.0;
    if( ignoreCaseEqual(evidenceType,"coding") ) {
      result = new GenePredEvidence(strId,intId,wght1);
    } else if( ignoreCaseEqual(evidenceType,"acc") ) {
      result = new GenePredEvidenceBndryAcc(strId,intId,wght1);
    } else if( ignoreCaseEqual(evidenceType,"don") ) {
      result = new GenePredEvidenceBndryDon(strId,intId,wght1);
    } else if( ignoreCaseEqual(evidenceType,"start") ) {
      result = new GenePredEvidenceBndryStart(strId,intId,wght1);
    } else if( ignoreCaseEqual(evidenceType,"overlap_start") ) {
      result = new GenePredEvidenceOverlapStart(strId,intId,wght1);
    } else if( ignoreCaseEqual(evidenceType,"stop") ) {
      result = new GenePredEvidenceBndryStop(strId,intId,wght1);
    } else if( ignoreCaseEqual(evidenceType,"intron") ) {
      result = new GenePredEvidenceIntron(strId,intId,wght1);
    } else {
      cerr<<"geneprediction Evidence type: ["<<evidenceType<<"] not supported, use: coding,acc,don,start,stop and intron "<<endl;
    }
  } else if( ignoreCaseEqual(predictionType,"homology") ) {
    const double wght2 = 1.0, wght5 = -1;
    if( ignoreCaseEqual(evidenceType,"coding") ) {
      result = new HomologyEvidence(strId,intId,wght1,wght2,wght5);
    } else if( ignoreCaseEqual(evidenceType,"acc") ) {
      result = new HomologyEvidenceBndryAcc(strId,intId,wght1,wght2,wght5,theSeq);
    } else if( ignoreCaseEqual(evidenceType,"don") ) {
      result = new HomologyEvidenceBndryDon(strId,intId,wght1,wght2,wght5,theSeq);
    } else if( ignoreCaseEqual(evidenceType,"start") ) {
      result = new HomologyEvidenceBndryStart(strId,intId,wght1,wght2,wght5,theSeq);
    } else if( ignoreCaseEqual(evidenceType,"stop") ) {
      result = new HomologyEvidenceBndryStop(strId,intId,wght1,wght2,wght5,theSeq);
    } else if( ignoreCaseEqual(evidenceType,"intron") ) {
      result = new HomologyEvidenceIntron(strId,intId,wght1,wght2,wght5,theSeq);
    } else {
      cerr<<"homology Evidence type: ["<<evidenceType<<"] not supported, use: coding,acc,don,start,stop and intron "<<endl;
    }
  } else if( ignoreCaseEqual(predictionType,"homgene") ) {
    //const double wght1 = 0.0;
    if( ignoreCaseEqual(evidenceType,"coding") ) {
      result = new HomGeneEvidence(strId,intId,wght1);
    } else if( ignoreCaseEqual(evidenceType,"acc") ) {
      result = new HomGeneEvidenceBndryAcc(strId,intId,wght1,theSeq);
    } else if( ignoreCaseEqual(evidenceType,"don") ) {
      result = new HomGeneEvidenceBndryDon(strId,intId,wght1,theSeq);
    } else if( ignoreCaseEqual(evidenceType,"start") ) {
      result = new HomGeneEvidenceBndryStart(strId,intId,wght1,theSeq);
    } else if( ignoreCaseEqual(evidenceType,"stop") ) {
      result = new HomGeneEvidenceBndryStop(strId,intId,wght1,theSeq);
    } else if( ignoreCaseEqual(evidenceType,"intron") ) {
      result = new HomGeneEvidenceIntron(strId,intId,wght1,theSeq);
    } else {
      cerr<<"homgene Evidence type: ["<<evidenceType<<"] not supported, use: coding,acc,don,start,stop and intron "<<endl;
    }
  } else if(ignoreCaseEqual(predictionType,"curation")) {
    // Quietly do nothing
  }  else {
    cerr<<"Don't recognize the "<<predictionType<<" data id evidenceType, skipping"<<endl;
  }
  return result;
}


EvidenceType* 
EvidenceTypeRegistry::createNewEvidenceTypeLinearBackwardCompatible(const string& evidenceType, const string& predictionType, double wght1,
					    const string& strId, int intId, const string& theSeq) {
    // determine data type by extension for now
  EvidenceType* result = NULL;
  if( ignoreCaseEqual(evidenceType,"spliceprediction") ) {
    if( ignoreCaseEqual(evidenceType,"acc") ) {
      result = new SpliceSiteAcc(strId,intId,1.0);
    } else if( ignoreCaseEqual(evidenceType,"acc") ) {
      result = new SpliceSiteDon(strId,intId,1.0);
    }
  } else if( ignoreCaseEqual(evidenceType,"geneprediction") ) {
    if( ignoreCaseEqual(evidenceType,"coding") ) {
      result = new GenePredEvidence(strId,intId,wght1);
    } else if( ignoreCaseEqual(evidenceType,"acc") ) {
      result = new GenePredEvidenceBndryAcc(strId,intId,wght1);
    } else if( ignoreCaseEqual(evidenceType,"don") ) {
      result = new GenePredEvidenceBndryDon(strId,intId,wght1);
    } else if( ignoreCaseEqual(evidenceType,"start") ) {
      result = new GenePredEvidenceBndryStart(strId,intId,wght1);
    } else if( ignoreCaseEqual(evidenceType,"stop") ) {
      result = new GenePredEvidenceBndryStop(strId,intId,wght1);
    }
  } else if( ignoreCaseEqual(evidenceType,"homology") ) {
    const double wght2 = 0.0;
    const double wght5 = -1;
    if( ignoreCaseEqual(evidenceType,"coding") ) {
      result = new HomologyEvidence(strId,intId,wght1,wght2,wght5);
    } else if( ignoreCaseEqual(evidenceType,"acc") ) {
      result = new HomologyEvidenceBndryAcc(strId,intId,wght1,wght2,wght5,theSeq);
    } else if( ignoreCaseEqual(evidenceType,"don") ) {
      result = new HomologyEvidenceBndryDon(strId,intId,wght1,wght2,wght5,theSeq);
    } else if( ignoreCaseEqual(evidenceType,"start") ) {
      result = new HomologyEvidenceBndryStart(strId,intId,wght1,wght2,wght5,theSeq);
    } else if( ignoreCaseEqual(evidenceType,"stop") ) {
      result = new HomologyEvidenceBndryStop(strId,intId,wght1,wght2,wght5,theSeq);
    }
  } 
  return result;
}
