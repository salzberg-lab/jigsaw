//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "dtree.h"
#include "moc1.h"
#include "DataDescrip.h"
#include "HomologyEvidence.h"
#include "GenePredEvidence.h"
#include "FeatVec.h"
#include <fstream>
#include <list>

/** this old tree testing code....**/

extern int no_of_dimensions;
extern int no_of_categories; 
int main(int argc, char** argv) {
  if( !argv[1] ) return -1;
  const int NT = 10;
  DTree dt(argv[1],NT);
  no_of_categories = 2;

  DataDescrip descrip;
  descrip.add( new GenePredEvidence("geneprediction4",3,1.0));
  descrip.add( new GenePredEvidence("geneprediction5",4,1.0));
  descrip.add( new GenePredEvidence("geneprediction6",5,1.0));
  descrip.add( new HomologyEvidence("homology1",0,1.0,1.0,1.0));
  descrip.add( new HomologyEvidence("homology2",1,1.0,1.0,1.0));
  descrip.add( new HomologyEvidence("homology3",2,1.0,1.0,1.0));
  
  if( !argv[2] ) return -1;
  ifstream ifs(argv[2]);
  if(!ifs) {
    cerr<<"Unable to open ["<<argv[2]<<"]\n";
    return -1;
  }
  FeatVec* vec = NULL;
  list<FeatVec*> lst;
  bool isFirst = false;
  while( (vec=descrip.read(ifs)) ) {
    lst.push_back(vec);
    double val;
    ifs>>val;
    if(isFirst) {
      no_of_dimensions = vec->getNumDim();
      isFirst = true;
    }
    POINT point;
    point.dimension = vec->toFloatArr();
    modify(&point,no_of_dimensions,val,dt.getRoot(),NT);
  }
  list<FeatVec*>::const_iterator iter = lst.begin(), stop = lst.end();
  int cnt = 0;
  while(iter != stop) {
    const FeatVec* vec = *iter;
    vec->val_print(cout);
    cout<<endl;
    POINT point;
    point.dimension = vec->toFloatArr();
    classify2(&point,no_of_dimensions,dt.getRoot(),NT);
    delete [] point.dimension;
    cout<<"val: "<<point.prob[0]<<endl;
    ++iter;
    ++cnt;
    
  }
  return 0;
}
