//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
#include "dtree.h"
#include "DataDescrip.h"
#include "HomologyEvidence.h"
#include "GenePredEvidence.h"
#include "moc1.h"
#include <stdexcept>
#include <sstream>
#include <list>
#if __GNUC__ < 4 || ( __GNUC__ == 4 && __GNUC_MINOR__ < 3 )
 #include <fstream.h>
#else 
 #include <fstream>
#endif

#if _RUNTIME_DEBUG >= 1
extern bool MYDEBUG;
#endif
//#define USE_DIST
using std::cerr;
using std::endl;
using std::cout;

int no_of_categories = 2;

static void tryread(const char* fstr, struct tree_node** root, int nt, const DataDescrip& descrip, int nv, int num_dim) {
  ifstream ifs(fstr);
  init_prob(root,nt);
  if(!ifs) {
    cerr<<"Unable to open ["<<fstr<<"]\n";
    return;
  }
  vector<FeatVec*> vec(nv);
  for(unsigned i = 0; i < vec.size(); ++i) {
    vec[i] = NULL;
  }
  bool doRead = true;
  unsigned debugCnt = 0;
  while(doRead) {
    for(unsigned i = 0; i < vec.size(); ++i) {
      if( !(vec[i] = descrip.read(ifs)) ) {
	doRead = false;
	if( i != 0 ) {
	  cerr<<"the expected feature vector size is not matching what we see in: "<<fstr<<endl;
	  cerr<<"i: "<<i<<" "<<vec.size()<<endl;
	  exit(1);
	}
	break;
      }
    }
    if( !doRead ) {
      break;
    }
    double val;
#ifdef USE_DIST
    int dist1,dist2;
    ifs>>dist1>>dist2>>val;
#else
    ifs>>val;
#endif
    POINT point;
    point.dimension = new float[num_dim];
    for(unsigned i = 0; i < vec.size(); ++i) {
      FeatVec* fv = vec[i];
      point.dimension = fv->toFloatArr(point.dimension,i*descrip.getNumDim());
      delete fv;
    }
#ifdef USE_DIST
    point.dimension[num_dim-2] = dist1;
    point.dimension[num_dim-1] = dist2;
#endif
#if _RUNTIME_DEBUG >= 4
    if(MYDEBUG) {
      int stop = 0;
      for(signed q = 0; q < num_dim; ++q) {
	cout<<point.dimension[q];
	if(q < num_dim-1)
	  cout<<" ";
      }
#ifdef USE_DIST
      cout<<dist1<<" "<<dist2<<" "<<val<<endl;
#endif
    }
#endif
    modify(&point,num_dim,val,root,nt);
    ++debugCnt;
  }
}


DTree::DTree(const string& fileName, 
	     const string& vfile, const DataDescrip& descrip, int nv ) : _root(NULL), _numDt(0), _numClus(0) {
  // 2 refers to number of nt's seperating the 3 feature vectors
#ifdef USE_DIST
  int num_dim = descrip.getNumDim()*nv + 2;
#else
  int num_dim = descrip.getNumDim()*nv;
#endif
  no_of_categories = 2;
  int nclus;
  list <tree_node*> lst;
  unsigned i = 0;
  while(1) {
    std::ostringstream strm;
    strm<<fileName<<"."<<(i+1)<<std::ends;
    try 
      {
	tree_node* tn = read_tree(strm.str().c_str(),num_dim,nclus);
	if(!tn) {
	  break;
	}
	lst.push_back(tn);
	++i;
      } 
    catch(std::runtime_error& msg) 
      {
	std::cerr<<"failed to read the following decision tree file: "<<strm.str()<<std::endl;
	std::cerr<<"I will try to move on but this should be fixed. "<<strm.str()<<std::endl;
	std::cerr<<"mesg: "<<msg.what()<<std::endl;
	break;
      }
  }
  _numDt = lst.size();
  _root = new tree_node * [ getNumDt() ];
  if( getNumDt() < 1 ) {
    cerr<<"Unable to get training data from specified directory: ["<<fileName<<"]"<<endl;
    cerr<<"Not able to continue.\n";
    exit(1);
  }
  list<tree_node*>::iterator iter = lst.begin(),stop = lst.end();
  unsigned cnt = 0;
  while(iter != stop) {
    _root[cnt] = *iter;
    ++iter;
    ++cnt;
  }
  
  tryread(vfile.c_str(),getRoot(),getNumDt(),descrip,nv, num_dim);
}

unsigned DTree::getNumDt() const {
  return _numDt;
}

double DTree::eval(const vector<FeatVec*>& vec) const {
  int no_of_dimensions = vec[0]->getNumDim()*vec.size();
  POINT point;
  point.dimension = new float[no_of_dimensions];
  for(unsigned i = 0; i < vec.size(); ++i) {
    const FeatVec* fv = vec[i];
    point.dimension = fv->toFloatArr(point.dimension,i*fv->getNumDim());
  }
#if _RUNTIME_DEBUG >= 1
  if(MYDEBUG) {
    for(unsigned l = 0; l < vec.size(); ++l ) {
      cout<<"[";
      vec[l]->val_print(cout);
      cout<<"]";
    }
    cout<<endl;
    for(int l = 0; l < no_of_dimensions; ++l) {
      cout<<point.dimension[l];
      if( l < no_of_dimensions-1 )
	cout<<" ";
    }
    cout<<endl;
  }
#endif
  classify2(&point,no_of_dimensions,_root,getNumDt());
  delete [] point.dimension;
#if _RUNTIME_DEBUG >= 1
  if(MYDEBUG) {
    cout<<"DT prediction: "<<point.prob[0]<<endl;
  }
#endif
 
  return point.prob[0];
}

int DTree::getClusterId(const vector<FeatVec*>& vec, const vector<double>& dvec) const {
  POINT point;
#ifdef USE_DIST
  int no_of_dimensions = vec[0]->getNumDim()*vec.size() + dvec.size();
#else
  int no_of_dimensions = vec[0]->getNumDim()*vec.size();
#endif
  point.dimension = new float[no_of_dimensions];
  for(unsigned i = 0; i < vec.size(); ++i) {
    const FeatVec* fv = vec[i];
    point.dimension = fv->toFloatArr(point.dimension,i*fv->getNumDim());
  }
#ifdef USE_DIST
  point.dimension[no_of_dimensions-2] = dvec[0];
  point.dimension[no_of_dimensions-1] = dvec[1];
#endif
#if _RUNTIME_DEBUG >= 4
  if(MYDEBUG) {
    for(unsigned l = 0; l < vec.size(); ++l ) {
      cout<<"[";
      vec[l]->val_print(cout);
      cout<<"]";
    }
    cout<<endl;
    for(int l = 0; l < no_of_dimensions; ++l) {
      cout<<point.dimension[l];
      if( l < no_of_dimensions-1 )
	cout<<" ";
    }
    cout<<endl;
  }
#endif
  int label = assignLabel(point.dimension,no_of_dimensions,_root[0]);
  assert(label != -1);
  delete [] point.dimension;
  return label;
}
