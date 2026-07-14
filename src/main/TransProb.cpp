#include "TransProb.h"
#include "Options.h"
#include "genemodels.h"
#include <assert.h>

TransProb TransProb::_instance;
static const int cDIM=6;

TransProb::TransProb() : _vec(2*cDIM) {
  for(unsigned i = 0; i < 2*cDIM; ++i) {
    _vec[i].resize(2*cDIM,0);
  }
  init();
}

const TransProb&
TransProb::getInstance() {
  return TransProb::_instance;
}

#if 0
void 
TransProb::write() const {
  const Options& options = (*Options::getInstancee());
  string pfile = "trans.txt";
  if( options._dtFile != "" ) {
    pfile = options._dtFile + string("/param.txt");
  }
  ofstream ofs(file.c_str());
  if( !ofs ) {
    ostringstream os;
    os<<"unable to create file: "<<file<<endl;
  }
  const int cStateTypes = 6;
  for(unsigned i = 0; i < cStateTypes; ++i) {
    for(unsigned j = 0; j < cStateTypes; ++j) {
      for(unsigned is = 0; is < 2; ++is) {
	dsu::Strand_t ist = (dsu::Strand_t)is;
	for(unsigned js = 0; js < 2; ++js) {
	  dsu::Strand_t jst = (dsu::Strand_t)js;
	  const double tprob = getProb(i,ist,jst);
	  if( tprob > 0 ) {
	    const string& pstate = genemodels::idx2str(i);
	    const string& pstrnd = dsu::strnd2str(ist);
	    const string& cstate = genemodels::idx2str(i);
	    const string& cstrnd = dsu::strnd2str(jst);
	    ofs<<pstate<<" "<<pstrnd<<" "<<cstate<<" "<<cstrnd<<" "<<tprob<<endl;
	  }
	}
      }
    }
  }
}
#endif

void
TransProb::init() {
  const dsu::Strand_t pos=dsu::ePos;
  const dsu::Strand_t neg=dsu::eNeg;
  const int itrn = Intron::getInstance().getScoreIdx();
  const int igre = Igr::getInstance().getScoreIdx();
  const int sngl = Single::getInstance().getScoreIdx();
  const int init = Initial::getInstance().getScoreIdx();
  const int term = Terminal::getInstance().getScoreIdx();
  const int intr = Internal::getInstance().getScoreIdx();

  setProb(sngl,pos,sngl,pos) = 1;
  setProb(sngl,pos,sngl,neg) = 1;
  setProb(sngl,neg,sngl,pos) = 1;
  setProb(sngl,neg,sngl,neg) = 1;
  setProb(sngl,pos,init,pos) = 1;
  setProb(sngl,neg,init,pos) = 1;
  setProb(sngl,pos,term,neg) = 1;
  setProb(sngl,neg,term,neg) = 1;
  setProb(sngl,neg,igre,neg) = 1;
  setProb(sngl,neg,igre,pos) = 1;
  setProb(sngl,pos,igre,pos) = 1;
  setProb(sngl,pos,igre,neg) = 1;

  setProb(init,neg,init,pos) = 1;
  setProb(init,neg,term,neg) = 1;
  setProb(init,neg,sngl,neg) = 1;
  setProb(init,neg,sngl,pos) = 1;
  setProb(init,pos,intr,pos) = 1;
  setProb(init,pos,term,pos) = 1;
  setProb(init,pos,itrn,pos) = 1;
  setProb(init,neg,igre,pos) = 1;
  setProb(init,neg,igre,neg) = 1;

  setProb(term,neg,intr,neg) = 1;
  setProb(term,neg,init,neg) = 1;
  setProb(term,pos,init,pos) = 1;
  setProb(term,pos,term,neg) = 1;
  setProb(term,pos,sngl,neg) = 1;
  setProb(term,pos,sngl,pos) = 1;
  setProb(term,pos,igre,pos) = 1;
  setProb(term,pos,igre,neg) = 1;
  setProb(term,neg,itrn,neg) = 1;

  setProb(intr,pos,intr,pos) = 1;
  setProb(intr,pos,term,pos) = 1;
  setProb(intr,neg,intr,neg) = 1;
  setProb(intr,neg,init,neg) = 1;
  setProb(intr,neg,itrn,neg) = 1;
  setProb(intr,pos,itrn,pos) = 1;
  setProb(intr,pos,intr,pos) = 1;
  setProb(intr,pos,term,pos) = 1;
  setProb(intr,neg,intr,neg) = 1;
  setProb(intr,neg,init,neg) = 1;

  setProb(itrn,pos,intr,pos) = 1;
  setProb(itrn,pos,term,pos) = 1;
  setProb(itrn,neg,init,neg) = 1;
  setProb(itrn,neg,intr,neg) = 1;

  setProb(igre,pos,sngl,pos) = 1;
  setProb(igre,neg,sngl,pos) = 1;
  setProb(igre,pos,sngl,neg) = 1;
  setProb(igre,neg,sngl,neg) = 1;
  setProb(igre,pos,init,pos) = 1;
  setProb(igre,neg,init,pos) = 1;
  setProb(igre,neg,term,neg) = 1;
  setProb(igre,pos,term,neg) = 1;
}

double&
TransProb::setProb(int ltype, dsu::Strand_t lstrnd, int rtype, dsu::Strand_t rstrnd) {
  return _vec[ ltype+(int)lstrnd*cDIM][rtype+(int)rstrnd*cDIM];
}
double
TransProb::getProb(int ltype, dsu::Strand_t lstrnd, int rtype, dsu::Strand_t rstrnd) const {
  assert(ltype != Igr::getInstance().getScoreIdx() || ltype != Intron::getInstance().getScoreIdx());
  assert(rtype != Igr::getInstance().getScoreIdx() || rtype != Intron::getInstance().getScoreIdx());
  if( ltype == -1 ) return 1.0;
  assert(rtype>=0&& rtype<cDIM);
  assert(ltype>=0&& ltype<cDIM);
  return _vec[ ltype+(int)lstrnd*cDIM][rtype+(int)rstrnd*cDIM];
}
