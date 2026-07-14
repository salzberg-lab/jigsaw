#include "MatrixIdx.h"
#include <assert.h>

static int
bIdx(dsu::Strand_t strnd, int frame, dsu::side_t side, dsu::gbnd_t type) {
  const int nfrm = ((int)(strnd)*3+frame);
  const int xdx = nfrm+((int)side*6);
  const int ddx = xdx +((int)type*12);
  assert(nfrm>=0&&nfrm<=5);
  assert(xdx>=0&&xdx<=11);
  assert(ddx>=0 && ddx<24);
  return ddx;
}

const MatrixIdx& 
MatrixIdx::getIdx() const {
  return *this;
}

int 
MatrixIdx::getIntIdx() const {
  return bIdx(getStrnd(),getFrame(),getSide(),getBnd());
}

void
MatrixIdx::init(const MatrixIdx& mi) {
  _frame = mi.getFrame();
  _strnd = mi.getStrnd();
  _side = mi.getSide();
  _bnd = mi.getBnd();
}

static const string plusStr = "+";
static const string minusStr = "-";
static const string leftStr = "left";
static const string rghtStr = "rght";
static const string midlStr = "midl";
static const string tendStr = "tend";

std::ostream&
operator<<(std::ostream& os, const MatrixIdx& mi) {
  const string& strnd=mi.getStrnd()==dsu::ePos?plusStr:minusStr;
  const string& side=mi.getSide()==dsu::eLeft?leftStr:rghtStr;
  const string& bnd=mi.getBnd()==dsu::eEndOfGene?tendStr:midlStr;
  os<<mi.getFrame()<<" "<<strnd<<" "<<side<<" "<<bnd;
  return os;
}
