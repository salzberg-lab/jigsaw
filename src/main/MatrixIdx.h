#ifndef _MATRIXIDX_H
#define _MATRIXIDX_H

#include "dsu.h"
#include <iostream>

class MatrixIdx {
 public:
  MatrixIdx() : _frame(-1), _strnd(dsu::eEither), _side(dsu::eSide), _bnd(dsu::eGbnd) { }
  MatrixIdx(const MatrixIdx& mi) : _frame(mi._frame), _strnd(mi._strnd), _side(mi._side), _bnd(mi._bnd)
    { }
  MatrixIdx(int frame, dsu::Strand_t strnd, dsu::side_t side, dsu::gbnd_t bnd) :
    _frame(frame), _strnd(strnd), _side(side), _bnd(bnd) 
    { }

  bool operator==(const MatrixIdx& mi) const {
    return _frame==mi._frame&&_strnd==mi._strnd&&_side==mi._side&&_bnd==mi._bnd;
  }
  void init(const MatrixIdx&);

  int getFrame() const { return _frame; }
  dsu::Strand_t getStrnd() const { return _strnd; }
  dsu::side_t getSide() const { return _side; }
  dsu::gbnd_t getBnd() const { return _bnd; }
  
  const MatrixIdx& getIdx() const;

  int getIntIdx() const;

  friend std::ostream&
  operator<<(std::ostream&,const MatrixIdx&);
 

 private:
  int _frame;
  dsu::Strand_t _strnd;
  dsu::side_t _side;
  dsu::gbnd_t _bnd;
};

#endif // _MATRIXIDX_H
