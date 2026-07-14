#ifndef _GFFFEATURE_H
#define _GFFFEATURE_H

#include <string>
#include <list>
#include <ostream>

class DataSetAbstractProduct;

class GffFeature {
 public:
  GffFeature() { }
  GffFeature(const GffFeature&);
  void fillFromBuffer(const std::string&);
  int getEnd5() const { return _end5; }
  int getEnd3() const { return _end3; }
  void setCol(int,const std::string&);
  void setEnd5(int end5) { _end5 = end5; }
  void setEnd3(int end3) { _end3 = end3; }

  const std::string& getType() const { return _c3; }
  char getStrnd() const { return _strnd; }
  int getFrame() const { return _frame; }
  const std::string& getCol3() const { return _c3; }
  const std::string& getCol8() const { return _c8; }
  void setCol5(const std::string& str) { _c5 = str; }
  void setCol8(const std::string& str) { _c8 = str; }
  void addComment(const std::string& str) { _c8 += str; }
  
  const std::string& getScore() const { return _c5; }
  static  DataSetAbstractProduct* makeDataSet(const std::list<GffFeature*>&);
  static void
    printLst(const std::list<GffFeature*>&);

  friend std::ostream& operator<<(std::ostream&, const GffFeature&);

 private:
  std::string _c1, _c2, _c3;
  int _end5, _end3;
  std::string _c5;
  char _strnd;
  int _frame;
  std::string _c8;
};


#endif // _GFFFEATURE_H
