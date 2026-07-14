#ifndef _PARAMFILE_H
#define _PARAMFILE_H

#include <string>
#include <list>

class ParamFile {
 public:
  ParamFile(const std::string&);

  static ParamFile* getInstance() {
    return _instance;
  }
  static void createInstance(const std::string& file) {
    _instance = new ParamFile(file);
  }

  void load(const std::string&);

  int getMaxSeqLen() const { return _maxSeqLen; }
  int getOverlapLen() const { return _overlapLen; }
  int getInternalExonLengthPenalty() const { return _internalExonLenPenalty; }
  int getIntronLengthPenalty() const { return _intronLenPenalty; }
  int getIgrLengthPenalty() const { return _igrLenPenalty; }
  int getMinIntronLength() const { return _minIntronLen; }
  int getMinIgrLength() const { return _minIgrLen; }
  int getAlignmentConnectionCutoff() const { return _minAlignmentConnectionCutoff; }

  double getIntronPenalty() const { return _intronPenalty; }
  double getIgrPenalty() const { return _igrPenalty; }
  double getInternalExonPenalty() const { return _internalExonPenalty; }
  
  const std::list<std::string>& getDonorList() const { return _donors; }
  const std::list<std::string>& getAcceptorList() const { return _acceptors; }
  
 private:
  void createDefault(const std::string&);

 private:
  int 
    _intronLenPenalty, _igrLenPenalty, _minIntronLen, _minIgrLen, _minAlignmentConnectionCutoff,
    _internalExonLenPenalty, _maxSeqLen, _overlapLen;

  double _intronPenalty, _igrPenalty, _internalExonPenalty;

  std::list<std::string> _donors;
  std::list<std::string> _acceptors;

  static ParamFile* _instance;
};

#endif // _PARAMFILE_H
