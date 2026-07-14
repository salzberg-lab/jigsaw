//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
EvScores& RegionPrediction::getScores(dsu::Strand_t strnd) { 
  const int idx = getIdx(strnd);
  return _score[idx];
}

EvScores& RegionPrediction::getScoresOpp(dsu::Strand_t strnd) { 
  const int idx = getIdx(strnd);
  return _scoreOpp[idx];
}

const EvScores& RegionPrediction::getScores(dsu::Strand_t strnd) const { 
  const int idx = getIdx(strnd);
  return _score[idx];
}

const EvScores& RegionPrediction::getScoresOpp(dsu::Strand_t strnd) const { 
  const int idx = getIdx(strnd);
  return _scoreOpp[idx];
}

const AnnotationScore& RegionPrediction::getScore(unsigned i, dsu::Strand_t strnd) const { 
  const int idx = getIdx(strnd);
  const EvScores& scores = _score[idx];
  assert(i < EvScores::eNumBnd);
  return scores[i];
}

const AnnotationScore& RegionPrediction::getScoreOpp(unsigned i, dsu::Strand_t strnd) const { 
  const int idx = getIdx(strnd);
  const EvScores& scores = _scoreOpp[idx];
  assert(i < EvScores::eNumBnd);
  return scores[i];
}

int RegionPrediction::getIdx(dsu::Strand_t strnd) const {
  assert(strnd == dsu::eNeg || strnd == dsu::ePos );
  int intStrnd = static_cast<int>(strnd);
  return intStrnd;
}
