#include "SingleClipped.h"
#include <assert.h>

SingleClipped::SingleClipped(const Locus& loc, const std::string& seq, int qual, int start, int clippedLen)
    : loc(loc),
      seq(seq),
      qual(qual),
      start(start),
      clippedLen(clippedLen) {
  assert(checkRep());
}

SingleClipped::~SingleClipped() {}

std::string SingleClipped::clippedSeq() const {
  return seq.substr(start, clippedLen);
}

bool SingleClipped::checkRep() {
  return start >= 0 && start + clippedLen <= seq.size();
}
