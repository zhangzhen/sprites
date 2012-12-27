#include "SingleClipped.h"
#include <assert.h>
#include <iostream>

SingleClipped::SingleClipped(const Locus& loc, const std::string& seq, int qual, int start, int clippedLen)
    : loc(loc),
      seq(seq),
      qual(qual),
      start(start),
      clippedLen(clippedLen) {
  assert(checkRep());
}

SingleClipped::~SingleClipped() {}

Locus SingleClipped::anchor() const { return loc; }

std::string SingleClipped::sequence() const { return seq; }

std::string SingleClipped::clippedSeq() const {
  return seq.substr(start, clippedLen);
}

int SingleClipped::getClippedLen() const {
  return clippedLen;
}

int SingleClipped::getMappedLen() const {
  return seq.size() - clippedLen;
}

bool SingleClipped::checkRep() {
  // std::cout << start << "\t" << clippedLen << "\t" << seq << std::endl;  
  return start >= 0 && start + clippedLen <= seq.size();
}
