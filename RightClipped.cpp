#include "RightClipped.h"
#include <assert.h>

RightClipped::RightClipped(const Locus& loc, const std::string& seq, int qual, int start, int clippedLen)
    : SingleClipped(loc, seq, qual, start, clippedLen) {
  assert(checkRep());
}

RightClipped::~RightClipped() {}

std::string RightClipped::mappedSeq() const {
  return seq.substr(0, start);
}

std::string RightClipped::type() const {
  return "RightClipped";
}

bool RightClipped::checkRep() {
  return start + clippedLen == seq.size();
}
