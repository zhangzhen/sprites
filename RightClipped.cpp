#include "RightClipped.h"
#include <assert.h>
#include <iostream>

RightClipped::RightClipped(const Locus& loc, const std::string& seq, const std::string& quals, int start, int clippedLen)
    : SingleClipped(loc, seq, quals, start, clippedLen) {
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
