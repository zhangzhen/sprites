#include "LeftClipped.h"
#include <assert.h>

LeftClipped::LeftClipped(const Locus& loc, const std::string& seq, const std::string& quals, int start, int clippedLen)
    : SingleClipped(loc, seq, quals, start, clippedLen) {
  assert(checkRep());
}

LeftClipped::~LeftClipped() {}

std::string LeftClipped::mappedSeq() const {
  return seq.substr(clippedLen);
}

std::string LeftClipped::type() const {
  return "LeftClipped";
}

bool LeftClipped::checkRep() {
  return start == 0;
}
