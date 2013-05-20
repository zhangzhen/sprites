#include "Deletion.h"
#include <cassert>

Deletion::Deletion(int referenceId, int start2, int end1, int offset) :
    referenceId(referenceId), start2(start2), end1(end1), offset(offset) {
  assert(start2 > 0 && length() > 0);
}

int Deletion::getStart1() const { return start2 - offset; }

int Deletion::getEnd1() const { return end1; }

int Deletion::getStart2() const { return start2; }

int Deletion::getEnd2() const { return end1 + offset; }

int Deletion::length() const {
  return end1 + offset - start2;
}

// double Deletion::score() const { return score; }
