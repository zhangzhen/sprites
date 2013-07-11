#include "Deletion.h"
#include <cassert>

Deletion::Deletion(int referenceId, int start2, int end1, int offset) :
    referenceId(referenceId), start2(start2), end1(end1), offset(offset) {
  assert(start2 > 0 && length() > 0);
}

int Deletion::getReferenceId() const { return referenceId; }

int Deletion::getStart1() const { return start2 - offset; }

int Deletion::getEnd1() const { return end1; }

int Deletion::getStart2() const { return start2; }

int Deletion::getEnd2() const { return end1 + offset; }

int Deletion::length() const {
  return end1 + offset - start2;
}

bool Deletion::overlaps(const Deletion& other) const {
  return referenceId == other.referenceId &&
      (getStart2() >= other.getStart2() && getStart2() < other.getEnd2() ||
       other.getStart2() >= getStart2() && other.getStart2() < getEnd2());
}

bool Deletion::operator <(const Deletion& other) const {
  if (referenceId != other.referenceId)
    return referenceId < other.referenceId;
  return getStart2() < other.getStart2();
}
// double Deletion::score() const { return score; }
