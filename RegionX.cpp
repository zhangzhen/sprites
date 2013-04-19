#include "RegionX.h"
#include <cassert>

RegionX::RegionX(int referenceId, int startPosition, int endPosition)
    : referenceId(referenceId), startPosition(startPosition), endPosition(endPosition) { assert(checkRep()); }

RegionX::~RegionX() {}

bool RegionX::checkRep() {
  return startPosition < endPosition;
}

int RegionX::getReferenceId() const { return referenceId; }

int RegionX::getStartPosition() const { return startPosition; }

int RegionX::getEndPosition() const { return endPosition; }

int RegionX::length() const {
  return endPosition - startPosition;
}

// bool Region::operator== (const Region& other) const {
//   return start == other.start &&
//       end == other.end;
// }

// bool Region::operator!= (const Region& other) const {
//   return !(*this == other);
// }
