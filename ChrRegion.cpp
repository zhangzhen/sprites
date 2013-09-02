#include "ChrRegion.h"
#include <sstream>
#include <cassert>

ChrRegion::ChrRegion(int id, int referenceId, int startPos, int endPos, int insertSize) :
    id(id), referenceId(referenceId), startPos(startPos), endPos(endPos), insertSize(insertSize) {}

int ChrRegion::getId() const { return id; }

int ChrRegion::getReferenceId() const { return referenceId; }

int ChrRegion::getStartPos() const { return startPos; }

int ChrRegion::getEndPos() const { return endPos; }

int ChrRegion::getInsertSize() const { return insertSize; }

size_t ChrRegion::length() const { return endPos - startPos; }

int ChrRegion::minDeletionLength(int mean, int std) const {
  assert(insertSize > mean);
  int delta = insertSize - mean;
  return delta < 150 ? std::max(delta - std, 0) : std::max(delta - 3 * std, 0);
}

int ChrRegion::maxDeletionLength(int mean, int std) const {
  assert(insertSize > mean);
  int delta = insertSize - mean;
  return delta < 150 ? delta + std : delta + 3 * std;
}

// bool ChrRegion::overlapsWith(const ChrRegion& other) const {
//   return false;
// }

// std::string ChrRegion::toString() const {
//   std::stringstream ss;
//   ss << chrom << ":" << startPos << "-" << endPos;
//   return ss.str();
// }

const EndPoint ChrRegion::getStart() const {
  return EndPoint(this, true);
}

const EndPoint ChrRegion::getEnd() const {
  return EndPoint(this, false);
}

// bool ChrRegion::compare(const ChrRegion& in1, const ChrRegion& in2) {
//   return in1.length() < in2.length();
// }

// bool ChrRegion::compareByStart(const ChrRegion& i1, const ChrRegion& i2) {
//   return i1.start < i2.start;
// }

// bool ChrRegion::compareByEnd(const ChrRegion& i1, const ChrRegion& i2) {
//   return i1.end < i2.end;
// }

std::ostream& operator <<(std::ostream& os, const ChrRegion& self) {
  return os << self.referenceId << ":" << self.startPos << "-" << self.endPos
            << "\t" << self.insertSize;
}

EndPoint::EndPoint(const ChrRegion* owner, bool start) :
    owner(owner), start(start) {
}

int EndPoint::ownerId() const {
  return owner->getId();
}

int EndPoint::position() const {
  if (start)
    return owner->getStartPos();
  return owner->getEndPos();
}

bool EndPoint::isStart() const {
  return start;
}

bool EndPoint::operator< (const EndPoint& other) const {
  assert(owner->getReferenceId() == other.owner->getReferenceId());
  if ((position() == other.position()) && (!start && other.start)) return true;
  return position() < other.position();
}

const ChrRegion* EndPoint::getOwner() const {
  return owner;
}
