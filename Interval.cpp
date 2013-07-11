#include "Interval.h"
#include <sstream>
#include <cassert>

Interval::Interval(int id, int referenceId, int startPos, int endPos, int insertSize) :
    id(id), referenceId(referenceId), startPos(startPos), endPos(endPos), insertSize(insertSize) {}

int Interval::getId() const { return id; }

int Interval::getReferenceId() const { return referenceId; }

int Interval::getStartPos() const { return startPos; }

int Interval::getEndPos() const { return endPos; }

int Interval::getInsertSize() const { return insertSize; }

size_t Interval::length() const { return endPos - startPos; }

int Interval::minDeletionLength(int mean, int std) const {
  assert(insertSize > mean);
  int delta = insertSize - mean;
  return delta < 150 ? std::max(delta - std, 0) : std::max(delta - 3 * std, 0);
}

int Interval::maxDeletionLength(int mean, int std) const {
  assert(insertSize > mean);
  int delta = insertSize - mean;
  return delta < 150 ? delta + std : delta + 3 * std;
}

// bool Interval::overlapsWith(const Interval& other) const {
//   return false;
// }

// std::string Interval::toString() const {
//   std::stringstream ss;
//   ss << chrom << ":" << startPos << "-" << endPos;
//   return ss.str();
// }

const EndPoint Interval::getStart() const {
  return EndPoint(this, true);
}

const EndPoint Interval::getEnd() const {
  return EndPoint(this, false);
}

// bool Interval::compare(const Interval& in1, const Interval& in2) {
//   return in1.length() < in2.length();
// }

// bool Interval::compareByStart(const Interval& i1, const Interval& i2) {
//   return i1.start < i2.start;
// }

// bool Interval::compareByEnd(const Interval& i1, const Interval& i2) {
//   return i1.end < i2.end;
// }

// std::ostream& operator <<(std::ostream& os, const Interval& self) {
//   return os << self.toString();
// }

EndPoint::EndPoint(const Interval* owner, bool start) :
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

const Interval* EndPoint::getOwner() const {
  return owner;
}
