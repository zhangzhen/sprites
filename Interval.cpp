#include "Interval.h"
#include <sstream>

Interval::Interval(unsigned id, std::string chrom, unsigned startPos, unsigned endPos) :
    id(id), chrom(chrom), startPos(startPos), endPos(endPos) {}

Interval::~Interval() {}

unsigned Interval::getId() const { return id; }

std::string Interval::getChrom() const { return chrom; }

unsigned Interval::getStartPos() const { return startPos; }

unsigned Interval::getEndPos() const { return endPos; }

size_t Interval::length() const { return endPos - startPos; }

bool Interval::overlapsWith(const Interval& other) const {
  return false;
}

std::string Interval::toString() const {
  std::stringstream ss;
  ss << chrom << ":" << startPos << "-" << endPos;
  return ss.str();
}

Point Interval::createStartPoint() {
  return Point(this, true);
}

Point Interval::createEndPoint() {
  return Point(this, false);
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

std::ostream& operator <<(std::ostream& os, const Interval& self) {
  return os << self.toString();
}
