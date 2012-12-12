#include "Region.h"
#include <assert.h>

Region::Region() {}

Region::Region(const Locus& start, const Locus& end)
    : start(start), end(end) { assert(checkRep()); }

Region::~Region() {}

bool Region::checkRep() {
  return start.chrom() == end.chrom() && start.position() <= end.position();
}

std::string Region::chrom() const { return start.chrom(); }

int Region::getStart() const { return start.position(); }

int Region::getEnd() const { return end.position(); }

int Region::length() const {
  return end.position() - start.position();
}

bool Region::operator== (const Region& other) const {
  return start == other.start &&
      end == other.end;
}

bool Region::operator!= (const Region& other) const {
  return !(*this == other);
}
