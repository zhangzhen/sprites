#include "Region.h"

Region::Region(const Locus& start, const Locus& end)
    : start(start), end(end) { assert(checkRep()); }

Region::~Region() {}

bool Region::checkRep() {
  return start.chrom() == end.chrom() && start.position <= end.position();
}

std::string chrom() { return start.chrom(); }

int start() { return start.position(); }

int end() { return end.position(); }
