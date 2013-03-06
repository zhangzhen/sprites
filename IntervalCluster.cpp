#include "IntervalCluster.h"
#include <sstream>

IntervalCluster::IntervalCluster() {}

IntervalCluster::~IntervalCluster() {}

void IntervalCluster::add(Interval *in) {
  elts.push_back(in);
}

bool IntervalCluster::empty() const {
  return elts.empty();
}

std::string IntervalCluster::toString() const {
  std::stringstream ss;
  ss << "#####################" << std::endl;
  for (size_t i = 0; i < elts.size(); ++i) {
    ss << *(elts[i]) << std::endl;
  }
  return ss.str();
}

Region2 IntervalCluster::focalRegion() const {
  Region2 r;
  r.low = elts.back()->getStartPos();
  r.high = elts.front()->getEndPos();
  return r;
}

std::ostream& operator <<(std::ostream& os, const IntervalCluster& self) {
  return os << self.toString();
}
