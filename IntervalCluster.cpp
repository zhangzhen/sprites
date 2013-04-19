#include "IntervalCluster.h"
#include <sstream>
#include <algorithm>

IntervalCluster::IntervalCluster() : valid(false) {}

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

Region2 IntervalCluster::focalRegion(int mean, int std) {
  if (!valid) {
    removeInvalidIntervals(mean + 4*std);
    valid = true;
  }
  int avgDeltaLength = avgInsertSize() - mean;
  Region2 r = { elts.back()->getStartPos(), elts.front()->getEndPos(), avgDeltaLength - 3*std, avgDeltaLength + 3*std };
  return r;
}

void IntervalCluster::removeInvalidIntervals(int threshold) {
  if (elts.size() <= 1) return;

  sort(elts.begin(), elts.end(), [](const Interval *in1, const Interval *in2) { return in1->length() > in2->length(); });
  int minLen = elts.back()->length();
  elts.erase(remove_if(elts.begin(), elts.end(), [minLen, threshold](Interval* in) { return in->length() - minLen > threshold; }), elts.end());
}

int IntervalCluster::avgInsertSize() const {
  return accumulate(elts.begin(), elts.end(), 0, [](int s, Interval* in) { return s + in->getInsertSize(); }) / elts.size();
}

std::ostream& operator <<(std::ostream& os, const IntervalCluster& self) {
  return os << self.toString();
}
