#include "IntervalCluster.h"
#include <iostream>
#include <algorithm>

IntervalCluster::IntervalCluster() : dirty(true) {}

void IntervalCluster::add(const Interval* in) {
  elts.push_back(in);
  if (!dirty) dirty = true;
}

// bool IntervalCluster::empty() const {
//   return elts.empty();
// }

// std::string IntervalCluster::toString() const {
//   std::stringstream ss;
//   ss << "#####################" << std::endl;
//   for (size_t i = 0; i < elts.size(); ++i) {
//     ss << *(elts[i]) << std::endl;
//   }
//   return ss.str();
// }

TargetRegion IntervalCluster::getTargetRegion(int mean, int std) {
  if (dirty) {
    removeAbnormalIntervals(mean + 4*std);
    dirty = false;
  }
  int avgDeltaLength = avgInsertSize() - mean;
  TargetRegion r = { elts.back()->getStartPos(), elts.front()->getEndPos(), avgDeltaLength - 3*std, avgDeltaLength + 3*std };
  std::cout << r.start << "\t" << r.end << "\t" << r.minDeletionLength << "\t" << r.maxDeletionLength << std::endl;
  return r;
}

void IntervalCluster::removeAbnormalIntervals(int threshold) {
  if (elts.size() <= 1) return;

  sort(elts.begin(), elts.end(), [](const Interval *in1, const Interval *in2) { return in1->length() > in2->length(); });
  int minLen = elts.back()->length();
  elts.erase(remove_if(elts.begin(), elts.end(), [minLen, threshold](const Interval* in) { return in->length() - minLen > threshold; }), elts.end());
  sort(elts.begin(), elts.end(), [](const Interval *in1, const Interval *in2) { return in1->getStartPos() < in2->getStartPos(); });
}

int IntervalCluster::avgInsertSize() const {
  return accumulate(elts.begin(), elts.end(), 0, [](int s, const Interval* in) { return s + in->getInsertSize(); }) / elts.size();
}

// std::ostream& operator <<(std::ostream& os, const IntervalCluster& self) {
//   return os << self.toString();
// }
