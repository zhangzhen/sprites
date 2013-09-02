#include "ChrRegionCluster.h"
// #include <iostream>
#include <algorithm>

ChrRegionCluster::ChrRegionCluster() : dirty(true) {}

void ChrRegionCluster::add(const ChrRegion* in) {
  elts.push_back(in);
  if (!dirty) dirty = true;
}

// bool ChrRegionCluster::empty() const {
//   return elts.empty();
// }

// std::string ChrRegionCluster::toString() const {
//   std::stringstream ss;
//   ss << "#####################" << std::endl;
//   for (size_t i = 0; i < elts.size(); ++i) {
//     ss << *(elts[i]) << std::endl;
//   }
//   return ss.str();
// }

TargetRegion ChrRegionCluster::getTargetRegion(int mean, int std) {
  if (dirty) {
    removeAbnormalChrRegions(mean + 4*std);
    dirty = false;
  }
  int avgDeltaLength = avgInsertSize() - mean;
  TargetRegion r = { elts.back()->getStartPos(), elts.front()->getEndPos(), avgDeltaLength - 3*std, avgDeltaLength + 3*std };
  // std::cout << r.start << "\t" << r.end << "\t" << r.minDeletionLength << "\t" << r.maxDeletionLength << std::endl;
  return r;
}

void ChrRegionCluster::removeAbnormalChrRegions(int threshold) {
  if (elts.size() <= 1) return;

  sort(elts.begin(), elts.end(), [](const ChrRegion *in1, const ChrRegion *in2) { return in1->length() > in2->length(); });
  int minLen = elts.back()->length();
  elts.erase(remove_if(elts.begin(), elts.end(), [minLen, threshold](const ChrRegion* in) { return in->length() - minLen > threshold; }), elts.end());
  sort(elts.begin(), elts.end(), [](const ChrRegion *in1, const ChrRegion *in2) { return in1->getStartPos() < in2->getStartPos(); });
}

int ChrRegionCluster::avgInsertSize() const {
  return accumulate(elts.begin(), elts.end(), 0, [](int s, const ChrRegion* in) { return s + in->getInsertSize(); }) / elts.size();
}

// std::ostream& operator <<(std::ostream& os, const ChrRegionCluster& self) {
//   return os << self.toString();
// }
