#include "ChrRegionCluster.h"
#include "DFinderHelper.h"
#include <iostream>
#include <iterator>
#include <algorithm>

ChrRegionCluster::ChrRegionCluster() : dirty(true) {}

void ChrRegionCluster::add(const ChrRegion* in) {
  elts.push_back(in);
  if (!dirty) dirty = true;
}

// bool ChrRegionCluster::empty() const {
//   return elts.empty();
// }

std::string ChrRegionCluster::toString() const {
  std::string str = ">>>>>>>>>>>>>>>>>>>>>>\n";
  for (auto itr = elts.begin(); itr != elts.end(); ++itr)
      str += (*itr)->toString() + "\n";
  return str;
}

void ChrRegionCluster::getOverlaps(int start, int end, std::vector<const ChrRegion*>& regions) const {
    copy_if(elts.begin(), elts.end(), back_inserter(regions),
	    [start, end](const ChrRegion *cr) {
		return ((start >= cr->getStartPos() && start < cr->getEndPos()) ||
			(cr->getStartPos() >= start && cr->getStartPos() < end));
	    });
}

bool ChrRegionCluster::getTargetRegion(int mean, int std, TargetRegion& regionOfInterest) {
  // if (dirty) {
  //   removeAbnormalChrRegions(mean + 4*std);
  //   dirty = false;
  // }
    // std::vector<const ChrRegion*> cleanRegions;
    // removeOutliersWithMad(cleanRegions);
    // if (cleanRegions.empty()) return false;
    // std::cout << "\n>>>>>>>>>>>>>>>>>>>>>>>> [" << cleanRegions.size() << "]" << std::endl;
    // std::transform(cleanRegions.begin(), cleanRegions.end(), std::ostream_iterator<const ChrRegion&>(std::cout, "\n"), [](const ChrRegion *cr) { return *cr; });
    // std::vector<int> lens;
    // toLengthList(lens);
    // const ChrRegion *cr = *max_element(elts.begin(), elts.end(), [](const ChrRegion *r1, const ChrRegion *r2) { return r1->length() < r2->length(); });
    const ChrRegion *cr = elts[elts.size() / 2];
    // const ChrRegion *cr = cleanRegions[0];
  int deltaLength = cr->getInsertSize() - mean;
  regionOfInterest = { cr->getStartPos(),
		       cr->getEndPos(),
		       std::max(0, deltaLength - 3*std),
		       deltaLength + 3*std,
		       elts.back()->getStartPos(),
		       elts.front()->getEndPos(),
		       elts.size() };
  // std::cout << r.start << "\t" << r.end << "\t" << r.minDeletionLength << "\t" << r.maxDeletionLength << std::endl;
  return true;
}

void ChrRegionCluster::removeAbnormalChrRegions(int threshold) {
  if (elts.size() <= 1) return;

  sort(elts.begin(), elts.end(), [](const ChrRegion *in1, const ChrRegion *in2) { return in1->length() > in2->length(); });
  int minLen = elts.back()->length();
  elts.erase(remove_if(elts.begin(), elts.end(), [minLen, threshold](const ChrRegion* in) { return in->length() - minLen > threshold; }), elts.end());
  sort(elts.begin(), elts.end(), [](const ChrRegion *in1, const ChrRegion *in2) { return in1->getStartPos() < in2->getStartPos(); });
}

void ChrRegionCluster::toLengthList(std::vector<int>& list) {
    list.resize(elts.size());
    std::transform(elts.begin(), elts.end(), list.begin(), [](const ChrRegion *cr) { return cr->length(); });
}

int ChrRegionCluster::median(std::vector<int>& v) {
    size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + n, v.end());
    return v[n];
}

int ChrRegionCluster::mad(std::vector<int>& v) {
    int m = median(v);
    std::vector<int> medians(v.size());
    std::transform(v.begin(), v.end(), medians.begin(), [m](int i) { return abs(i - m); });
    return median(medians);
}

void ChrRegionCluster::removeOutliersWithMad(std::vector<const ChrRegion*>& cleanRegions) {
    if (elts.size() == 0) return;
    if (elts.size() == 1) {
	cleanRegions.push_back(elts[0]);
	return;
    }
    if (elts.size() == 2 && isInconsistent(elts[0]->length(), elts[1]->length())) {
	cleanRegions.push_back(elts[0]);
	cleanRegions.push_back(elts[1]);
	return;
    }
    std::vector<int> lens;
    toLengthList(lens);
    int m1 = median(lens);
    int m2 = mad(lens);
    for (auto itr = elts.begin(); itr != elts.end(); ++itr) {
	int l = (*itr)->length();
	if (l >= m1 - 3 * m2 && l <= m1 + 3 * m2) cleanRegions.push_back(*itr);
    }
}

int ChrRegionCluster::average(const std::vector<int>& v) {
  return accumulate(v.begin(), v.end(), 0) / v.size();
}

std::ostream& operator <<(std::ostream& os, const ChrRegionCluster& crc) {
  return os << crc.toString();
}
