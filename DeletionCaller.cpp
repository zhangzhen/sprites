#include "DeletionCaller.h"
#include <algorithm>

const int MinDelLen = 50;

void DeletionCaller::callAll(const std::vector<Region2>& in,
                             std::vector<Contig>& cons1,
                             std::vector<Contig>& cons2,
                             std::vector<Region>& calls,
                             int minOverlapLen,
                             int maxMismatches) {
  for (size_t i = 0; i < in.size(); i++) {
    std::vector<Contig>::iterator first1 = lower_bound(cons1.begin(), cons1.end(), in[i].low, Contig::compare);
    std::vector<Contig>::iterator last1 = upper_bound(cons1.begin(), cons1.end(), in[i].high, Contig::compare2);
    std::vector<Contig>::iterator first2 = lower_bound(cons2.begin(), cons2.end(), in[i].low, Contig::compare);
    std::vector<Contig>::iterator last2 = upper_bound(cons2.begin(), cons2.end(), in[i].high, Contig::compare2);
    callOne(first1, last1, first2, last2, calls, in[i], minOverlapLen, maxMismatches);
  }
}

void DeletionCaller::callOne(std::vector<Contig>::iterator first1,
                             std::vector<Contig>::iterator last1,
                             std::vector<Contig>::iterator first2,
                             std::vector<Contig>::iterator last2,
                             std::vector<Region>& calls,
                             const Region2& region,
                             int minOverlapLen,
                             int maxMismatches) {
  for (std::vector<Contig>::iterator it = first1; it != last1; it++) {
    first2 = upper_bound(first2, last2, *it);
    for (std::vector<Contig>::iterator it2 = first2; it2 != last2; it2++) {
      // if ((*it).upperBoundOfVarlength(*it2) < region.minDeltaLength) continue;
      // if ((*it).lowerBoundOfVarlength(*it2) > region.maxDeltaLength) break;
      int offset = 0;
      int overlapLen = (*it).overlaps(*it2, minOverlapLen, maxMismatches, offset);
      unsigned s = (*it).getAnchor().position();
      unsigned t = (*it2).getAnchor().position() + offset;
      if (overlapLen > 0 && t - s >= MinDelLen) {
        Locus end = (*it2).getAnchor();
        calls.push_back(Region((*it).getAnchor(), Locus(end.chrom(), end.position() + offset), overlapLen));
        return;
      }
    }
  }
}
