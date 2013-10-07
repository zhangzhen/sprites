#ifndef _CHRREGIONCLUSTER_H_
#define _CHRREGIONCLUSTER_H_

#include "ChrRegion.h"
#include "TargetRegion.h"
#include <string>
#include <vector>

class ChrRegionCluster {
 private:
  std::vector<const ChrRegion*> elts;
  bool dirty;
 public:
  ChrRegionCluster();
  void add(const ChrRegion *in);
  // bool empty() const;
  std::string toString() const;
  void getOverlaps(int start, int end, std::vector<const ChrRegion*>& regions) const;
  bool getTargetRegion(int mean, int std, TargetRegion& regionOfInterest);
  friend std::ostream& operator <<(std::ostream& os, const ChrRegionCluster& crc)
      ;
 private:
  void removeAbnormalChrRegions(int threshold);
  void removeOutliersWithMad(std::vector<const ChrRegion*>& cleanRegions);
  void toLengthList(std::vector<int>& list);
  static int average(const std::vector<int>& v);
  static int median(std::vector<int>& v);
  static int mad(std::vector<int>& v);
};

#endif /* _CHRREGIONCLUSTER_H_ */
