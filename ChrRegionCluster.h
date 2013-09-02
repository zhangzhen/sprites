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
  // std::string toString() const;
  TargetRegion getTargetRegion(int mean, int std);
  // friend std::ostream& operator <<(std::ostream& os, const ChrRegionCluster& self)
      ;
 private:
  void removeAbnormalChrRegions(int threshold);
  int avgInsertSize() const;
};

#endif /* _CHRREGIONCLUSTER_H_ */
