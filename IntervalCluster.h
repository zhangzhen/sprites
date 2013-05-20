#ifndef _INTERVALCLUSTER_H_
#define _INTERVALCLUSTER_H_

#include "Interval.h"
#include "TargetRegion.h"
#include <string>
#include <vector>

class IntervalCluster {
 private:
  std::vector<const Interval*> elts;
  bool dirty;
 public:
  IntervalCluster();
  void add(const Interval *in);
  // bool empty() const;
  // std::string toString() const;
  TargetRegion getTargetRegion(int mean, int std);
  // friend std::ostream& operator <<(std::ostream& os, const IntervalCluster& self)
      ;
 private:
  void removeAbnormalIntervals(int threshold);
  int avgInsertSize() const;
};

#endif /* _INTERVALCLUSTER_H_ */
