#ifndef _INTERVALCLUSTER_H_
#define _INTERVALCLUSTER_H_

#include "Interval.h"
#include <string>
#include <vector>


struct Region2 {
  int low;
  int high;
  int minDeltaLength;
  int maxDeltaLength;
};

class IntervalCluster {
 private:
  std::vector<Interval*> elts;
  bool valid;
 public:
  IntervalCluster();
  virtual ~IntervalCluster();
  void add(Interval *in);
  bool empty() const;
  std::string toString() const;
  Region2 focalRegion(int mean, int std);
  friend std::ostream& operator <<(std::ostream& os, const IntervalCluster& self);
 private:
  void removeInvalidIntervals(int threshold);
  int avgInsertSize() const;
};

#endif /* _INTERVALCLUSTER_H_ */
