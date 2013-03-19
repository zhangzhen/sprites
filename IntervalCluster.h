#ifndef _INTERVALCLUSTER_H_
#define _INTERVALCLUSTER_H_

#include "Interval.h"
#include <string>
#include <vector>


struct Region2 {
  unsigned low;
  unsigned high;
};

class IntervalCluster {
 private:
  std::vector<Interval*> elts;
 public:
  IntervalCluster();
  virtual ~IntervalCluster();
  void add(Interval *in);
  bool empty() const;
  std::string toString() const;
  Region2 focalRegion() const;
  friend std::ostream& operator <<(std::ostream& os, const IntervalCluster& self);
  void removeInvalidIntervals(unsigned threshold);
};

#endif /* _INTERVALCLUSTER_H_ */
