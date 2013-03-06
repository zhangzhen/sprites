#ifndef INTERVAL_H
#define INTERVAL_H

#include <string>
#include "Point.h"


class Interval {
 private:
  unsigned id;
  std::string chrom;
  unsigned startPos;
  unsigned endPos;
 public:
  Interval(unsigned id, std::string chrom, unsigned startPos, unsigned endPos);
  virtual ~Interval();
  unsigned getId() const;
  std::string getChrom() const;
  unsigned getStartPos() const;
  unsigned getEndPos() const;
  size_t length() const;
  bool overlapsWith(const Interval& other) const;
  std::string toString() const;
  Point createStartPoint();
  Point createEndPoint();
  
  // static bool compareByStart(const Interval& i1, const Interval& i2);
  // static bool compareByEnd(const Interval& i1, const Interval& i2);
  friend std::ostream& operator <<(std::ostream& os, const Interval& self);  
};

#endif /* INTERVAL_H */
