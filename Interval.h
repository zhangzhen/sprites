#ifndef INTERVAL_H
#define INTERVAL_H

#include <string>
#include "Point.h"


class Interval {
 private:
  int id;
  std::string chrom;
  int startPos;
  int endPos;
  int insertSize;
 public:
  Interval(int id, std::string chrom, int startPos, int endPos, int insertSize);
  virtual ~Interval();
  int getId() const;
  std::string getChrom() const;
  int getStartPos() const;
  int getEndPos() const;
  int getInsertSize() const;
  size_t length() const;
  // bool overlapsWith(const Interval& other) const;
  std::string toString() const;
  Point createStartPoint();
  Point createEndPoint();
  // static bool compare(const Interval& in1, const Interval& in2);
  
  // static bool compareByStart(const Interval& i1, const Interval& i2);
  // static bool compareByEnd(const Interval& i1, const Interval& i2);
  friend std::ostream& operator <<(std::ostream& os, const Interval& self);  
};

#endif /* INTERVAL_H */
