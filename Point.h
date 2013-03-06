#ifndef _POINT_H_
#define _POINT_H_

#include <string>

class Interval;

class Point {
 private:
  Interval *interval;
  bool start;
 public:
  Point(Interval *interval, bool start);
  virtual ~Point();
  unsigned intervalId() const;
  std::string chrom() const;
  unsigned position() const;
  bool isStart() const;
  bool operator< (const Point& other) const;
  Interval* getInterval() const;
};

#endif /* _POINT_H_ */
