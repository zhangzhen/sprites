#include "Point.h"
#include "Interval.h"

Point::Point(Interval *interval, bool start) : interval(interval), start(start) {}

Point::~Point() {}

unsigned Point::intervalId() const { return interval->getId(); }

std::string Point::chrom() const { return interval->getChrom(); }

unsigned Point::position() const {
  if (start)
    return interval->getStartPos();
  return interval->getEndPos();
}

bool Point::isStart() const { return start; }

bool Point::operator< (const Point& other) const {
  if ((position() == other.position()) && (!start && other.start)) return true;
  return position() < other.position();
}

Interval* Point::getInterval() const { return interval; }
