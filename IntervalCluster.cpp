#include "IntervalCluster.h"
#include <sstream>
#include <algorithm>

IntervalCluster::IntervalCluster() {}

IntervalCluster::~IntervalCluster() {}

void IntervalCluster::add(Interval *in) {
  elts.push_back(in);
}

bool IntervalCluster::empty() const {
  return elts.empty();
}

std::string IntervalCluster::toString() const {
  std::stringstream ss;
  ss << "#####################" << std::endl;
  for (size_t i = 0; i < elts.size(); ++i) {
    ss << *(elts[i]) << std::endl;
  }
  return ss.str();
}

Region2 IntervalCluster::focalRegion() const {
  Region2 r = { elts.back()->getStartPos(), elts.front()->getEndPos() };
  return r;
}

void IntervalCluster::removeInvalidIntervals(unsigned threshold) {
  if (elts.size() <= 1) return;
  std::vector<size_t> idx(elts.size());
  for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
  sort(idx.begin(), idx.end(),
       [&elts](size_t i1, size_t i2) {return elts[i1]->length() > elts[i2]->length(); });
  size_t i = 0;
  while (elts[idx[i]]->length() - elts[idx[idx.size() - 1]]->length() > threshold
         && i < idx.size()) {
    elts.erase(elts.begin() + i);
    i++;
  }
}

std::ostream& operator <<(std::ostream& os, const IntervalCluster& self) {
  return os << self.toString();
}
