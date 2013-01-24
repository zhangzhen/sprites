#include "SingleClippedCluster.h"
#include <algorithm>
#include <sstream>

std::ostream& operator <<(std::ostream& stream, SingleClippedCluster& self) {
  stream << self.str();
  return stream;
}

bool comp(const std::string& s1, const std::string& s2) {
  return s1.size() < s2.size();
}

bool SingleClippedCluster::comp2(SingleClipped* sc1, SingleClipped* sc2) {
  return sc1->clippedSeq().size() < sc2->clippedSeq().size();
}

SingleClippedCluster::SingleClippedCluster(const Locus& anchor)
    : anchor(anchor) {}

SingleClippedCluster::~SingleClippedCluster() {}

void SingleClippedCluster::add(SingleClipped* cl) {
  cls.push_back(cl);
}

size_t SingleClippedCluster::size() { return cls.size(); }

Locus SingleClippedCluster::getAnchor() { return anchor; }

std::string SingleClippedCluster::assembleClipped() {
  std::vector<std::string> seqs;
  for (int i = 0; i < cls.size(); ++i) {
    seqs.push_back(cls[i]->clippedSeq());
  }
  return assemble(seqs);
}

std::string SingleClippedCluster::assembleMapped() {
  std::vector<std::string> seqs;
  for (int i = 0; i < cls.size(); ++i) {
    seqs.push_back(cls[i]->mappedSeq());
  }
  return assemble(seqs);
}

std::string SingleClippedCluster::assemble(const std::vector<std::string>& seqs) {
  return *max_element(seqs.begin(), seqs.end(), comp);
}

std::string SingleClippedCluster::str() {
  std::stringstream sstream;
  sstream << anchor << "\t[" << size() << "]" << std::endl;
  sstream << contig().sequence();
  return sstream.str();
}
