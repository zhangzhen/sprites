#include "SingleClippedCluster.h"
#include <algorithm>
#include <sstream>
#include <cassert>

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

std::string SingleClippedCluster::consensus() {
  assert(size() > 0);
  if (size() == 1) return cls[0]->sequence();
  
  std::vector<int> lens1(size());
  for (int i = 0; i < lens1.size(); i++) lens1[i] = cls[i]->lengthOfLeftPart();
  int N1 = secondLargest(lens1);

  std::vector<int> diffs(size());
  for (int i = 0; i < diffs.size(); i++)
    diffs[i] = cls[i]->lengthOfLeftPart() - N1;
  std::vector<int> lens2(size());
  for (int i = 0; i < lens2.size(); i++) lens2[i] = cls[i]->lengthOfRightPart();
  
  int N = N1 + secondLargest(lens2);
  std::string seq;
  for (int i = 0; i < N; i++) {
    std::map<char, int> bases;
    std::map<char, int> quals;
    for (int j = 0; j < size(); j++) {
      int ind = diffs[j] + i;
      if (ind < 0 || ind >= cls[j]->length()) continue;
      ++bases[cls[j]->at(ind)];
      quals[cls[j]->at(ind)] += cls[j]->qual(ind) - '!';
    }
    for (std::map<char, int>::iterator it = bases.begin(); it != bases.end(); ++it)
      quals[it->first] /= it->second;
    seq += correctBase(bases, quals);
  }
  
}

char SingleClippedCluster::correctBase(std::map<char, int>& bases, std::map<char, int>& quals) {
  std::map<char, int>::iterator it;
  int maxnum = 0;
  for (std::map<char, int>::iterator it = bases.begin(); it != bases.end(); ++it)
    if (maxnum < it->second) maxnum = it->second;
  int maxqual= 0;
  char ch;
  for (std::map<char, int>::iterator it = bases.begin(); it != bases.end(); ++it) {
    if (it->second == maxnum && maxqual < quals[it->first]) {
      maxqual = quals[it->first];
      ch = it->first;
    }
  }
  return ch;
}

int SingleClippedCluster::secondLargest(const std::vector<int>& lens) {
  assert(lens.size() > 1);
  int max = lens[0];
  int secondMax = lens[1];
  if (max < secondMax) std::swap(max, secondMax);

  for (int i = 2; i < lens.size(); i++) {
    if (max <= lens[i]) {
      secondMax = max;
      max = lens[i];
    } else if (secondMax < lens[i]) {
      secondMax = lens[i];
    }
  }
  return secondMax;
}
