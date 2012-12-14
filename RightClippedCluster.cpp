#include "RightClippedCluster.h"
#include <sstream>
#include <iostream>

RightClippedCluster::RightClippedCluster(const Locus& anchor)
    : SingleClippedCluster(anchor) {}

RightClippedCluster::~RightClippedCluster() {}

Contig RightClippedCluster::contig() {
  int marker = assembleMapped().size();
  std::string seq = assembleMapped() + assembleClipped();
  return Contig(seq, anchor, marker, cls.size());
}

std::string RightClippedCluster::str() {
  int maxMappedLen = assembleMapped().size();
  std::stringstream sstream;
  sstream << SingleClippedCluster::str();
  for (std::vector<SingleClipped*>::iterator itr = cls.begin();
       itr != cls.end();
       ++itr)
    sstream << std::string(maxMappedLen - (*itr)->getMappedLen(), ' ')
            << (*itr)->sequence()
            << std::endl;
  return sstream.str();  
}
