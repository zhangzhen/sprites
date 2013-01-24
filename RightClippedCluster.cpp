#include "RightClippedCluster.h"
#include <sstream>
#include <algorithm>
#include <iostream>

RightClippedCluster::RightClippedCluster(const Locus& anchor)
    : SingleClippedCluster(anchor) {}

RightClippedCluster::~RightClippedCluster() {}

Contig RightClippedCluster::contig() {
  SingleClipped *sc = *min_element(cls.begin(), cls.end(), SingleClippedCluster::comp2);
  int marker = sc->mappedSeq().size();
  return Contig(sc->sequence(), anchor, marker, cls.size(), true);
}

std::string RightClippedCluster::str() {
  int maxMappedLen = assembleMapped().size();
  std::stringstream sstream;
  sstream << SingleClippedCluster::str() << std::endl;
  for (std::vector<SingleClipped*>::iterator itr = cls.begin();
       itr != cls.end();
       ++itr)
    sstream << std::string(maxMappedLen - (*itr)->getMappedLen(), ' ')
            << (*itr)->sequence()
            << std::endl;
  sstream << std::string(maxMappedLen, ' ') << '^';
  return sstream.str();  
}
