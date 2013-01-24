#include "LeftClippedCluster.h"
#include <algorithm>
#include <sstream>

LeftClippedCluster::LeftClippedCluster(const Locus& anchor)
    : SingleClippedCluster(anchor) {}

LeftClippedCluster::~LeftClippedCluster() {}

Contig LeftClippedCluster::contig() {
  SingleClipped *sc = *max_element(cls.begin(), cls.end(), SingleClippedCluster::comp2);
  int marker = sc->clippedSeq().size();
  return Contig(sc->sequence(), anchor, marker, cls.size(), false);
}

std::string LeftClippedCluster::str() {
  int maxClippedLen = assembleClipped().size();
  std::stringstream sstream;
  sstream << SingleClippedCluster::str() << std::endl;
  for (std::vector<SingleClipped*>::iterator itr = cls.begin();
       itr != cls.end();
       ++itr)
    sstream << std::string(maxClippedLen - (*itr)->getClippedLen(), ' ')
            << (*itr)->sequence()
            << std::endl;
  sstream << std::string(maxClippedLen, ' ') << '^';
  return sstream.str();
}
