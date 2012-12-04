#include "RightClippedCluster.h"

RightClippedCluster::RightClippedCluster(const Locus& anchor)
    : SingleClippedCluster(anchor) {}

RightClippedCluster::~RightClippedCluster() {}

Contig RightClippedCluster::contig() {
  int marker = assembleMapped().size();
  std::string seq = assembleMapped() + assembleClipped();
  return Contig(seq, anchor, marker, cls.size());
}
