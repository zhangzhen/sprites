#include "LeftClippedCluster.h"

LeftClippedCluster::LeftClippedCluster(const Locus& anchor)
    : SingleClippedCluster(anchor) {}

LeftClippedCluster::~LeftClippedCluster() {}

Contig LeftClippedCluster::contig() {
  int marker = assembleClipped().size();
  std::string seq = assembleClipped() + assembleMapped();
  return Contig(seq, anchor, marker, cls.size());
}
