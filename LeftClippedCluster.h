#ifndef LEFT_CLIPPED_CLUSTER_H
#define LEFT_CLIPPED_CLUSTER_H

#include "SingleClippedCluster.h"

class LeftClippedCluster : public SingleClippedCluster {
 public:
  LeftClippedCluster(const Locus& anchor);
  virtual ~LeftClippedCluster();
  Contig contig();
};

#endif /* LEFT_CLIPPED_CLUSTER_H */
