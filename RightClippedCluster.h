#ifndef RIGHT_CLIPPED_CLUSTER_H
#define RIGHT_CLIPPED_CLUSTER_H

#include "SingleClippedCluster.h"

class RightClippedCluster : public SingleClippedCluster {
 public:
  RightClippedCluster(const Locus& anchor);
  virtual ~RightClippedCluster();
  Contig contig();
  std::string str();
};

#endif /* RIGHT_CLIPPED_CLUSTER_H */
