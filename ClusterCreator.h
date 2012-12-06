#ifndef CLUSTER_CREATOR_H
#define CLUSTER_CREATOR_H

#include "SingleClippedCluster.h"

class ClusterCreator {
 public:
  virtual SingleClippedCluster* createCluster(const Locus& anchor) = 0;
};

template <class TheCluster>
class StandardClusterCreator : public ClusterCreator {
 public:
  virtual SingleClippedCluster* createCluster(const Locus& anchor);
};

template <class TheCluster>
SingleClippedCluster* StandardClusterCreator<TheCluster>::createCluster(const Locus& anchor) {
  return new TheCluster(anchor);
}

#endif /* CLUSTER_CREATOR_H */
