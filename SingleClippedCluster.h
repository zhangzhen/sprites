#ifndef SINGLE_CLIPPED_CLUSTER_H
#define SINGLE_CLIPPED_CLUSTER_H

#include "SingleClipped.h"
#include "Contig.h"
#include <vector>

class SingleClippedCluster {
 protected:
  Locus anchor;
  std::vector<SingleClipped*> cls;
 public:
  SingleClippedCluster(const Locus& anchor);
  void add(SingleClipped* cl);
  size_t size();
  virtual ~SingleClippedCluster();
  virtual Contig contig() = 0;
  virtual std::string str();
  friend std::ostream& operator <<(std::ostream& stream, SingleClippedCluster& self);
 protected:
  std::string assembleClipped();
  std::string assembleMapped();
 private:
  std::string assemble(const std::vector<std::string>& seqs);
};

#endif /* SINGLE_CLIPPED_CLUSTER_H */
