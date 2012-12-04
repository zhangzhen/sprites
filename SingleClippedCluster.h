#ifndef SINGLE_CLIPPED_CLUSTER_H
#define SINGLE_CLIPPED_CLUSTER_H

#include "SingleClipped.h"
#include "Contig.h"
#include <vector>

bool comp(const std::string& s1, const std::string& s2);

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
 protected:
  std::string assembleClipped();
  std::string assembleMapped();
 private:
  std::string assemble(const std::vector<std::string>& seqs);
};

#endif /* SINGLE_CLIPPED_CLUSTER_H */
