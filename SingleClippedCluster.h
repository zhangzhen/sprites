#ifndef SINGLE_CLIPPED_CLUSTER_H
#define SINGLE_CLIPPED_CLUSTER_H

#include "SingleClipped.h"
#include "Contig.h"
#include <vector>
#include <map>

class SingleClippedCluster {
 protected:
  Locus anchor;
  std::vector<SingleClipped*> cls;
 public:
  SingleClippedCluster(const Locus& anchor);
  void add(SingleClipped* cl);
  size_t size();
  Locus getAnchor();
  std::string consensus();
  virtual ~SingleClippedCluster();
  virtual Contig contig() = 0;
  virtual std::string str();
  friend std::ostream& operator <<(std::ostream& stream, SingleClippedCluster& self);
 protected:
  std::string assembleClipped();
  std::string assembleMapped();
  static bool comp2(SingleClipped* sc1, SingleClipped* sc2);
 private:
  std::string assemble(const std::vector<std::string>& seqs);
  static char correctBase(std::map<char, int>& bases, std::map<char, int>& quals);
  static int secondLargest(const std::vector<int>& lens);
};

#endif /* SINGLE_CLIPPED_CLUSTER_H */
