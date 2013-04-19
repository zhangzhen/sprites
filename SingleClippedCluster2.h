#ifndef _SINGLECLIPPEDCLUSTER2_H_
#define _SINGLECLIPPEDCLUSTER2_H_

#include "SingleClipped2.h"
#include "Contig2.h"
#include <vector>
#include <map>

class SingleClippedCluster2 {
 public:
  SingleClippedCluster2(const SingleClipped2* clip);
  virtual ~SingleClippedCluster2();
  void add(const SingleClipped2* clip);
  size_t size() const;
  int referenceId() const;
  int clipPosition() const;
  std::string consensus() const;
  Contig2* contig() const;
  // friend std::ostream& operator <<(std::ostream& stream, SingleClippedCluster& self);
 private:
  std::vector<const SingleClipped2*> clips;
  int localClipPosition() const;
  static char correctBase(const std::map<char, int>& bases, const std::map<char, int>& quals);
  static int secondLargest(const std::vector<int>& lens);
};

#endif /* _SINGLECLIPPEDCLUSTER2_H_ */
