#ifndef _SINGLECLIPPEDCLUSTER2_H_
#define _SINGLECLIPPEDCLUSTER2_H_


#include "SingleClipped2.h"
#include "Contig.h"
#include <vector>
#include <map>

class SingleClippedCluster2 {
 public:
  SingleClippedCluster2(const SingleClipped2* clip);
  virtual ~SingleClippedCluster2();
  void add(const SingleClipped2* const clip);
  size_t size() const;
  int referenceId() const;
  int position() const;
  std::string consensus() const;
  // friend std::ostream& operator <<(std::ostream& stream, SingleClippedCluster& self);
 private:
  std::vector<const SingleClipped2*> clips;
  static char correctBase(const std::map<char, int>& bases, const std::map<char, int>& quals);
  static int secondLargest(const std::vector<int>& lens);
};

#endif /* _SINGLECLIPPEDCLUSTER2_H_ */
