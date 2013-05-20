#ifndef _SOFTCLIPCLUSTER_H_
#define _SOFTCLIPCLUSTER_H_

#include "SoftClip.h"
#include "Consensus.h"
#include <vector>
#include <map>

class SoftClipCluster {
 public:
  SoftClipCluster(const SoftClip* clip);
  void add(const SoftClip* clip);
  size_t size() const;
  int referenceId() const;
  int clipPosition() const;
  Consensus getConsensus() const;
  // friend std::ostream& operator <<(std::ostream& stream, SingleClippedCluster& self);
 private:
  std::vector<const SoftClip*> clips;

  std::string consensusSequence() const;
  int localClipPosition() const;
  static char correctBase(const std::map<char, int>& bases, const std::map<char, int>& quals);
  static int secondLargest(const std::vector<int>& lens);
};

#endif /* _SOFTCLIPCLUSTER_H_ */
