#ifndef _CONSENSUS_H_
#define _CONSENSUS_H_

#include <string>
#include "Overlap.h"

class Consensus
{
 private:
  int referenceId;
  int clipPosition;
  int localClipPosition;
  std::string sequence;
  int numClips;
  
 public:
  Consensus(int referenceId, int clipPosition, int localClipPosition, const std::string& sequence, int numClips);
  
  int minDeletionLength(const Consensus& other) const;
  int maxDeletionLength(const Consensus& other) const;
  
  bool overlaps(const Consensus& other, int minOverlapLength, double maxMismatchRate, Overlap& overlap) const;
  // bool operator== (const Consensus& other) const;
  bool operator< (const Consensus& other) const;
  // friend std::ostream& operator <<(std::ostream& stream, const Consensus& self);
  // const std::string getReferenceName() const;
  int getReferenceId() const;
  int getClipPosition() const;
  // int leftmostPosition() const;
  // int rightmostPosition() const;
  // const std::string& getSequence() const;
  size_t length() const;
  int getLocalClipPosition() const;
  static bool compare(const Consensus& cons, const int pos);
  static bool compare2(const int pos, const Consensus& cons);
 private:
  // static bool equals(const std::string& s1, const std::string& s2, int maxMismatches, int& numMismatches);
};

#endif /* _CONSENSUS_H_ */

