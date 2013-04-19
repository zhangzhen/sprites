#ifndef _CONTIG2_H_
#define _CONTIG2_H_

#include <string>
#include "Overlap.h"

class Contig2
{
 private:
  int referenceId;
  int clipPosition;
  int localClipPosition;
  std::string sequence;
  int numClips;
  
 public:
  Contig2(int referenceId, int clipPosition, int localClipPosition, const std::string& sequence, int numClips);
  virtual ~Contig2();
  int lowerBoundOfVarLength(const Contig2& other) const;
  int upperBoundOfVarLength(const Contig2& other) const;
  
  bool overlaps(const Contig2& other, int minOverlapLength, double maxMismatchRate, Overlap* const overlap) const;
  // bool operator== (const Contig2& other) const;
  // bool operator< (const Contig2& other) const;
  // friend std::ostream& operator <<(std::ostream& stream, const Contig2& self);
  // const std::string getReferenceName() const;
  int getReferenceId() const;
  int getClipPosition() const;
  int leftmostPosition() const;
  int rightmostPosition() const;
  // const std::string& getSequence() const;
  size_t length() const;
  int getLocalClipPosition() const;
  // static bool compare(const Contig& c, const Locus2& locus);
  // static bool compare2(const Locus2& locus, const Contig& c);
 private:
  static bool equals(const std::string& s1, const std::string& s2, int maxMismatches, int& numMismatches);
};

#endif /* _CONTIG2_H_ */
