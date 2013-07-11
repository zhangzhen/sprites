#ifndef _OVERLAP_H_
#define _OVERLAP_H_

#include "Deletion.h"
#include <string>
#include <iostream>
#include <vector>

class SoftClip;

class Overlap
{
 public:
  Overlap();
  Overlap(const SoftClip *first, const SoftClip *second, int length, int numMismatches, int offset);
  // int getLength() const;
  // int getNumMismatches() const;
  double score() const;
  Deletion getDeletion() const;
  static bool equals(const std::string& s1, const std::string& s2, int maxMismatches, int& numMismatches);
  static bool getHighScoreOverlap(std::vector<Overlap> overlaps, Overlap& ov);

  friend std::ostream& operator <<(std::ostream& stream, const Overlap& o);

  int deletionLength() const;
  
 private:

  const SoftClip *first;
  const SoftClip *second;
  int length;
  int numMismatches;
  int offset;
};

#endif /* _OVERLAP_H_ */
