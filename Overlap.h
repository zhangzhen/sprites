#ifndef _OVERLAP_H_
#define _OVERLAP_H_

#include "Deletion.h"
#include <string>

class Overlap
{
 public:
  Overlap();
  Overlap(int referenceId, int clipPosition1, int clipPosition2, int numClips1, int numClips2, int length, int numMismatches, int offset);
  // int getLength() const;
  // int getNumMismatches() const;
  double score() const;
  Deletion getDeletion() const;
  static bool equals(const std::string& s1, const std::string& s2, int maxMismatches, int& numMismatches);
  
 private:
  int deletionLength() const;

  int referenceId;
  int clipPosition1;
  int clipPosition2;
  int numClips1;
  int numClips2;
  int length;
  int numMismatches;
  int offset;
};

#endif /* _OVERLAP_H_ */
