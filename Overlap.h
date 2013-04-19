#ifndef _OVERLAP_H_
#define _OVERLAP_H_

#include "RegionX.h"

class Overlap
{
 public:
  Overlap(int referenceId, int clipPosition1, int clipPosition2, int numClips1, int numClips2, int length, int numMismatches, int offset);
  virtual ~Overlap();
  // int getLength() const;
  // int getNumMismatches() const;
  int regionLength() const;
  double score() const;
  const RegionX* primaryRegion() const;
  const RegionX* secondaryRegion() const;
  
 private:
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
