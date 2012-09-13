#ifndef OVERLAP_INCLUDED
#define OVERLAP_INCLUDED

#include "Clip.h"

class Overlap
{
 public:
  Overlap(Clip* left, Clip* right);
  virtual ~Overlap();
  
 private:
  Clip* left;
  Clip* right;
  int len;
  int mismatches;
};

#endif // OVERLAP_INCLUDED
