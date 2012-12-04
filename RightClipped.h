#ifndef RIGHT_CLIPPED_H
#define RIGHT_CLIPPED_H

#include "SingleClipped.h"

class RightClipped : public SingleClipped {
 public:
  RightClipped(const Locus& loc, const std::string& seq, int qual, int start, int clippedLen);
  virtual ~RightClipped();
  std::string mappedSeq() const;
  std::string type() const;

 private:
  bool checkRep();
};

#endif /* RIGHT_CLIPPED_H */
