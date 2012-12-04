#ifndef LEFT_CLIPPED_H
#define LEFT_CLIPPED_H

#include "SingleClipped.h"

class LeftClipped : public SingleClipped {
 public:
  LeftClipped(const Locus& loc, const std::string& seq, int qual, int start, int clippedLen);
  virtual ~LeftClipped();
  std::string mappedSeq() const;
  std::string type() const;
 private:
  bool checkRep();
};

#endif /* LEFT_CLIPPED_H */
