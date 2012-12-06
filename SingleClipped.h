#ifndef SINGLE_CLIPPED_H
#define SINGLE_CLIPPED_H

#include "Locus.h"

class SingleClipped {
 protected:
  Locus loc;
  std::string seq;
  int qual;
  int start;
  int clippedLen;
  
 public:
  SingleClipped(const Locus& loc, const std::string& seq, int qual, int start, int clippedLen);
  virtual ~SingleClipped();
  Locus anchor() const;
  // std::string chr() const;
  // int position() const;
  std::string clippedSeq() const;
  virtual std::string mappedSeq() const = 0;
  virtual std::string type() const = 0;

 private:
  bool checkRep();
};

#endif
