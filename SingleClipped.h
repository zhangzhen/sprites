#ifndef SINGLE_CLIPPED_H
#define SINGLE_CLIPPED_H

#include "Locus.h"

class SingleClipped {
 protected:
  Locus loc;
  std::string seq;
  std::string quals;
  int start;
  int clippedLen;
  
 public:
  SingleClipped(const Locus& loc, const std::string& seq, const std::string& quals, int start, int clippedLen);
  virtual ~SingleClipped();
  Locus anchor() const;
  // std::string chr() const;
  // int position() const;
  int length() const;
  int lengthOfLeftPart() const;
  int lengthOfRightPart() const;
  char at(int i) const;
  char qual(int i) const;
  std::string sequence() const;
  std::string clippedSeq() const;
  int getClippedLen() const;
  int getMappedLen() const;
  virtual std::string mappedSeq() const = 0;
  virtual std::string type() const = 0;

 private:
  bool checkRep();
};

#endif
