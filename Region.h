#ifndef REGION_H
#define REGION_H

#include "Locus.h"

class Region
{
 private:
  Locus start;
  Locus end;
  int overlapLen;
  
 public:
  Region();
  Region(const Locus& start, const Locus& end, int overlapLen);
  virtual ~Region();

  std::string chrom() const;
  int getStart() const;
  int getEnd() const;
  int length() const;
  int overlapLength() const;

  bool operator== (const Region& other) const;
  bool operator!= (const Region& other) const;

 private:
  bool checkRep();
};

#endif /* REGION_H */
