#ifndef REGION_H
#define REGION_H

#include "Locus.h"

class Region
{
 private:
  Locus start;
  Locus end;
  
 public:
  Region(const Locus& start, const Locus& end);
  virtual ~Region();

  std::string chrom() const;
  int start() const;
  int end() const;

 private:
  bool checkRep();
};

#endif /* REGION_H */
