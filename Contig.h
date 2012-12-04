#ifndef CONTIG_H
#define CONTIG_H

#include "Locus.h"
#include <string>

bool equals2(const std::string s1, const std::string s2, int mismatches=0);

class Contig
{
 private:
  std::string seq;
  Locus anchor;
  int marker;
  int num;
  
 public:
  Contig(const std::string& seq, const Locus& anchor, int marker, int num);
  virtual ~Contig();
  bool overlaps(const Contig& other, int mismatches = 0) const;
  bool operator==(const Contig& other) const;
  std::string sequence() const;
};

#endif /* CONTIG_H */
