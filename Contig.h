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
  bool overlaps(const Contig& other, int minSupportSize, int minOverlapLen, double mismatchRate) const;
  bool operator== (const Contig& other) const;
  bool operator< (const Contig& other) const;
  friend std::ostream& operator <<(std::ostream& stream, const Contig& self);
  std::string sequence() const;
  Locus getAnchor() const;
};

#endif /* CONTIG_H */
