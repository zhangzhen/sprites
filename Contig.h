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
  bool proximal;
  
 public:
  Contig(const std::string& seq, const Locus& anchor, int marker, int num, bool proximal);
  virtual ~Contig();
  bool overlaps(const Contig& other, int minSupportSize, int minOverlapLen, double mismatchRate, int& offset) const;
  bool operator== (const Contig& other) const;
  bool operator< (const Contig& other) const;
  friend std::ostream& operator <<(std::ostream& stream, const Contig& self);
  std::string sequence() const;
  Locus getAnchor() const;
  int getMarker() const { return marker; }
  bool getProximal() const { return proximal; }
 private:
  static bool overlaps2(const Contig& c1, const Contig& c2, int minSupportSize, int minOverlapLen, double mismatchRate, int& offset);
};

#endif /* CONTIG_H */
