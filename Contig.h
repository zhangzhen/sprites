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
  Contig();
  Contig(const std::string& seq, const Locus& anchor, int marker, int num, bool proximal);
  virtual ~Contig();
  int overlaps(const Contig& other, int minSupportSize, int minOverlapLen, double mismatchRate, int& offset) const;
  bool operator== (const Contig& other) const;
  bool operator< (const Contig& other) const;
  friend std::ostream& operator <<(std::ostream& stream, const Contig& self);
  std::string sequence() const;
  unsigned position() const;
  Locus getAnchor() const;
  int getMarker() const { return marker; }
  bool getProximal() const { return proximal; }
  static bool compare(const Contig& c, unsigned v);
  static bool compare2(unsigned v, const Contig& c);
 private:
  static int overlaps2(const Contig& c1, const Contig& c2, int minSupportSize, int minOverlapLen, double mismatchRate, int& offset);
};

#endif /* CONTIG_H */
