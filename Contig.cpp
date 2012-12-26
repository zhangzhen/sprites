#include "Contig.h"
#include <math.h>
#include <assert.h>
#include <iostream>

bool equals2(const std::string s1, const std::string s2, int mismatches) {
  int cnt = 0;
  if (s1.size() != s2.size()) return false;
  // assert(s1.size() == s2.size());
  for (size_t i = 0; i < s1.size(); ++i) {
    if (s1[i] != s2[i]) ++cnt;
    if (cnt > mismatches) return false;
  }
  return true;
}

std::ostream& operator <<(std::ostream& stream, const Contig& self) {
  stream << self.anchor;
  // stream << self.marker << std::endl;
  // stream << self.seq << std::endl;
  return stream;
}

Contig::Contig(const std::string& seq, const Locus& anchor, int marker, int num, bool proximal)
    : seq(seq), anchor(anchor), marker(marker), num(num), proximal(proximal) {}

Contig::~Contig() {}

bool Contig::overlaps(const Contig& other, int minSupportSize, int minOverlapLen, double mismatchRate) const {
  assert(proximal ^ other.proximal);
  if (proximal) return Contig::overlaps2(*this, other, minSupportSize, minOverlapLen, mismatchRate);
  return Contig::overlaps2(other, *this, minSupportSize, minOverlapLen, mismatchRate);
}

bool Contig::overlaps2(const Contig& c1, const Contig& c2, int minSupportSize, int minOverlapLen, double mismatchRate) {
  if (c1.anchor.chrom() != c2.anchor.chrom() ||
      c1.num + c2.num < minSupportSize) return false;
  int m = c2.marker;
  int s1, s2;
  while (m <= c2.seq.size()) {
    if (c1.marker <= m) {
      s1 = 0;
      s2 = m - c1.marker;
    } else {
      s1 = c1.marker - m;
      s2 = 0;
    }
    int len = std::min(c1.marker, m) + std::min(c1.seq.size() - c1.marker, c2.seq.size() - m);
    m++;
    if (len < minOverlapLen) continue;
    int mismatches = (int)ceil(mismatchRate * len);
    std::cout << c1.seq.substr(s1, len) << std::endl;
    std::cout << c2.seq.substr(s2, len) << std::endl << std::endl;    
    if (equals2(c1.seq.substr(s1, len), c2.seq.substr(s2, len), mismatches)) return true;
  }
  return false;
}

bool Contig::operator== (const Contig& other) const {
  return seq == other.seq &&
      anchor == other.anchor &&
      marker == other.marker &&
      num == other.num;
}

bool Contig::operator< (const Contig& other) const {
  return anchor < other.anchor;
}

std::string Contig::sequence() const {
  return seq;
}

Locus Contig::getAnchor() const {
  return anchor;
}
