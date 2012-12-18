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

Contig::Contig(const std::string& seq, const Locus& anchor, int marker, int num)
    : seq(seq), anchor(anchor), marker(marker), num(num) {}

Contig::~Contig() {}

bool Contig::overlaps(const Contig& other, int minSupportSize, int minOverlapLen, double mismatchRate) const {
  if (anchor.chrom() != other.anchor.chrom() ||
      num + other.num < minSupportSize)
    return false;
  int s1, s2;
  if (marker <= other.marker) {
    s1 = 0;
    s2 = other.marker - marker;
  } else {
    s1 = marker - other.marker;
    s2 = 0;
  }
  int len = std::min(marker, other.marker) + std::min(seq.size() - marker, other.seq.size() - other.marker);
  if (len < minOverlapLen) return false;
  int mismatches = (int)ceil(mismatchRate*len);
  return equals2(seq.substr(s1, len), other.seq.substr(s2, len), mismatches);
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
