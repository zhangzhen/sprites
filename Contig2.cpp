#include "Contig2.h"
// #include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>

bool Contig2::equals(const std::string& s1, const std::string& s2, int maxMismatches, int& numMismatches) {
  int cnt = 0;
  if (s1.size() != s2.size()) return false;
  // assert(s1.size() == s2.size());
  for (size_t i = 0; i < s1.size(); ++i) {
    if (s1[i] != s2[i]) ++cnt;
    if (cnt > maxMismatches) return false;
  }
  numMismatches = cnt;
  return true;
}

// std::ostream& operator <<(std::ostream& stream, const Contig& self) {
//   stream << self.anchor;
//   // stream << self.marker << std::endl;
//   // stream << self.seq << std::endl;
//   return stream;
// }

// Contig2::Contig2() {}

Contig2::Contig2(int referenceId, int clipPosition, int localClipPosition, const std::string& sequence, int numClips)
    : referenceId(referenceId), clipPosition(clipPosition), localClipPosition(localClipPosition), sequence(sequence), numClips(numClips) {}

Contig2::~Contig2() {}

int Contig2::lowerBoundOfVarLength(const Contig2& other) const {
  assert(clipPosition < other.clipPosition);
  return other.clipPosition - clipPosition;
}

int Contig2::upperBoundOfVarLength(const Contig2& other) const {
  assert(clipPosition < other.clipPosition);
  return other.clipPosition + std::min(length(), other.length()) - localClipPosition - other.localClipPosition - clipPosition;
}

bool Contig2::overlaps(const Contig2& other, int minOverlapLength, double maxMismatchRate, Overlap* overlap) const {
  assert(referenceId == other.referenceId);
  int l1 = length();
  int l2 = other.length();
  int minLen = std::min(l1, l2);
  int minNumClips = std::min(numClips, other.numClips);
  int maxNumClips = std::max(numClips, other.numClips);

  for (int i = minOverlapLength; i < minLen; i++) {
    int maxMismatches = ceil(i * maxMismatchRate / (0.8*minNumClips + 0.2*maxNumClips));
    const std::string& suffix = sequence.substr(l1 - i, i);
    const std::string& prefix = other.sequence.substr(0, i);
    int numMismatches;
    if (equals(suffix, prefix, maxMismatches, numMismatches)) {
      overlap = new Overlap(referenceId, clipPosition, other.clipPosition, numClips, other.numClips, i, numMismatches, i - localClipPosition - other.localClipPosition);
      return true;
    }
  }
  return false;
}

// bool Contig::operator== (const Contig& other) const {
//   return seq == other.seq &&
//       anchor == other.anchor &&
//       marker == other.marker &&
//       num == other.num;
// }

// bool Contig::operator< (const Contig& other) const {
//   return anchor < other.anchor;
// }

// const std::string& Contig2::getSequence() const {
//   return sequence;
// }

size_t Contig2::length() const {
  return sequence.length();
}

int Contig2::getReferenceId() const {
  return referenceId;
}

int Contig2::getClipPosition() const {
  return clipPosition;
}

int Contig2::getLocalClipPosition() const {
  return localClipPosition;
}

// bool Contig::compare(const Contig& c, unsigned v) {
//   return c.position() < v;
// }

// bool Contig::compare2(unsigned v, const Contig& c) {
//   return v < c.position();
// }
