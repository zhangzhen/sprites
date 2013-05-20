#include "Consensus.h"
// #include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>

bool Consensus::equals(const std::string& s1, const std::string& s2, int maxMismatches, int& numMismatches) {
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

Consensus::Consensus(int referenceId, int clipPosition, int localClipPosition, const std::string& sequence, int numClips)
    : referenceId(referenceId), clipPosition(clipPosition), localClipPosition(localClipPosition), sequence(sequence), numClips(numClips) {
}

int Consensus::minDeletionLength(const Consensus& other) const {
  assert(clipPosition < other.clipPosition);
  return other.clipPosition - clipPosition;
}

int Consensus::maxDeletionLength(const Consensus& other) const {
  assert(clipPosition < other.clipPosition);
  return other.clipPosition + std::min(length(), other.length()) - localClipPosition - other.localClipPosition - clipPosition;
}

bool Consensus::overlaps(const Consensus& other, int minOverlapLength, double maxMismatchRate, Overlap& overlap) const {
  assert(referenceId == other.referenceId);
  int l1 = length();
  int l2 = other.length();
  int initClippedLength = l1 - localClipPosition + other.localClipPosition;
  int start = std::max(minOverlapLength, initClippedLength);
  int end = std::min(l1, l2);

  for (int i = start; i <= end; i++) {
    int maxMismatches = ceil(i * maxMismatchRate / (0.8*std::min(numClips, other.numClips) + 0.2*std::max(numClips, other.numClips)));
    const std::string& suffix = sequence.substr(l1 - i);
    const std::string& prefix = other.sequence.substr(0, i);
    int numMismatches;
    if (equals(suffix, prefix, maxMismatches, numMismatches)) {
      overlap = Overlap(referenceId, clipPosition, other.clipPosition, numClips, other.numClips, i, numMismatches, i - initClippedLength);
      return true;
    }
  }
  return false;
}

// bool Consensus::operator== (const Consensus& other) const {
//   return seq == other.seq &&
//       anchor == other.anchor &&
//       marker == other.marker &&
//       num == other.num;
// }

bool Consensus::operator< (const Consensus& other) const {
  return clipPosition < other.clipPosition;
}

// const std::string& Consensus::getSequence() const {
//   return sequence;
// }

size_t Consensus::length() const {
  return sequence.length();
}

int Consensus::getReferenceId() const {
  return referenceId;
}

int Consensus::getClipPosition() const {
  return clipPosition;
}

int Consensus::getLocalClipPosition() const {
  return localClipPosition;
}

bool Consensus::compare(const Consensus& cons, const int pos) {
  return cons.clipPosition < pos;
}

bool Consensus::compare2(const int pos, const Consensus& cons) {
  return pos < cons.clipPosition;
}
