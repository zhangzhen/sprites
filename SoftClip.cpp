#include "SoftClip.h"
#include <cassert>
#include <algorithm>
#include <iostream>

SoftClip::SoftClip(int refId, int pos, int clipPos, const std::string& seq, const std::string& quals)
    : refId(refId),
      pos(pos),
      clipPos(clipPos),
      seq(seq),
      quals(quals) {
}

int SoftClip::referenceId() const { return refId; }

int SoftClip::position() const { return pos; }

const std::string& SoftClip::sequence() const { return seq; }

int SoftClip::length() const { return seq.length(); }

int SoftClip::lengthOfLeftPart() const { return clipPos; }

int SoftClip::lengthOfRightPart() const { return length() - clipPos; }

char SoftClip::at(int i) const { return seq[i]; }

char SoftClip::qual(int i) const { return quals[i]; }

bool SoftClip::compare(SoftClip* s1, SoftClip* s2) {
  return (*s1).position() < (*s2).position();
}

int SoftClip::minDeletionLength(const SoftClip& other) const {
  assert(pos <= other.pos);
  return other.pos - pos;
}

int SoftClip::maxDeletionLength(const SoftClip& other) const {
  assert(pos < other.pos);
  return other.pos + std::min(length(), other.length()) - clipPos - other.clipPos - pos;
}

bool SoftClip::overlaps(const SoftClip& other, int minOverlapLength, double maxMismatchRate, Overlap& overlap) const {
  assert(refId == other.refId);
  int l1 = length();
  int l2 = other.length();

  if (lengthOfLeftPart() >= other.lengthOfLeftPart()) {
    int initClippedLength = lengthOfRightPart() + other.lengthOfLeftPart();
    int start = std::max(minOverlapLength, initClippedLength);
    int end = std::min(l1, l2);
    for (int i = start; i <= end; i++) {
      int maxMismatches = ceil(i * maxMismatchRate);
      const std::string& suffix = seq.substr(l1 - i);
      const std::string& prefix = other.seq.substr(0, i);
      int numMismatches;
      if (Overlap::equals(suffix, prefix, maxMismatches, numMismatches)) {
        overlap = Overlap(this, &other, i, numMismatches, i - initClippedLength);
        return true;
      }
    }
  } else {
    int start = lengthOfLeftPart() + other.lengthOfRightPart();
    int end = std::max(minOverlapLength,std::max(lengthOfLeftPart(), other.lengthOfRightPart()));
    for (int i = start; i >= end; i--) {
      int maxMismatches = ceil(i * maxMismatchRate);
      const std::string& prefix = seq.substr(0, i);
      const std::string& suffix = other.seq.substr(l2 - i);
      int numMismatches;
      if (Overlap::equals(suffix, prefix, maxMismatches, numMismatches)) {
        overlap = Overlap(this, &other, i, numMismatches, start - i);
        return true;
      }
    }
  }
  
  return false;
}

bool SoftClip::compare1(SoftClip* o, const int pos) {
  return o->pos < pos;
}

bool SoftClip::compare2(const int pos, SoftClip* o) {
  return pos < o->pos;
}
