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

bool SoftClip::compareL(SoftClip* s1, SoftClip* s2) {
  if (s1->pos != s2->pos)
    return s1->pos < s2->pos;
  return s1->lengthOfLeftPart() < s2->lengthOfLeftPart();
}

bool SoftClip::compareR(SoftClip* s1, SoftClip* s2) {
  if (s1->pos != s2->pos)
    return s1->pos < s2->pos;
  return s1->lengthOfRightPart() > s2->lengthOfRightPart();
}

int SoftClip::minDeletionLength(const SoftClip& other) const {
  return std::max(other.pos - pos, 0);
}

int SoftClip::maxDeletionLength(const SoftClip& other) const {
  return other.pos - pos + std::min(lengthOfLeftPart(), other.lengthOfRightPart());
}

bool SoftClip::overlaps(const SoftClip& other, int minOverlapLength, double maxMismatchRate, Overlap& overlap) const {
  assert(refId == other.refId);
  if (maxDeletionLength(other) < 0) return false;

  int l1 = length();
  int l2 = other.length();

  // if (pos == 33115831) std::cout << lengthOfLeftPart() << '\t' << other.lengthOfLeftPart() << std::endl;

  for (int i = std::max(0, minOverlapLength - lengthOfRightPart() - other.lengthOfLeftPart()); i < std::min(lengthOfLeftPart(), other.lengthOfRightPart()); ++i) {
    int n1 = std::min(lengthOfLeftPart(), other.lengthOfLeftPart() + i);
    int n2 = std::min(lengthOfRightPart(), other.lengthOfRightPart() - i);
    int n = n1 + n2;
    int st1 = (lengthOfLeftPart() > other.lengthOfLeftPart() + i) ? lengthOfLeftPart() - other.lengthOfLeftPart() - i : 0;
    int st2 = (lengthOfLeftPart() > other.lengthOfLeftPart() + i) ? 0 : other.lengthOfLeftPart() + i - lengthOfLeftPart();
    int maxMismatches = ceil(n * maxMismatchRate);
    int numMismatches;
    if (Overlap::equals(seq.substr(st1, n), other.seq.substr(st2, n), maxMismatches, numMismatches)) {
      overlap = Overlap(this, &other, n, numMismatches, i);
      // if (pos == 30255998) std::cout << overlap << std::endl;
      return true;
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

std::ostream& operator <<(std::ostream& stream, const SoftClip& o) {
  stream << o.refId << ":" << o.pos << std::endl;
  stream << std::string(o.lengthOfLeftPart(), ' ') << '+' << std::endl;
  stream << o.seq << std::endl;
  return stream;
}
