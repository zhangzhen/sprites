#include "SoftClip.h"

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
