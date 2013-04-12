#include "SingleClipped2.h"

SingleClipped2::SingleClipped2(int refId, int pos, int clipPos, const std::string& seq, const std::string& quals)
    : refId(refId),
      pos(pos),
      clipPos(clipPos),
      seq(seq),
      quals(quals) {
}

SingleClipped2::~SingleClipped2() {}

int SingleClipped2::referenceId() const { return refId; }

int SingleClipped2::position() const { return pos; }

const std::string& SingleClipped2::sequence() const { return seq; }

int SingleClipped2::length() const { return seq.length(); }

int SingleClipped2::lengthOfLeftPart() const { return clipPos; }

int SingleClipped2::lengthOfRightPart() const { return length() - clipPos; }

char SingleClipped2::at(int i) const { return seq[i]; }

char SingleClipped2::qual(int i) const { return quals[i]; }
