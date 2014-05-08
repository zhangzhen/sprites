#include "SoftClip.h"
#include <cassert>
#include <algorithm>
#include <iostream>

SoftClip::SoftClip() {
}

SoftClip::SoftClip(int referenceId,
                   int position,
                   int leftmostPosition,
                   int clipPosition,
                   int matePosition,
                   bool is_reverse,
                   bool is_mate_reverse,
                   int clippedSize,
                   int offset,
                   const std::string &sequence) :
    referenceId(referenceId),
    position(position),
    leftmostPosition(leftmostPosition),
    clipPosition(clipPosition),
    matePosition(matePosition),
    is_reverse(is_reverse),
    is_mate_reverse(is_mate_reverse),
    clippedSize(clippedSize),
    offset(offset),
    sequence(sequence)
{
    assert(clippedSize < sequence.size());
}

// todo
int SoftClip::getAdjustedPosition() const
{
    return clipPosition + offset;
}

int SoftClip::getLeftmostPosition() const
{
    return leftmostPosition;
}

bool SoftClip::isReverse() const
{
    return is_reverse;
}

bool SoftClip::isLeft() const
{
    return position == clipPosition;
}

bool SoftClip::isTypeIForLeftBp() const
{
    return position != clipPosition && is_reverse && !is_mate_reverse && position > matePosition;
}

bool SoftClip::isTypeIIForLeftBp() const
{
    return position != clipPosition && !is_reverse && is_mate_reverse && position < matePosition;
}

bool SoftClip::isTypeIForRightBp() const
{
    return position == clipPosition && !is_reverse && is_mate_reverse && position < matePosition;
}

bool SoftClip::isTypeIIForRightBp() const
{
    return position == clipPosition && is_reverse && !is_mate_reverse && position > matePosition;
}

bool SoftClip::isForLeftBp() const
{
    return isTypeIForLeftBp() || isTypeIIForLeftBp();
}

bool SoftClip::isForRightBp() const
{
    return isTypeIForRightBp() || isTypeIIForRightBp();
}

int SoftClip::getReferenceId() const { return referenceId; }

int SoftClip::getPosition() const { return position; }

int SoftClip::getClipPosition() const
{
    return clipPosition;
}

const std::string& SoftClip::getSequence() const { return sequence; }

int SoftClip::size() const { return sequence.size(); }

int SoftClip::getClippedSize() const { return clippedSize; }

std::ostream& operator <<(std::ostream& stream, const SoftClip& o) {
  stream << o.referenceId << ":" << o.clipPosition << std::endl;
  stream << std::string(o.clippedSize, ' ') << '+' << std::endl;
  stream << o.sequence << std::endl;
  return stream;
}
