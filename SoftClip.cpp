#include "SoftClip.h"
#include <cassert>
#include <algorithm>
#include <iostream>


size_t SoftClip::getClippedSize() const
{
    return clippedSize;
}

SoftClip::SoftClip() {
}

SoftClip::SoftClip(int referenceId,
                   int position,
                   int leftmostPosition,
                   int clipPosition,
                   int matePosition,
                   bool is_reverse,
                   bool is_mate_reverse,
                   size_t clippedSize,
                   const std::string &sequence) :
    referenceId(referenceId),
    position(position),
    leftmostPosition(leftmostPosition),
    clipPosition(clipPosition),
    matePosition(matePosition),
    is_reverse(is_reverse),
    is_mate_reverse(is_mate_reverse),
    clippedSize(clippedSize),
    sequence(sequence)
{
    assert(clippedSize < sequence.size());
}

int SoftClip::getLeftmostPosition() const
{
    return leftmostPosition;
}

bool SoftClip::isReverse() const
{
    return is_reverse;
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

int SoftClip::getClipPositionInRead() const {
    if (isForRightBp()) return clippedSize;
    return size() - clippedSize;
}

const std::string& SoftClip::getSequence() const { return sequence; }

size_t SoftClip::size() const { return sequence.size(); }

std::ostream& operator <<(std::ostream& stream, const SoftClip& o) {
  stream << o.referenceId << ":" << o.clipPosition << std::endl;
  stream << std::string(o.clippedSize, ' ') << '+' << std::endl;
  stream << o.sequence << std::endl;
  return stream;
}
