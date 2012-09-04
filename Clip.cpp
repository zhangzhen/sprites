#include "Clip.h"

Clip::Clip(int refId, int refPosition, int readPosition, int size, const std::string& readSeq)
    :refId(refId),
     refPosition(refPosition),
     readPosition(readPosition),
     size(size),
     readSeq(readSeq) {
}

Clip::~Clip() {
}

ClipType LeftClip::getType() {
  return Left;
}

ClipType RightClip::getType() {
  return Right;
}
