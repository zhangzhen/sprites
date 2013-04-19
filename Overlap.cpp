#include "Overlap.h"

Overlap::Overlap(int referenceId, int clipPosition1, int clipPosition2, int numClips1, int numClips2, int length, int numMismatches, int offset):
    referenceId(referenceId),
    clipPosition1(clipPosition1),
    clipPosition2(clipPosition2),
    numClips1(numClips1),
    numClips2(numClips2),
    length(length),
    numMismatches(numMismatches),
    offset(offset) {}

Overlap::~Overlap() {}

// int Overlap::getLength() const { return length; }

// int Overlap::getNumMismatches() const { return numMismatches; }

int Overlap::regionLength() const {
  return clipPosition2 + offset - clipPosition1;
}

double Overlap::score() const {
  return 0;
}

const RegionX* Overlap::primaryRegion() const {
  return new RegionX(referenceId, clipPosition1, clipPosition2 + offset);
}

const RegionX* Overlap::secondaryRegion() const {
  return new RegionX(referenceId,clipPosition1 - offset, clipPosition2);
}
