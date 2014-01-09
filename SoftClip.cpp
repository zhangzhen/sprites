#include "SoftClip.h"
#include "DFinderHelper.h"
#include <cassert>
#include <algorithm>
#include <iostream>

using namespace std;

int SoftClip::getReferenceId() const
{
    return referenceId;
}


int SoftClip::getClipPosition() const
{
    return clipPosition;
}


std::string SoftClip::getSequence() const
{
    return sequence;
}

bool SoftClip::searchRangeForSpanningPair(const Library& library, int &start, int &end) const
{
    int readLength = library.getReadLength();
    int mean = library.getMeanOfInsertSize();
    int sd = library.getSdOfInsertSize();

    if (isPartiallyMapped() && isLeftPartClipped() && orientation)
    {
        start = clipPosition;
        end = clipPosition - 2 * readLength + mean + 3 * sd;
        return true;
    }
    if (isPartiallyMapped() && !isLeftPartClipped() && !orientation)
    {
        start = clipPosition + readLength - mean - 3 * sd;
        end = clipPosition - readLength;
        return true;
    }
    return false;
}


void SoftClip::searchRangeForSoftClip(const SoftClip &orig, int minSize, int& start, int& end) const
{
    start = mateOrientation ? mateGenomePosition + length() - orig.lengthOfClippedPart() :
                              max(clipPosition + minSize,
                                  mateGenomePosition - getMeanOfInsertSize() - orig.lengthOfClippedPart() + 2 * length() - 3 *getSdOfInsertSize());
    end = mateOrientation ? min(mateGenomePosition + getMeanOfInsertSize() - orig.lengthOfClippedPart() - length() + 3 * getSdOfInsertSize(),
                                clipPosition - minSize - orig.lengthOfClippedPart()) :
                            mateGenomePosition;
}

int SoftClip::getInsertSize() const
{
    return insertSize;
}

bool SoftClip::isLeftPartClipped() const
{
    return genomePosition == clipPosition;
}


int SoftClip::getOrientation() const
{
    return orientation;
}


int SoftClip::getMateGenomePosition() const
{
    return mateGenomePosition;
}

bool SoftClip::isPartiallyMapped() const
{
    return lengthOfClippedPart() != 0;
}

int SoftClip::lengthOfClippedPart() const
{
    return isLeftPartClipped() ? localClipPosition : length() - localClipPosition;
}

int SoftClip::getMeanOfInsertSize() const
{
    return readgroup->getMeanOfInsertSize();
}

int SoftClip::getSdOfInsertSize() const
{
    return readgroup->getSdOfInsertSize();
}

bool SoftClip::isPartOfSpanningPair() const
{
    return isDiscordantInsertSize() &&
            ((isLeftPartClipped() && !orientation) ||
             (!isLeftPartClipped() && orientation));
}

bool SoftClip::isDiscordantInsertSize() const
{
    return insertSize > getMeanOfInsertSize() + DiscordantScalar * getSdOfInsertSize();
}

SoftClip::SoftClip(int referenceId, int genomePosition, int clipPosition, int localClipPosition, int orientation,
                   int mateGenomePosition, int mateOrientation, int insertSize, const std::string &sequence, ReadGroup *readgroup)
    : referenceId(referenceId),
      genomePosition(genomePosition),
      clipPosition(clipPosition),
      localClipPosition(localClipPosition),
      orientation(orientation),
      mateGenomePosition(mateGenomePosition),
      mateOrientation(mateOrientation),
      insertSize(abs(insertSize)),
      sequence(sequence),
      readgroup(readgroup)
{}

int SoftClip::length() const { return sequence.length(); }

int SoftClip::lengthOfLeftPart() const { return localClipPosition; }

int SoftClip::lengthOfRightPart() const { return length() - localClipPosition; }

char SoftClip::at(int i) const { return sequence[i]; }

bool SoftClip::overlaps(const SoftClip &other, int minOverlapLength, double maxMismatchRate, Overlap &overlap) const
{
    if (isLeftPartClipped()) return other.overlaps2(*this, minOverlapLength, maxMismatchRate, overlap);
    return overlaps2(other, minOverlapLength, maxMismatchRate, overlap);
}

bool SoftClip::compareL(SoftClip* s1, SoftClip* s2) {
    if (s1->clipPosition != s2->clipPosition)
        return s1->clipPosition < s2->clipPosition;
    return s1->lengthOfLeftPart() < s2->lengthOfLeftPart();
}

bool SoftClip::compareR(SoftClip* s1, SoftClip* s2) {
    if (s1->clipPosition != s2->clipPosition)
        return s1->clipPosition < s2->clipPosition;
    return s1->lengthOfRightPart() > s2->lengthOfRightPart();
}

bool SoftClip::overlaps2(const SoftClip& other, int minOverlapLength, double maxMismatchRate, Overlap& overlap) const {
    assert(referenceId == other.referenceId);

    int l1 = length();
    int l2 = other.length();

    int mn = max(max(lengthOfClippedPart(), other.lengthOfClippedPart()), minOverlapLength);
    int mx = min(l1, l2);

    for (int i = mx; i >= mn; i--) {
        int maxMismatches = ceil(i * maxMismatchRate);
        int numMismatches;
        if (Overlap::equals(sequence.substr(l1 - i, i), other.sequence.substr(0, i), maxMismatches, numMismatches)) {
            overlap = Overlap(this, &other, i, numMismatches, i - lengthOfClippedPart() - other.lengthOfClippedPart());
            return true;
        }
    }

    return false;
}

bool SoftClip::overlapWith(const SoftClip& other, int minOverlapLength, double maxMismatchRate, Overlap& overlap) const {
    int mismatches, len;
    if (::overlaps2(sequence, other.sequence, maxMismatchRate, len, mismatches) &&
            len >= minOverlapLength) {
        int offset = len - lengthOfRightPart() - other.lengthOfLeftPart();
        overlap = Overlap(this, &other, len, mismatches, offset);
        return true;
    }
    return false;
}

bool SoftClip::compare1(SoftClip* o, const int pos) {
    return o->clipPosition < pos;
}

bool SoftClip::compare2(const int pos, SoftClip* o) {
    return pos < o->clipPosition;
}

std::ostream& operator <<(std::ostream& stream, const SoftClip& o) {
    stream << o.referenceId << ":" << o.clipPosition << std::endl;
    stream << std::string(o.lengthOfLeftPart(), ' ') << '+' << std::endl;
    stream << o.sequence << std::endl;
    return stream;
}
