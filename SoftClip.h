#ifndef _SOFTCLIP_H_
#define _SOFTCLIP_H_

#include <string>
#include "Overlap.h"
#include "readgroup.h"
#include "library.h"

class SoftClip {
private:
    int referenceId;

    int genomePosition;

    int clipPosition;
    int orientation;
    int localClipPosition;

    int mateGenomePosition;
    int mateOrientation;
    int insertSize;
    std::string sequence;

    bool isPartiallyMapped() const;

    bool overlaps2(const SoftClip& other, int minOverlapLength, double maxMismatchRate, Overlap& overlap) const;

    ReadGroup *readgroup;

    static const int DiscordantScalar = 5;

public:

    SoftClip(int referenceId, int genomePosition, int clipPosition, int localClipPosition, int orientation, int mateGenomePosition,
             int mateOrientation, int insertSize, const std::string& sequence, ReadGroup *readgroup);

    int length() const;
    int lengthOfLeftPart() const;
    int lengthOfRightPart() const;
    char at(int i) const;

    bool overlaps(const SoftClip& other, int minOverlapLength, double maxMismatchRate, Overlap& overlap) const;
    bool overlapWith(const SoftClip& other, int minOverlapLength, double maxMismatchRate, Overlap& overlap) const;

    static bool compareL(SoftClip* s1, SoftClip* s2);
    static bool compareR(SoftClip* s1, SoftClip* s2);
    static bool compare1(SoftClip* o, const int clipPosition);
    static bool compare2(const int clipPosition, SoftClip* o);

    friend std::ostream& operator <<(std::ostream& stream, const SoftClip& o);
    int getReferenceId() const;
    int getClipPosition() const;
    std::string getSequence() const;

    bool isLeftPartClipped() const;

    bool searchRangeForSpanningPair(const Library& library, int& start, int& end) const;
    void searchRangeForSoftClip(const SoftClip& orig, int minSize, int& start, int& end) const;

    int getInsertSize() const;
    int lengthOfClippedPart() const;

    int getMeanOfInsertSize() const;
    int getSdOfInsertSize() const;

    bool isPartOfSpanningPair() const;
    int getOrientation() const;
    int getMateGenomePosition() const;

    bool isDiscordantInsertSize() const;

};

#endif /* _SOFTCLIP_H_ */
