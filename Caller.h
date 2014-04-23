#ifndef CALLER_H
#define CALLER_H

#include "api/BamReader.h"
#include "SoftClip.h"
#include "Deletion.h"
#include "Thirdparty/multiple_alignment.h"

#include <string>
#include <vector>

struct TargetRegion
{
    int referenceId;
    int start;
    int end;
};

struct SequenceOverlapPair
{
    std::string sequence[2];
    SequenceOverlap overlap;

    static bool sortByOverlapLengthDesc(const SequenceOverlapPair& a, const SequenceOverlapPair& b) { return a.overlap.getOverlapLength() > b.overlap.getOverlapLength(); }

};

typedef std::vector<SequenceOverlapPair> SequenceOverlapPairVector;

class Caller
{
public:
    Caller(const std::string& filename,
           int minOverlap,
           double minIdentity,
           int insertMean,
           int insertStd);
    virtual ~Caller();

    bool call(const SoftClip& clip, Deletion& del);

private:
    bool call(const SoftClip& clip, const TargetRegion& region, Deletion& del);

    bool getSuppMatePositions(const SoftClip& clip, std::vector<int>& matePositions);
    TargetRegion getTargetRegion(const SoftClip &clip,
                         int matePosition);
    void getTargetRegions(const SoftClip &clip,
                          std::vector<int>& matePositions,
                          std::vector<TargetRegion>& regions);

    MultipleAlignment buildMultipleAlignment(const std::string& query, const TargetRegion& region);
    void retrieveMatches(const std::string& query, const TargetRegion& region, SequenceOverlapPairVector& result);

    BamTools::BamReader reader;
    int minOverlap;
    double minIdentity;
    int insertMean;
    int insertStd;
};

#endif // CALLER_H
