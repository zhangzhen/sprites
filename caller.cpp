#include "caller.h"
#include <cassert>
#include <numeric>

#include "Thirdparty/overlapper.h"

using namespace std;
using namespace BamTools;

Caller::Caller(const string &filename,
               int minOverlap,
               double minIdentity,
               int insertMean,
               int insertStd) :
    minOverlap(minOverlap),
    minIdentity(minIdentity),
    insertMean(insertMean),
    insertStd(insertStd)
{
    if (!reader.Open(filename))
        error("Could not open the input BAM file.");
}

Caller::~Caller()
{
    reader.Close();
}

Deletion Caller::call(const SoftClip &clip)
{
    vector<Deletion> deletions;

    if ((!clip.isReverse() && clip.isLeft()) || (clip.isReverse() && !clip.isLeft()))
    {
        vector<int> matePositions;
        getSuppMatePositions(clip, matePositions);
        vector<TargetRegion> regions;
        getTargetRegions(clip, matePositions, regions);
        for (auto itr = regions.begin(); itr != regions.end(); ++itr)
        {
            deletions.push_back(call(clip, *itr));
        }
    }
}

Deletion Caller::call(const SoftClip &clip, const TargetRegion &region)
{
    MultipleAlignment ma = buildMultipleAlignment(query, region);
    return getDeletion(ma);
}

void Caller::getSuppMatePositions(const SoftClip &clip, vector<int> &matePositions)
{
    int start, end;
    if (!clip.isReverse())
    {
        start = clip.getPosition() + clip.getClippedSize();
        end = start + insertMean + 3 * insertStd + clip.size();
    }
    else
    {
        end = clip.getPosition() + clip.size();
        start = end > insertMean + 3 * insertStd + clip.size() ? end - insertMean - 3 * insertStd - clip.size() :
                                                                 0;
    }
    assert(start < end);
    reader.SetRegion(clip.getReferenceId(), start, clip.getReferenceId(), end);

    BamAlignment al;
    while(reader.GetNextAlignmentCore(al))
    {
        if ((!clip.isReverse() && al.IsReverseStrand() && al.Position > al.MatePosition) ||
                (clip.isReverse() && !al.IsMateReverseStrand() && al.Position < al.MatePosition))
            matePositions.push_back(al.MatePosition);
    }
}

TargetRegion Caller::getTargetRegion(const SoftClip& clip, int matePosition)
{
    int start = !clip.isReverse() ? matePosition :
                                    matePosition - 3 * insertStd;
    int end = !clip.isReverse() ? matePosition + clip.size() + insertMean + 3 * insertStd :
                                  matePosition + clip.size();
    return {clip.getReferenceId(), start, end};
}

void Caller::getTargetRegions(const SoftClip& clip,
                              std::vector<int> &matePositions,
                              std::vector<TargetRegion> &regions)
{
    sort(matePositions.begin(), matePositions.end());

    if (matePositions.size() == 1)
    {
        regions.push_back(getTargetRegion(clip, matePositions[0]));
        return;
    }

    std::vector<int> diffs(matePositions.size());

    adjacent_difference(matePositions.begin(), matePositions.end(), diffs);

    for (int i = 1; i < diffs.size(); ++i)
    {
        if (abs(diffs[i]) <= clip.size() + insertMean + 3 * insertStd)
            regions.pop_back();
        regions.push_back(getTargetRegion(clip, matePositions[i]));
    }
}

MultipleAlignment Caller::buildMultipleAlignment(const string &query, const TargetRegion &region)
{
    SequenceOverlapPairVector overlap_vector;
    retrieveMatches(query, region, overlap_vector);

    MultipleAlignment multiple_alignment;
    multiple_alignment.addBaseSequence("query", query, "");
    for(size_t i = 0; i < overlap_vector.size(); ++i)
        multiple_alignment.addOverlap("null", overlap_vector[i].sequence[1], "", overlap_vector[i].overlap);
    return multiple_alignment;
}

void Caller::retrieveMatches(const string &query, const TargetRegion &region, SequenceOverlapPairVector &result)
{
    reader.SetRegion(region.referenceId, region.start, region.referenceId, region.end);
    BamAlignment al;
    while (reader.GetNextAlignment(al))
    {
        SequenceOverlap overlap = Overlapper::computeOverlap(query, al.QueryBases, ungapped_params);
        if (overlap.getOverlapLength() >= minOverlap &&
                overlap.getPercentIdentity() / 100 >= minIdentity)
        {
            SequenceOverlapPair op;
            op.sequence[0] = query;
            op.sequence[1] = al.QueryBases;
            op.overlap = overlap;
            result.push_back(op);
        }
    }

}

