#include "Caller.h"
#include "error.h"

#include <cassert>
#include <numeric>
#include <algorithm>

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
    if (!reader.LocateIndex())
        error("Could not locate index.");
}

Caller::~Caller()
{
    reader.Close();
}

bool Caller::call(const SoftClip &clip, Deletion &del)
{
    bool result = false;
    vector<Deletion> deletions;

    if ((!clip.isReverse() && clip.isLeft()) || (clip.isReverse() && !clip.isLeft()))
    {
        vector<int> matePositions;
        if (!getSuppMatePositions(clip, matePositions))
            return false;
        vector<TargetRegion> regions;
        getTargetRegions(clip, matePositions, regions);
        for (auto itr = regions.begin(); itr != regions.end(); ++itr)
        {
            Deletion d;
            if (call(clip, *itr, d)) {
                deletions.push_back(d);
                result = true;
            }
        }
    }
    return result;
}

bool Caller::call(const SoftClip &clip, const TargetRegion &region, Deletion &del)
{
    MultipleAlignment ma = buildMultipleAlignment(clip, region);
    if (ma.getNumRows() > 1)
        ma.print(110);
    return false;
}

bool Caller::getSuppMatePositions(const SoftClip &clip, vector<int> &matePositions)
{
    bool result = false;

    int start, end;
    if (!clip.isReverse())
    {
        start = clip.getLeftmostPosition();
        end = start + insertMean + 3 * insertStd + clip.size();
    }
    else
    {
        end = clip.getLeftmostPosition() + clip.size();
        start = end > insertMean + 3 * insertStd + clip.size() ? end - insertMean - 3 * insertStd - clip.size() :
                                                                 0;
    }
    assert(start < end);
    if (!reader.SetRegion(clip.getReferenceId(), start, clip.getReferenceId(), end)) {
        cerr << "Could not set the region.";
        return false;
    }

    BamAlignment al;
    while(reader.GetNextAlignmentCore(al))
    {
        if ((!clip.isReverse() && al.IsReverseStrand() && al.Position > al.MatePosition) ||
                (clip.isReverse() && !al.IsMateReverseStrand() && al.Position < al.MatePosition)) {
            matePositions.push_back(al.MatePosition);
            result = true;
        }
    }
    return result;
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

    regions.push_back(getTargetRegion(clip, matePositions[0]));
    if (matePositions.size() == 1)
    {
        return;
    }

    std::vector<int> diffs(matePositions.size());

    adjacent_difference(matePositions.begin(), matePositions.end(), diffs.begin());

    for (size_t i = 1; i < diffs.size(); ++i) {
        int start;
        TargetRegion tr = getTargetRegion(clip, matePositions[i]);
        if (abs(diffs[i]) <= clip.size() + insertMean + 3 * insertStd) {
            start = regions.back().start;
            regions.pop_back();
        } else {
            start = tr.start;
        }
        regions.push_back({tr.referenceId, start, tr.end});
    }
}

MultipleAlignment Caller::buildMultipleAlignment(const SoftClip &clip, const TargetRegion &region)
{
    SequenceOverlapPairVector overlap_vector;
    retrieveMatches(clip, region, overlap_vector);

    MultipleAlignment multiple_alignment;
    multiple_alignment.addBaseSequence("query", clip.getSequence(), "", clip.getClipPosition());
    for(size_t i = 0; i < overlap_vector.size(); ++i)
        multiple_alignment.addOverlap("null", overlap_vector[i].sequence[1], "", overlap_vector[i].position[1], overlap_vector[i].overlap);
    return multiple_alignment;
}

void Caller::retrieveMatches(const SoftClip &clip, const TargetRegion &region, SequenceOverlapPairVector &result)
{
    assert(region.start < region.end);
    if(!reader.SetRegion(region.referenceId, region.start, region.referenceId, region.end)) {
        cerr << "Could not set the region.";
        return;
    }

    BamAlignment al;
    string query = clip.getSequence();
    while (reader.GetNextAlignment(al))
    {
        if ((clip.isLeft() && !clip.isReverse() && al.Position + al.Length >= clip.getClipPosition()) ||
                (!clip.isLeft() && clip.isReverse() && al.Position <= clip.getClipPosition()))
            break;
        SequenceOverlap overlap = Overlapper::computeOverlap(query, al.QueryBases, ungapped_params);
        if (overlap.getOverlapLength() >= minOverlap &&
                overlap.getPercentIdentity() / 100 >= minIdentity)
        {
            SequenceOverlapPair op;
            op.sequence[0] = query;
            op.sequence[1] = al.QueryBases;
            op.position[0] = clip.getClipPosition();
            op.position[1] = al.Position;
            op.overlap = overlap;
            result.push_back(op);
        }
    }

}

