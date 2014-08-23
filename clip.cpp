#include "clip.h"
#include "error.h"
#include "Thirdparty/overlapper.h"
#include "Helper.h"

#include <algorithm>

using namespace std;
using namespace BamTools;

AbstractClip::AbstractClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const string &sequence, const vector<CigarOp>& cigar)
    : referenceId(referenceId),
      mapPosition(mapPosition),
      clipPosition(clipPosition),
      matePosition(matePosition),
      sequence(sequence),
      cigar(cigar) {
}

int AbstractClip::length() const {
    return sequence.length();
}

int AbstractClip::leftmostPosition() const {
    if (cigar[0].Type == 'S') return mapPosition - cigar[0].Length;
    return mapPosition;
}

AbstractClip::~AbstractClip() {
}

Deletion AbstractClip::call(BamReader &reader, FaidxWrapper &faidx, int insLength, int minOverlap, double minIdentity)
{
    return call(faidx, tRegions(reader, insLength), minOverlap, minIdentity);
}

void AbstractClip::tRegions(BamReader &reader, int insLength, std::vector<TargetRegion> &regions) {

    vector<int> anchors;
    fetchAnchors(reader, insLength, anchors);

    if (anchors.empty()) return;

    sort(anchors.begin(), anchors.end());

    if (anchors.size() == 1) {
        regions.push_back(tRegion(anchors[0], insLength));
        return;
    }

    std::vector<int> diffs(anchors.size());

    adjacent_difference(anchors.begin(), anchors.end(), diffs.begin());

    for (size_t i = 1; i < diffs.size(); ++i) {
        int start;
        TargetRegion tr = tRegion(Helper::getReferenceName(referenceId), anchors[i], insLength);
        if (abs(diffs[i]) <= length() + insLength) {
            start = regions.back().start;
            regions.pop_back();
        } else {
            start = tr.start;
        }
        regions.push_back({tr.referenceName, start, tr.end});
    }

}


ForwardBClip::ForwardBClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const string &sequence, const vector<CigarOp>& cigar)
    : AbstractClip(referenceId, mapPosition, clipPosition, matePosition, sequence, cigar) {
}

Deletion ForwardBClip::call(FaidxWrapper &faidx, const std::vector<TargetRegion> &regions, int minOverlap, double minIdentity)
{
    for (auto it = regions.begin(); it != regions.end(); ++it) {
        SequenceOverlap overlap = Overlapper::computeOverlapSW((*it).sequence(faidx), sequence, ungapped_params);
        if (overlap.getOverlapLength() > minOverlap &&
                overlap.getPercentIdentity() >= minIdentity &&
                overlap.match[1].start == 0) {
            int rightBp = clipPosition - 1;
            int leftBp = (overlap.getOverlapLength()  > cigar[0].Length) ? (*it).start + overlap.match[0].start + cigar[0].Length -1
                    : (*it).start + overlap.match[0].end;
            int len = leftBp - rightBp;
            if (overlap.getOverlapLength() < cigar[0].Length) len += cigar[0].Length - overlap.getOverlapLength();
            return Deletion(referenceId,leftBp, rightBp, len);
        }
    }
}

TargetRegion ForwardBClip::tRegion(const string &referenceName, int anchor, int insLength) {
    return { referenceName,
              anchor,
              anchor + length() + insLength };
}

void ForwardBClip::fetchAnchors(BamReader &reader, int insLength, std::vector<int> &anchors) {

    int start = leftmostPosition();
    int end = start + insLength + length();

    if (!reader.SetRegion(referenceId, start, referenceId, end))
        error("Could not set the region.");

    BamAlignment al;
    while(reader.GetNextAlignmentCore(al)) {
        if (al.IsReverseStrand() && !al.IsMateReverseStrand() && al.Position > al.MatePosition) {
            anchors.push_back(al.MatePosition);
        }
    }

}


ReverseBClip::ReverseBClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const string &sequence, const std::vector<CigarOp> &cigar)
    : AbstractClip(referenceId, mapPosition, clipPosition, matePosition, sequence, cigar) {

}

Deletion ReverseBClip::call(const std::vector<Reference> &references, int minOverlap, double minIdentity)
{

}

void ReverseBClip::fetchAnchors(BamReader &reader, int insLength, std::vector<int> &anchors) {
    anchors.push_back(matePosition);
}

TargetRegion ReverseBClip::tRegion(const string &referenceName, int anchor, int insLength) {
    return { referenceName,
             anchor,
             anchor + insLength
           };
}


ForwardEClip::ForwardEClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const string &sequence, const std::vector<CigarOp> &cigar)
    : AbstractClip(referenceId, mapPosition, clipPosition, matePosition, sequence, cigar) {
}

Deletion ForwardEClip::call(FaidxWrapper &faidx, const std::vector<TargetRegion> &regions, int minOverlap, double minIdentity)
{

}

void ForwardEClip::fetchAnchors(BamReader &reader, int insLength, std::vector<int> &anchors) {
    anchors.push_back(matePosition);
}

TargetRegion ForwardEClip::tRegion(const string &referenceName, int anchor, int insLength) {
    return {referenceName,
             anchor + length(),
             anchor - insLength };
}


ReverseEClip::ReverseEClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const string &sequence, const std::vector<CigarOp> &cigar)
    : AbstractClip(referenceId, mapPosition, clipPosition, matePosition, sequence, cigar) {
}

Deletion ReverseEClip::call(FaidxWrapper &faidx, const std::vector<TargetRegion> &regions, int minOverlap, double minIdentity)
{
    for (auto it = regions.begin(); it != regions.end(); ++it) {
        SequenceOverlap overlap = Overlapper::computeOverlapSW((*it).sequence(faidx), sequence, ungapped_params);
        if (overlap.getOverlapLength() > minOverlap &&
                overlap.getPercentIdentity() >= minIdentity &&
                overlap.match[1].end == overlap.length[1] - 1) {
            int rightBp = (*it).start + overlap.match[0].start - 1;
            int leftBp = (overlap.getOverlapLength()  > cigar[cigar.size() - 1].Length) ? clipPosition - overlap.getOverlapLength() + cigar[cigar.size() - 1].Length
                    : clipPosition;
            leftBp--;   // left breakpoint refers the position of the last base prior to the clipped part conforming to the VCF format.
            int len = leftBp - rightBp;
            if (overlap.getOverlapLength() < cigar[cigar.size() - 1].Length) len += cigar[cigar.size() - 1].Length - overlap.getOverlapLength();
            return Deletion(referenceId,leftBp, rightBp, len);
        }
    }
    error("No deletion is found.");
}

void ReverseEClip::fetchAnchors(BamReader &reader, int insLength, std::vector<int> &anchors) {
    int end = leftmostPosition() + length();
    int start = end - insLength + length();

    if (!reader.SetRegion(referenceId, start, referenceId, end))
        error("Could not set the region.");

    BamAlignment al;
    while(reader.GetNextAlignmentCore(al)) {
        if (!al.IsReverseStrand() && al.IsMateReverseStrand() && al.Position < al.MatePosition) {
            anchors.push_back(al.MatePosition);
        }
    }
}

TargetRegion ReverseEClip::tRegion(const string &referenceName, int anchor, int insLength) {
    return { referenceName,
             anchor - insLength,
             anchor + length() };
}

string TargetRegion::sequence(FaidxWrapper& faidx)
{
    return faidx.fetch(referenceName, start, end);
}
