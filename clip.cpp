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
      cigar(cigar),
      conflictFlag(false) {
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
    string refName = Helper::getReferenceName(reader, referenceId);

    vector<IRange> ranges;
    fetchSpanningRanges(reader, insLength, ranges);

    if (ranges.empty()) error("No deletion is found");

    vector<TargetRegion> regions;
    toTargetRegions(refName, insLength, ranges, regions);

    return call(faidx, regions, minOverlap, minIdentity);
}

bool AbstractClip::hasConflictWith(AbstractClip *other) {
    if (getType() == other->getType()) return false;
    return abs(clipPosition - other->clipPosition) < Helper::CONFLICT_THRESHOLD;
}

bool AbstractClip::getConflictFlag() const
{
    return conflictFlag;
}

void AbstractClip::setConflictFlag(bool value)
{
    conflictFlag = value;
}



ForwardBClip::ForwardBClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const string &sequence, const vector<CigarOp>& cigar)
    : AbstractClip(referenceId, mapPosition, clipPosition, matePosition, sequence, cigar) {
}

Deletion ForwardBClip::call(FaidxWrapper &faidx, const std::vector<TargetRegion> &regions, int minOverlap, double minIdentity)
{
//    error("No deletion is found.");
    for (auto it = regions.begin(); it != regions.end(); ++it) {
        string s1 = (*it).sequence(faidx);
        reverse(s1.begin(), s1.end());
        string s2 = sequence;
        reverse(s2.begin(), s2.end());
        SequenceOverlap overlap = Overlapper::computeOverlapSW(s1, s2, minOverlap, minIdentity);

        for (size_t i = 0; i < 2; ++i)
            overlap.match[i].flipStrand(overlap.length[i]);

        int delta = overlap.getOverlapLength() - cigar[0].Length;
        int offset = 0;
        for (auto &ci: cigar) {
            if (ci.Type == 'D') offset += ci.Length;
            else if (ci.Type == 'I') offset -= ci.Length;
        }
        int rightBp = clipPosition + offset;
        int leftBp = (*it).start + overlap.match[0].start + cigar[0].Length - 1;
        int len = leftBp - rightBp + 1;
        int start1 = delta > 0 ? leftBp : leftBp + delta;
        int start2 = delta > 0 ? leftBp + delta : leftBp;
        int end1 = delta > 0 ? rightBp : rightBp + delta;
        int end2 = delta > 0 ? rightBp + delta : rightBp;
        if (len > Helper::SVLEN_THRESHOLD) continue;
        return Deletion((*it).referenceName, start1, start2, end1, end2, len);
    }
    error("No deletion is found.");
}

string ForwardBClip::getType()
{
    return "FH";
}

void ForwardBClip::fetchSpanningRanges(BamReader &reader, int insLength, std::vector<IRange> &ranges)
{
    int start = clipPosition;
    int end = clipPosition + insLength - 2 * length();

    if (!reader.SetRegion(referenceId, start - 1, referenceId, end))
        error("Could not set the region.");

    BamAlignment al;
    while(reader.GetNextAlignment(al)) {
        string xt;
        al.GetTag("XT", xt);
        xt = xt.substr(0,1);
        if (al.IsReverseStrand() && !al.IsMateReverseStrand() && al.RefID == al.MateRefID
                && al.MapQuality > 0 && xt == "U"
                && al.Position > al.MatePosition && al.MatePosition + length() - Helper::SVLEN_THRESHOLD <= clipPosition) {
            ranges.push_back({al.MatePosition + 1, al.Position + 1});
        }
    }

}

void ForwardBClip::toTargetRegions(const string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions)
{
    int rightmostPos = clipPosition + length();

    std::vector<IRange> newRanges(ranges.size());
    transform(ranges.begin(), ranges.end(), newRanges.begin(), [=](const IRange &ran) { IRange r = {ran.start, ran.start + insLength - length()}; return r; });
    std::vector<IdCluster> idClusters;
    clusterRanges(newRanges, idClusters);
    transform(idClusters.begin(), idClusters.end(), back_inserter(regions)
              , [=](const IdCluster &id) { int e = newRanges[id[id.size() - 1]].end; TargetRegion tr = { referenceName, newRanges[id[0]].start, (e > rightmostPos) ? rightmostPos : e}; return tr;});
}


/*
ReverseBClip::ReverseBClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const string &sequence, const std::vector<CigarOp> &cigar)
    : AbstractClip(referenceId, mapPosition, clipPosition, matePosition, sequence, cigar) {

}

void ReverseBClip::fetchSpanningRanges(BamReader &reader, int insLength, std::vector<IRange> &ranges)
{

}

void ReverseBClip::toTargetRegions(const string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions)
{

}

Deletion ReverseBClip::call(FaidxWrapper &faidx, const std::vector<TargetRegion> &regions, int minOverlap, double minIdentity)
{

}


ForwardEClip::ForwardEClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const string &sequence, const std::vector<CigarOp> &cigar)
    : AbstractClip(referenceId, mapPosition, clipPosition, matePosition, sequence, cigar) {
}

void ForwardEClip::fetchSpanningRanges(BamReader &reader, int insLength, std::vector<IRange> &ranges)
{
    ranges.push_back({clipPosition, matePosition});
}

void ForwardEClip::toTargetRegions(const string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions)
{
    int pe = ranges[0].end + length();
    int leftmostPos = clipPosition;
    int len1 = length() - cigar[cigar.size() - 1].Length;
    int cPrime = ranges[0].end  + len1 - insLength;
    if (cPrime < leftmostPos) cPrime = leftmostPos;
    regions.push_back({referenceName, cPrime, pe});
}

Deletion ForwardEClip::call(FaidxWrapper &faidx, const std::vector<TargetRegion> &regions, int minOverlap, double minIdentity)
{
    string s1 = regions[0].sequence(faidx);
    SequenceOverlap overlap = Overlapper::computeOverlapSW(s1, sequence, minOverlap, minIdentity, ungapped_params);
    if (overlap.getOverlapLength() >= minOverlap &&
            overlap.getPercentIdentity() >= minIdentity * 100) {
        int rightBp = regions[0].start + overlap.match[0].start - 1;
        int leftBp = (overlap.getOverlapLength()  > cigar[cigar.size() - 1].Length) ? clipPosition - overlap.getOverlapLength() + cigar[cigar.size() - 1].Length
                : clipPosition;
        leftBp--;   // left breakpoint refers the position of the last base prior to the clipped part conforming to the VCF format.
        int len = leftBp - rightBp;
        if (overlap.getOverlapLength() < cigar[cigar.size() - 1].Length) len += cigar[cigar.size() - 1].Length - overlap.getOverlapLength();
        if (len > Helper::SVLEN_THRESHOLD) error("No deletion was found.");
        return Deletion(regions[0].referenceName, leftBp, leftBp, rightBp, rightBp, len);
    }
    error("No deletion was found.");
}
*/


ReverseEClip::ReverseEClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const string &sequence, const std::vector<CigarOp> &cigar)
    : AbstractClip(referenceId, mapPosition, clipPosition, matePosition, sequence, cigar) {
}

void ReverseEClip::fetchSpanningRanges(BamReader &reader, int insLength, std::vector<IRange> &ranges)
{
    int start = clipPosition - insLength + length();
    int end = clipPosition - length();

    if (!reader.SetRegion(referenceId, start - 1, referenceId, end))
        error("Could not set the region.");

    BamAlignment al;
    while(reader.GetNextAlignment(al)) {
        string xt;
        al.GetTag("XT", xt);
        xt = xt.substr(0,1);
        if (al.Position < start - 1) continue;
        if (!al.IsReverseStrand() && al.IsMateReverseStrand() && al.RefID == al.MateRefID
                && al.MapQuality > 0 && xt == "U"
                && al.Position < al.MatePosition && al.MatePosition >= clipPosition - Helper::SVLEN_THRESHOLD) {
            ranges.push_back({al.Position + 1, al.MatePosition + 1});
        }
    }

}

void ReverseEClip::toTargetRegions(const string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions)
{
    int leftmostPos = clipPosition - length();

    std::vector<IRange> newRanges(ranges.size());
    transform(ranges.begin(), ranges.end(), newRanges.begin(), [=](const IRange &ran) { IRange r = {ran.end - insLength + 2 * length(), ran.end + length()}; return r; });
    std::vector<IdCluster> idClusters;
    clusterRanges(newRanges, idClusters);
    transform(idClusters.begin(), idClusters.end(), back_inserter(regions)
              , [=](const IdCluster &id) { int s = newRanges[id[0]].start; TargetRegion tr = { referenceName, (s < leftmostPos) ? leftmostPos : s, newRanges[id[id.size() - 1]].end}; return tr;});

}

Deletion ReverseEClip::call(FaidxWrapper &faidx, const std::vector<TargetRegion> &regions, int minOverlap, double minIdentity)
{
//    error("No deletion is found.");
    for (auto it = regions.rbegin(); it != regions.rend(); ++it) {
        string s1 = (*it).sequence(faidx);
        SequenceOverlap overlap = Overlapper::computeOverlapSW(s1, sequence, minOverlap, minIdentity);

        int delta = overlap.getOverlapLength() - cigar[cigar.size() - 1].Length;
        int offset = 0;
        for (auto &ci: cigar) {
            if (ci.Type == 'D') offset += ci.Length;
            else if (ci.Type == 'I') offset -= ci.Length;
        }
        int leftBp = clipPosition - 1 - offset;
        int rightBp = (*it).start + overlap.match[0].end - cigar[cigar.size() - 1].Length + 1;
        int len = leftBp - rightBp + 1;
        int start1 = delta > 0 ? leftBp - delta : leftBp;
        int start2 = delta > 0 ? leftBp : leftBp - delta;
        int end1 = delta > 0 ? rightBp - delta : rightBp;
        int end2 = delta > 0 ? rightBp : rightBp - delta;

        if (len > Helper::SVLEN_THRESHOLD) continue;
        return Deletion((*it).referenceName,start1, start2, end1, end2, len);
    }
    error("No deletion is found.");
}

string ReverseEClip::getType()
{
    return "RH";
}


string TargetRegion::sequence(FaidxWrapper& faidx) const
{
    return faidx.fetch(referenceName, start - 1, end - 1);
}
