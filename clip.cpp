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
    string refName = Helper::getReferenceName(reader, referenceId);

    vector<IRange> ranges;
    fetchSpanningRanges(reader, insLength, ranges);

    if (ranges.empty()) error("No deletion is found");

    vector<TargetRegion> regions;
    toTargetRegions(refName, insLength, ranges, regions);

    return call(faidx, regions, minOverlap, minIdentity);
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
        SequenceOverlap overlap = Overlapper::computeOverlapSW(s1, s2, ungapped_params);
        for (size_t i = 0; i < 2; ++i)
            overlap.match[i].flipStrand(overlap.length[i]);

        if (mapPosition == 60735213) {
//            cout << overlap.getOverlapLength() << endl;
//            overlap.printAlignment((*it).sequence(faidx), sequence);
        }
        if (overlap.getOverlapLength() >= minOverlap &&
                overlap.getPercentIdentity() >= minIdentity * 100) {
            int delta = cigar[0].Length - overlap.getOverlapLength();
            int rightBp = clipPosition - 1;
            int leftBp = (*it).start + overlap.match[0].end;
            if (delta < 0) leftBp += delta;
            int len = leftBp - rightBp;
            if (delta > 0) len += delta;
            if (len > Helper::SVLEN_THRESHOLD) break;
//            overlap.printAlignment((*it).sequence(faidx), sequence);
            return Deletion((*it).referenceName, leftBp, rightBp, len);
        }
    }
    error("No deletion is found.");
}

void ForwardBClip::fetchSpanningRanges(BamReader &reader, int insLength, std::vector<IRange> &ranges)
{
    int start = clipPosition;
    int end = clipPosition + insLength - 2 * length();

    if (!reader.SetRegion(referenceId, start - 1, referenceId, end))
        error("Could not set the region.");

    BamAlignment al;
    while(reader.GetNextAlignment(al)) {
//        if (al.Position < start - 1) continue;
        if (al.IsReverseStrand() && !al.IsMateReverseStrand() && al.RefID == al.MateRefID
                && al.Position > al.MatePosition && al.MatePosition + length() - Helper::SVLEN_THRESHOLD <= clipPosition) {
            ranges.push_back({al.MatePosition + 1, al.Position + 1});
        }
    }
//    if (clipPosition == 58927037) {
//        cout << "Debug: placeholder" << endl;
//    }

}

void ForwardBClip::toTargetRegions(const string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions)
{
    int rightmostPos = clipPosition + Helper::SVLEN_THRESHOLD;

    std::vector<IRange> newRanges(ranges.size());
    transform(ranges.begin(), ranges.end(), newRanges.begin(), [=](const IRange &ran) { IRange r = {ran.start, ran.start + insLength - length()}; return r; });
    std::vector<IdCluster> idClusters;
    clusterRanges(newRanges, idClusters);
    transform(idClusters.begin(), idClusters.end(), back_inserter(regions)
              , [=](const IdCluster &id) { int e = newRanges[id[id.size() - 1]].end; TargetRegion tr = { referenceName, newRanges[id[0]].start, (e > rightmostPos) ? rightmostPos : e}; return tr;});

//    sort(ranges.begin(), ranges.end());

//    for (auto it = ranges.begin(); it != ranges.end(); ++it) {
//        if ((*it).start > rightmostPos + length()) break;
//        int s = (*it).start;
//        int e = insLength - length() - ((*it).end - clipPosition) + (*it).start;
//        if (e > rightmostPos) {
//            e = rightmostPos;
//            regions.push_back({referenceName, s, e});
//            break;
//        }
//        regions.push_back({referenceName, s, e});
//    }
}


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
    SequenceOverlap overlap = Overlapper::computeOverlapSW(s1, sequence);
//    overlap.printAlignment(s1, sequence);
//    if (mapPosition == 2863866) {
//        cout << "Debug: placehold" << endl;
//        overlap.printAlignment(s1, sequence);
//    }
    if (overlap.getOverlapLength() >= minOverlap &&
            overlap.getPercentIdentity() >= minIdentity * 100) {
        int rightBp = regions[0].start + overlap.match[0].start - 1;
        int leftBp = (overlap.getOverlapLength()  > cigar[cigar.size() - 1].Length) ? clipPosition - overlap.getOverlapLength() + cigar[cigar.size() - 1].Length
                : clipPosition;
        leftBp--;   // left breakpoint refers the position of the last base prior to the clipped part conforming to the VCF format.
        int len = leftBp - rightBp;
        if (overlap.getOverlapLength() < cigar[cigar.size() - 1].Length) len += cigar[cigar.size() - 1].Length - overlap.getOverlapLength();
        if (len > Helper::SVLEN_THRESHOLD) error("No deletion was found.");
//        overlap.printAlignment(s1, sequence);
        return Deletion(regions[0].referenceName, leftBp, rightBp, len);
    }
    error("No deletion was found.");
}



ReverseEClip::ReverseEClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const string &sequence, const std::vector<CigarOp> &cigar)
    : AbstractClip(referenceId, mapPosition, clipPosition, matePosition, sequence, cigar) {
}

void ReverseEClip::fetchSpanningRanges(BamReader &reader, int insLength, std::vector<IRange> &ranges)
{
    int start = clipPosition - insLength + length();
    int end = clipPosition - length();

//    if (mapPosition == 33621625) {
//        cout << start << "\t" << end << endl;
//        cout << "-------------------------" << endl;
//    }

    if (!reader.SetRegion(referenceId, start - 1, referenceId, end))
        error("Could not set the region.");

    BamAlignment al;
    while(reader.GetNextAlignmentCore(al)) {
        if (al.Position < start - 1) continue;
        if (!al.IsReverseStrand() && al.IsMateReverseStrand() && al.RefID == al.MateRefID
                && al.Position < al.MatePosition && al.MatePosition >= clipPosition - Helper::SVLEN_THRESHOLD) {
            ranges.push_back({al.Position + 1, al.MatePosition + 1});
        }
    }
}

void ReverseEClip::toTargetRegions(const string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions)
{
    int leftmostPos = clipPosition;

    std::vector<IRange> newRanges(ranges.size());
    transform(ranges.begin(), ranges.end(), newRanges.begin(), [=](const IRange &ran) { IRange r = {ran.end - insLength + 2 * length(), ran.end + length()}; return r; });
    std::vector<IdCluster> idClusters;
    clusterRanges(newRanges, idClusters);
    transform(idClusters.begin(), idClusters.end(), back_inserter(regions)
              , [=](const IdCluster &id) { int s = newRanges[id[0]].start; TargetRegion tr = { referenceName, (s < leftmostPos) ? leftmostPos : s, newRanges[id[id.size() - 1]].end}; return tr;});

//    sort(ranges.begin(), ranges.end(),
//         [](const IRange &lhs, const IRange &rhs) { if (lhs.end != rhs.end) return lhs.end > rhs.end; return lhs.start > rhs.end; });

//    for (auto it = ranges.begin(); it != ranges.end(); ++it) {
//        if ((*it).end < leftmostPos) break;
////        if ((*it).start > clipPosition) continue;
//        int e = (*it).end + length();
//        int s = clipPosition - (*it).start + (*it).end - (insLength - length());
//        if (s < leftmostPos) {
//            s = leftmostPos;
//            regions.push_back({referenceName, s, e});
//            break;
//        }
//        regions.push_back({referenceName, s, e});
//    }
}

Deletion ReverseEClip::call(FaidxWrapper &faidx, const std::vector<TargetRegion> &regions, int minOverlap, double minIdentity)
{
//    error("No deletion is found.");
    for (auto it = regions.rbegin(); it != regions.rend(); ++it) {
        string s1 = (*it).sequence(faidx);
        SequenceOverlap overlap = Overlapper::computeOverlapSW(s1, sequence, ungapped_params);
        if (mapPosition == 19341369) {
//            cout << overlap.getOverlapLength() << "\t" << (*it).start << endl;
//            overlap.printAlignment((*it).sequence(faidx), sequence);
        }
        if (overlap.getOverlapLength() >= minOverlap &&
                overlap.getPercentIdentity() >= minIdentity * 100) {
            int delta = cigar[cigar.size() - 1].Length - overlap.getOverlapLength();
            int rightBp = (*it).start + overlap.match[0].start - 1;
            int leftBp = clipPosition - 1;
            if (delta < 0) leftBp += delta;
            int len = leftBp - rightBp;
            if (delta > 0) len += delta;
            if (len > Helper::SVLEN_THRESHOLD) break;
//            overlap.printAlignment((*it).sequence(faidx), sequence);
            return Deletion((*it).referenceName,leftBp, rightBp, len);
        }
    }
    error("No deletion is found.");
}


string TargetRegion::sequence(FaidxWrapper& faidx) const
{
    return faidx.fetch(referenceName, start - 1, end - 1);
}
