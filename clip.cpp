#include "clip.h"
#include "error.h"
#include "Helper.h"
#include <iterator>
#include <algorithm>
#include <sstream>
#include <utility>

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

Deletion AbstractClip::call(BamReader &reader, FaidxWrapper &faidx, int insLength, int minOverlap, double minIdentity, int minMapQual)
{
    string refName = Helper::getReferenceName(reader, referenceId);

    vector<IRange> ranges;
    fetchSpanningRanges(reader, insLength, ranges, minMapQual);
//    vector<int> sizes;
//    fecthSizesForSpanningPairs(reader, insLength, sizes);

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

int AbstractClip::maxEditDistanceForSoftclippedPart()
{
    if (lengthOfSoftclippedPart() >= 20) return 2;
    return 1;
}



ForwardBClip::ForwardBClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const string &sequence, const vector<CigarOp>& cigar)
    : AbstractClip(referenceId, mapPosition, clipPosition, matePosition, sequence, cigar) {
}

/*
Deletion ForwardBClip::call(FaidxWrapper &faidx, const std::vector<TargetRegion> &regions, int minOverlap, double minIdentity) {
    for (auto it = regions.rbegin(); it != regions.rend(); ++it) {
        if ((*it).length() < lengthOfSoftclippedPart()) continue;
        string s1 = (*it).sequence(faidx);
        SequenceOverlap overlap = Overlapper::computeOverlapSG(s1, softclippedPart());
        if (overlap.edit_distance > maxEditDistanceForSoftclippedPart()) continue;
        int leftEnd = (*it).start + overlap.match[0].end;
        int offsetToRight = offsetFromThatEnd((*it).referenceName, faidx, leftEnd);
        int rightEnd = clipPosition;
        int offsetToLeft = offsetFromThisEnd((*it).referenceName, faidx);
        int start1 = leftEnd - offsetToLeft;
        int start2 = leftEnd + offsetToRight;
        int end1 = rightEnd - offsetToLeft;
        int end2 = rightEnd + offsetToRight;
        int len = start1 - end1 + 1;
        if (len > Helper::SVLEN_THRESHOLD) continue;
        return Deletion((*it).referenceName, start1, start2, end1, end2, len, getType());
    }
    error("No deletion is found.");
}
*/

Deletion ForwardBClip::call(FaidxWrapper &faidx, const std::vector<TargetRegion> &regions, int minOverlap, double minIdentity)
{
//    error("No deletion is found.");
    ScoreParam score_param(1, -1, 2, 4);
    for (auto it = regions.begin(); it != regions.end(); ++it) {
        string s1 = (*it).sequence(faidx);
        reverse(s1.begin(), s1.end());
        string s2 = sequence;
        reverse(s2.begin(), s2.end());

        SequenceOverlap overlap;

        try {
            overlap = Overlapper::computeOverlapSW2(s1, s2, minOverlap, minIdentity, ungapped_params);
        } catch (ErrorException& ex) {
            continue;
        }

        for (size_t i = 0; i < 2; ++i)
            overlap.match[i].flipStrand(overlap.length[i]);

//        overlap = Overlapper::alignPrefix(s1, s2, ungapped_params);

//        if (s2 == "GCCTACAGAGTGCAGAGCCAGCCCAGGACAGGGGACAATTACACAGGCGATGGTCCTAAGAACCGAACCTTCCAATCCCAAAACTCTAGACAGGTATCCAA")
//            cout << s1 << endl;
//        overlap = Overlapper::ageAlignPrefix(s1, s2, score_param);
//        if (!overlap.isQualified(minOverlap, minIdentity))
//            continue;

        int delta = overlap.getOverlapLength() - lengthOfSoftclippedPart();
        int offset = 0;
        for (auto &ci: cigar) {
            if (ci.Type == 'D') offset += ci.Length;
            else if (ci.Type == 'I') offset -= ci.Length;
        }
        int rightBp = clipPosition + offset;
        int leftBp = (*it).start + overlap.match[0].start + lengthOfSoftclippedPart() - 1;

//        int delta = overlap.match[1].length() - lengthOfSoftclippedPart();
//        int rightBp = clipPosition;
//        // leftBp might need to be adjusted.
//        int leftBp = (*it).start + overlap.match[0].start + lengthOfSoftclippedPart() - 1;

        int len = leftBp - rightBp + 1;
        int start1 = delta > 0 ? leftBp : leftBp + delta;
        int start2 = delta > 0 ? leftBp + delta : leftBp;
        int end1 = delta > 0 ? rightBp : rightBp + delta;
        int end2 = delta > 0 ? rightBp + delta : rightBp;

//        if (start2 == 23483811) {
//            cout << overlap << endl;
//            cout << s2 << endl;
//            cout << s1 << endl;
//        }

        if (len > Helper::SVLEN_THRESHOLD) continue;
        return Deletion((*it).referenceName, start1, start2, end1, end2, len, getType());
    }
    error("No deletion is found.");
}

string ForwardBClip::getType()
{
    return "5F";
}

void ForwardBClip::fetchSpanningRanges(BamReader &reader, int insLength, std::vector<IRange> &ranges, int minMapQual)
{
    // SVSeq2.length
//    int start = leftmostPosition();
    int start = clipPosition;
//    int end = start + insLength + length();
    int end = start + insLength - 2 * length();

    if (start > end) error("the region is invalid.");

    if (!reader.SetRegion(referenceId, start - 1, referenceId, end))
        error("Could not set the region.");

    BamAlignment al;
    while(reader.GetNextAlignment(al)) {
//        string xt;
//        al.GetTag("XT", xt);
//        xt = xt.substr(0,1);
        if (al.IsReverseStrand() && !al.IsMateReverseStrand() && al.RefID == al.MateRefID
                && al.MapQuality >= minMapQual //&& xt == "U"
                && al.Position > al.MatePosition && al.MatePosition + length() - Helper::SVLEN_THRESHOLD <= clipPosition) {
            ranges.push_back({al.MatePosition + 1, al.Position + 1});
        }
    }

}

void ForwardBClip::fecthSizesForSpanningPairs(BamReader &reader, int insLength, std::vector<int> &sizes)
{
    int start = clipPosition;
    int end = clipPosition + insLength + length();

    if (!reader.SetRegion(referenceId, start - 1, referenceId, end))
        error("Could not set the region.");

    vector<pair<int, int> > records;
    BamAlignment al;
    while(reader.GetNextAlignment(al)) {
        if (al.IsReverseStrand() && !al.IsMateReverseStrand() && al.RefID == al.MateRefID
                && al.MapQuality > 0 && al.Position > al.MatePosition) {
            records.push_back(make_pair(abs(al.InsertSize), al.Position - clipPosition));
        }
    }
    sort(records.begin(), records.end(), [](const pair<int,int>& r1, const pair<int,int>& r2){ return r1.first < r2.first; });
    cout << ">" << clipPosition << "," << mapPosition << endl;
    transform(records.begin(), records.end(), ostream_iterator<string>(cout, " "), [](const pair<int,int>& r){
        stringstream ss;
        ss << "(" << r.first << "," << r.second << ")";
        return ss.str();
    });
    cout << endl;
}

void ForwardBClip::toTargetRegions(const string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions)
{
    int rightmostPos = clipPosition + length();

    std::vector<IRange> newRanges(ranges.size());
//    transform(ranges.begin(), ranges.end(), newRanges.begin(), [=](const IRange &ran) { IRange r = {ran.start, ran.start + insLength + length()}; return r; });
    transform(ranges.begin(), ranges.end(), newRanges.begin(), [=](const IRange &ran) { IRange r = {ran.start, ran.start + insLength - length()}; return r; });
    std::vector<IdCluster> idClusters;
    clusterRanges(newRanges, idClusters);
// Replace with the merging method used by SVSeq2
//    sort(std::begin(newRanges), std::end(newRanges));
//    clusterRanges2(newRanges, idClusters);
    for (auto &elt : idClusters) {
        int s = newRanges[elt.front()].start;
        if (s > rightmostPos) break;
        int e = newRanges[elt.back()].end;
        if (e > rightmostPos) e = rightmostPos;
        if (s > e) break;
        regions.push_back({referenceName, s, e});
    }
}


/*

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

void ReverseEClip::fetchSpanningRanges(BamReader &reader, int insLength, std::vector<IRange> &ranges, int minMapQual)
{
    // Experiment ID: SVSeq2.length
    int start = clipPosition - insLength + length();
//    int end = leftmostPosition() + length();
//    int start = end - insLength - length();
    if (start < 0) start = 0;
    int end = clipPosition - length();

    if (start > end) error("the region is invalid.");

    if (!reader.SetRegion(referenceId, start - 1, referenceId, end))
        error("Could not set the region.");

    BamAlignment al;
    while(reader.GetNextAlignment(al)) {
//        string xt;
//        al.GetTag("XT", xt);
//        xt = xt.substr(0,1);
        if (al.Position < start - 1) continue;
        if (!al.IsReverseStrand() && al.IsMateReverseStrand() && al.RefID == al.MateRefID
                && al.MapQuality >= minMapQual //&& xt == "U"
                && al.Position < al.MatePosition && al.MatePosition >= clipPosition - Helper::SVLEN_THRESHOLD) {
            ranges.push_back({al.Position + 1, al.MatePosition + 1});
        }
    }

}

void ReverseEClip::fecthSizesForSpanningPairs(BamReader &reader, int insLength, std::vector<int> &sizes)
{

}

void ReverseEClip::toTargetRegions(const string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions)
{
    int leftmostPos = clipPosition - length();

    std::vector<IRange> newRanges(ranges.size());
//    transform(ranges.begin(), ranges.end(), newRanges.begin(), [=](const IRange &ran) { IRange r = {ran.end - insLength - length(), ran.end}; return r; });
    transform(ranges.begin(), ranges.end(), newRanges.begin(), [=](const IRange &ran) { IRange r = {ran.end - insLength + 2 * length(), ran.end + length()}; return r; });
    std::vector<IdCluster> idClusters;
    clusterRanges(newRanges, idClusters);
// Replace with the merging method used by SVSeq2
//    sort(std::begin(newRanges), std::end(newRanges));
//    clusterRanges2(newRanges, idClusters);
    for (auto &elt : idClusters) {
        int e = newRanges[elt.back()].end;
        if (e < leftmostPos) continue;
        int s = newRanges[elt.front()].start;
        if (s < leftmostPos) s = leftmostPos;
        if (s > e) continue;
        regions.push_back({referenceName, s, e});
    }

}

/*
Deletion ReverseEClip::call(FaidxWrapper &faidx, const std::vector<TargetRegion> &regions, int minOverlap, double minIdentity) {
    for (auto it = regions.rbegin(); it != regions.rend(); ++it) {
        if ((*it).length() < lengthOfSoftclippedPart()) continue;
        string s1 = (*it).sequence(faidx);
        SequenceOverlap overlap = Overlapper::computeOverlapSG(s1, softclippedPart());
        if (overlap.edit_distance > maxEditDistanceForSoftclippedPart()) continue;
        int leftEnd = clipPosition - 1;
        int offsetToRight = offsetFromThisEnd((*it).referenceName, faidx);
        int rightEnd = (*it).start + overlap.match[0].start;
        int offsetToLeft = offsetFromThatEnd((*it).referenceName, faidx, rightEnd);
        int start1 = leftEnd - offsetToLeft;
        int start2 = leftEnd + offsetToRight;
        int end1 = rightEnd - offsetToLeft;
        int end2 = rightEnd + offsetToRight;
        int len = start1 - end1 + 1;
        if (len > Helper::SVLEN_THRESHOLD) continue;
        return Deletion((*it).referenceName, start1, start2, end1, end2, len, getType());
    }
    error("No deletion is found.");
}
*/

Deletion ReverseEClip::call(FaidxWrapper &faidx, const std::vector<TargetRegion> &regions, int minOverlap, double minIdentity)
{
//    error("No deletion is found.");
    ScoreParam score_param(1, -1, 2, 4);

    for (auto it = regions.rbegin(); it != regions.rend(); ++it) {
        string s1 = (*it).sequence(faidx);
        string s2 = sequence;
        SequenceOverlap overlap;

        try {
            overlap = Overlapper::computeOverlapSW2(s1, s2, minOverlap, minIdentity, ungapped_params);
        } catch (ErrorException& ex) {
            continue;
        }

//        overlap = Overlapper::alignSuffix(s1, s2, ungapped_params);
//        overlap = Overlapper::ageAlignSuffix(s1, s2, score_param);
//        if (!overlap.isQualified(minOverlap, minIdentity))
//            continue;

        int delta = overlap.getOverlapLength() - lengthOfSoftclippedPart();
        int offset = 0;
        for (auto &ci: cigar) {
            if (ci.Type == 'D') offset += ci.Length;
            else if (ci.Type == 'I') offset -= ci.Length;
        }
        int leftBp = clipPosition - 1 - offset;
        int rightBp = (*it).start + overlap.match[0].end - lengthOfSoftclippedPart() + 1;

//        int delta = overlap.match[1].length() - lengthOfSoftclippedPart();
//        int leftBp = clipPosition - 1;
//        // rightBp might need to be adjusted
//        int rightBp = (*it).start + overlap.match[0].end - lengthOfSoftclippedPart() + 1;

        int len = leftBp - rightBp + 1;
        int start1 = delta > 0 ? leftBp - delta : leftBp;
        int start2 = delta > 0 ? leftBp : leftBp - delta;
        int end1 = delta > 0 ? rightBp - delta : rightBp;
        int end2 = delta > 0 ? rightBp : rightBp - delta;

//        if (start2 == 54151129) {
//            cout << overlap << endl;
//            cout << s2 << endl;
//            cout << s1 << endl;
//        }

        if (len > Helper::SVLEN_THRESHOLD) continue;
        return Deletion((*it).referenceName, start1, start2, end1, end2, len, getType());
    }
    error("No deletion is found.");
}

string ReverseEClip::getType()
{
    return "5R";
}


int ForwardBClip::offsetFromThisEnd(string referenceName, FaidxWrapper &faidx)
{
    return numOfThelongestSuffix(softclippedPart(),
                                 faidx.fetch(referenceName,
                                             clipPosition - lengthOfSoftclippedPart(),
                                             clipPosition - 1));
}

int ForwardBClip::offsetFromThatEnd(string referenceName, FaidxWrapper &faidx, int orignal)
{
    return numOfTheLongestPrefix(mappedPart(),
                                 faidx.fetch(referenceName,
                                             orignal + 1,
                                             orignal + lengthOfMappedPart()));
}

int ReverseEClip::offsetFromThisEnd(string referenceName, FaidxWrapper &faidx)
{
    return numOfTheLongestPrefix(softclippedPart(),
                                 faidx.fetch(referenceName,
                                             clipPosition,
                                             clipPosition + lengthOfSoftclippedPart() - 1));
}

int ReverseEClip::offsetFromThatEnd(string referenceName, FaidxWrapper &faidx, int orignal)
{
    return numOfThelongestSuffix(mappedPart(),
                                 faidx.fetch(referenceName,
                                             orignal - lengthOfMappedPart(),
                                             orignal - 1));
}


ReverseBClip::ReverseBClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const string &sequence, const std::vector<CigarOp> &cigar)
    : AbstractClip(referenceId, mapPosition, clipPosition, matePosition, sequence, cigar) {
}

string ReverseBClip::getType()
{
    return "3R";
}

Deletion ReverseBClip::call(FaidxWrapper &faidx, const std::vector<TargetRegion> &regions, int minOverlap, double minIdentity)
{
    string s1 = regions[0].sequence(faidx);
    reverse(s1.begin(), s1.end());
    string s2 = sequence;
    reverse(s2.begin(), s2.end());
    SequenceOverlap overlap = Overlapper::computeOverlapSW2(s1, s2, minOverlap, minIdentity, ungapped_params);

    for (size_t i = 0; i < 2; ++i)
        overlap.match[i].flipStrand(overlap.length[i]);

    int delta = overlap.getOverlapLength() - lengthOfSoftclippedPart();
    int offset = 0;
    for (auto &ci: cigar) {
        if (ci.Type == 'D') offset += ci.Length;
        else if (ci.Type == 'I') offset -= ci.Length;
    }
    int rightBp = clipPosition + offset;
    int leftBp = regions[0].start + overlap.match[0].start + lengthOfSoftclippedPart() - 1;

    int len = leftBp - rightBp + 1;
    int start1 = delta > 0 ? leftBp : leftBp + delta;
    int start2 = delta > 0 ? leftBp + delta : leftBp;
    int end1 = delta > 0 ? rightBp : rightBp + delta;
    int end2 = delta > 0 ? rightBp + delta : rightBp;
    if (len > Helper::SVLEN_THRESHOLD) error("No deletion is found.");
    return Deletion(regions[0].referenceName, start1, start2, end1, end2, len, getType());

    error("No deletion is found.");
}

void ReverseBClip::fetchSpanningRanges(BamReader &reader, int insLength, std::vector<IRange> &ranges, int minMapQual)
{
    ranges.push_back({matePosition + 1, clipPosition + 1});
}

void ReverseBClip::fecthSizesForSpanningPairs(BamReader &reader, int inslength, std::vector<int> &sizes)
{
}

void ReverseBClip::toTargetRegions(const string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions)
{
    // ran.start, ran.start + insLength - length()
    int rightmostPos = clipPosition + length();
    int s = ranges[0].start;
    int e = s + insLength - length();
    if (e > rightmostPos) e = rightmostPos;
    if (s > e) return;
    regions.push_back({referenceName, s, e});
}

int ReverseBClip::lengthOfSoftclippedPart()
{
    return cigar[0].Length;
}

string ReverseBClip::softclippedPart()
{
    return sequence.substr(0, lengthOfSoftclippedPart());
}

string ReverseBClip::mappedPart()
{
    return sequence.substr(lengthOfSoftclippedPart());
}

int ReverseBClip::offsetFromThisEnd(string referenceName, FaidxWrapper &faidx)
{
}

int ReverseBClip::offsetFromThatEnd(string referenceName, FaidxWrapper &faidx, int orignal)
{
}

ForwardEClip::ForwardEClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const string &sequence, const std::vector<CigarOp> &cigar)
    : AbstractClip(referenceId, mapPosition, clipPosition, matePosition, sequence, cigar) {
}

string ForwardEClip::getType()
{
    return "3F";
}

Deletion ForwardEClip::call(FaidxWrapper &faidx, const std::vector<TargetRegion> &regions, int minOverlap, double minIdentity)
{
    string s1 = regions[0].sequence(faidx);
    SequenceOverlap overlap = Overlapper::computeOverlapSW2(s1, sequence, minOverlap, minIdentity, ungapped_params);

    int delta = overlap.getOverlapLength() - lengthOfSoftclippedPart();
    int offset = 0;
    for (auto &ci: cigar) {
        if (ci.Type == 'D') offset += ci.Length;
        else if (ci.Type == 'I') offset -= ci.Length;
    }
    int leftBp = clipPosition - 1 - offset;
    int rightBp = regions[0].start + overlap.match[0].end - lengthOfSoftclippedPart() + 1;

    int len = leftBp - rightBp + 1;
    int start1 = delta > 0 ? leftBp - delta : leftBp;
    int start2 = delta > 0 ? leftBp : leftBp - delta;
    int end1 = delta > 0 ? rightBp - delta : rightBp;
    int end2 = delta > 0 ? rightBp : rightBp - delta;

    if (len <= Helper::SVLEN_THRESHOLD) {
        return Deletion(regions[0].referenceName, start1, start2, end1, end2, len, getType());
    }

    error("No deletion is found.");
}

void ForwardEClip::fetchSpanningRanges(BamReader &reader, int insLength, std::vector<IRange> &ranges, int minMapQual)
{
    ranges.push_back({clipPosition + 1, matePosition + 1});
}

void ForwardEClip::fecthSizesForSpanningPairs(BamReader &reader, int inslength, std::vector<int> &sizes)
{
}

void ForwardEClip::toTargetRegions(const string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions)
{
    // ran.end - insLength + 2 * length(), ran.end + length()
    int leftmostPos = clipPosition;
    int s = ranges[0].end - insLength + 2 * length();
    if (s < leftmostPos) s = leftmostPos;
    int e = ranges[0].end + length();
    if (s > e) return;
    regions.push_back({referenceName, s, e});
}

int ForwardEClip::lengthOfSoftclippedPart()
{
    return cigar[cigar.size() - 1].Length;
}

string ForwardEClip::softclippedPart()
{
    return sequence.substr(lengthOfMappedPart());
}

string ForwardEClip::mappedPart()
{
    return sequence.substr(0, lengthOfMappedPart());
}

int ForwardEClip::offsetFromThisEnd(string referenceName, FaidxWrapper &faidx)
{
}

int ForwardEClip::offsetFromThatEnd(string referenceName, FaidxWrapper &faidx, int orignal)
{
}
