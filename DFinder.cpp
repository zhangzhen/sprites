#include <algorithm>
#include <queue>
#include <set>
#include <iterator>

#include "DFinder.h"
#include "ChrRegionCluster.h"
#include "statisticscalculator.h"
#include "error.h"

using namespace std;
using namespace BamTools;

DFinder::DFinder(const std::string& filename, const string& referenceName, int minOverlapLength, double maxMismatchRate) :
    referenceName(referenceName), minOverlapLength(minOverlapLength), maxMismatchRate(maxMismatchRate)
{
    if (!r1.Open(filename)) {
        error("could not open BAM file");
    }

    if (!r2.Open(filename)) {
        error("could not open BAM file");
    }
    if (!r2.LocateIndex()) {
        if (!r2.CreateIndex())
            error("could not find or create index");
    }

    referenceId = r1.GetReferenceID(referenceName);
    if (referenceId == -1)
    {
        error("could find reference id");
    }

    SamHeader header = r1.GetHeader();
    for (auto itr = header.ReadGroups.Begin(); itr != header.ReadGroups.End(); ++itr) {
        if (!libraries.count((*itr).Library)) libraries[(*itr).Library] = new Library((*itr).Library);
        readgroups[(*itr).ID] = new ReadGroup((*itr).ID, libraries[(*itr).Library]);
        libraries[(*itr).Library]->add((*itr).ID);
    }

    for (auto itr = libraries.begin(); itr != libraries.end(); ++itr) {
        itr->second->setReadLength(100);
        StatisticsCalculator sta(filename, itr->second);
        cout << "Libary: " << itr->second->getName() << endl;
        sta.calcMean();
        cout << "Mean of insert size: " << itr->second->getMeanOfInsertSize() << endl;
        sta.calcSd();
        cout << "SD of insert size: " << itr->second->getSdOfInsertSize() << endl;
    }

}

DFinder::~DFinder() {
    //    for (int i = 0; i < size; ++i) {
    //        for_each(leftClips[i].begin(), leftClips[i].end(), DeletePtr<SoftClip>());
    //        for_each(leftParts[i].begin(), leftParts[i].end(), DeletePtr<SoftClip>());
    //        for_each(rightClips[i].begin(), rightClips[i].end(), DeletePtr<SoftClip>());
    //        for_each(rightParts[i].begin(), rightParts[i].end(), DeletePtr<SoftClip>());
    //        for_each(intervals[i].begin(), intervals[i].end(), DeletePtr<ChrRegion>());
    //    }
    r1.Close();
    r2.Close();
}

bool DFinder::isValidAlignment(const BamTools::BamAlignment &al)
{
    return !al.IsDuplicate() &&
            al.IsPaired() &&
            al.IsMapped() &&
            al.IsMateMapped() &&
            al.RefID == referenceId &&
            al.RefID == al.MateRefID &&
            al.MapQuality >= MapQualityThreshold;
}

bool DFinder::isCorrectOrientation(const BamTools::BamAlignment &al)
{
    return (!al.IsReverseStrand() && al.IsMateReverseStrand() && al.MatePosition > al.Position) ||
            (al.IsReverseStrand() && !al.IsMateReverseStrand() && al.Position > al.MatePosition);
}

bool DFinder::isValidPartnerCandidateFor(const SoftClip &orig, const SoftClip& partner)
{
    return (orig.isLeftPartClipped() && !partner.isLeftPartClipped() && !partner.getOrientation()) ||
            (orig.isLeftPartClipped() && !partner.isLeftPartClipped() && partner.getOrientation() && partner.getMateGenomePosition() > orig.getClipPosition()) ||
            (!orig.isLeftPartClipped() && partner.isLeftPartClipped() && partner.getOrientation()) ||
            (!orig.isLeftPartClipped() && partner.isLeftPartClipped() && !partner.getOrientation() && partner.getMateGenomePosition() < orig.getClipPosition());
}

void DFinder::callToFile(const string& filename) {
    vector<Deletion> calls;
    call(calls);

    ofstream out(filename.c_str());
    out << "Chromosome\tType\tStart\tEnd\tLength" << endl;
    for(auto itr = calls.begin(); itr != calls.end(); ++itr) {
        out << referenceName << "\tDEL\t" << (*itr).getStart2() << "\t" << (*itr).getEnd2() << "\t" << (*itr).length()
                                                           << endl;
    }
}

void DFinder::call(vector<Deletion>& calls)
{
    BamAlignment al;
    string rg;
    while(r1.GetNextAlignment(al))
    {
        Deletion del;
        vector<int> clipSizes, readPositions, genomePositions;
        if (al.IsProperPair() &&
                al.GetSoftClips(clipSizes, readPositions, genomePositions) &&
                clipSizes.size() == 1 &&
                al.GetTag("RG", rg) &&
                callFor(SoftClip(al.RefID, al.Position, genomePositions[0],
                                 al.Position == genomePositions[0] ? readPositions[0] : al.Length - clipSizes[0],
                                 !al.IsReverseStrand(), al.MatePosition, !al.IsMateReverseStrand(),
                                 al.InsertSize, al.QueryBases, readgroups[rg]), del))
        {
            calls.push_back(del);
        }
    }
}

bool DFinder::callFor(const SoftClip& softclip, Deletion& del)
{
    if (softclip.isPartOfSpanningPair())
    {
        return callFor(softclip, softclip, del);
    }

    vector<SoftClip> spanners;
    getSpanningPairsFor(softclip, spanners);
    for (auto itr = spanners.begin(); itr != spanners.end(); ++itr)
    {
        if (callFor(softclip, *itr, del)) return true;
    }
    return false;
}

void DFinder::getSpanningPairsFor(const SoftClip &softclip, vector<SoftClip> &results)
{
    for (auto itr = libraries.begin(); itr != libraries.end(); ++itr)
    {
        int start;
        int end;
        Library lib = *itr->second;
        if (!softclip.searchRangeForSpanningPair(lib, start, end)) continue;

        if(!r2.SetRegion(referenceId, start, referenceId, end))
            error("could not set region");

        BamAlignment al;
        while (r2.GetNextAlignment(al))
        {
            string rg;
            vector<int> clipSizes, readPositions, genomePositions;
            if (isValidAlignment(al) &&
                    al.GetTag("RG", rg) &&
                    lib.contains(rg) &&
                    (!softclip.getOrientation() && !al.IsReverseStrand() && al.IsMateReverseStrand() && al.MatePosition > al.Position ||
                     softclip.getOrientation() && al.IsReverseStrand() && !al.IsMateReverseStrand() && al.Position > al.MatePosition) &&
                    !al.GetSoftClips(clipSizes, readPositions, genomePositions))
            {
                SoftClip cl(al.RefID, al.Position, al.Position, 0, !al.IsReverseStrand(), al.MatePosition,
                            !al.IsMateReverseStrand(), al.InsertSize, al.QueryBases, readgroups[rg]);
                if (cl.isDiscordantInsertSize())
                    results.push_back(cl);
            }
        }
    }
}

void DFinder::findPartnerCandidates(const SoftClip &softclip, int start, int end, vector<SoftClip>& softclips, vector<SoftClip>& reads)
{
    BamAlignment al;
    if (!r2.SetRegion(referenceId, start, referenceId, end))
    {
        error("could not set region");
    }
    while (r2.GetNextAlignment(al))
    {
        vector<int> clipSizes, readPositions, genomePositions;
        string rg;

        if (isValidAlignment(al) &&
                isCorrectOrientation(al) &&
                al.GetTag("RG", rg))
        {
            bool isSoftClip = al.GetSoftClips(clipSizes, readPositions, genomePositions);
            if (clipSizes.size() > 1) continue;

            int clipPosition, localClipPosition;

            if (!isSoftClip)
            {
                clipPosition = softclip.isLeftPartClipped() ? al.Position : al.GetEndPosition();
                localClipPosition = softclip.isLeftPartClipped() ? 0 : al.Length;
            }
            else
            {
                clipPosition = genomePositions[0];
                localClipPosition = softclip.isLeftPartClipped() ? al.Length - clipSizes[0] : readPositions[0];
            }

            SoftClip candidate(al.RefID, al.Position, clipPosition, localClipPosition, !al.IsReverseStrand(), al.MatePosition,
                               !al.IsMateReverseStrand(), al.InsertSize, al.QueryBases, readgroups[rg]);
            if (!isValidPartnerCandidateFor(softclip, candidate)) continue;
            if (!isSoftClip)
            {
                reads.push_back(candidate);
            }
            else
            {
                softclips.push_back(candidate);
            }
        }
    }
}

bool DFinder::callFor(const SoftClip &softclip, const SoftClip& spanner, Deletion &del)
{
    int start, end;
    spanner.searchRangeForSoftClip(softclip, lengthThreshold, start, end);
    vector<SoftClip> softclips, reads;
    findPartnerCandidates(softclip, start, end, softclips, reads);

    return callFor(softclip, softclips, del) || callFor(softclip, reads, del);
}

bool DFinder::callFor(const SoftClip &softclip, const vector<SoftClip> &candidates, Deletion &del)
{
    for (auto itr = candidates.begin(); itr != candidates.end(); ++itr)
    {
        Overlap overlap;
        if (softclip.overlaps(*itr, minOverlapLength, maxMismatchRate, overlap) && overlap.deletionLength() >=  lengthThreshold)
        {
            del = overlap.getDeletion();
            return true;
        }
    }
    return false;
}
