#include "Caller.h"
#include "error.h"
#include "SoftClipReader.h"

#include <cassert>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <iterator>

#include "Thirdparty/overlapper.h"

using namespace std;
using namespace BamTools;

Caller::Caller(const string& filename, const Parameters &params) :
    params(params) {
    idCount = 0;
    pReader = new SoftClipReader(filename, params.minClip, params.mode);

    if (!bamReader.Open(filename))
        error("Could not open the input BAM file.");
    if (!bamReader.LocateIndex())
        error("Could not locate index.");
}

Caller::~Caller() {
    delete pReader;
}

void Caller::readClipsForRightBp(std::vector<SoftClip> &clips) {
    SequenceOverlap overlap = Overlapper::computeOverlapSW(
                                  "TTATCCACCTGTTTCCACCACCTCTCCCTAGAGATGCAGCTGCTCTCATGCCCCTGAATT"
                                  "AGACTCTGAGCCTACAGAGTGCAGAGCCAGCCCAGGACAGGGGACAATTACACAGGCGAA"
                                  "GGAGGATGGAGCAGGTTGTCCATGGCCTCCCATCTGCTCCATCCTCCTTCGCCTCATTCA"
                                  "GGCGATGGTCCTAAGAACCGAACCTTCCAATCCCAAAACTCTAGACAGGTATCCAATACC"
                                  "TACTGTGTTTTTGTAGAAGAAGTACAGCACCATGTTGGCAAGTCGGGAGTAGCACCAATG"
                                  "C",
                                  "GCCTACAGAGTGCAGAGCCAGCCCAGGACAGGGGACAATTACACAGGCGATGGTCCTAAG"
                                  "AACCGAACCTTCCAATCCCAAAACTCTAGACAGGTATCCAA",
                                  ungapped_params);

    cout << endl;

    overlap.printAlignment("TTATCCACCTGTTTCCACCACCTCTCCCTAGAGATGCAGCTGCTCTCATGCCCCTGAATT"
                           "AGACTCTGAGCCTACAGAGTGCAGAGCCAGCCCAGGACAGGGGACAATTACACAGGCGAA"
                           "GGAGGATGGAGCAGGTTGTCCATGGCCTCCCATCTGCTCCATCCTCCTTCGCCTCATTCA"
                           "GGCGATGGTCCTAAGAACCGAACCTTCCAATCCCAAAACTCTAGACAGGTATCCAATACC"
                           "TACTGTGTTTTTGTAGAAGAAGTACAGCACCATGTTGGCAAGTCGGGAGTAGCACCAATG"
                           "C",
                           "GCCTACAGAGTGCAGAGCCAGCCCAGGACAGGGGACAATTACACAGGCGATGGTCCTAAGAACCGAACCTTCCAATCCCAAAACTCTAGACAGGTATCCAA");

    std::vector<SoftClip> buffer;

    SoftClip clip;
    while (pReader->getSoftClip(clip)) {
        if (clip.isForRightBp()) buffer.push_back(clip);
    }

    refineClips(buffer, clips);
}

void Caller::call(const std::vector<SoftClip> &clips, std::vector<Deletion> &dels) {
    for (auto it = clips.begin(); it != clips.end(); ++it) {
        Deletion del;
        if (call(*it, del)) {
            dels.push_back(del);
        }
    }
}

void Caller::output(const string &filename, const std::vector<Deletion> &dels) {
    ofstream out(filename.c_str());
    out << "ID\tCHROM\tSTART\tEND\tTYPE\tLENGTH\tALT\tHOMSEQ\tGENOTYPE" << endl;
    copy(dels.begin(), dels.end(), ostream_iterator<Deletion>(out, "\n"));
}

// require orig be ordered by clipping position
void Caller::refineClips(const std::vector<SoftClip> &orig, std::vector<SoftClip> &result) {
    std::vector<SoftClip> buffer;

    for (auto it = orig.begin(); it != orig.end(); ++it) {
        if (buffer.empty()) {
            buffer.push_back(*it);
            continue;
        }
        if (buffer[0].getClipPosition() == (*it).getClipPosition()) {
            buffer.push_back(*it);
            continue;
        }
        result.push_back(chooseBestClipFrom(buffer));
        displayBuffer(buffer);
        buffer.clear();
        buffer.push_back(*it);
    }
    if (!buffer.empty()) {
        result.push_back(chooseBestClipFrom(buffer));
        displayBuffer(buffer);
        buffer.clear();
    }
}

SoftClip Caller::chooseBestClipFrom(std::vector<SoftClip> &buffer) {
    assert(!buffer.empty());
    sort(buffer.begin(), buffer.end(), [](const SoftClip& c1, const SoftClip& c2) {
        return c1.getClippedSize() < c2.getClippedSize();
    });
    return buffer[buffer.size() - 1];
}

void Caller::displayBuffer(const std::vector<SoftClip> &buffer) {
    if (buffer.size() == 1) {
        //cout << buffer[0].getSequence() << endl;
        return;
    }
    int min = buffer[0].getClippedSize();
    int max = buffer[buffer.size() - 1].getClippedSize();

    cout << endl << ">>>>>>>>>>>>>>>\t" << buffer[0].getClipPosition() << "\t" << buffer.size() << endl;
    for (auto itr = buffer.begin(); itr != buffer.end(); ++itr) {
        string s;
        if ((*itr).isForRightBp()) {
            s.append(max - (*itr).getClippedSize(), ' ');
        } else {
            s.append((*itr).getClippedSize() - min, ' ');
        }
        s.append((*itr).getSequence());
        cout << s << endl;
    }

}

bool Caller::call(const SoftClip &clip, Deletion &del) {
    if (clip.isForRightBp()) {
        vector<int> matePositions;
        if (!getSuppMatePositions(clip, matePositions))
            return false;
        vector<TargetRegion> regions;
        getTargetRegions(clip, matePositions, regions);
        for (auto itr = regions.begin(); itr != regions.end(); ++itr) {
            if (call(clip, *itr, del)) {
                return true;
            }
        }
    }
    return false;
}

bool Caller::call(const SoftClip &clip, const TargetRegion &region, Deletion &del) {
    assert(clip.isForRightBp());

    vector<SoftClip> buffer;

    if (!pReader->setRegion(region.referenceId, region.start, region.referenceId, region.end))
        error("could not set region");

    SoftClip other;
    while (pReader->getSoftClip(other)) {
        if (other.isForLeftBp()) {
            buffer.push_back(other);
        }
    }

//    sort(buffer.begin(), buffer.end(), [](const SoftClip& c1, const SoftClip& c2) { return c1.getClipPosition() < c2.getClipPosition(); });
//    vector<SoftClip> refined;
//    refineClips(buffer, refined);

    for (auto it = buffer.begin(); it != buffer.end(); ++it) {
        SequenceOverlap overlap = Overlapper::computeOverlap(clip.getSequence(), (*it).getSequence(), ungapped_params);
        if (overlap.getOverlapLength() >= params.minOverlap &&
                overlap.getPercentIdentity() / 100 >= params.minIdentity) {
            if (createDeletionFrom(overlap, clip, (*it), del)) {
                cout << endl << ">>>>>>>>>>>>>>>" << endl;
                overlap.printAlignment(clip.getSequence(), (*it).getSequence());
                return true;
            }
        }
    }

    return false;
}

bool Caller::getSuppMatePositions(const SoftClip &clip, vector<int> &matePositions) {
    bool result = false;

    int start, end;
    if (!clip.isReverse()) {
        start = clip.getLeftmostPosition();
        end = start + params.insertMean + 3 * params.insertSd + clip.size();
    } else {
        end = clip.getLeftmostPosition() + clip.size();
        start = end > params.insertMean + 3 * params.insertSd + clip.size() ? end - params.insertMean - 3 * params.insertSd - clip.size() :
                0;
    }
    assert(start < end);

    if (!bamReader.SetRegion(clip.getReferenceId(), start, clip.getReferenceId(), end))
        error("Could not set the region.");

    BamAlignment al;
    while(bamReader.GetNextAlignmentCore(al)) {
        if ((!clip.isReverse() && al.IsReverseStrand() && al.Position > al.MatePosition) ||
                (clip.isReverse() && !al.IsMateReverseStrand() && al.Position < al.MatePosition)) {
            matePositions.push_back(al.MatePosition);
            result = true;
        }
    }
    return result;
}

TargetRegion Caller::getTargetRegion(const SoftClip& clip, int matePosition) {
    int start = !clip.isReverse() ? matePosition :
                matePosition - 3 * params.insertSd;
    int end = !clip.isReverse() ? matePosition + clip.size() + params.insertMean + 3 * params.insertSd :
              matePosition + clip.size();
    return {clip.getReferenceId(), start, end};
}

void Caller::getTargetRegions(const SoftClip& clip,
                              std::vector<int> &matePositions,
                              std::vector<TargetRegion> &regions) {
    sort(matePositions.begin(), matePositions.end());

    regions.push_back(getTargetRegion(clip, matePositions[0]));
    if (matePositions.size() == 1) {
        return;
    }

    std::vector<int> diffs(matePositions.size());

    adjacent_difference(matePositions.begin(), matePositions.end(), diffs.begin());

    for (size_t i = 1; i < diffs.size(); ++i) {
        int start;
        TargetRegion tr = getTargetRegion(clip, matePositions[i]);
        if (abs(diffs[i]) <= clip.size() + params.insertMean + 3 * params.insertSd) {
            start = regions.back().start;
            regions.pop_back();
        } else {
            start = tr.start;
        }
        regions.push_back({tr.referenceId, start, tr.end});
    }
}

bool Caller::createDeletionFrom(const SequenceOverlap &overlap, const SoftClip& c1, const SoftClip& c2, Deletion& del) {
    int diff = (c1.getClipPositionInRead() - overlap.match[0].start) -
               (c2.getClipPositionInRead() - overlap.match[1].start);
    int length = c2.getClipPosition() - c1.getClipPosition() + diff;
    if (length > SVLEN_THRESHOLD) return false;
    int leftBp = (diff < 0) ? c2.getClipPosition() + diff : c2.getClipPosition();
    del = Deletion(getId(), getReferenceName(c1.getReferenceId()), leftBp - 1, c1.getClipPosition() - 1, length);
    return true;
}

string Caller::getReferenceName(int referenceId) const {
    assert(referenceId >= 0 && referenceId < bamReader.GetReferenceCount());
    return bamReader.GetReferenceData()[referenceId].RefName;
}

int Caller::getId() {
    return idCount++;
}
