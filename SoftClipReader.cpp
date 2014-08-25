#include "SoftClipReader.h"
#include "error.h"

#include <vector>

using namespace std;
using namespace BamTools;

SoftClipReader::SoftClipReader(const string &filename, int minClip, int mode) :
    minClip(minClip), mode(mode) {
    if (!reader.Open(filename))
        error("Could not open the input BAM file.");
    if (!reader.LocateIndex())
        error("Could not locate the index file");
}

SoftClipReader::~SoftClipReader() {
    reader.Close();
}

int SoftClipReader::getReferenceId(const string &referenceName) {
    return reader.GetReferenceID(referenceName);
}

bool SoftClipReader::getSoftClip(SoftClip &clip) {
    BamAlignment al;
    while (reader.GetNextAlignment(al)) {
        vector<int> clipSizes, readPositions, genomePositions;
        if (!al.GetSoftClips(clipSizes, readPositions, genomePositions)) continue;
        int size = clipSizes.size();

        if (inEnhancedMode()) {
            if (!al.IsReverseStrand() && al.IsMateReverseStrand() && al.Position < al.MatePosition &&
                    al.Position != genomePositions[size - 1] && clipSizes[size - 1] > minClip &&
                    (size == 1 || (size == 2 && clipSizes[0] <= minClip))) {
                clip = SoftClip(al.RefID,
                                al.Position + 1,
                                al.Position + 1 - ((size == 2) ? clipSizes[0] : 0),
                                genomePositions[size - 1] + 1,
                                al.MatePosition + 1,
                                al.IsReverseStrand(),
                                al.IsMateReverseStrand(),
                                clipSizes[size - 1],
                                al.QueryBases);
                return true;
            }
            if (al.IsReverseStrand() && !al.IsMateReverseStrand() && al.Position > al.MatePosition &&
                    al.Position == genomePositions[0] && clipSizes[0] > minClip &&
                    (size == 1 || (size == 2 && clipSizes[1] <= minClip))) {
                clip = SoftClip(al.RefID,
                                al.Position + 1,
                                al.Position + 1 - clipSizes[0],
                                genomePositions[0] + 1,
                                al.MatePosition + 1,
                                al.IsReverseStrand(),
                                al.IsMateReverseStrand(),
                                clipSizes[0],
                                al.QueryBases);
                return true;
            }
        } else if (al.IsProperPair()) {
            if (!al.IsReverseStrand() && al.Position == genomePositions[0] &&
                    clipSizes[0] > minClip &&
                    (size == 1 ||
                     (size == 2 && clipSizes[1] <= minClip))) {
                clip = SoftClip(al.RefID,
                                al.Position + 1,
                                al.Position - clipSizes[0] + 1,
                                genomePositions[0] + 1,
                                al.MatePosition + 1,
                                al.IsReverseStrand(),
                                al.IsMateReverseStrand(),
                                clipSizes[0],
                                al.QueryBases);
                return true;
            }
            if (al.IsReverseStrand() && al.Position != genomePositions[size - 1] &&
                    clipSizes[size - 1] > minClip &&
                    (size == 1 ||
                     (size == 2 && clipSizes[0] <= minClip))) {
                clip = SoftClip(al.RefID,
                                al.Position + 1,
                                al.Position + 1 - ((size == 2) ? clipSizes[0] : 0),
                                genomePositions[size - 1] + 1,
                                al.MatePosition + 1,
                                al.IsReverseStrand(),
                                al.IsMateReverseStrand(),
                                clipSizes[size - 1],
                                al.QueryBases);
                return true;
            }
        }

    }

    return false;
}

bool SoftClipReader::setRegion(int leftRefId, int leftPosition, int rightRefId, int rightPosition) {
    return reader.SetRegion(leftRefId, leftPosition, rightRefId, rightPosition);
}

bool SoftClipReader::inEnhancedMode() const {
    return mode == 1;
}
