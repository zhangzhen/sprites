#include "SoftClipReader.h"
#include "error.h"

#include <vector>

using namespace std;
using namespace BamTools;

SoftClipReader::SoftClipReader(const string &filename, int minClip) :
    minClip(minClip)
{
    if (!reader.Open(filename))
        error("Could not open the input BAM file.");
    if (!reader.LocateIndex())
        error("Could not locate the index file");
}

SoftClipReader::~SoftClipReader()
{
    reader.Close();
}

int SoftClipReader::getReferenceId(const string &referenceName)
{
    return reader.GetReferenceID(referenceName);
}

bool SoftClipReader::getSoftClip(SoftClip &clip)
{
    BamAlignment al;
    while (reader.GetNextAlignment(al)) {
        vector<int> clipSizes, readPositions, genomePositions;
        if (al.IsProperPair() &&
                al.GetSoftClips(clipSizes, readPositions, genomePositions))
        {
            int size = clipSizes.size();
            if (!al.IsReverseStrand() && al.Position == genomePositions[0] &&
                    clipSizes[0] > minClip &&
                    (size == 1 ||
                     (size == 2 && clipSizes[1] <= minClip)))
            {
                clip = SoftClip(al.RefID,
                                al.Position,
                                al.Position - clipSizes[0],
                                genomePositions[0],
                                al.MatePosition,
                                al.IsReverseStrand(),
                                clipSizes[0],
                                0,
                                al.QueryBases);
                return true;
            }
            if (al.IsReverseStrand() && al.Position != genomePositions[size - 1] &&
                    clipSizes[size - 1] > minClip &&
                    (size == 1 ||
                     (size == 2 && clipSizes[0] <= minClip)))
            {
                clip = SoftClip(al.RefID,
                                al.Position,
                                al.Position - (size == 2) ? clipSizes[0] : 0,
                                genomePositions[size - 1],
                                al.MatePosition,
                                al.IsReverseStrand(),
                                clipSizes[size - 1],
                                0,
                                al.QueryBases);
                return true;
            }
        }
    }

    return false;
}

bool SoftClipReader::setRegion(int leftRefId, int leftPosition, int rightRefId, int rightPosition)
{
    return reader.SetRegion(leftRefId, leftPosition, rightRefId, rightPosition);
}
