#include "softclipreader.h"
#include "error.h"

#include <vector>

using namespace std;
using namespace BamTools;

SoftClipReader::SoftClipReader(const string &filename)
{
    if (!reader.Open(filename))
        error("Could not open the input BAM file.");
}

SoftClipReader::~SoftClipReader()
{
    reader.Close();
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
                    clipSizes[0] > 5 &&
                    (size == 1 ||
                     (size == 2 && clipSizes[1] <= 5)))
            {
                clip = SoftClip(al.RefID,
                                al.Position,
                                genomePositions[0],
                                al.MatePosition,
                                al.IsReverseStrand(),
                                clipSizes[0],
                                0,
                                al.QueryBases);
                return true;
            }
            if (al.IsReverseStrand() && al.Position != genomePositions[size - 1] &&
                    clipSizes[size - 1] > 5 &&
                    (size == 1 ||
                     (size == 2 && clipSizes[0] <= 5)))
            {
                clip = SoftClip(al.RefID,
                                al.Position,
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
