#include "ClipReader.h"
#include "error.h"
#include "api/BamAlgorithms.h"

using namespace std;
using namespace BamTools;

ClipReader::ClipReader(const string &filename, int allowedNum, int mode, int minMapQual, int isizeCutoff)
    : allowedNum(allowedNum), mode(mode), minMapQual(minMapQual), isizeCutoff(isizeCutoff)
{
    if (!reader.Open(filename))
        error("Could not open the input BAM file.");
    if (!reader.LocateIndex())
        error("Could not locate the index file");
}

ClipReader::~ClipReader()
{
    reader.Close();
}

bool ClipReader::setRegion(int leftRefId, int leftPosition, int rightRefId, int rightPosition)
{
    return reader.SetRegion(leftRefId, leftPosition, rightRefId, rightPosition);
}

int ClipReader::getReferenceId(const string &referenceName)
{
    return reader.GetReferenceID(referenceName);
}

string ClipReader::getReferenceName(int referenceId)
{
    assert(referenceId >= 0 && referenceId < reader.GetReferenceCount());
    return reader.GetReferenceData()[referenceId].RefName;
}

AbstractClip *ClipReader::nextClip() {
    BamAlignment al;
    while (reader.GetNextAlignment(al)) {
        vector<int> clipSizes, readPositions, genomePositions;
//        if (!al.GetSoftClips(clipSizes, readPositions, genomePositions)) continue;
        if (al.MapQuality < minMapQual || !al.GetSoftClips(clipSizes, readPositions, genomePositions)) continue;
        int size = clipSizes.size();

        if (al.IsProperPair()) {
            if (!al.IsReverseStrand() && al.Position == genomePositions[0] &&
                    clipSizes[0] >= allowedNum &&
                    (size == 1 ||
                     (size == 2 && clipSizes[1] <= 5))) {
                return new ForwardBClip(al.RefID,
                                        al.Position + 1,
                                        genomePositions[0] + 1,
                                        al.MatePosition + 1,
                                        al.QueryBases,
                                        al.CigarData);
            }
            if (al.IsReverseStrand() && al.Position != genomePositions[size - 1] &&
                    clipSizes[size - 1] >= allowedNum &&
                    (size == 1 ||
                     (size == 2 && clipSizes[0] <= 5))) {
                return new ReverseEClip(al.RefID,
                                        al.Position + 1,
                                        genomePositions[size - 1] + 1,
                                        al.MatePosition + 1,
                                        al.QueryBases,
                                        al.CigarData);
            }
        }

        if (inEnhancedMode()) {
            if (al.RefID != al.MateRefID || abs(al.InsertSize) <= isizeCutoff)
                continue;
            if ((al.AlignmentFlag == 161 || al.AlignmentFlag == 97) && al.Position < al.MatePosition &&
                    clipSizes[size - 1] >= allowedNum &&
                    (size == 1 || (size == 2 && clipSizes[0] <= 5))) {
                return new ForwardEClip(al.RefID,
                                        al.Position + 1,
                                        genomePositions[size - 1] + 1,
                                        al.MatePosition + 1,
                                        al.QueryBases,
                                        al.CigarData);
            }
            if ((al.AlignmentFlag == 81 || al.AlignmentFlag == 145) && al.Position > al.MatePosition &&
                    clipSizes[0] >= allowedNum &&
                    (size == 1 || (size == 2 && clipSizes[1] <= 5))) {
                return new ReverseBClip(al.RefID,
                                        al.Position + 1,
                                        genomePositions[0] + 1,
                                        al.MatePosition + 1,
                                        al.QueryBases,
                                        al.CigarData);
            }
        }

    }
    return NULL;
}

bool ClipReader::inEnhancedMode() const
{
    return mode == 1;
}
