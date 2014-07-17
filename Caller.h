#ifndef CALLER_H
#define CALLER_H

#include "SoftClip.h"
#include "Deletion.h"
#include "Thirdparty/overlapper.h"
#include "Parameters.h"
#include "SoftClipReader.h"

#include <string>
#include <vector>

struct TargetRegion
{
    int referenceId;
    int start;
    int end;
};

class Caller
{
public:
    Caller(const std::string &filename, const Parameters& params);

    virtual ~Caller();

    void readClipsForRightBp(std::vector<SoftClip>& clips);

    void call(const std::vector<SoftClip>& clips, std::vector<Deletion>& dels);

    void output(const std::string& filename, const std::vector<Deletion>& dels);

private:
    void refineClips(const std::vector<SoftClip>& orig, std::vector<SoftClip>& result);
    SoftClip chooseBestClipFrom(std::vector<SoftClip> &buffer);
    void displayBuffer(const std::vector<SoftClip>& buffer);

    bool call(const SoftClip& clip, Deletion& del);

    bool call(const SoftClip& clip, const TargetRegion& region, Deletion& del);

    bool getSuppMatePositions(const SoftClip& clip, std::vector<int>& matePositions);
    TargetRegion getTargetRegion(const SoftClip &clip,
                         int matePosition);
    void getTargetRegions(const SoftClip &clip,
                          std::vector<int>& matePositions,
                          std::vector<TargetRegion>& regions);

    bool createDeletionFrom(const SequenceOverlap& overlap, const SoftClip &c1, const SoftClip &c2, Deletion& del);

    std::string getReferenceName(int referenceId) const;

    int getId();

    SoftClipReader *pReader;
    BamTools::BamReader bamReader;

    Parameters params;

    int idCount;

};

#endif // CALLER_H
