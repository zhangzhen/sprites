#ifndef CLIP_H
#define CLIP_H

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "Deletion.h"
#include "FaidxWrapper.h"
#include "range.h"

#include <string>
#include <vector>

struct TargetRegion
{
    std::string referenceName;
    int start;
    int end;

    std::string sequence(FaidxWrapper &faidx) const;
};


class AbstractClip {
public:
    AbstractClip(int referenceId, int mapPosition, int clipPosition,
                 int matePosition, const std::string& sequence,
                 const std::vector<BamTools::CigarOp>& cigar);

    int length() const;

    int leftmostPosition() const;

    virtual ~AbstractClip();

    Deletion call(BamTools::BamReader& reader, FaidxWrapper &faidx, int insLength, int minOverlap, double minIdentity);

protected:

    virtual Deletion call(FaidxWrapper &faidx, const std::vector<TargetRegion>& regions, int minOverlap, double minIdentity) = 0;
    virtual void fetchSpanningRanges(BamTools::BamReader &reader, int insLength, std::vector<IRange> &ranges) = 0;
    virtual void toTargetRegions(const std::string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions) = 0;

    int referenceId;
    int mapPosition;
    int clipPosition;
    int matePosition;
    std::string sequence;
    std::vector<BamTools::CigarOp> cigar;
};

class ForwardBClip : public AbstractClip {
public:
    ForwardBClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const std::string& sequence, const std::vector<BamTools::CigarOp>& cigar);

private:
    virtual void fetchSpanningRanges(BamTools::BamReader &reader, int insLength, std::vector<IRange> &ranges);
    virtual void toTargetRegions(const std::string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions);

    virtual Deletion call(FaidxWrapper &faidx, const std::vector<TargetRegion>& regions, int minOverlap, double minIdentity);
};

class ReverseBClip : public AbstractClip {
public:
    ReverseBClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const std::string& sequence, const std::vector<BamTools::CigarOp>& cigar);

private:
    virtual void fetchSpanningRanges(BamTools::BamReader &reader, int insLength, std::vector<IRange> &ranges);
    virtual void toTargetRegions(const std::string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions);

    virtual Deletion call(FaidxWrapper &faidx, const std::vector<TargetRegion>& regions, int minOverlap, double minIdentity);
};

class ForwardEClip : public AbstractClip {
public:
    ForwardEClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const std::string& sequence, const std::vector<BamTools::CigarOp>& cigar);

private:
    virtual void fetchSpanningRanges(BamTools::BamReader &reader, int insLength, std::vector<IRange> &ranges);
    virtual void toTargetRegions(const std::string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions);

    virtual Deletion call(FaidxWrapper &faidx, const std::vector<TargetRegion>& regions, int minOverlap, double minIdentity);
};

class ReverseEClip : public AbstractClip {
public:
    ReverseEClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const std::string& sequence, const std::vector<BamTools::CigarOp>& cigar);

private:
    virtual void fetchSpanningRanges(BamTools::BamReader &reader, int insLength, std::vector<IRange> &ranges);
    virtual void toTargetRegions(const std::string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions);

    virtual Deletion call(FaidxWrapper &faidx, const std::vector<TargetRegion>& regions, int minOverlap, double minIdentity);
};

#endif // CLIP_H
