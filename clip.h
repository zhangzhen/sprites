#ifndef CLIP_H
#define CLIP_H

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "Deletion.h"
#include "FaidxWrapper.h"
#include "range.h"
#include "Thirdparty/overlapper.h"

#include <string>
#include <vector>
#include <sstream>

struct TargetRegion
{
    std::string referenceName;
    int start;
    int end;

    std::string sequence(FaidxWrapper &faidx) const {
        return faidx.fetch(referenceName, start, end);
    }

    int length() const {
        return end - start + 1;
    }
};


class AbstractClip {
public:
    AbstractClip(int referenceId, int mapPosition, int clipPosition,
                 int matePosition, const std::string& sequence,
                 const std::vector<BamTools::CigarOp>& cigar);

    int length() const;

    int leftmostPosition() const;
    int getClipPosition() const {
        return clipPosition;
    }

    virtual ~AbstractClip();

    Deletion call(BamTools::BamReader& reader, FaidxWrapper &faidx, int insLength, int minOverlap, double minIdentity, int minMapQual);

    bool hasConflictWith(AbstractClip *other);
    virtual std::string getType() = 0;
    bool getConflictFlag() const;
    void setConflictFlag(bool value);
    std::string toString() {
        std::stringstream ss;
        ss << getClipPosition() << "\t" << getType();
        return ss.str();
    }

protected:

    virtual Deletion call(FaidxWrapper &faidx, const std::vector<TargetRegion>& regions, int minOverlap, double minIdentity) = 0;
    virtual void fetchSpanningRanges(BamTools::BamReader &reader, int insLength, std::vector<IRange> &ranges, int minMapQual) = 0;
    virtual void fecthSizesForSpanningPairs(BamTools::BamReader &reader, int inslength, std::vector<int>& sizes) = 0;
    virtual void toTargetRegions(const std::string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions) = 0;

    virtual int lengthOfSoftclippedPart() = 0;
    int lengthOfMappedPart() {
        return sequence.size() - lengthOfSoftclippedPart();
    }

    int maxEditDistanceForSoftclippedPart();

    virtual std::string softclippedPart() = 0;
    virtual std::string mappedPart() = 0;
    virtual int offsetFromThisEnd(std::string referenceName, FaidxWrapper& faidx) = 0;
    virtual int offsetFromThatEnd(std::string referenceName, FaidxWrapper& faidx, int orignal) = 0;

    int referenceId;
    int mapPosition;
    int clipPosition;
    int matePosition;
    std::string sequence;
    std::vector<BamTools::CigarOp> cigar;

    bool conflictFlag;
};

class ForwardBClip : public AbstractClip {
public:
    ForwardBClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const std::string& sequence, const std::vector<BamTools::CigarOp>& cigar);

private:
    virtual void fetchSpanningRanges(BamTools::BamReader &reader, int insLength, std::vector<IRange> &ranges, int minMapQual);
    virtual void fecthSizesForSpanningPairs(BamTools::BamReader& reader, int insLength, std::vector<int>& sizes);
    virtual void toTargetRegions(const std::string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions);

    virtual Deletion call(FaidxWrapper &faidx, const std::vector<TargetRegion>& regions, int minOverlap, double minIdentity);    

    // AbstractClip interface
public:
    std::string getType();

    // AbstractClip interface
protected:
    int lengthOfSoftclippedPart() {
        return cigar[0].Length;
    }

    std::string softclippedPart() {
        return sequence.substr(0, lengthOfSoftclippedPart());
    }

    std::string mappedPart() {
        return sequence.substr(lengthOfSoftclippedPart());
    }

    int offsetFromThisEnd(std::string referenceName, FaidxWrapper &faidx);
    int offsetFromThatEnd(std::string referenceName, FaidxWrapper &faidx, int orignal);
};

class ReverseBClip : public AbstractClip {

    // AbstractClip interface
public:
    ReverseBClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const std::string& sequence, const std::vector<BamTools::CigarOp>& cigar);
    std::string getType();

protected:
    Deletion call(FaidxWrapper &faidx, const std::vector<TargetRegion> &regions, int minOverlap, double minIdentity);
    void fetchSpanningRanges(BamTools::BamReader &reader, int insLength, std::vector<IRange> &ranges, int minMapQual);
    void fecthSizesForSpanningPairs(BamTools::BamReader &reader, int inslength, std::vector<int> &sizes);
    void toTargetRegions(const std::string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions);
    int lengthOfSoftclippedPart();
    std::string softclippedPart();
    std::string mappedPart();
    int offsetFromThisEnd(std::string referenceName, FaidxWrapper &faidx);
    int offsetFromThatEnd(std::string referenceName, FaidxWrapper &faidx, int orignal);
};

class ForwardEClip : public AbstractClip {

    // AbstractClip interface
public:
    ForwardEClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const std::string& sequence, const std::vector<BamTools::CigarOp>& cigar);
    std::string getType();

protected:
    Deletion call(FaidxWrapper &faidx, const std::vector<TargetRegion> &regions, int minOverlap, double minIdentity);
    void fetchSpanningRanges(BamTools::BamReader &reader, int insLength, std::vector<IRange> &ranges, int minMapQual);
    void fecthSizesForSpanningPairs(BamTools::BamReader &reader, int inslength, std::vector<int> &sizes);
    void toTargetRegions(const std::string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions);
    int lengthOfSoftclippedPart();
    std::string softclippedPart();
    std::string mappedPart();
    int offsetFromThisEnd(std::string referenceName, FaidxWrapper &faidx);
    int offsetFromThatEnd(std::string referenceName, FaidxWrapper &faidx, int orignal);
};

class ReverseEClip : public AbstractClip {
public:
    ReverseEClip(int referenceId, int mapPosition, int clipPosition, int matePosition, const std::string& sequence, const std::vector<BamTools::CigarOp>& cigar);

private:
    virtual void fetchSpanningRanges(BamTools::BamReader &reader, int insLength, std::vector<IRange> &ranges, int minMapQual);
    virtual void fecthSizesForSpanningPairs(BamTools::BamReader& reader, int insLength, std::vector<int>& sizes);
    virtual void toTargetRegions(const std::string &referenceName, int insLength, std::vector<IRange> &ranges, std::vector<TargetRegion> &regions);

    virtual Deletion call(FaidxWrapper &faidx, const std::vector<TargetRegion>& regions, int minOverlap, double minIdentity);

    // AbstractClip interface
public:
    std::string getType();

    // AbstractClip interface
protected:
    int lengthOfSoftclippedPart() {
        return cigar[cigar.size() - 1].Length;
    }

    std::string softclippedPart() {
        return sequence.substr(lengthOfMappedPart());
    }

    std::string mappedPart() {
        return sequence.substr(0, lengthOfMappedPart());
    }

    int offsetFromThisEnd(std::string referenceName, FaidxWrapper &faidx);
    int offsetFromThatEnd(std::string referenceName, FaidxWrapper &faidx, int orignal);
};

#endif // CLIP_H
