#ifndef SOFTCLIPREADER_H
#define SOFTCLIPREADER_H

#include "SoftClip.h"
#include "api/BamReader.h"
#include "clip.h"

#include <string>

class SoftClipReader
{
public:
    // 0 indicates the standard mode and 1 indicates the enhanced mode, which reads reads of type 2 besides type 1
    SoftClipReader(const std::string& filename, int minClip, int mode);
    virtual ~SoftClipReader();

    int getReferenceId(const std::string& referenceName);

    bool getSoftClip(SoftClip& clip);
    bool setRegion(int leftRefId, int leftPosition, int rightRefId, int rightPosition);

    int getMinClip() const;

private:
    BamTools::BamReader reader;
    int minClip;
    int mode;

    bool inEnhancedMode() const;
};

class ClipReader {
public:
    // 0 indicates the standard mode and 1 indicates the enhanced mode, which reads reads of type 2 besides type 1
    ClipReader(const std::string& filename, int allowedNum, int mode);
    virtual ~SoftClipReader();

    bool setRegion(int leftRefId, int leftPosition, int rightRefId, int rightPosition);

    int getReferenceId(const std::string& referenceName);
    std::string getReferenceName(int referenceId);

    int getAllowedNum() const;

    AbstractClip* nextClip();

private:
    BamTools::BamReader reader;
    int allowedNum;
    int mode;

    bool inEnhancedMode() const;
};

#endif // SOFTCLIPREADER_H
