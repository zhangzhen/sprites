#ifndef SOFTCLIPREADER_H
#define SOFTCLIPREADER_H

#include "SoftClip.h"
#include "api/BamReader.h"

#include <string>

class SoftClipReader
{
public:
    SoftClipReader(const std::string& filename, int minClip);
    virtual ~SoftClipReader();

    int getReferenceId(const std::string& referenceName);

    bool getSoftClip(SoftClip& clip, bool plus=false);
    bool setRegion(int leftRefId, int leftPosition, int rightRefId, int rightPosition);

private:
    BamTools::BamReader reader;
    int minClip;
};

#endif // SOFTCLIPREADER_H
