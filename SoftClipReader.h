#ifndef SOFTCLIPREADER_H
#define SOFTCLIPREADER_H

#include "SoftClip.h"
#include "api/BamReader.h"

class SoftClipReader
{
public:
    SoftClipReader(const std::string& filename);
    virtual ~SoftClipReader();

    bool getSoftClip(SoftClip& clip);

private:
    BamTools::BamReader reader;
};

#endif // SOFTCLIPREADER_H
