#ifndef CLIPREADER_H
#define CLIPREADER_H

#include "clip.h"

class ClipReader
{
public:
    // 0 indicates the standard mode and 1 indicates the enhanced mode, which reads reads of type 2 besides type 1
    ClipReader(const std::string& filename, int allowedNum, int mode, int minMapQual, int isizeCutoff);
    virtual ~ClipReader();

    bool setRegion(int leftRefId, int leftPosition, int rightRefId, int rightPosition);

    int getReferenceId(const std::string& referenceName);
    std::string getReferenceName(int referenceId);

    int getAllowedNum() const;

    AbstractClip* nextClip();

private:    
    BamTools::BamReader reader;
    int allowedNum;
    int mode;
    int minMapQual;
    int isizeCutoff;

    bool inEnhancedMode() const;
};

#endif // CLIPREADER_H
