#ifndef SOFTCLIPCOUNTER_H
#define SOFTCLIPCOUNTER_H

#include "SoftClipReader.h"

class SoftClipCounter
{
public:
    SoftClipCounter(SoftClipReader *reader, int radius);
    int countLeftBp(int refId, int position, int plus=false);
    int countRightBp(int refId, int position, int plus=false);

private:
    SoftClipReader *reader;
    int radius;
};

#endif // SOFTCLIPCOUNTER_H
