#ifndef SOFTCLIPCOUNTER_H
#define SOFTCLIPCOUNTER_H

#include "SoftClipReader.h"

class SoftClipCounter
{
public:
    SoftClipCounter(SoftClipReader *reader, int radius);
    int count(int refId, int position);

private:
    SoftClipReader *reader;
    int radius;
};

#endif // SOFTCLIPCOUNTER_H
