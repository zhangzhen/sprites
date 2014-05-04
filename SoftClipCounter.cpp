#include "SoftClipCounter.h"
#include "error.h"

SoftClipCounter::SoftClipCounter(SoftClipReader *reader, int radius) :
    reader(reader), radius(radius)
{
}

int SoftClipCounter::count(int refId, int position)
{
    int cnt = 0;
    if (!reader->setRegion(refId, position - radius, refId, position + radius))
        error("could not set region");
    SoftClip clip;
    while (reader->getSoftClip(clip)) {
        cnt++;
    }
    return cnt;
}
