#include "SoftClipCounter.h"
#include "error.h"

using namespace std;

SoftClipCounter::SoftClipCounter(SoftClipReader *reader, int radius) :
    reader(reader), radius(radius)
{
}

int SoftClipCounter::countLeftBp(int refId, int position, int plus)
{
    int cnt = 0;
    if (!reader->setRegion(refId, position - radius, refId, position + radius / 2))
        error("could not set region");
    SoftClip clip;
    while (reader->getSoftClip(clip, plus)) {
        if (clip.isForLeftBp()) {
//            cout << "\t" << clip.getPosition() << "\t" << clip.getClipPosition() << endl;
            cnt++;
        }
    }
    return cnt;
}

int SoftClipCounter::countRightBp(int refId, int position, int plus)
{
    int cnt = 0;
    if (!reader->setRegion(refId, position - radius / 2, refId, position + radius))
        error("could not set region");
    SoftClip clip;
    while (reader->getSoftClip(clip, plus)) {
        if (clip.isForLeftBp()) {
//            cout << "\t" << clip.getPosition() << "\t" << clip.getClipPosition() << endl;
            cnt++;
        }
    }
    return cnt;
}
