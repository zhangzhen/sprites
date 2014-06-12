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
    if (!reader->setRegion(refId, position - radius, refId, position))
        error("could not set region");
    SoftClip clip;
    while (reader->getSoftClip(clip, plus)) {
        int diff = position - clip.getLeftmostPosition();
//        if (clip.isForLeftBp() && diff > 0 && diff < clip.getSequence().length()) {
        if (clip.isForLeftBp()) {
            cout << "\t" << clip.getClipPosition() - 1 << "\t" << clip.getClipPosition() - 1 - position << "\t" << clip.getSequence() << endl;
            cnt++;
        }
    }
    return cnt;
}

int SoftClipCounter::countRightBp(int refId, int position, int plus)
{
    int cnt = 0;
    if (!reader->setRegion(refId, position, refId, position + radius))
        error("could not set region");
    SoftClip clip;
    while (reader->getSoftClip(clip, plus)) {
        if (clip.isForRightBp()) {
            cout << "\t" << clip.getClipPosition() - 1 << "\t" << clip.getClipPosition() - 1 - position << "\t" << clip.getSequence() << endl;
            cnt++;
        }
    }
    return cnt;
}
