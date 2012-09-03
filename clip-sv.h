#ifndef CLIPSV_INCLUDED
#define CLIPSV_INCLUDED

#include <string>
#include <vector>
#include "api/BamReader.h"
#include "Clip.h"

void CountAlignments(BamTools::BamReader& reader);
void getClips(BamTools::BamReader& reader, std::vector<Clip>& leftClips, std::vector<Clip>& rightClips);
void countClipLength(std::vector<Clip>& clips, int lenValues[], int binWidth);
void writeData(ofstream& output, int lenValues[], int nBins);
bool findFirstClipInRange(const std::vector<Clip>& clips, int min, int max, Clip& cl);
int countMismatches(const std::string& s1, const std::string& s2);
int countOverlappedReads(const std::vector<Clip>& clips, const Clip cl);

#endif // CLIPSV_INCLUDED
