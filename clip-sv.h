#ifndef CLIPSV_INCLUDED
#define CLIPSV_INCLUDED

#include <string>
#include <vector>
#include <set>
#include "api/BamReader.h"
#include "Clip.h"

const int MAX_BINS = 50;
const int READ_LENGTH = 70;

void countAlignments(BamTools::BamReader& reader);
void getClips(BamTools::BamReader& reader, std::vector<Clip*>& leftClips, std::vector<Clip*>& rightClips);
void countClipsInLength(const std::vector<Clip*>& clips, int stats[], int binWidth);
void outputData(std::ostream& output, int lenValues[], int nBins);
Clip* findFirstClipInRange(const std::vector<Clip*>& clips, int min, int max);
int countMismatches(const std::string& s1, const std::string& s2);
int countOverlappedReads(const std::vector<Clip*>& clips, Clip* cl);
void countClipsInLengthOneToFive(const std::vector<Clip*>& clips, int stats[]);
void matchClips(const std::vector<Clip*>& cls1, const std::vector<Clip*>& cls2, std::map<Clip*, std::vector<Clip*> >& matches);
void clusterClips(const std::vector<Clip*>& clips, std::map<int, std::set<Clip*> >& clusters);

#endif // CLIPSV_INCLUDED
