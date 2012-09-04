#include <math.h>
#include "clip-sv.h"
#include "error.h"

const int MAX_BINS = 50;
const int READ_LENGTH = 75;

void CountAlignments(BamTools::BamReader& reader) {
  int count = 0;
  int count2 = 0;

  BamTools::BamAlignment al;
  while (reader.GetNextAlignmentCore(al)) {
    std::vector<int> clipSizes, readPositions, genomePositions;
    count++;
    if (al.GetSoftClips(clipSizes, readPositions, genomePositions))
      count2++;
  }

  std::cout << "You got " << count << " alignments." << std::endl;
  std::cout << "You got " << count2 << " alignments with soft clipping." << std::endl;

}

void getClips(BamTools::BamReader& reader, std::vector<Clip*>& leftClips, std::vector<Clip*>& rightClips) {
  BamTools::BamAlignment al;
  while (reader.GetNextAlignment(al)) {
    std::vector<int> clipSizes, readPositions, genomePositions;
    if (!al.GetSoftClips(clipSizes, readPositions, genomePositions)) {
      continue;
    }
    if (clipSizes.size() > 1) {
      continue;
    }
    if (al.Position == genomePositions[0]) { // left clip
      leftClips.push_back(new LeftClip(al.RefID, genomePositions[0], readPositions[0], clipSizes[0], al.QueryBases));
    } else if (al.Length == readPositions[0]+clipSizes[0]) { // right clip
      rightClips.push_back(new RightClip(al.RefID, genomePositions[0], readPositions[0], clipSizes[0], al.QueryBases));
    }
  }
}

void countClipLength(std::vector<Clip>& clips, int lenValues[], int binWidth) {
  for (int i = 0; i < MAX_BINS; i++) {
    lenValues[i] = 0;
  }
  for (std::vector<Clip>::iterator itr = clips.begin(); itr != clips.end(); ++itr) {
    int bin = (int)((*itr).getSize() - 1) / binWidth;
    lenValues[bin]++;
  }
}

void writeData(std::ofstream& output, int lenValues[], int nBins) {
  for (int i = 0; i < nBins; i++) {
    output << lenValues[i] << "\t";
  }
  output << std::endl;
}

Clip* findFirstClipInRange(const std::vector<Clip*>& clips, int min, int max) {
  for (std::vector<Clip*>::const_iterator itr = clips.begin(); itr != clips.end(); ++itr) {
    int s = (*itr)->getSize();
    if (s >= min && s <= max) {
      return *itr;
    }
  }
  
  error("No such clip available.");
}

int countMismatches(const std::string& s1, const std::string& s2) {
  int count = 0;
  for (int i = 0; i < s1.size(); i++) {
    if (s1[i] != s2[i]) {
      count++;
    }
  }
  return count;
}

int countOverlappedReads(const std::vector<Clip*>& clips, const Clip* cl) {
  int count = 0;
  for (std::vector<Clip*>::const_iterator itr = clips.begin(); itr != clips.end(); ++itr) {
    int n = (*itr)->getSize();
    if (cl->getSize() + n > READ_LENGTH) {
      continue;
    }
    if (countMismatches(cl->getReadSeq(), (*itr)->getReadSeq()) < 2) {
      count++;
    }
  }
  return count;
}
