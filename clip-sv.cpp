#include <math.h>
#include "clip-sv.h"

using namespace std;
using namespace BamTools;

const int MAX_BINS = 50;
const int READ_LENGTH = 75;

void CountAlignments(BamReader& reader) {
  int count = 0;
  int count2 = 0;

  BamAlignment al;
  while (reader.GetNextAlignmentCore(al)) {
    vector<int> clipSizes, readPositions, genomePositions;
    count++;
    if (al.GetSoftClips(clipSizes, readPositions, genomePositions))
      count2++;
  }

  cout << "You got " << count << " alignments." << endl;
  cout << "You got " << count2 << " alignments with soft clipping." << endl;

}

void getClips(BamReader& reader, vector<Clip>& leftClips, vector<Clip>& rightClips) {
  BamAlignment al;
  while (reader.GetNextAlignment(al)) {
    vector<int> clipSizes, readPositions, genomePositions;
    if (!al.GetSoftClips(clipSizes, readPositions, genomePositions)) {
      continue;
    }
    if (clipSizes.size() > 1) {
      continue;
    }
    Clip cl(al.RefID, genomePositions[0], readPositions[0], clipSizes[0], al.QueryBases);
    if (al.Position == genomePositions[0]) { // left clip
      leftClips.push_back(cl);
    } else if (al.Length == readPositions[0]+clipSizes[0]) { // right clip
      rightClips.push_back(cl);
    }
  }
}

void countClipLength(vector<Clip>& clips, int lenValues[], int binWidth) {
  for (int i = 0; i < MAX_BINS; i++) {
    lenValues[i] = 0;
  }
  for (vector<Clip>::iterator itr = clips.begin(); itr != clips.end(); ++itr) {
    int bin = (int)((*itr).getSize() - 1) / binWidth;
    lenValues[bin]++;
  }
}

void writeData(ofstream& output, int lenValues[], int nBins) {
  for (int i = 0; i < nBins; i++) {
    output << lenValues[i] << "\t";
  }
  output << endl;
}

bool findFirstClipInRange(const vector<Clip>& clips, int min, int max, Clip& cl) {
  for (vector<Clip>::const_iterator itr = clips.begin(); itr != clips.end(); ++itr) {
    int s = (*itr).getSize();
    if (s >= min && s <= max) {
      cl = *itr;
      return true;
    }
  }
  
  return false;
}

int countMismatches(const string& s1, const string& s2) {
  int count = 0;
  for (int i = 0; i < s1.size(); i++) {
    if (s1[i] != s2[i]) {
      count++;
    }
  }
  return count;
}

int countOverlappedReads(const vector<Clip>& clips, const Clip cl) {
  int count = 0;
  for (vector<Clip>::const_iterator itr = clips.begin(); itr != clips.end(); ++itr) {
    int n = (*itr).getSize();
    if (cl.getSize() + n > READ_LENGTH) {
      continue;
    }
    if (countMismatches(cl.getReadSeq(), (*itr).getReadSeq()) < 2) {
      count++;
    }
  }
  return count;
}
