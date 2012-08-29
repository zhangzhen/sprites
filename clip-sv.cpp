#include <string>
#include <vector>
#include <math.h>
#include "api/BamReader.h"
#include "Clip.h"

using namespace std;
using namespace BamTools;

const int MAX_BINS = 50;

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

int main(int argc, char *argv[]) {
  string filename(argv[1]);
  vector<Clip> leftClips, rightClips;
  BamReader reader;
  if (!reader.Open(filename)) {
    cerr << "Could not open input BAM file." << endl;
    exit(1);
  }
  getClips(reader, leftClips, rightClips);
  reader.Close();
  cout << "#left breakpoints: " << leftClips.size() << endl;
  cout << "#right breakpoints: " << rightClips.size() << endl;

  int readLength = 75;
  int binWidth = 5;
  int nBins = ceil(readLength/binWidth);
  int lenValues[MAX_BINS];
  countClipLength(leftClips, lenValues, binWidth);
  ofstream output;
  output.open("results.txt");
  writeData(output, lenValues, nBins);

  //CountAlignments(filename);
  return 0;
}
