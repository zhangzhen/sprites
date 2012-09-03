#include "clip-sv.h"

int main(int argc, char *argv[]) {
  string filename(argv[1]);
  vector<Clip> leftClips, rightClips;
  BamTools::BamReader reader;
  if (!reader.Open(filename)) {
    cerr << "Could not open input BAM file." << endl;
    exit(1);
  }
  getClips(reader, leftClips, rightClips);
  reader.Close();
  // cout << "#left breakpoints: " << leftClips.size() << endl;
  // cout << "#right breakpoints: " << rightClips.size() << endl;

  Clip cl;
  if (findFirstClipInRange(leftClips, 6, 10, cl)) {
    cout << "You got " << countOverlappedReads(rightClips, cl) << " clips matched." << endl;
  }
  // int readLength = 75;
  // int binWidth = 5;
  // int nBins = ceil(readLength/binWidth);
  // int lenValues[MAX_BINS];
  // ofstream output;
  // output.open("results.txt");

  // countClipLength(leftClips, lenValues, binWidth);
  // writeData(output, lenValues, nBins);
  // countClipLength(rightClips, lenValues, binWidth);
  // writeData(output, lenValues, nBins);

  //CountAlignments(filename);
  return 0;
}
