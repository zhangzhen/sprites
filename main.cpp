#include <math.h>
#include "clip-sv.h"
#include "error.h"


int main(int argc, char *argv[]) {
  std::string filename(argv[1]);
  std::vector<Clip*> leftClips;
  std::vector<Clip*> rightClips;
  BamTools::BamReader reader;
  if (!reader.Open(filename)) {
    std::cerr << "Could not open input BAM file." << std::endl;
    exit(1);
  }
  getClips(reader, leftClips, rightClips);
  std::cout << "#left breakpoints: " << leftClips.size() << std::endl;
  std::cout << "#right breakpoints: " << rightClips.size() << std::endl;
  // try {
  //   Clip* cl = findFirstClipInRange(leftClips, 6, 10);
  //   std::cout << "You got " << countOverlappedReads(rightClips, cl) << " clips matched." << std::endl;
  // } catch (ErrorException & ex) {
  //   std::cerr << ex.getMessage() << std::endl;
  // }
  int stats[5];
  std::ofstream output;
  output.open("results.txt");
  countClipsInLengthOneToFive(leftClips, stats);
  outputData(output, stats, 5);
  countClipsInLengthOneToFive(rightClips, stats);
  outputData(output, stats, 5);
  output.close();
  
  int binWidth = 5;
  int nBins = ceil(READ_LENGTH/binWidth);
  int lenValues[MAX_BINS];
  std::ofstream output2;
  output2.open("results2.txt");

  countClipsInLength(leftClips, lenValues, binWidth);
  outputData(output2, lenValues, nBins);
  countClipsInLength(rightClips, lenValues, binWidth);
  outputData(output2, lenValues, nBins);
  output2.close();
  
  // reader.Rewind();
  // countAlignments(reader);
  reader.Close();
  return 0;
}
