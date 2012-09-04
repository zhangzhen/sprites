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
  reader.Close();
  // std::cout << "#left breakpoints: " << leftClips.size() << std::endl;
  // std::cout << "#right breakpoints: " << rightClips.size() << std::endl;
  try {
    Clip* cl = findFirstClipInRange(leftClips, 6, 10);
    std::cout << "You got " << countOverlappedReads(rightClips, cl) << " clips matched." << std::endl;
  } catch (ErrorException & ex) {
    std::cerr << ex.getMessage() << std::endl;
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
