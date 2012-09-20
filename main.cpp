#include <math.h>
#include <ctime>
#include "clip-sv.h"
#include "error.h"

void callSVs(std::string filename);
void callSVs2(std::string filename);

int main(int argc, char *argv[]) {
  std::string filename(argv[1]);
  callSVs(filename);
  return 0;
}

void callSVs(std::string filename) {
  std::vector<Clip*> leftClips;
  std::vector<Clip*> rightClips;
  BamTools::BamReader reader;
  if (!reader.Open(filename)) {
    std::cerr << "Could not open input BAM file." << std::endl;
    return;
  }
  clock_t startTime, stopTime;
  double elapsedTime;
  
  startTime = clock();
  getClips(reader, leftClips, rightClips);
  stopTime = clock();
  elapsedTime = double(stopTime - startTime);
  std::cout << "getClips() elapsed execution time: " << elapsedTime << " (sec)" << std::endl;
  std::cout << "#left breakpoints: " << leftClips.size() << std::endl;
  std::cout << "#right breakpoints: " << rightClips.size() << std::endl;
  
  std::vector<Clip*> LCs, RCs;
  startTime = clock();
  extractClipsForDels(leftClips, rightClips, LCs, RCs);
  stopTime = clock();
  elapsedTime = double(stopTime - startTime);
  std::cout << "extractClipsForDels() elapsed execution time: " << elapsedTime << " (sec)" << std::endl;
  
  std::vector<Breakpoint> bps;
  startTime = clock();
  buildBreakpoints(LCs, RCs, bps);
  stopTime = clock();
  elapsedTime = double(stopTime - startTime);
  std::cout << "buildBreakpoints() elapsed execution time: " << elapsedTime << " (sec)" << std::endl;
  std::cout << "#Breakpoints: " << bps.size() << std::endl;

  std::vector<std::vector<Breakpoint> > clusters;
  startTime = clock();
  clusterBreakpoints(bps, clusters);
  stopTime = clock();
  elapsedTime = double(stopTime - startTime);
  std::cout << "clusterBreakpoints() elapsed execution time: " << elapsedTime << " (sec)" << std::endl;
  std::cout << "#Deletions: " << clusters.size() << std::endl;

  freeClips(leftClips);
  freeClips(rightClips);
}

void callSVs2(std::string filename) {
  std::vector<Clip*> leftClips;
  std::vector<Clip*> rightClips;
  BamTools::BamReader reader;
  if (!reader.Open(filename)) {
    std::cerr << "Could not open input BAM file." << std::endl;
    return;
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
}
