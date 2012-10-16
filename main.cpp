#include <math.h>
#include <ctime>
#include "clip-sv.h"
#include "error.h"

void callSVs(BamTools::BamReader& reader);
void outputClips(BamTools::BamReader& reader);

int main(int argc, char *argv[]) {
  BamTools::BamReader reader;
  std::string filename(argv[1]);
  if (!reader.Open(filename)) {
    std::cerr << "Could not open input BAM file." << std::endl;
    return -1;
  }
  // outputClips(reader);
  callSVs(reader);
  reader.Close();
  return 0;
}

void outputClips(BamTools::BamReader& reader) {
  std::vector<Clip*> leftClips;
  std::vector<Clip*> rightClips;
  getClips(reader, leftClips, rightClips, 5);
  tofile("left_clips.txt", leftClips);
  tofile("right_clips.txt", rightClips);

  freeClips(leftClips);
  freeClips(rightClips);
}

void callSVs(BamTools::BamReader& reader) {
  std::vector<Clip*> leftClips;
  std::vector<Clip*> rightClips;
  time_t startTime;
  double elapsedTime;
  
  // startTime = time(NULL);
  getClips(reader, leftClips, rightClips);
  // elapsedTime = difftime(time(NULL), startTime);
  // std::cout << "getClips() elapsed execution time: " << elapsedTime << " (sec)" << std::endl;
  std::cout << "#left clips: " << leftClips.size() << std::endl;
  std::cout << "#right rights: " << rightClips.size() << std::endl;
  
  std::vector<Clip*> LCs, RCs;
  startTime = time(NULL);
  extractClipsForDels(leftClips, rightClips, LCs, RCs);
  elapsedTime = difftime(time(NULL), startTime);
  std::cout << "extractClipsForDels() elapsed execution time: " << elapsedTime << " (sec)" << std::endl;
  
  std::vector<Breakpoint> bps;
  startTime = time(NULL);
  buildBreakpoints(LCs, RCs, bps);
  elapsedTime = difftime(time(NULL), startTime);
  std::cout << "buildBreakpoints() elapsed execution time: " << elapsedTime << " (sec)" << std::endl;
  std::cout << "#Breakpoints: " << bps.size() << std::endl;

  std::vector<std::vector<Breakpoint> > clusters;
  // startTime = time(NULL);
  clusterBreakpoints(bps, clusters);
  // elapsedTime = difftime(time(NULL), startTime);
  // std::cout << "clusterBreakpoints() elapsed execution time: " << elapsedTime << " (sec)" << std::endl;
  std::cout << "#Deletions: " << clusters.size() << std::endl;

  freeClips(leftClips);
  freeClips(rightClips);
}
