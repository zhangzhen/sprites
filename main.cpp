#include <math.h>
#include <ctime>
#include <algorithm>
#include "clip-sv.h"
#include "error.h"
#include "ClusterCreator.h"
#include "LeftClippedCluster.h"
#include "RightClippedCluster.h"

void callDelsFromBam(BamTools::BamReader& reader,
                     std::string output,
                     double mismatchRate,
                     int minClusterSize,
                     int minLen,
                     int maxLen);
void callSVs(BamTools::BamReader& reader, std::string sv_filename, int minlen);
void outputClips(BamTools::BamReader& reader);

int main(int argc, char *argv[]) {
  char *progname;
  int minClusterSize = 2;
  int minCallLen = 0;
  int maxCallLen = 5000000;
  double mismatchRate = 0.0;
  std::string outFilename;
  int c, status = 0;

  progname = argv[0];
  while ((c = getopt(argc, argv, "c:l:m:x:o:")) != -1)
    switch (c) {
      case 'c':
        minClusterSize = atoi(optarg);
        break;
      case 'l':
        minCallLen = atoi(optarg);
        break;
      case 'm':
        maxCallLen = atoi(optarg);
        break;
      case 'x':
        mismatchRate = atof(optarg);
        break;
      case 'o':
        outFilename = std::string(optarg);
        break;
      case '?':
      default:
        status = 1;
        break;
    }

  if (optind == argc || outFilename.empty()) { status = 1; }
  if (status) {
    std::cerr << "Usage: "
              << progname
              << " [-c minimalClusterSize] [-l minimalCallLength] [-m maximalCallLength] [-x mismatchRate] -o outputFilename bamFilename"
              << std::endl;
    return status;
  }
  
  BamTools::BamReader reader;
  std::string bamFilename(argv[optind]);
  if (!reader.Open(bamFilename)) {
    std::cerr << "Could not open input BAM file." << std::endl;
    return 1;
  }

  callDelsFromBam(reader, outFilename, mismatchRate, minClusterSize, minCallLen, maxCallLen);
  reader.Close();
  return 0;
}

void outputClips(BamTools::BamReader& reader) {
  std::vector<Clip*> leftClips;
  std::vector<Clip*> rightClips;
  getClips(reader, leftClips, rightClips);
  tofile("left_clips.txt", leftClips);
  tofile("right_clips.txt", rightClips);

  freeClips(leftClips);
  freeClips(rightClips);
}

bool compareBreakpoints(Breakpoint one, Breakpoint two) {
  if(one.getX() != two.getX())
    return one.getX() < two.getX();
  return one.getY() < two.getY();
}

class ShorterThan
{
 public:
  /* Accept and store an int parameter */
  explicit ShorterThan(size_t maxLength) : length(maxLength) {}
  /* Return whether the string length is less than the stored int. */
  bool operator() (const StructVar& sv) const
  {
    return sv.length() < length;
  }
 private:
  const size_t length;
};

void callSVs(BamTools::BamReader& reader, std::string sv_filename, int minlen) {
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
  extractClipsForDels(leftClips, rightClips, LCs, RCs, 12, 75);
  elapsedTime = difftime(time(NULL), startTime);
  std::cout << "extractClipsForDels() elapsed execution time: " << elapsedTime << " (sec)" << std::endl;
  std::cout << "#left clips for deletions: " << LCs.size() << std::endl;
  std::cout << "#right rights for deletions: " << RCs.size() << std::endl;
  
  std::vector<Breakpoint> bps;
  startTime = time(NULL);
  buildBreakpoints(LCs, RCs, bps);
  elapsedTime = difftime(time(NULL), startTime);
  std::cout << "buildBreakpoints() elapsed execution time: " << elapsedTime << " (sec)" << std::endl;
  std::cout << "#Breakpoints: " << bps.size() << std::endl;

  sort(bps.begin(), bps.end(), compareBreakpoints);

  std::vector<std::vector<Breakpoint> > groups;
  // startTime = time(NULL);
  groupBreakpoints(bps, groups);
  // elapsedTime = difftime(time(NULL), startTime);
  // std::cout << "groupBreakpoints() elapsed execution time: " << elapsedTime << " (sec)" << std::endl;
  std::vector<StructVar> calls;
  makeCalls(groups, calls, minlen);
  
  std::cout << "#Calls: " << calls.size() << std::endl;

  std::vector<StructVar> trueSvs;
  getTrueSvs(sv_filename, trueSvs);
  ShorterThan st(minlen);
  trueSvs.erase(remove_if(trueSvs.begin(), trueSvs.end(), st), trueSvs.end());
  std::cout << "#True SVs: " << trueSvs.size() << std::endl;
  evaluateCalls(calls, trueSvs);

  freeClips(leftClips);
  freeClips(rightClips);
}

void callDelsFromBam(BamTools::BamReader& reader,
                     std::string output,
                     double mismatchRate,
                     int minClusterSize,
                     int minLen,
                     int maxLen) {
  time_t startTime;
  double elapsedTime;

  // Step 1: Loading clipped reads
  std::vector<SingleClipped*> lefts, rights;
  startTime = time(NULL);
  loadClippeds(reader, lefts, rights);
  elapsedTime = difftime(time(NULL), startTime);
  std::cout << "Execution time of Step 1: "
            << elapsedTime
            << " (sec)"
            << std::endl;

  // Step 2: Clustering two collections of clipped reads respectively
  StandardClusterCreator<RightClippedCluster> cluCreator1;
  StandardClusterCreator<LeftClippedCluster> cluCreator2;
  std::vector<SingleClippedCluster*> clus1, clus2;
  startTime = time(NULL);
  clusterClippeds(rights, clus1, cluCreator1, minClusterSize);
  clusterClippeds(lefts, clus2, cluCreator2, minClusterSize);
  elapsedTime = difftime(time(NULL), startTime);
  std::cout << "Execution time of Step 2: "
            << elapsedTime
            << " (sec)"
            << std::endl;

  // Step 3: Obtaining contigs from two collections of clusters
  std::vector<Contig> cons1, cons2;
  startTime = time(NULL);
  obtainContigs(clus1, cons1);
  obtainContigs(clus2, cons2);
  elapsedTime = difftime(time(NULL), startTime);
  std::cout << "Execution time of Step 3: "
            << elapsedTime
            << " (sec)"
            << std::endl;

  // Step 4: Calling SVs and selecting them by length
  std::vector<Region> calls;
  startTime = time(NULL);
  callDeletions(cons1, cons2, calls, mismatchRate);
  selectCallsByLength(calls, minLen, maxLen);
  elapsedTime = difftime(time(NULL), startTime);
  std::cout << "Execution time of Step 4: "
            << elapsedTime
            << " (sec)"
            << std::endl;

  // Step 5: Outputing results to a file
  startTime = time(NULL);
  outputCalls(output, calls);
  elapsedTime = difftime(time(NULL), startTime);
  std::cout << "Execution time of Step 5: "
            << elapsedTime
            << " (sec)"
            << std::endl;
}
