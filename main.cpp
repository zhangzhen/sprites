#include <math.h>
#include <ctime>
#include <algorithm>
#include "clip-sv.h"
#include "error.h"
#include "ClusterCreator.h"
#include "LeftClippedCluster.h"
#include "RightClippedCluster.h"

const std::string ControlFilename = "../sv/chr22_report.txt";
const int MinDelLen = 50;

void callDelsFromBam(BamTools::BamReader& reader,
                     std::string output,
                     double mismatchRate,
                     int minSupportSize,
                     int minOverlapLen);
void callSVs(BamTools::BamReader& reader, std::string sv_filename, int minlen);
void outputClips(BamTools::BamReader& reader);

int main(int argc, char *argv[]) {
  char *progname;
  int minSupportSize = 2;
  int minOverlapLen = 10;
  double mismatchRate = 0.0;
  std::string outFilename;
  int c, status = 0;

  progname = argv[0];
  while ((c = getopt(argc, argv, "k:l:x:o:")) != -1)
    switch (c) {
      case 'k':
        minSupportSize = atoi(optarg);
        break;
      case 'l':
        minOverlapLen = atoi(optarg);
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
              << " [-k minimalSupportSize] [-l minimalOverlapLength] [-x mismatchRate] -o outputFilename bamFilename"
              << std::endl;
    return status;
  }

  BamTools::BamReader r1, r2;  
  BamTools::BamReader reader;
  std::string bamFilename(argv[optind]);
  std::vector<Window> windows;
  r1.Open(bamFilename);
  r2.Open(bamFilename);
  r1.LocateIndex();
  r2.LocateIndex();
  loadWindowsFromBam(r1, r2, windows, 240);
  std::cout << windows.size() << std::endl;
  return 0;
  
  if (!reader.Open(bamFilename)) {
    std::cerr << "Could not open input BAM file." << std::endl;
    return 1;
  }

  callDelsFromBam(reader, outFilename, mismatchRate, minSupportSize, minOverlapLen);
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
                     int minSupportSize,
                     int minOverlapLen) {
  time_t startTime;
  double elapsedTime;

  // Step 1: Loading clipped reads
  std::vector<SingleClipped*> lefts, rights;
  // startTime = time(NULL);
  loadClippeds(reader, lefts, rights);
  std::cout << "#rights: " << rights.size() << std::endl;  
  std::cout << "#lefts " << lefts.size() << std::endl;  
  
  // elapsedTime = difftime(time(NULL), startTime);
  // std::cout << "Execution time of Step 1: "
  //           << elapsedTime
  //           << " (sec)"
  //           << std::endl;

  // Step 2: Clustering two collections of clipped reads respectively
  StandardClusterCreator<RightClippedCluster> cluCreator1;
  StandardClusterCreator<LeftClippedCluster> cluCreator2;
  std::vector<SingleClippedCluster*> clus1, clus2;
  // startTime = time(NULL);
  sort(rights.begin(), rights.end(), compSC);
  clusterClippeds(rights, clus1, cluCreator1);
  std::cout << "#clusters1: " << clus1.size() << std::endl;
  // for (size_t i = 0; i < clus1.size(); ++i)
  //   std::cout << *clus1[i] << std::endl;
  // std::cout << *clus1[0];
  sort(lefts.begin(), lefts.end(), compSC);
  clusterClippeds(lefts, clus2, cluCreator2);
  std::cout << "#clusters2: " << clus2.size() << std::endl;
  std::vector<Region> controls;
  loadControls(ControlFilename, controls, MinDelLen);
  // std::cout << "Minimal Distance between two adjacent variants: "
  //           << minDistance(controls)
  //           << std::endl;
  showControlContexts(controls, clus1, clus2);
  return;
  // Locus anc1("22", 15000726);
  // SingleClippedCluster* clu1 = cluCreator1.createCluster(anc1);
  // std::cout << **lower_bound(clus1.begin(), clus1.end(), clu1, comp) << std::endl;
  // Locus anc2("22", 15001350);
  // SingleClippedCluster* clu2 = cluCreator2.createCluster(anc2);
  // std::cout << **lower_bound(clus2.begin(), clus2.end(), clu2, comp) << std::endl;
  // std::cout << *clus2[0];
  // elapsedTime = difftime(time(NULL), startTime);
  // std::cout << "Execution time of Step 2: "
  //           << elapsedTime
  //           << " (sec)"
  //           << std::endl;

  // Step 3: Obtaining contigs from two collections of clusters
  std::vector<Contig> cons1, cons2;
  // startTime = time(NULL);
  obtainContigs(clus1, cons1);
  std::cout << "#contigs1: " << cons1.size() << std::endl;    
  obtainContigs(clus2, cons2);
  std::cout << "#contigs2: " << cons2.size() << std::endl;
  // elapsedTime = difftime(time(NULL), startTime);
  // std::cout << "Execution time of Step 3: "
  //           << elapsedTime
  //           << " (sec)"
  //           << std::endl;

  // Step 4: Calling SVs
  std::vector<Region> calls;
  startTime = time(NULL);
  callDeletions(cons1, cons2, calls, minSupportSize, minOverlapLen, mismatchRate);
  std::cout << "#calls: " << calls.size() << std::endl;
  elapsedTime = difftime(time(NULL), startTime);
  std::cout << "Execution time of Step 4: "
            << elapsedTime
            << " (sec)"
            << std::endl;

  // Step 5: Outputing results to a file
  // startTime = time(NULL);
  outputCalls(output, calls);
  // elapsedTime = difftime(time(NULL), startTime);
  // std::cout << "Execution time of Step 5: "
  //           << elapsedTime
  //           << " (sec)"
  //           << std::endl;
}
