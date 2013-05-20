#include <cmath>
#include <ctime>
#include <algorithm>
#include "error.h"
#include "DFinder.h"

const std::string ControlFilename = "../sv/chr22_report.txt";

int main(int argc, char *argv[]) {
  char *progname;
  int mean = 200;
  int std = 10;
  int minOverlapLen = 10;
  double maxMismatchRate = 0.1;
  std::string outFilename;
  int c, status = 0;

  progname = argv[0];
  while ((c = getopt(argc, argv, "m:s:l:x:o:")) != -1)
    switch (c) {
      case 'm':
        mean = atoi(optarg);
        break;
      case 's':
        std = atoi(optarg);
        break;
      case 'l':
        minOverlapLen = atoi(optarg);
        break;
      case 'x':
        maxMismatchRate = atof(optarg);
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
    std::cerr << "Usage: " << progname << "[options...] <BAM>" << std::endl
              << "Options: "
              << " -o FILE\tWrite output to <file>, required"
              << " -m NUM\tSpecify the mean of insert size, required"
              << " -s\tStd of insert size [-l minimal length of ] [-x mismatchRate] -o outputFilename bamFilename"
              << std::endl;
    return status;
  }

  std::string filename(argv[optind]);
  DFinder dfinder(filename, mean, std, minOverlapLen, maxMismatchRate);
  dfinder.callTo(outFilename);
  return 0;
}

// void callDelsFromBam(BamTools::BamReader& reader,
//                      const std::vector<Region2>& regs,
//                      std::string output,
//                      int maxMismatches,
//                      int minOverlapLen) {
//   time_t startTime;
//   double elapsedTime;

//   // Step 1: Loading clipped reads
//   std::vector<SingleClipped*> lefts, rights;
//   // startTime = time(NULL);
//   loadClippeds(reader, lefts, rights);
//   std::cout << "#rights: " << rights.size() << std::endl;  
//   std::cout << "#lefts " << lefts.size() << std::endl;  
  
//   // elapsedTime = difftime(time(NULL), startTime);
//   // std::cout << "Execution time of Step 1: "
//   //           << elapsedTime
//   //           << " (sec)"
//   //           << std::endl;

//   // Step 2: Clustering two collections of clipped reads respectively
//   StandardClusterCreator<RightClippedCluster> cluCreator1;
//   StandardClusterCreator<LeftClippedCluster> cluCreator2;
//   std::vector<SingleClippedCluster*> clus1, clus2;
//   // startTime = time(NULL);
//   sort(rights.begin(), rights.end(), compSC);
//   clusterClippeds(rights, clus1, cluCreator1);
//   std::cout << "#clusters1: " << clus1.size() << std::endl;
//   // for (size_t i = 0; i < clus1.size(); ++i)
//   //   std::cout << *clus1[i] << std::endl;
//   // std::cout << *clus1[0];
//   sort(lefts.begin(), lefts.end(), compSC);
//   clusterClippeds(lefts, clus2, cluCreator2);
//   std::cout << "#clusters2: " << clus2.size() << std::endl;
//   // std::vector<Region> controls;
//   // loadControls(ControlFilename, controls, MinDelLen);
//   // std::cout << "Minimal Distance between two adjacent variants: "
//   //           << minDistance(controls)
//   //           << std::endl;
//   // showControlContexts(controls, clus1, clus2);
//   // return;
//   // Locus anc1("22", 15000726);
//   // SingleClippedCluster* clu1 = cluCreator1.createCluster(anc1);
//   // std::cout << **lower_bound(clus1.begin(), clus1.end(), clu1, comp) << std::endl;
//   // Locus anc2("22", 15001350);
//   // SingleClippedCluster* clu2 = cluCreator2.createCluster(anc2);
//   // std::cout << **lower_bound(clus2.begin(), clus2.end(), clu2, comp) << std::endl;
//   // std::cout << *clus2[0];
//   // elapsedTime = difftime(time(NULL), startTime);
//   // std::cout << "Execution time of Step 2: "
//   //           << elapsedTime
//   //           << " (sec)"
//   //           << std::endl;

//   // Step 3: Obtaining contigs from two collections of clusters
//   std::vector<Contig> cons1, cons2;
//   // startTime = time(NULL);
//   obtainContigs(clus1, cons1);
//   std::cout << "#contigs1: " << cons1.size() << std::endl;    
//   obtainContigs(clus2, cons2);
//   std::cout << "#contigs2: " << cons2.size() << std::endl;
//   // elapsedTime = difftime(time(NULL), startTime);
//   // std::cout << "Execution time of Step 3: "
//   //           << elapsedTime
//   //           << " (sec)"
//   //           << std::endl;

//   // Step 4: Calling SVs
//   std::vector<Region> calls;
//   startTime = time(NULL);
//   DeletionCaller::callAll(regs, cons1, cons2, calls, minOverlapLen, maxMismatches);
//   std::cout << "#calls: " << calls.size() << std::endl;
//   elapsedTime = difftime(time(NULL), startTime);
//   std::cout << "Execution time of Step 4: "
//             << elapsedTime
//             << " (sec)"
//             << std::endl;

//   // Step 5: Outputing results to a file
//   // startTime = time(NULL);
//   outputCalls(output, calls);
//   // elapsedTime = difftime(time(NULL), startTime);
//   // std::cout << "Execution time of Step 5: "
//   //           << elapsedTime
//   //           << " (sec)"
//   //           << std::endl;
// }
