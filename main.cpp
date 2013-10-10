#include <cmath>
#include <ctime>
#include <algorithm>
#include "error.h"
#include "DFinder.h"

// const std::string ControlFilename = "../sv/chr22_report.txt";

int main(int argc, char *argv[]) {
  char *progname;
  int mean = 200;
  int std = 10;
  int minOverlapLen = 15;
  double maxMismatchRate = 0.1;
  double discordant = 5.0;
  std::string outFilename;
  int c, status = 0;

  progname = argv[0];
  while ((c = getopt(argc, argv, "m:s:l:x:d:o:")) != -1)
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
      case 'd':
        discordant = atof(optarg);
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
  DFinder dfinder(filename, mean, std, minOverlapLen, maxMismatchRate, discordant);
  dfinder.callToFile(outFilename);
  // dfinder.checkAgainstGoldStandard("../gold-standard-svs/venter-chr22-known-dels.txt");
  // dfinder.checkAgainstGoldStandard("../gold-standard-svs/NA19312-chr20-known-dels.txt");
  // dfinder.printOverlaps("../gold-standard-svs/venter-chr22-known-dels.txt", 75);
  return 0;
}
