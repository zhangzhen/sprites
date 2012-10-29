#include <math.h>
#include <ctime>
#include <algorithm>
#include "clip-sv.h"
#include "error.h"

void callSVs(BamTools::BamReader& reader, std::string sv_filename, int minlen);
void outputClips(BamTools::BamReader& reader);

int main(int argc, char *argv[]) {
  char *sv_filename = NULL;
  char *sv_minlen = NULL;
  int c;

  opterr = 0;
  while ((c = getopt(argc, argv, "f:l:")) != -1)
    switch (c) {
      case 'f':
        sv_filename = optarg;
        break;
      case 'l':
        sv_minlen = optarg;
        break;
      case '?':
        if (optopt == 'f' || optopt =='l')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint(optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort();        
    }

  if (optind == argc)
  {
    fprintf(stderr, "A bam file is required.\n");
    return 1;
  }
  std::string bam_name(argv[optind]);
  BamTools::BamReader reader;
  if (!reader.Open(bam_name)) {
    std::cerr << "Could not open input BAM file." << std::endl;
    return 1;
  }
  // outputClips(reader);
  callSVs(reader, std::string(sv_filename), atoi(sv_minlen));
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
  std::cout << "#Calls: " << groups.size() << std::endl;

  std::vector<StructVar> trueSvs;
  getTrueSvs(sv_filename, trueSvs);
  ShorterThan st(minlen);
  trueSvs.erase(remove_if(trueSvs.begin(), trueSvs.end(), st), trueSvs.end());
  std::vector<StructVar> calls;
  evaluateCalls(calls, trueSvs);

  freeClips(leftClips);
  freeClips(rightClips);
}
