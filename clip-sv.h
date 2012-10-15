#ifndef CLIPSV_INCLUDED
#define CLIPSV_INCLUDED

#include <string>
#include <vector>
#include <set>
#include "api/BamReader.h"
#include "Clip.h"

const int MAX_BINS = 50;
const int READ_LENGTH = 70;

class Breakpoint;

void countAlignments(BamTools::BamReader& reader);
bool compareClips(Clip* one, Clip* two);
void getClips(BamTools::BamReader& reader, std::vector<Clip*>& leftClips, std::vector<Clip*>& rightClips);
void countClipsInLength(const std::vector<Clip*>& clips, int stats[], int binWidth);
void outputData(std::ostream& output, int lenValues[], int nBins);
Clip* findFirstClipInRange(const std::vector<Clip*>& clips, int min, int max);
bool isOverlapped(Clip* c1, Clip* c2);
int countMismatches(const std::string& s1, const std::string& s2);
int countOverlappedReads(const std::vector<Clip*>& clips, Clip* cl);
void countClipsInLengthOneToFive(const std::vector<Clip*>& clips, int stats[]);
bool findMateIndex(Clip* lc, std::vector<Clip*>& RCs, int& index);
bool binarySearch2(int key, std::vector<int>& vec, int& index);
void extractClipsForDels(std::vector<Clip*>& inLCs, std::vector<Clip*>& inRCs, std::vector<Clip*>& outLCs, std::vector<Clip*>& outRCs);
void createOverlapGraph(const std::vector<Clip*>& LCs, const std::vector<Clip*>& RCs, std::vector<std::vector<int> >& g);
// void matchClips(std::vector<Clip*>& cls1, std::vector<Clip*>& cls2, std::map<Clip*, std::vector<Clip*> >& matches);
// void clusterClips(const std::vector<Clip*>& clips, std::map<int, std::set<Clip*> >& clusters);
void buildBreakpoints(const std::vector<Clip*>& LCs, const std::vector<Clip*>& RCs, std::vector<Breakpoint>& bps);
void clusterBreakpoints(const std::vector<Breakpoint>& bps, std::vector<std::vector<Breakpoint> >& clusters);
void freeClips(std::vector<Clip*> cl);

class Matching {
 public:
  int match() {
    int cnt = 0;
    for (int i = 0; i < m; i++) {
      seen.assign(n, 0);
      if (bpm(i)) cnt++;
    }
    return cnt;
  }
  
  Matching(const std::vector<std::vector<int> >& graph, int m, int n) : graph(graph), m(m), n(n) {
    matchL.assign(m, -1);
    matchR.assign(n, -1);
  }

  int getMateL(int v) {
    return matchL[v];
  }

  int getMateR(int v) {
    return matchR[v];
  }
  
 private:
  bool bpm(int u) {
    for (std::vector<int>::iterator itr = graph[u].begin(); itr != graph[u].end(); ++itr) {
      int v = *itr;
      if (seen[v]) continue;
      seen[v] = true;
      if (matchR[v] < 0 || bpm(matchR[v])) {
        matchL[u] = v;
        matchR[v] = u;
        return true;
      }
    }
    return false;
  }
  int m, n;
  std::vector<std::vector<int> > graph;
  std::vector<int> matchL, matchR;
  std::vector<bool> seen;
};

class Breakpoint {
 public:
  Breakpoint(int x, int y) : x(x), y(y) {
  }
  
  int center() const {
    return x + (y - x)/2;
  }

  int length() const {
    return y - x;
  }

  bool inSameCluster(const Breakpoint& other) const {
    if (center() == other.center() && length() == other.length()) return true;
    return false;
  }

 private:
  int x, y;
};
#endif // CLIPSV_INCLUDED
