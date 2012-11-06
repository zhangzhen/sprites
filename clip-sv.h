#ifndef CLIPSV_INCLUDED
#define CLIPSV_INCLUDED

#include <string>
#include <vector>
#include <set>
#include "api/BamReader.h"
#include "Clip.h"


class Breakpoint;
struct StructVar {
  std::string chr;
  int left;                             // 1-based, included
  int right;                            // 1-based, included

  StructVar(std::string chr, int left, int right) : chr(chr), left(left), right(right) {
  }

  StructVar(const StructVar& other) {
    chr = other.chr;
    left = other.left;
    right = other.right;
  }

  int length() const {
    return right - left + 1;
  }

  bool operator== (const StructVar& other) {
    return chr == other.chr && left == other.left && right == other.right;
  }

  bool operator< (const StructVar& other) const {
    if (chr != other.chr)
      return chr < other.chr;
    if (left != other.left)
      return left < other.left;
    return right < other.right;
  }
};


void countAlignments(BamTools::BamReader& reader);
bool compareClips(Clip* one, Clip* two);
void getClips(BamTools::BamReader& reader, std::vector<Clip*>& leftClips, std::vector<Clip*>& rightClips);
void tofile(std::string filename, const std::vector<Clip*>& clips);
bool equals(const std::string s1, const std::string s2, int mismatches=0);
bool isOverlapped(Clip* c1, Clip* c2, int mismatches=0);
bool findMateIndex(Clip* lc, std::vector<Clip*>& RCs, int& index);
void extractClipsForDels(std::vector<Clip*>& inLCs, std::vector<Clip*>& inRCs, std::vector<Clip*>& outLCs, std::vector<Clip*>& outRCs, int min, int max);
void createOverlapGraph(const std::vector<Clip*>& LCs, const std::vector<Clip*>& RCs, std::vector<std::vector<int> >& g);
// void clusterClips(const std::vector<Clip*>& clips, std::map<int, std::set<Clip*> >& clusters);
void buildBreakpoints(const std::vector<Clip*>& LCs, const std::vector<Clip*>& RCs, std::vector<Breakpoint>& bps);
void groupBreakpoints(const std::vector<Breakpoint>& bps, std::vector<std::vector<Breakpoint> >& groups);
void makeCalls(const std::vector<std::vector<Breakpoint> >& groups, std::vector<StructVar>& calls, int minlen);
// void clusterBreakpoints(const std::vector<Breakpoint>& bps, std::vector<std::vector<Breakpoint> >& clusters);
void getTrueSvs(std::string filename, std::vector<StructVar>& trueSvs);
std::set<StructVar> findOverlaps(StructVar t, const std::vector<StructVar>& cs1, const std::vector<StructVar>& cs2);
bool cmp1(StructVar sv1, StructVar sv2);
bool cmp2(StructVar sv1, StructVar sv2);
void evaluateCalls(const std::vector<StructVar>& calls, const std::vector<StructVar>& trueSvs);
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

  int getX() const {                    // 1-based, included
    return x;
  }

  int getY() const {                    // 1-based, not included
    return y;
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
