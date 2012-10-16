#include <algorithm>
#include "clip-sv.h"
#include "error.h"
const int READ_LENGTH = 75;

void countAlignments(BamTools::BamReader& reader) {
  int count = 0;
  int count2 = 0;
  int count3 = 0;

  BamTools::BamAlignment al;
  while (reader.GetNextAlignmentCore(al)) {
    std::vector<int> clipSizes, readPositions, genomePositions;
    count++;
    if (!al.IsMapped()) {
      count3++;
      continue;
    }
    if (al.GetSoftClips(clipSizes, readPositions, genomePositions))
      count2++;
    if (count <= 10) {
      std::cout << al.Position << std::endl;
    }
  }

  std::cout << "#alignments: " << count << std::endl;
  std::cout << "#alignments-soft clipping: " << count2 << " " << std::endl;
  std::cout << "#alignments-unmapped: " << count3 << " " << std::endl;

}

bool compareClips(Clip* one, Clip* two) {
  return one->getPosition() < two->getPosition();
}

void getClips(BamTools::BamReader& reader, std::vector<Clip*>& leftClips, std::vector<Clip*>& rightClips, int cutoff) {
  BamTools::BamAlignment al;
  // int report[5] = {0};
  while (reader.GetNextAlignment(al)) {
    std::vector<int> clipSizes, readPositions, genomePositions;
    if (!al.IsMapped()) {
      continue;
    }
    if (!al.GetSoftClips(clipSizes, readPositions, genomePositions)) {
      continue;
    }
    // if (clipSizes.size() > 0 && clipSizes.size() <=5)
    //   report[clipSizes.size()-1]++;
    if (clipSizes.size() > 1 or clipSizes[0] < cutoff) {
      continue;
    }
    if (al.Position == genomePositions[0]) { // left clip - or readPositions[i] == clipSizes[i]
      leftClips.push_back(new LeftClip(al.RefID, genomePositions[0], readPositions[0], clipSizes[0], al.QueryBases));
    } else if (al.Length == readPositions[0]+clipSizes[0]) { // right clip
      rightClips.push_back(new RightClip(al.RefID, genomePositions[0], readPositions[0], clipSizes[0], al.QueryBases));
    }
  }
  sort(leftClips.begin(), leftClips.end(), compareClips);
  sort(rightClips.begin(), rightClips.end(), compareClips);
}

std::map<int, size_t> generateClipReport(std::vector<Clip*> clips) {
  std::map<int, size_t> result;
  for (size_t i = 0; i < clips.size(); ++i)
    ++result[clips[i]->getSize()];
  return result;
}

void outputClipReport(std::string filename, std::map<int,size_t> report) {
  std::ofstream out(filename.c_str());
  std::cout << "length\tcount" << std::endl;
  for (std::map<int, size_t>::iterator itr = report.begin(); itr != report.end(); ++itr)
    std::cout << itr->first << itr->second << std::endl;
  out.close();
}

int countMismatches(const std::string& s1, const std::string& s2) {
  if (s1.size() != s2.size()) {
    error("Two strings must have the same size");
  }
  int count = 0;
  for (int i = 0; i < s1.size(); i++) {
    if (s1[i] != s2[i]) {
      count++;
    }
  }
  return count;
}

int countOverlappedReads(const std::vector<Clip*>& clips, Clip* cl) {
  int count = 0;
  for (std::vector<Clip*>::const_iterator itr = clips.begin(); itr != clips.end(); ++itr) {
    if ((*itr)->getType() == cl->getType()) {
      continue;
    }
    int len = cl->getSize() + (*itr)->getSize();
    if (len > READ_LENGTH) {
      continue;
    }
    if (isOverlapped(cl, *itr)) {
      count++;
    }
  }
  return count;
}

bool isOverlapped(Clip* c1, Clip* c2) {
  std::string ls, rs;
  int len = c1->getSize() + c2->getSize();
  if (c1->getType() == Left) {
    ls = c1->getReadSeq().substr(0, len);
    rs = c2->getReadSeq().substr(c2->getReadPosition()-c1->getSize());
  } else {
    ls = c2->getReadSeq().substr(0, len);
    rs = c1->getReadSeq().substr(c1->getReadPosition()-c2->getSize());
  }
  if (ls == rs) {                       // might consider mismatches between two strings
    return true;
  }
  return false;
}

bool findMateIndex(Clip* lc, std::vector<Clip*>& RCs, int& index) {
  int low = 0;
  int high = RCs.size()-1;
  int key = lc->getPosition();
  if ((RCs[low])->getPosition() >= key) return false;
  if ((RCs[high])->getPosition() < key) {
    index = high;
    return true;
  }
  while (high >= low) {
    int mid = low + (high - low)/2;
    int mval = (RCs[mid])->getPosition();
    if (mval < key && (RCs[mid+1])->getPosition() >= key) {
      index = mid;
      return true;
    }
    if (mval < key) low = mid + 1;
    else high = mid - 1;
  }
}

bool binarySearch2(int key, std::vector<int>& vec, int& index) {
  int low = 0;
  int high = vec.size()-1;
  if (vec[low] >= key) return false;
  if (vec[high] < key) {
    index = high;
    return true;
  }
  while (high >= low) {
    int mid = low + (high-low)/2;
    if (vec[mid] < key && vec[mid+1] >= key) {
      index = mid;
      return true;
    }
    if (vec[mid] < key) low = mid + 1;
    else high = mid - 1;
  }
}

void extractClipsForDels(std::vector<Clip*>& inLCs, std::vector<Clip*>& inRCs, std::vector<Clip*>& outLCs, std::vector<Clip*>& outRCs) {
  std::set<Clip*> LCs, RCs;
  int index;
  for (size_t i=0; i<inLCs.size(); ++i) {
    if (!findMateIndex(inLCs[i], inRCs, index)) continue;
    for (size_t j = 0; j <= index; ++j) {
      int len = inLCs[i]->getSize() + inRCs[j]->getSize();
      if (len > READ_LENGTH) {
        continue;
      }
      if (isOverlapped(inLCs[i], inRCs[j])) {
        if (!LCs.count(inLCs[i])) {
          LCs.insert(inLCs[i]);
        }
        if (!RCs.count(inRCs[j])) {
          RCs.insert(inRCs[j]);
        }
      }
    }
  }
  outLCs.assign(LCs.begin(), LCs.end());
  outRCs.assign(RCs.begin(), RCs.end());
}

void createOverlapGraph(const std::vector<Clip*>& LCs, const std::vector<Clip*>& RCs, std::vector<std::vector<int> >& g) {
  for (size_t i = 0; i < LCs.size(); ++i) {
    std::vector<int> nbs;
    for (size_t j = 0; j < RCs.size(); ++j) {
      if (isOverlapped(LCs[i], RCs[j]))
        nbs.push_back(j);
    }
    g.push_back(nbs);
  }
}

void countClipsInLengthOneToFive(const std::vector<Clip*>& clips, int stats[]) {
  for (int i=0; i<5; i++) {
    stats[i] = 0;
  }
  for (std::vector<Clip*>::const_iterator itr = clips.begin(); itr != clips.end(); ++itr) {
    if ((*itr)->getSize() <= 5) {
      stats[(*itr)->getSize()-1]++;
    }
  }
}

// void clusterClips(const std::vector<Clip*>& clips, std::map<int, std::set<Clip*> >& clusters) {
//   for (std::vector<Clip*>::const_iterator itr = clips.begin(); itr != clips.end(); ++itr) {
//     int pos = (*itr)->getPosition();
//     if (!clusters.count(pos))
//       clusters[pos] = std::set<Clip*>();
//     clusters[pos].insert(*itr);
//   }
// }

void buildBreakpoints(const std::vector<Clip*>& LCs, const std::vector<Clip*>& RCs, std::vector<Breakpoint>& bps) {
  int m, n;
  m = LCs.size();
  n = RCs.size();
  std::vector<std::vector<int> > g;
  createOverlapGraph(LCs, RCs, g);
  Matching ma(g, m, n);
  int cnt = ma.match();
  for (int i = 0; i < m; ++i) {
    int v = ma.getMateL(i);
    if (v == -1) continue;
    bps.push_back(Breakpoint(LCs[i]->getPosition(), RCs[v]->getPosition()));
  }
}

void clusterBreakpoints(const std::vector<Breakpoint>& bps, std::vector<std::vector<Breakpoint> >& clusters) {
  std::vector<int> labels(bps.size(), -1);
  int cnt = 0;
  for (int i = 0; i < bps.size(); ++i) {
    if (labels[i] != -1) continue;
    labels[i] = cnt;
    for (int j = i+1; j < bps.size(); ++j) {
      if(labels[j] == -1 && bps[i].inSameCluster(bps[j]))
        labels[j] = cnt;
    }
    cnt++;
  }
  for (int i = 0; i < cnt; ++i) {
    std::vector<Breakpoint> clu;
    clusters.push_back(clu);
  }
  for (int i = 0; i < bps.size(); ++i) {
    clusters[labels[i]].push_back(bps[i]);
  }
}

void freeClips(std::vector<Clip*> cl) {
  while (!cl.empty()) delete cl.back(), cl.pop_back();
}
