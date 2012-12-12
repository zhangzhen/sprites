#include "LeftClipped.h"
#include "RightClipped.h"
#include "ClippedCreator.h"
#include <algorithm>
#include "clip-sv.h"
#include "error.h"

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

void loadClippeds(BamTools::BamReader& reader,
                  std::vector<SingleClipped*>& lefts,
                  std::vector<SingleClipped*>& rights) {
  BamTools::BamAlignment al;
  StandardClippedCreator<LeftClipped> lCreator;
  StandardClippedCreator<RightClipped> rCreator;
  
  while (reader.GetNextAlignment(al)) {
    std::vector<int> lens, cutpoints, anchors;
    if (!al.IsMapped()) {
      continue;
    }
    if (!al.GetSoftClips(lens, cutpoints, anchors)) {
      continue;
    }
    if (lens.size() > 1) {
      continue;
    }
    std::stringstream ss;
    ss << al.RefID;
    Locus anchor(ss.str(), anchors[0]);
    if (cutpoints[0] == lens[0]) { // left clipped
      lefts.push_back(lCreator.createClipped(anchor, al.QueryBases, al.MapQuality, 0, lens[0]));
    } else if (al.Length == cutpoints[0] + lens[0]) { // right clipped
      rights.push_back(rCreator.createClipped(anchor, al.QueryBases, al.MapQuality, cutpoints[0], lens[0]));
    }
  }  
}

void clusterClippeds(std::vector<SingleClipped*>& clis,
                     std::vector<SingleClippedCluster*>& clus,
                     ClusterCreator& creator,
                     int cutoff) {
  assert(clis.size() > 0);
  SingleClippedCluster *clu = creator.createCluster(clis[0]->anchor());
  for (int i = 0; i < clis.size() - 1; ++i) {
    clu->add(clis[i]);
    if (clis[i]->anchor() != clis[i + 1]->anchor()) {
      if (clu->size() >= cutoff) clus.push_back(clu);
      clu = creator.createCluster(clis[i + 1]->anchor());
    }
  }
  clu->add(clis.back());
  if (clu->size() >= cutoff) clus.push_back(clu);
}

void obtainContigs(const std::vector<SingleClippedCluster*>& clus,
                   std::vector<Contig> contigs) {
  for (int i = 0; i < clus.size(); ++i)
    contigs.push_back(clus[i]->contig());
}

bool findFirstRegion(std::vector<Contig>::iterator first,
                     std::vector<Contig>::iterator last,
                     const Contig& con,
                     double mismatchRate,
                     Region& region) {
  for (std::vector<Contig>::iterator itr = first; itr != last; ++itr) {
    if (con.overlaps(*itr, mismatchRate)) {
      region = Region(con.getAnchor(), (*itr).getAnchor());
      return true;
    }
  }
  return false;
}

void callDeletions(std::vector<Contig>& cons1,
                   std::vector<Contig>& cons2,
                   std::vector<Region>& calls,
                   double mismatchRate) {
  std::vector<Contig>::iterator last = cons2.end();
  for (std::vector<Contig>::reverse_iterator ritr = cons1.rbegin();
       ritr != cons1.rend();
       ++ritr) {
    std::vector<Contig>::iterator first = upper_bound(cons2.begin(), last, *ritr);
    if (first == last) return;
    Region reg;
    if (findFirstRegion(first, last, *ritr, mismatchRate, reg)) {
      calls.push_back(reg);
      if (first == cons2.begin()) return;
      last = first;
    }
  }
}

void selectCallsByLength(std::vector<Region>& calls,
                         int minLen,
                         int maxLen) {
  LengthSelector ls(minLen, maxLen);
  calls.erase(remove_if(calls.begin(), calls.end(), ls), calls.end());  
}

void outputCalls(std::string filename,
                 const std::vector<Region>& calls) {
  std::ofstream out(filename.c_str());
  out << "chromosome\ttype\tstart\tend" << std::endl;
  for (size_t i = calls.size() - 1; i >= 0; --i)
    out << calls[i].chrom() << "\t"
        << "deletion" << "\t"
        << calls[i].getStart() << "\t"
        << calls[i].getEnd() << std::endl;
  out.close();
}

void getClips(BamTools::BamReader& reader, std::vector<Clip*>& leftClips, std::vector<Clip*>& rightClips) {
  BamTools::BamAlignment al;
  while (reader.GetNextAlignment(al)) {
    std::vector<int> clipSizes, readPositions, genomePositions;
    if (!al.IsMapped()) {
      continue;
    }
    if (!al.GetSoftClips(clipSizes, readPositions, genomePositions)) {
      continue;
    }
    if (clipSizes.size() > 1) {
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

void tofile(std::string filename, const std::vector<Clip*>& clips) {
  std::ofstream out(filename.c_str());
  out << "length\tread_pos\tref_pos" << std::endl;
  for (size_t i = 0; i < clips.size(); ++i)
    out << clips[i]->getSize() << "\t" << clips[i]->getReadPosition() << "\t" << clips[i]->getPosition() << std::endl;
  out.close();
}

bool equals(const std::string s1, const std::string s2, int mismatches) {
  int cnt = 0;
  if (s1.size() != s2.size()) return false;
  for (int i = 0; i < s1.size(); ++i) {
    if (s1[i] != s2[i]) ++cnt;
    if (cnt > mismatches) return false;
  }
  return true;
}

bool isOverlapped(Clip* c1, Clip* c2, int mismatches) {
  std::string ls, rs;
  int len = c1->getSize() + c2->getSize();
  if (c1->getType() == Left) {
    ls = c1->getReadSeq().substr(0, len);
    rs = c2->getReadSeq().substr(c2->getReadPosition()-c1->getSize());
  } else {
    ls = c2->getReadSeq().substr(0, len);
    rs = c1->getReadSeq().substr(c1->getReadPosition()-c2->getSize());
  }
  if (equals(ls, rs, mismatches)) {                       // might consider mismatches between two strings
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

void extractClipsForDels(std::vector<Clip*>& inLCs, std::vector<Clip*>& inRCs, std::vector<Clip*>& outLCs, std::vector<Clip*>& outRCs, int min, int max) {
  std::set<Clip*> LCs, RCs;
  int index;
  for (size_t i=0; i<inLCs.size(); ++i) {
    if (!findMateIndex(inLCs[i], inRCs, index)) continue;
    for (size_t j = 0; j <= index; ++j) {
      int len = inLCs[i]->getSize() + inRCs[j]->getSize();
      if (len < min || len >= max) {
        continue;
      }
      if (isOverlapped(inLCs[i], inRCs[j], 1)) {
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
    bps.push_back(Breakpoint(LCs[i]->getPosition()+1, RCs[v]->getPosition()+1));
  }
}

// void buildBreakpoints(const std::vector<Clip*>& LCs, const std::vector<Clip*>& RCs, std::vector<Breakpoint>& bps) {
//   for (size_t i = 0; i < LCs.size(); ++i) {
//     for (size_t j = 0; j < RCs.size(); ++j)
//       if (isOverlapped(LCs[i], RCs[j]))
//         bps.push_back(Breakpoint(LCs[i]->getPosition()+1, RCs[j]->getPosition()+1));
//   }
// }

void groupBreakpoints(const std::vector<Breakpoint>& bps, std::vector<std::vector<Breakpoint> >& groups) {
  groups.push_back(std::vector<Breakpoint>(1, bps[0]));
  for (size_t i = 0; i < bps.size()-1; ++i) {
    if (bps[i].getY() > bps[i+1].getX())
      groups.back().push_back(bps[i+1]);
    else {
      groups.push_back(std::vector<Breakpoint>(1, bps[i+1]));
    }
  }
}

void makeCalls(const std::vector<std::vector<Breakpoint> >& groups, std::vector<StructVar>& calls, int minlen) {
  for (size_t i = 0; i < groups.size(); ++i) {
    if (groups[i][0].getY() - groups[i][0].getX() < minlen) continue;
    StructVar sv("22", groups[i][0].getX(), groups[i][0].getY()-1);
    calls.push_back(sv);
  }
}

// void clusterBreakpoints(const std::vector<Breakpoint>& bps, std::vector<std::vector<Breakpoint> >& clusters) {
//   std::vector<int> labels(bps.size(), -1);
//   int cnt = 0;
//   for (int i = 0; i < bps.size(); ++i) {
//     if (labels[i] != -1) continue;
//     labels[i] = cnt;
//     for (int j = i+1; j < bps.size(); ++j) {
//       if(labels[j] == -1 && bps[i].inSameCluster(bps[j]))
//         labels[j] = cnt;
//     }
//     cnt++;
//   }
//   for (int i = 0; i < cnt; ++i) {
//     std::vector<Breakpoint> clu;
//     clusters.push_back(clu);
//   }
//   for (int i = 0; i < bps.size(); ++i) {
//     clusters[labels[i]].push_back(bps[i]);
//   }
// }

void getTrueSvs(std::string filename, std::vector<StructVar>& trueSvs) {
  std::ifstream input(filename.c_str());  
  std::string line;
  getline(input, line);                 // skip the header line

  std::string chr, cls, seq;
  int begin, end;
  while (input >> chr >> begin >> end >> cls >> seq) {
    if (cls != "DEL") continue;
    StructVar sv(chr, begin, end-1);
    trueSvs.push_back(sv);
  }
}

std::set<StructVar> findOverlaps(StructVar t, const std::vector<StructVar>& cs1, const std::vector<StructVar>& cs2) {
  std::set<StructVar> s;
  for (size_t i = 0; i < cs1.size(); ++i) {
    if (cs1[i].right < t.left) continue;
    if (cs1[i].left > t.left) break;
    s.insert(cs1[i]);
  }
  for (size_t j = 0; j < cs2.size(); ++j) {
    if (cs2[j].left > t.right) continue;
    if (cs2[j].right < t.right) break;
    s.insert(cs2[j]);
  }
  return s;
}

bool cmp1(StructVar sv1, StructVar sv2) {
  if (sv1.left != sv2.left)
    return sv1.left < sv2.left;
  return sv1.right < sv2.right;
}

bool cmp2(StructVar sv1, StructVar sv2) {
  if (sv1.right != sv2.right)
    return sv1.right > sv2.right;
  return sv1.left > sv2.left;
}

void evaluateCalls(const std::vector<StructVar>& calls, const std::vector<StructVar>& trueSvs) {
  std::vector<std::set<StructVar> > result;
  std::vector<StructVar> cs1 = calls;
  std::vector<StructVar> cs2 = calls;
  sort(cs1.begin(), cs1.end(), cmp1);
  sort(cs2.begin(), cs2.end(), cmp2);  
  for (size_t i = 0; i < trueSvs.size(); ++i) {
    std::set<StructVar> s = findOverlaps(trueSvs[i], cs1, cs2);
    if (s.size() > 0) result.push_back(s);
  }
  std::cout << "#Identified calls: " << result.size() << std::endl;
}

void freeClips(std::vector<Clip*> cl) {
  while (!cl.empty()) delete cl.back(), cl.pop_back();
}
