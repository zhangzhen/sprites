#include "LeftClipped.h"
#include "RightClipped.h"
#include "ClippedCreator.h"
#include "RightClippedCluster.h"
#include "LeftClippedCluster.h"
#include <algorithm>
#include <queue>
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
  int numOfMismatches, editDist;
  // int n = 0;

  const BamTools::RefVector references = reader.GetReferenceData();
  while (reader.GetNextAlignment(al)) {
    std::vector<int> lens, cutpoints, anchors;
    if (!al.IsMapped() ||
        !al.GetSoftClips(lens, cutpoints, anchors) ||
        lens.size() > 2) {
      continue;
    }
    // if (al.GetTag("XM", numOfMismatches) &&
    //     al.GetTag("NM", editDist) &&
    //     numOfMismatches < editDist) continue;
    std::string refname = references[al.RefID].RefName;
    if (lens.size() == 2) {
      // std::cout << lens[0] << "\t" << cutpoints[0] << "\t" << anchors[0] << std::endl;
      // std::cout << lens[1] << "\t" << cutpoints[1] << "\t" << anchors[1] << std::endl;
      // std::cout << al.QueryBases << std::endl;
      Locus anc0(refname, anchors[0]);
      // std::cout << al.QueryBases.substr(0, cutpoints[1]) << std::endl;
      lefts.push_back(lCreator.createClipped(anc0, al.QueryBases.substr(0, al.Length - lens[1]), al.Qualities, 0, lens[0]));
      Locus anc1(refname, anchors[1]);
      // std::cout << std::string(lens[0], ' ') << al.QueryBases.substr(lens[0]) << std::endl;
      rights.push_back(rCreator.createClipped(anc1, al.QueryBases.substr(lens[0]), al.Qualities, al.Length - lens[0] - lens[1], lens[1]));
      continue;
    }
    Locus anchor(refname, anchors[0]);
    if (cutpoints[0] == lens[0]) { // left clipped
      lefts.push_back(lCreator.createClipped(anchor, al.QueryBases, al.Qualities, 0, lens[0]));
    } else { // right clipped
      rights.push_back(rCreator.createClipped(anchor, al.QueryBases, al.Qualities, al.Length - lens[0], lens[0]));
    }
  }

}

bool compSC(SingleClipped* sc1, SingleClipped* sc2) {
  return sc1->anchor() < sc2->anchor();
}

void clusterClippeds(std::vector<SingleClipped*>& clis,
                     std::vector<SingleClippedCluster*>& clus,
                     ClusterCreator& creator) {
  assert(clis.size() > 0);
  assert(is_sorted(clis.begin(), clis.end(), compSC));
  SingleClippedCluster *clu = creator.createCluster(clis[0]->anchor());
  for (int i = 0; i < clis.size() - 1; ++i) {
    clu->add(clis[i]);
    if (clu->getAnchor() < clis[i + 1]->anchor()) {
      clus.push_back(clu);
      clu = creator.createCluster(clis[i + 1]->anchor());
    }
  }
  clu->add(clis.back());
  clus.push_back(clu);
}

void loadControls(const std::string& filename,
                  std::vector<Region>& controls,
                  int minLen) {
  std::ifstream input(filename.c_str());  
  std::string line;
  getline(input, line);                 // skip the header line

  std::string chr, cls, seq;
  int start, end;
  while (input >> chr >> start >> end >> cls >> seq) {
    if (cls != "DEL" ||
        end - start < minLen) continue;
    controls.push_back(Region(Locus(chr, start), Locus(chr, end), 10));
  }
}

int minDistance(const std::vector<Region>& controls) {
  std::vector<Region>::const_iterator first = controls.begin();
  std::vector<Region>::const_iterator last = controls.end();
  assert(first + 1 < last);
  std::vector<Region>::const_iterator next = first;
  std::vector<int> distances;
  while (++next != last) {
    distances.push_back((*next).getStart() - (*first).getEnd());
    first = next;
  }
  return *min_element(distances.begin(), distances.end());
}

bool comp(SingleClippedCluster* clu1, SingleClippedCluster* clu2) {
  return clu1->getAnchor() < clu2->getAnchor();
}

void showBreakpointNeighborhood(SingleClippedCluster* clu,
                               std::vector<SingleClippedCluster*>& clus,
                               int range) {
  Locus l = clu->getAnchor();
  std::cout << l << std::endl;
  int cnt = 0;
  int center = clu->getAnchor().position();
  std::vector<SingleClippedCluster*>::iterator init, curr;
  init = lower_bound(clus.begin(), clus.end(), clu, comp);
  Contig c = (*init)->contig();
  if (c.getAnchor().position() > center && !c.getProximal()) curr = init - 1;
  else curr = init;
  while (abs(center - (*curr)->getAnchor().position()) <= range) {
    std::cout << cnt << ")" << std::endl;
    std::cout << **curr << std::endl;
    if (c.getProximal()) curr++;
    else curr--;
    cnt++;
  }
  
  curr = init - 1;
  Contig c2 = (*curr)->contig();
  while (c.getProximal() &&
         center - c2.getAnchor().position() < c2.sequence().size() - c2.getMarker()) {
    std::cout << cnt << ")" << std::endl;
    std::cout << **curr << std::endl;
    cnt++;
    curr--;
    c2 = (*curr)->contig();
  }
  
  if (c.getAnchor().position() > center) curr = init;
  else curr = init + 1;
  c2 = (*curr)->contig();
  while (!c.getProximal() &&
         c2.getAnchor().position() - l.position() < c2.getMarker()) {
    std::cout << cnt << ")" << std::endl;
    std::cout << **curr << std::endl;
    cnt++;
    curr++;
    c2 = (*curr)->contig();
  }
}

bool showSingleAnchorContext(SingleClippedCluster* clu,
                             std::vector<SingleClippedCluster*>& clus,
                             Contig& cont) {
  Locus l = clu->getAnchor();
  std::cout << l << std::endl;
  std::vector<SingleClippedCluster*>::iterator low, up;
  low = lower_bound(clus.begin(), clus.end(), clu, comp);
  up = upper_bound(clus.begin(), clus.end(), clu, comp);
  if (low != up) {
    std::cout << **low << std::endl;
    cont = (*low)->contig();
    return true;
  }
  Contig c = (*up)->contig();
  Contig c2 = (*(up-1))->contig();
  if ((!c2.getProximal() &&
       l.position() > c2.getAnchor().position() &&
       l.position() - c2.getAnchor().position() + c2.getMarker() + 1 <= c2.sequence().size()) ||
      (c2.getProximal() &&
       l.position() > c2.getAnchor().position() &&
       l.position() - c2.getAnchor().position() < c2.sequence().size() - c2.getMarker())) {
    std::cout << **(up-1) << std::endl;
    cont = c2;
    return true;
  }
  if ((c.getProximal() &&
       c.getAnchor().position() > l.position() &&
       c.getAnchor().position() - l.position() <= c.getMarker()) ||
      (!c.getProximal() &&
       c.getAnchor().position() > l.position() &&
       c.getAnchor().position() - l.position() < c.getMarker())) {
    std::cout << **up << std::endl;
    cont = c;
    return true;
  }
  return false;
}

void showControlContexts(const std::vector<Region>& controls,
                         std::vector<SingleClippedCluster*>& clus1,
                         std::vector<SingleClippedCluster*>& clus2) {
  assert(is_sorted(clus1.begin(), clus1.end(), comp));
  assert(is_sorted(clus2.begin(), clus2.end(), comp));
  StandardClusterCreator<RightClippedCluster> cluCreator1;
  StandardClusterCreator<LeftClippedCluster> cluCreator2;
  int cnt = 0;
  for (size_t i = 0; i < controls.size(); ++i) {
    std::cout << "> Control: " << i << "\tLength: " << controls[i].length() << std::endl;
    Locus start(controls[i].chrom(), controls[i].getStart());
    SingleClippedCluster* clu1 = cluCreator1.createCluster(start);
    Contig c1;
    bool b1 = showSingleAnchorContext(clu1, clus1, c1);
    // showBreakpointNeighborhood(clu1, clus1, 30);
    std::cout << "---------------------------" << std::endl;
    Locus end(controls[i].chrom(), controls[i].getEnd());
    SingleClippedCluster* clu2 = cluCreator2.createCluster(end);
    Contig c2;
    bool b2 = showSingleAnchorContext(clu2, clus2, c2);
    // showBreakpointNeighborhood(clu2, clus2, 30);
    if (b1 && b2) {
      int marker1 = c1.getMarker() + start.position() - c1.getAnchor().position();
      int marker2 = c2.getMarker() + end.position() - c2.getAnchor().position();
      std::cout << "******************************" << std::endl;
      int delta = abs(marker1-marker2);
      if (marker1 > marker2) {
        std::cout << std::string(marker1, ' ') << "+" << std::endl;
        std::cout << c1.sequence() << std::endl;
        std::cout << std::string(c1.getMarker(), ' ') << "^" << std::endl;
        std::cout << std::string(delta, ' ') << c2.sequence() << std::endl;
        std::cout << std::string(c2.getMarker() + delta, ' ') << "^" << std::endl;
      } else {
        std::cout << std::string(marker2, ' ') << "+" << std::endl;
        std::cout << std::string(delta, ' ') << c1.sequence() << std::endl;
        std::cout << std::string(c1.getMarker() + delta, ' ') << "^" << std::endl;
        std::cout << c2.sequence() << std::endl;
        std::cout << std::string(c2.getMarker(), ' ') << "^" << std::endl;
      }
    }
    if (!b1 || !b2) cnt++;
    std::cout << std::endl;
  }
  std::cout << cnt << std::endl;
}

void obtainContigs(const std::vector<SingleClippedCluster*>& clus,
                   std::vector<Contig>& contigs) {
  for (int i = 0; i < clus.size(); ++i)
    contigs.push_back(clus[i]->contig());
}

// bool findFirstRegion(std::vector<Contig>::iterator first,
//                      std::vector<Contig>::iterator last,
//                      const Contig& con,
//                      int minSupportSize,
//                      int minOverlapLen,
//                      double mismatchRate,
//                      Region& region) {
//   for (std::vector<Contig>::iterator itr = first; itr != last; ++itr) {
//     int offset = 0;
//     int overlapLen = con.overlaps(*itr, minSupportSize, minOverlapLen, mismatchRate, offset);
//     if (overlapLen > 0) {
//       Locus end = (*itr).getAnchor();
//       region = Region(con.getAnchor(), Locus(end.chrom(), end.position() + offset), overlapLen);
//       return true;
//     }
//   }
//   return false;
// }

// void callDeletions(std::vector<Contig>& cons1,
//                    std::vector<Contig>& cons2,
//                    std::vector<Region>& calls,
//                    int minSupportSize,
//                    int minOverlapLen,
//                    double mismatchRate) {
//   std::vector<Contig>::iterator last = cons2.end();
//   for (std::vector<Contig>::reverse_iterator ritr = cons1.rbegin();
//        ritr != cons1.rend();
//        ++ritr) {
//     std::vector<Contig>::iterator first = upper_bound(cons2.begin(), last, *ritr);
//     if (first == cons2.end()) continue;
//     Region reg;
//     if (findFirstRegion(first, last, *ritr, minSupportSize, minOverlapLen, mismatchRate, reg)) {
//       calls.push_back(reg);
//       last = first;
//     }
//   }
// }

// void callDeletion(std::vector<Contig>::iterator first1,
//                   std::vector<Contig>::iterator last1,
//                   std::vector<Contig>::iterator first2,
//                   std::vector<Contig>::iterator last2,
//                   std::vector<Region>& calls,
//                   int minSupportSize,
//                   int minOverlapLen,
//                   double mismatchRate) {
//   for (std::vector<Contig>::iterator it = first1; it != last1; it++) {
//     first2 = upper_bound(first2, last2, *it);
//     for (std::vector<Contig>::iterator it2 = first2; it2 != last2; it2++) {
//       int offset = 0;
//       int overlapLen = (*it).overlaps(*it2, minSupportSize, minOverlapLen, mismatchRate, offset);
//       unsigned s = (*it).getAnchor().position();
//       unsigned t = (*it2).getAnchor().position() + offset;
//       if (overlapLen > 0 && t - s >= MinDelLen) {
//         Locus end = (*it2).getAnchor();
//         calls.push_back(Region((*it).getAnchor(), Locus(end.chrom(), end.position() + offset), overlapLen));
//         return;
//       }
//     }
//   }
// }

// void callDeletions2(const std::vector<Region2>& in,
//                     std::vector<Contig>& cons1,
//                     std::vector<Contig>& cons2,
//                     std::vector<Region>& calls,
//                     int minSupportSize,
//                     int minOverlapLen,
//                     double mismatchRate) {
//   for (size_t i = 0; i < in.size(); i++) {
//     std::vector<Contig>::iterator first1 = lower_bound(cons1.begin(), cons1.end(), in[i].low, Contig::compare);
//     std::vector<Contig>::iterator last1 = upper_bound(cons1.begin(), cons1.end(), in[i].high, Contig::compare2);
//     std::vector<Contig>::iterator first2 = lower_bound(cons2.begin(), cons2.end(), in[i].low, Contig::compare);
//     std::vector<Contig>::iterator last2 = upper_bound(cons2.begin(), cons2.end(), in[i].high, Contig::compare2);
//     callDeletion(first1, last1, first2, last2, calls, minSupportSize, minOverlapLen, mismatchRate);
//   }
// }

void outputCalls(std::string filename,
                 const std::vector<Region>& calls) {
  std::ofstream out(filename.c_str());
  out << "chromosome\ttype\tstart\tend\toverlap length" << std::endl;
  for (std::vector<Region>::const_iterator it = calls.begin();
       it != calls.end();
       it++)
    out << (*it).chrom()
        << "\t"
        << "deletion"
        <<"\t"
        << (*it).getStart()
        << "\t"
        << (*it).getEnd()
        << "\t"
        << (*it).overlapLength()
        << std::endl;
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
