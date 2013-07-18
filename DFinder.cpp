#include <algorithm>
#include <queue>
#include "DFinder.h"
#include "IntervalCluster.h"

DFinder::DFinder(const std::string& filename, int meanInsertSize, int stdInsertSize, int minOverlapLength, double maxMismatchRate, double discordant) :
    meanInsertSize(meanInsertSize), stdInsertSize(stdInsertSize), minOverlapLength(minOverlapLength), maxMismatchRate(maxMismatchRate), discordant(discordant) {
  loadFrom(filename);
  // assert(leftClips[19].size() > 0);
  // assert(rightClips[19].size() > 0);
  for (int i = 0; i < size; ++i) {
    sort(leftClips[i].begin(), leftClips[i].end(), SoftClip::compareL);
    sort(rightClips[i].begin(), rightClips[i].end(), SoftClip::compareR);
  }
}

DFinder::~DFinder() {
  for (int i = 0; i < size; ++i) {
    for_each(leftClips[i].begin(), leftClips[i].end(), DeletePtr<SoftClip>());
    for_each(rightClips[i].begin(), rightClips[i].end(), DeletePtr<SoftClip>());
    for_each(intervals[i].begin(), intervals[i].end(), DeletePtr<Interval>());
  }
}

void DFinder::loadFrom(const std::string& filename) {
  BamTools::BamReader r1, r2;
  r1.Open(filename);
  r2.Open(filename);
  r1.LocateIndex();
  r2.LocateIndex();

  references = r1.GetReferenceData();
  size = r1.GetReferenceCount();
  // for (auto itr = references.begin(); itr != references.end(); ++itr)
  //   std::cout << (*itr).RefName << "\t" << (*itr).RefLength << std::endl;
  leftClips.resize(size);
  rightClips.resize(size);
  intervals.resize(size);

  int threshold = meanInsertSize + round(discordant * stdInsertSize);
  BamTools::BamAlignment ba1;
  int cnt = 0;
  int n_pairs = 0, found_pairs = 0;
  while (r1.GetNextAlignment(ba1)) {
    // load a soft-clip
    std::vector<int> clipSizes, readPositions, genomePositions;
    // int xt;
    // if (!ba1.GetTag("XT", xt)) xt = 'U';
    
    if (ba1.IsDuplicate()) continue;
    if (ba1.GetSoftClips(clipSizes, readPositions, genomePositions) &&
        // (xt == 'M' || xt =='U') &&
        clipSizes.size() <= 2) {
      // if (ba1.Position == 15395140) {
      //   std::cout << ba1.QueryBases << std::endl;
      //   std::cout << ba1.AlignedBases << std::endl;
      //   std::cout << clipSizes[0] << "\t" << readPositions[0] << "\t"
      //             << genomePositions[0] << "\t" << ba1.GetEndPosition() << std::endl;
      // }
      if (clipSizes.size() == 2) {
        leftClips[ba1.RefID].push_back(new SoftClip(ba1.RefID,
                                                    genomePositions[0],
                                                    readPositions[0],
                                                    ba1.QueryBases.substr(0, ba1.Length - clipSizes[1]),
                                                    ba1.Qualities.substr(0, ba1.Length - clipSizes[1])));
        rightClips[ba1.RefID].push_back(new SoftClip(ba1.RefID,
                                                     genomePositions[1],
                                                     ba1.Length - clipSizes[0] - clipSizes[1],
                                                     ba1.QueryBases.substr(clipSizes[0]),
                                                     ba1.Qualities.substr(clipSizes[0])));      
      } else if (readPositions[0] == clipSizes[0]) { // left clip
        leftClips[ba1.RefID].push_back(new SoftClip(ba1.RefID,
                                                    genomePositions[0],
                                                    readPositions[0],
                                                    ba1.QueryBases,
                                                    ba1.Qualities));
      } else {
        rightClips[ba1.RefID].push_back(new SoftClip(ba1.RefID,
                                                    genomePositions[0],
                                                    ba1.Length - clipSizes[0],
                                                    ba1.QueryBases,
                                                    ba1.Qualities));
      }
      // continue;
    }
    // load an interval of paired-ends
    if (!ba1.IsPaired() ||
        // ba1.IsProperPair() ||
        !ba1.IsMapped() ||
        !ba1.IsMateMapped() ||
        ba1.RefID != ba1.MateRefID ||
        // ba1.InsertSize < meanInsertSize ||
        ba1.InsertSize < threshold ||
        ba1.Position >= ba1.MatePosition) continue;
        
    if (!ba1.IsReverseStrand() &&
        // r2.Jump(ba1.RefID, ba1.MatePosition) &&
        ba1.IsMateReverseStrand()) {
      // bool found = false;
      // BamTools::BamAlignment ba2;
      // n_pairs++;
      // while (r2.GetNextAlignment(ba2)) {
      //   if (ba2.Position > ba1.MatePosition) break;
      //   if (ba1.Name == ba2.Name &&
      //       ba1.IsFirstMate() != ba2.IsFirstMate()) {
      //     found = true;
      //     break;
      //   }
      // }
      // if (!found) continue;
      // found_pairs++;
      // int xt2;
      // if (!ba2.GetTag("XT", xt2)) xt2 = 'U';
      // if ((xt == 'U' || xt == 'M') &&
      //     (xt2 == 'U' || xt2 == 'M')) {
      intervals[ba1.RefID].push_back(new Interval(cnt,
                                                  ba1.RefID,
                                                  ba1.GetEndPosition(),
                                                  ba1.MatePosition,
                                                  ba1.InsertSize));
      cnt++;
      // }
    }
  }
  // std::cout << n_pairs << "\t" << found_pairs << std::endl;
}

void DFinder::callToFile(const std::string& filename) {
  std::vector<Deletion> calls;
  call(filename, calls);

  std::ofstream out(filename.c_str());
  out << "Chromosome\tType\tStart\tEnd\tLength" << std::endl;
  for(auto itr = calls.begin(); itr != calls.end(); ++itr) {
    out << references[(*itr).getReferenceId()].RefName << "\tDEL\t" << (*itr).getStart2() << "\t" << (*itr).getEnd2() << "\t" << (*itr).length()
        << std::endl;
  }
}

void DFinder::callToVcf(const std::string& filename) {
  std::vector<Deletion> calls;
  call(filename, calls);
  
  std::ofstream out(filename.c_str());
  out << "##fileformat=VCFv4.1" << std::endl;
  out << "##reference=human_g1k_v37" << std::endl;
  out << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">" << std::endl;
  out << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << std::endl;
  out << "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">" << std::endl;
  out << "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">" << std::endl;
  out << "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << std::endl;
  out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << std::endl;
  out << "##INFO=<ID=SAMPLES,Number=.,Type=String,Description=\"List of samples\">" << std::endl;
  out << "##ALT=<ID=DEL,Description=\"Deletion\">" << std::endl;
  out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;

  for(auto itr = calls.begin(); itr != calls.end(); ++itr) {
    out << references[(*itr).getReferenceId()].RefName << "\t" << (*itr).getStart2() << "\t.\tN\t<DEL>\t.\tPASS\t"
        << "IMPRECISE;SVTYPE=DEL;END=" << (*itr).getEnd2() << ";SVLEN=-" << (*itr).length()
        << std::endl;
  }
}

void DFinder::call(const std::string& filename, std::vector<Deletion>& calls) {
  std::vector<Deletion> res;
  for (int i = 0; i < size; i++) {
    std::vector<TargetRegion> regions;
    // for (auto itr = intervals[i].begin(); itr != intervals[i].end(); ++itr) {
    //   regions.push_back({(*itr)->getStartPos(), (*itr)->getEndPos(),
    //           (*itr)->minDeletionLength(meanInsertSize, stdInsertSize),
    //           (*itr)->maxDeletionLength(meanInsertSize, stdInsertSize)});
    // }
    identifyTargetRegions(i, regions);
    // for (auto itr = regions.begin(); itr != regions.end(); ++itr)
    //   std::cout << references[i].RefName << "\t"
    //             << (*itr).start << "\t"
    //             << (*itr).end << "\t"
    //             << (*itr).minDeletionLength << "\t"
    //             << (*itr).maxDeletionLength << "\t"
    //             << std::endl;

    // for (auto itr1 = lRegions[i].begin(); itr1 != lRegions[i].end(); ++itr1)
    //   std::cout << references[i].RefName << "\t"
    //             << (*itr1).start << "\t"
    //             << (*itr1).end << "\t"
    //             << (*itr1).minDeletionLength << "\t"
    //             << (*itr1).maxDeletionLength << "\t"
    //             << std::endl;

    callAllDeletions(regions, leftClips[i], rightClips[i], SoftClip::compare1, SoftClip::compare2, res);
  }
  mergeCalls(res, calls);
}

void DFinder::mergeCalls(std::vector<Deletion>& in, std::vector<Deletion>& out) {
  auto first = in.begin();
  auto last = in.end();
  sort(first, last);
  if (first == last) return;
  auto next = first;
  Deletion d = *first;
  while (++next != last) {
    if ((*first).overlaps(*next)) {
      if (d.length() < (*next).length()) d = *next;
    } else {
      out.push_back(d);
      d = *next;
    }
    ++first;
  }
  out.push_back(d);
}

void removeSuperIntervals(const std::vector<Interval*>& input, std::vector<const Interval*>& remainder) {
  std::vector<EndPoint> endpoints;
  for (auto itr = input.begin(); itr != input.end(); ++itr) {
    endpoints.push_back((*itr)->getStart());
    endpoints.push_back((*itr)->getEnd());
  }
  sort(endpoints.begin(), endpoints.end());
  int prev_id = -1;
  for (auto itr = endpoints.begin(); itr != endpoints.end(); ++itr) {
    int curr_id = (*itr).ownerId();
    if (!(*itr).isStart() && prev_id < curr_id) {
      remainder.push_back((*itr).getOwner());
      prev_id = curr_id;
    }
  }
}

void clusterIntervals(const std::vector<const Interval*>& remainder, std::vector<IntervalCluster>& clusters) {
  std::vector<EndPoint> endpoints;
  for (auto itr = remainder.begin(); itr != remainder.end(); ++itr) {
    endpoints.push_back((*itr)->getStart());
    endpoints.push_back((*itr)->getEnd());
  }
  sort(endpoints.begin(), endpoints.end());
  std::queue<const Interval*> q;
  for (auto itr = endpoints.begin(); itr != endpoints.end(); ++itr) {
    if ((*itr).isStart())
      q.push((*itr).getOwner());
    else {
      if (q.empty()) continue;
      IntervalCluster clu;
      while (!q.empty()) {
        clu.add(q.front());
        q.pop();
      }
      clusters.push_back(clu);
    }
  }
}

void DFinder::identifyTargetRegions(int referenceId, std::vector<TargetRegion>& regions) {
  // gets rid of intervals that contain sub-intervals
  std::vector<const Interval*> remainder;
  removeSuperIntervals(intervals[referenceId], remainder);
  // clusters pairwise overlapping intervals 
  // std::vector<IntervalCluster> clusters;
  // clusterIntervals(remainder, clusters);

  for (auto itr = remainder.begin(); itr != remainder.end(); ++itr) {
    regions.push_back({(*itr)->getStartPos(), (*itr)->getEndPos(),
            (*itr)->minDeletionLength(meanInsertSize, stdInsertSize),
            (*itr)->maxDeletionLength(meanInsertSize, stdInsertSize)});
  }

  // for (auto itr = regions.begin(); itr != regions.end(); ++itr) {
  //   std::cout << (*itr).start << "\t" << (*itr).end << "\t" << (*itr).minDeletionLength << "\t" << (*itr).maxDeletionLength << std::endl;
  // }
}

void DFinder::printOverlaps(const std::string& filename, int readlength) {
  std::ifstream in(filename.c_str());
  
  int start, end, length;
  std::string refname;
  int cnt = 0, matched = 0;

  while(in >> refname >> start >> end >> length) {
  
    bool found = false;
    int id;
    for (int i = 0; i < size; ++i) {
      if (references[i].RefName == refname) {
        id = i;
        found = true;
        break;
      }
    }

    std::cout << "#######################################" << std::endl;
    std::cout << refname << ":" << start << "-" << end << "\t" << length << std::endl;
    cnt++;

    if (!found) continue;
    auto first1 = upper_bound(leftClips[id].begin(), leftClips[id].end(), end - readlength, SoftClip::compare2) - 1;
    auto last1 = upper_bound(leftClips[id].begin(), leftClips[id].end(), end, SoftClip::compare2) - 1;
    auto first2 = lower_bound(rightClips[id].begin(), rightClips[id].end(), start, SoftClip::compare1);
    auto last2 =  lower_bound(rightClips[id].begin(), rightClips[id].end(), start + readlength, SoftClip::compare1);
    Overlap overlap;
    if (overlaps(first1, last1, first2, last2, (length < 150) ? std::max(0, length - stdInsertSize) : std::max(0, length - 3 * stdInsertSize), (length < 150) ? length + stdInsertSize : length + 3 * stdInsertSize, overlap)) {
      std::cout << overlap << std::endl;
      matched++;
    }
  }
  std::cout << "#Matched: " << matched << std::endl;
  std::cout << "Matching rate: " << float(matched) / cnt << std::endl;
}
