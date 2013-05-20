#include <algorithm>
#include <queue>
#include "DFinder.h"
#include "IntervalCluster.h"
#include "SoftClipCluster.h"

const int MinDeleltionLength = 50;

DFinder::DFinder(const std::string& filename, int meanInsertSize, int stdInsertSize, int minOverlapLength, double maxMismatchRate) :
    meanInsertSize(meanInsertSize), stdInsertSize(stdInsertSize), minOverlapLength(minOverlapLength), maxMismatchRate(maxMismatchRate) {
  loadFrom(filename);
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
  size = references.size();
  leftClips.resize(size);
  rightClips.resize(size);
  intervals.resize(size);
  
  int threshold = meanInsertSize + 4 * stdInsertSize;
  BamTools::BamAlignment ba1;
  int cnt = 0;
  while (r1.GetNextAlignment(ba1)) {
    // load a soft-clip
    std::vector<int> clipSizes, readPositions, genomePositions;
    if (ba1.GetSoftClips(clipSizes, readPositions, genomePositions) &&
        clipSizes.size() < 3) {
      if (clipSizes.size() == 2) {
        leftClips[ba1.RefID].push_back(new SoftClip(ba1.RefID,
                                                    genomePositions[0],
                                                    readPositions[0],
                                                    ba1.QueryBases.substr(0, ba1.Length - clipSizes[1]),
                                                    ba1.Qualities.substr(0, ba1.Length - clipSizes[1])));
        rightClips[ba1.RefID].push_back(new SoftClip(ba1.RefID,
                                                     genomePositions[1],
                                                     readPositions[1] - readPositions[0],
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
                                                    readPositions[0],
                                                    ba1.QueryBases,
                                                    ba1.Qualities));
      }
      continue;
    }
    // load an interval of paired-ends
    if (!ba1.IsPaired() ||
        !ba1.IsMapped() ||
        !ba1.IsMateMapped() ||
        ba1.RefID != ba1.MateRefID ||
        ba1.Position >= ba1.MatePosition ||
        ba1.InsertSize < threshold) continue;
        
    if (!ba1.IsReverseStrand() &&
        ba1.IsMateReverseStrand() &&
        r2.Jump(ba1.MateRefID, ba1.Position + ba1.InsertSize - 1)) {
      bool found = false;
      BamTools::BamAlignment ba2;
      while (r2.GetNextAlignment(ba2)) {
        if (ba2.Position > ba1.MatePosition) break;
        if (ba1.Name == ba2.Name &&
            ba1.IsFirstMate() != ba2.IsFirstMate()) {
          found = true;
          break;
        }
      }
      if (!found) continue;
      uint8_t t1, t2;
      if (!ba1.GetTag("XT", t1) || !ba2.GetTag("XT", t2)) continue;
      if ((t1 == 'U' && t2 == 'U')) {
        intervals[ba1.RefID].push_back(new Interval(cnt,
                                         ba1.RefID,
                                         ba1.GetEndPosition() + 1,
                                         ba1.MatePosition,
                                         ba1.InsertSize));
        cnt++;
      }
    }
  }

}

void DFinder::callTo(const std::string& filename) {
  std::ofstream out(filename.c_str());
  out << "chromosome\ttype\tstart1\tstart2\tend1\tend2" << std::endl;
  for (int i = 0; i < size; i++) {
    std::vector<TargetRegion> regions;
    identifyTargetRegions(i, regions);
    
    std::vector<Consensus> consensuses1;
    std::vector<Consensus> consensuses2;
    computeConsensuses(i, consensuses1, consensuses2);
    std::vector<Deletion> calls;
    callAllDeletions(regions, consensuses1, consensuses2, calls);

    for(auto itr = calls.begin(); itr != calls.end(); ++itr) {
      out << references[i].RefName << "\tdeletion\t"
          << (*itr).getStart1() << "\t" << (*itr).getStart2() << "\t"
          << (*itr).getEnd1() << "\t" << (*itr).getEnd2() << std::endl;
    }
  }
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
  std::vector<IntervalCluster> clusters;
  clusterIntervals(remainder, clusters);

  for (auto itr = clusters.begin(); itr != clusters.end(); ++itr) {
    regions.push_back((*itr).getTargetRegion(meanInsertSize, stdInsertSize));
  }
  
}

void clusterSoftClips(const std::vector<SoftClip*>& clips,
                     std::vector<SoftClipCluster>& clusters) {
  assert(clips.size() > 0);
  SoftClipCluster clu(clips[0]);
  for (int i = 0; i < clips.size() - 1; ++i) {
    clu.add(clips[i]);
    if (clu.clipPosition() < clips[i + 1]->position()) {
      clusters.push_back(clu);
      clu = SoftClipCluster(clips[i + 1]);
    }
  }
  clu.add(clips.back());
  clusters.push_back(clu);
}

void acquireConsensuses(const std::vector<SoftClip*>& clips,
                        std::vector<Consensus>& consensuses) {
  std::vector<SoftClipCluster> clusters;
  clusterSoftClips(clips, clusters);
  for (auto itr = clusters.begin(); itr != clusters.end(); ++itr) {
    consensuses.push_back((*itr).getConsensus());
  }
}

void DFinder::computeConsensuses(int referenceId, std::vector<Consensus>& consensuses1, std::vector<Consensus>& consensuses2) {
  acquireConsensuses(leftClips[referenceId], consensuses1);
  acquireConsensuses(rightClips[referenceId], consensuses2);
}

void DFinder::callAllDeletions(const std::vector<TargetRegion>& regions,
                               const std::vector<Consensus>& consensuses1,
                               const std::vector<Consensus>& consensuses2,
                               std::vector<Deletion>& calls) {
  for (auto itr = regions.begin(); itr != regions.end(); ++itr) {
    auto first1 = lower_bound(consensuses1.begin(), consensuses1.end(), (*itr).start, Consensus::compare);
    auto last1 = upper_bound(consensuses1.begin(), consensuses1.end(), (*itr).end, Consensus::compare2);
    auto first2 = lower_bound(consensuses2.begin(), consensuses2.end(), (*itr).start, Consensus::compare);
    auto last2 = upper_bound(consensuses2.begin(), consensuses2.end(), (*itr).end, Consensus::compare2);
    callOneDeletion(first1, last1, first2, last2, *itr, calls);
  }
}

void DFinder::callOneDeletion(std::vector<Consensus>::const_iterator first1,
                              std::vector<Consensus>::const_iterator last1,
                              std::vector<Consensus>::const_iterator first2,
                              std::vector<Consensus>::const_iterator last2,
                              const TargetRegion& region,
                              std::vector<Deletion>& calls) {
  for (auto itr1 = first2; itr1 != last2; itr1++) {
    first1 = upper_bound(first1, last1, *itr1);
    for (auto itr2 = first1; itr2 != last1; itr2++) {
      if ((*itr1).minDeletionLength(*itr2) < region.minDeletionLength) continue;
      if ((*itr1).maxDeletionLength(*itr2) > region.maxDeletionLength) break;
      Overlap overlap;
      if((*itr1).overlaps(*itr2, minOverlapLength, maxMismatchRate, overlap)) {
        calls.push_back(overlap.getDeletion());
        return;
      }
    }
  }
}
