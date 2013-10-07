#ifndef _DFINDER_H_
#define _DFINDER_H_

#include <iostream>
#include <iterator>
#include <map>
#include <utility>
#include "SoftClip.h"
#include "ChrRegion.h"
#include "ChrRegionCluster.h"
#include "TargetRegion.h"
#include "api/BamReader.h"
#include "DFinderHelper.h"


template<typename T>
class DeletePtr {
 public:
  void operator() (T* ptr) {
    delete ptr;
  }
};

class DFinder
{
 public:
  DFinder(const std::string& filename, int meanInsertSize, int stdInsertSize, int minOverlapLength, double maxMismatchRate, double discordant);
  ~DFinder();
  void callToFile(const std::string& filename);
  void callToVcf(const std::string& filename);
  void printOverlaps(const std::string& filename, int readlength);
  void checkAgainstGoldStandard(const std::string& filename);

 private:
  void call(const std::string& filename, std::vector<Deletion>& calls);

  void loadFrom();

  bool isLargeInsertSize(int insertSize);

  bool getMateOf(const BamTools::BamAlignment& it, BamTools::BamAlignment& itsMate);

  void identifyTargetRegions(int referenceId, std::vector<TargetRegion>& regions);

  static void mergeCalls(std::vector<Deletion>& in, std::vector<Deletion>& out);

  // void computeConsensuses(int referenceId, std::vector<Consensus>& consensuses1, std::vector<Consensus>& consensuses2);
  int numOfClipsIn(const TargetRegion& region, const std::vector<SoftClip*>& clips);

  template <typename T, typename Compare1, typename Compare2>
  void callAllDeletions(const std::vector<TargetRegion>& regions,
                        const T& consensuses1,
                        const T& consensuses2,
                        Compare1 comp1,
                        Compare2 comp2,
                        std::vector<Deletion>& calls);

  template <typename ForwardIterator>
  void callOneDeletion(ForwardIterator first1,
                       ForwardIterator last1,
                       ForwardIterator first2,
                       ForwardIterator last2,
                       const TargetRegion& region,
                       std::vector<Deletion>& calls);

  template <typename ForwardIterator>
  bool overlaps(ForwardIterator first1,
                ForwardIterator last1,
                ForwardIterator first2,
                ForwardIterator last2,
		const TargetRegion& region,
                Overlap& overlap);

  bool findReferenceId(const std::string& name, int& id);

  int meanInsertSize;
  int stdInsertSize;
  int minOverlapLength;
  double maxMismatchRate;
  double discordant;
  BamTools::RefVector references;
  int size;
  static const int lengthThreshold = 50;
  BamTools::BamReader r1, r2;

  std::vector<std::vector<SoftClip*> > leftClips;
  std::vector<std::vector<SoftClip*> > rightClips;
  std::vector<std::vector<ChrRegion*> > intervals;

  struct MyInterval {
    std::string refname;
    int start;
    int end;
    int length;
  };

  bool loadMyIntervals(const std::string& filename, std::vector<MyInterval>& out);
  bool checkMyInterval(const MyInterval& myInterval, int refId, const std::vector<ChrRegionCluster>& regions);
};

template <typename T, typename Compare1, typename Compare2>
void DFinder::callAllDeletions(const std::vector<TargetRegion>& regions,
                               const T& consensuses1,
                               const T& consensuses2,
                               Compare1 comp1,
                               Compare2 comp2,
                               std::vector<Deletion>& calls) {
  for (auto itr = regions.begin(); itr != regions.end(); ++itr) {
      /* std::cout << ">>>>>>>>>>>>>>>>>>>>>>>" << std::endl; */
      /* std::cout << *itr << std::endl; */
    auto first1 = lower_bound(consensuses1.begin(), consensuses1.end(), (*itr).start, comp1);
    auto last1 = upper_bound(consensuses1.begin(), consensuses1.end(), (*itr).end, comp2);
    if (last1 == consensuses1.end()) continue;
    auto first2 = lower_bound(consensuses2.begin(), consensuses2.end(), (*itr).start, comp1);
    if (first2 == consensuses2.end()) continue;
    auto last2 = upper_bound(consensuses2.begin(), consensuses2.end(), (*itr).end, comp2);
    /* if ((*itr).start == 30255871) { */
    /*   for (auto itr2 = first1; itr2 != last1 + 1; ++itr2) */
    /*     std::cout << **itr2 << std::endl; */
    /*   std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl; */
    /*   for (auto itr2 = first2 - 1; itr2 != last2; ++itr2) { */
    /*     std::cout << **itr2 << std::endl; */
    /*   } */
    /* } */
    callOneDeletion(first1, last1, first2, last2, *itr, calls);
  }
}

template <typename ForwardIterator>
void DFinder::callOneDeletion(ForwardIterator first1,
                              ForwardIterator last1,
                              ForwardIterator first2,
                              ForwardIterator last2,
                              const TargetRegion& region,
                              std::vector<Deletion>& calls) {
  Overlap ov;
  if (overlaps(first1, last1, first2, last2, region, ov))
    calls.push_back(ov.getDeletion());
}

template <typename ForwardIterator>
bool DFinder::overlaps(ForwardIterator first1,
                       ForwardIterator last1,
                       ForwardIterator first2,
                       ForwardIterator last2,
		       const TargetRegion& region,
                       Overlap& overlap) {
    std::map<std::pair<int,int>,std::vector<Overlap> > overlaps;
  /* std::vector<Overlap> ovs; */
  for (auto itr1 = first2; itr1 != last2; ++itr1) {
    for (auto itr2 = last1 - 1; itr2 != first1 - 1; --itr2) {
      if ((*itr1)->maxDeletionLength(**itr2) < std::max(lengthThreshold, region.minDeletionLength)) break;
      if ((*itr1)->minDeletionLength(**itr2) > region.maxDeletionLength) continue;
      Overlap ov;
      if((*itr1)->overlaps(**itr2, minOverlapLength, maxMismatchRate, ov) &&
         ov.deletionLength() >= std::max(lengthThreshold, region.minDeletionLength) &&
         ov.deletionLength() <= region.maxDeletionLength &&
	 ov.score() < maxMismatchRate) {
	/*   if (region.start == 16969499) { */
	/*   std::cout << ">>>>>>>>>>>>>>>>>>>>>>" << std::endl; */
	/*   std::cout << ov << std::endl; */
	/* } */
	  overlaps[std::make_pair(ov.start(), ov.end())].push_back(ov);
        /* ovs.push_back(ov); */
      }
    }
  }
  if (overlaps.empty()) { return false; }
  if (overlaps.size() == 1) {
      return Overlap::getHighScoreOverlap(overlaps.begin()->second, overlap);
  }
  auto itr = overlaps.begin();
  auto res = itr++;
  for (; itr != overlaps.end(); ++itr) {
      if (!contains(res->first, itr->first)) {
	  return Overlap::getHighScoreOverlap(res->second, overlap);
      }
      if ((res->first.second - res->first.first) < (itr->first.second - itr->first.first))
	  res = itr;
      /* std::cout << "(" << itr->first.first << ", " << itr->first.second << ")" << std::endl; */
  }
  /* copy(ovs.begin(), ovs.end(), std::ostream_iterator<Overlap>(std::cout, "\n")); */
  return Overlap::getHighScoreOverlap(res->second, overlap);
  /* return Overlap::getHighScoreOverlap(ovs, overlap); */
}


#endif /* _DFINDER_H_ */
