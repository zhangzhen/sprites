#ifndef _DFINDER_H_
#define _DFINDER_H_

#include "SoftClip.h"
#include "Interval.h"
#include "TargetRegion.h"
#include "api/BamReader.h"


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
  DFinder(const std::string& filename, int meanInsertSize, int stdInsertSize, int minOverlapLength, double maxMismatchRate);
  ~DFinder();
  void callToFile(const std::string& filename);
  void callToVcf(const std::string& filename);
  void printOverlaps(const std::string& filename, int readlength);
  
 private: 
  void call(const std::string& filename, std::vector<Deletion>& calls);
  
  void loadFrom(const std::string& filename);
  
  void identifyTargetRegions(int referenceId, std::vector<TargetRegion>& regions);

  // void computeConsensuses(int referenceId, std::vector<Consensus>& consensuses1, std::vector<Consensus>& consensuses2);

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
                int minDeletionLength,
                int maxDeletionLength,
                Overlap& overlap);
  
  int meanInsertSize;
  int stdInsertSize;
  int minOverlapLength;
  double maxMismatchRate;
  BamTools::RefVector references;
  int size;
  
  std::vector<std::vector<SoftClip*> > leftClips;
  std::vector<std::vector<SoftClip*> > rightClips;
  std::vector<std::vector<Interval*> > intervals;
};

template <typename T, typename Compare1, typename Compare2>
void DFinder::callAllDeletions(const std::vector<TargetRegion>& regions,
                               const T& consensuses1,
                               const T& consensuses2,
                               Compare1 comp1,
                               Compare2 comp2,
                               std::vector<Deletion>& calls) {
  for (auto itr = regions.begin(); itr != regions.end(); ++itr) {
    auto first1 = upper_bound(consensuses1.begin(), consensuses1.end(), (*itr).start, comp2);
    auto last1 = upper_bound(consensuses1.begin(), consensuses1.end(), (*itr).end, comp2);
    auto first2 = lower_bound(consensuses2.begin(), consensuses2.end(), (*itr).start, comp1);
    auto last2 = lower_bound(consensuses2.begin(), consensuses2.end(), (*itr).end, comp1);
    // if ((*itr).start == 5942980) {
    //   for (auto itr2 = first1; itr2 != last1 + 1; ++itr2)
    //     std::cout << **itr2 << std::endl;
    //   std::cout << "2222222222222" << std::endl;
    //   for (auto itr2 = first2 - 1; itr2 != last2; ++itr2) {
    //     std::cout << **itr2 << std::endl;
    //   }
    // }
    callOneDeletion(first1, last1 + 1, first2 - 1, last2, *itr, calls);
  }
}

template <typename ForwardIterator>
void DFinder::callOneDeletion(ForwardIterator first1,
                              ForwardIterator last1,
                              ForwardIterator first2,
                              ForwardIterator last2,
                              const TargetRegion& region,
                              std::vector<Deletion>& calls) {
  for (auto itr1 = first2; itr1 != last2; ++itr1) {
    for (auto itr2 = last1; itr2 != first1 - 1; --itr2) {
      // if (region.start == 5942980) {
      //   std::cout << (*itr1)->position() << '-' << (*itr2)->position() << std::endl;
      //   std::cout << ((*itr1)->minDeletionLength(**itr2) < region.minDeletionLength) << '\t' << ((*itr1)->maxDeletionLength(**itr2) > region.maxDeletionLength) << std::endl;
      // }
      if ((*itr1)->maxDeletionLength(**itr2) < region.minDeletionLength) break;
      if ((*itr1)->minDeletionLength(**itr2) > region.maxDeletionLength) continue;
      Overlap overlap;
      if((*itr1)->overlaps(**itr2, minOverlapLength, maxMismatchRate, overlap) &&
         overlap.deletionLength() >= region.minDeletionLength &&
         overlap.deletionLength() <= region.maxDeletionLength) {
        calls.push_back(overlap.getDeletion());
        std::cout << overlap << std::endl;
        return;
      }
    }
  }
}

template <typename ForwardIterator>
bool DFinder::overlaps(ForwardIterator first1,
                       ForwardIterator last1,
                       ForwardIterator first2,
                       ForwardIterator last2,
                       int minDeletionLength,
                       int maxDeletionLength,
                       Overlap& overlap) {
  // bool found = false;
  // std::vector<Overlap> ovs;
  for (auto itr1 = first2; itr1 != last2; ++itr1) {
    for (auto itr2 = last1; itr2 != first1 - 1; --itr2) {
      // if ((*itr1)->maxDeletionLength(**itr2) < minDeletionLength) {
      //   if (found) return Overlap::getBestOverlap(ovs, overlap);
      //   return false;
      // }
      // if ((*itr1)->minDeletionLength(**itr2) > maxDeletionLength) continue;
      Overlap ov;
      if((*itr1)->overlaps(**itr2, minOverlapLength, maxMismatchRate, ov) &&
         ov.deletionLength() >= minDeletionLength &&
         ov.deletionLength() <= maxDeletionLength) {
        overlap = ov;
        return true;
        // ovs.push_back(ov);
        // found = true;
        // break;
      }
    }
  }
  // return Overlap::getBestOverlap(ovs, overlap);
  return false;
}


#endif /* _DFINDER_H_ */
