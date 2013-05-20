#ifndef _DFINDER_H_
#define _DFINDER_H_

#include "SoftClip.h"
#include "Interval.h"
#include "TargetRegion.h"
#include "Consensus.h"
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
  void callTo(const std::string& filename);
  
 private: 
  void loadFrom(const std::string& filename);
  
  void identifyTargetRegions(int referenceId, std::vector<TargetRegion>& regions);

  void computeConsensuses(int referenceId, std::vector<Consensus>& consensuses1, std::vector<Consensus>& consensuses2);

  void callAllDeletions(const std::vector<TargetRegion>& regions,
                        const std::vector<Consensus>& consensuses1,
                        const std::vector<Consensus>& consensuses2,
                        std::vector<Deletion>& calls);
  void callOneDeletion(std::vector<Consensus>::const_iterator first1,
                       std::vector<Consensus>::const_iterator last1,
                       std::vector<Consensus>::const_iterator first2,
                       std::vector<Consensus>::const_iterator last2,
                       const TargetRegion& region,
                       std::vector<Deletion>& calls);
  
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

#endif /* _DFINDER_H_ */
