#ifndef _DFINDER_H_
#define _DFINDER_H_

#include <iostream>
#include <iterator>
#include <map>
#include <utility>
#include "SoftClip.h"
#include "Interval.h"
#include "ChrRegion.h"
#include "ChrRegionCluster.h"
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
  DFinder(const std::string& filename, int meanInsertSize, int stdInsertSize, int minOverlapLength, double maxMismatchRate, int numOfStd);
  ~DFinder();
  void callToFile(const std::string& filename);
  void callToVcf(const std::string& filename);
  void printOverlaps(const std::string& filename, int readlength);
  void checkAgainstGoldStandard(const std::string& filename);

 private:
  void call(const std::string& filename, std::vector<Deletion>& calls);

  void loadFrom();

  bool isLargeInsertSize(int insertSize);
  bool isProperInsertSize(int insertSize);

  /* bool getMateOf(const BamTools::BamAlignment& it, BamTools::BamAlignment& itsMate); */

  /* void identifyTargetRegions(int referenceId, std::vector<TargetRegion>& regions); */

  /* static void mergeCalls(std::vector<Deletion>& in, std::vector<Deletion>& out); */

  void removeLargeChrRegions(std::vector<ChrRegion*>& regions);

  void clusterChrRegions(const std::vector<ChrRegion*>& remainder, std::vector<ChrRegionCluster>& clusters);

  bool callDeletionInCluster(const ChrRegionCluster& cluster, Deletion& deletion);

  bool getOverlapInRegion(const ChrRegion& region, Overlap& overlap);

  /* int numOfClipsIn(const TargetRegion& region, const std::vector<SoftClip*>& clips); */

  bool findReferenceId(const std::string& name, int& id);
  Interval getIntervalOfLeftClips(const ChrRegion& regionOfInterest);
  Interval getIntervalOfRightClips(const ChrRegion& regionOfInterest);

  static void getSoftClipsIn(const Interval& interval, const std::vector<SoftClip*>& input, std::vector<SoftClip*>& output);

  int meanInsertSize;
  int stdInsertSize;
  int minOverlapLength;
  double maxMismatchRate;
  int numOfStd;
  BamTools::RefVector references;
  int size;
  static const int lengthThreshold = 50;
  static const int MapQualityThreshold = 10;
  BamTools::BamReader r1;

  std::vector<std::vector<SoftClip*> > leftClips;
  std::vector<std::vector<SoftClip*> > leftParts;

  std::vector<std::vector<SoftClip*> > rightClips;
  std::vector<std::vector<SoftClip*> > rightParts;

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

#endif /* _DFINDER_H_ */
