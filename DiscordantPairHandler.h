#include <string>
#include "IntervalCluster.h"
#include "api/BamReader.h"

class DiscordantPairHandler {
 public:
  static void identifyFocalRegions(const std::string& filename,
                                   std::vector<Region2>& regions,
                                   int mu,
                                   int sigma);

 private:
  static void loadIntervalsFromBam(BamTools::BamReader& r1,
                                   BamTools::BamReader& r2,
                                   std::vector<Interval*>& intervals,
                                   int threshold);
  static void removeNestingIntervals(const std::vector<Interval*>& in,
                                     std::vector<Interval*>& out);
  static void clusterIntervals(const std::vector<Interval*>& in,
                               std::vector<IntervalCluster>& clus,
                               int threshold);

};

