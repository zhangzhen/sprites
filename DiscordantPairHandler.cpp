#include "DiscordantPairHandler.h"
#include <algorithm>
#include <queue>

void DiscordantPairHandler::identifyFocalRegions(const std::string& filename,
                                                 std::vector<Region2>& regions,
                                                 int mu,
                                                 int sigma) {
  BamTools::BamReader r1, r2;
  r1.Open(filename);
  r2.Open(filename);
  r1.LocateIndex();
  r2.LocateIndex();

  std::vector<Interval*> in;
  int th1 = mu + 4 * sigma;
  loadIntervalsFromBam(r1, r2, in, th1);
  // for (size_t i = 0; i < in.size(); i++)
  //   std::cout << i[i]->getId() << "\t"
  //             << in[i]->getStartPos() << "\t"
  //             << in[i]->getEndPos() << std::endl;
  // return 0;
  std::vector<Interval*> out;
  removeNestingIntervals(in, out);
  // for (size_t i = 0; i < out.size(); i++)
  //   std::cout << out[i]->getStartPos() << "\t"
  //             << out[i]->getEndPos() << "\t"
  //             << out[i]->length() << "\t"
  //             << std::endl;
  std::vector<IntervalCluster> clus;
  int th2 = 6 * sigma;
  clusterIntervals(out, clus, th2);
  // for (size_t i = 0; i < clus.size(); i++) {
  //   std::cout << clus[i] << std::endl;
  // }
  // return 0;

  for (size_t i = 0; i < clus.size(); i++) {
    Region2 r2 = clus[i].focalRegion();
    regions.push_back(r2);
    // std::cout << r2.low << "\t"
    //           << r2.high << std::endl;
  }
  r1.Close();
  r2.Close();
}

void DiscordantPairHandler::loadIntervalsFromBam(BamTools::BamReader& r1,
                                                 BamTools::BamReader& r2,
                                                 std::vector<Interval*>& intervals,
                                                 int threshold) {
  const BamTools::RefVector references = r1.GetReferenceData();
  BamTools::BamAlignment ba1;
  unsigned cnt = 0;
  while (r1.GetNextAlignment(ba1)) {
    if (!ba1.IsPaired() ||
        !ba1.IsMapped() ||
        !ba1.IsMateMapped() ||
        ba1.RefID != ba1.MateRefID ||
        ba1.Position >= ba1.MatePosition ||
        ba1.InsertSize <= threshold) continue;
        
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
        intervals.push_back(new Interval(cnt,
                                         references[ba1.RefID].RefName,
                                         ba1.GetEndPosition() + 1,
                                         ba1.MatePosition));
        cnt++;
      }
    }
  }
}

/*
  @param in should be sorted by start position
 */
void DiscordantPairHandler::removeNestingIntervals(const std::vector<Interval*>& in,
                                                   std::vector<Interval*>& out) {
  std::vector<Point> pts;
  for (size_t i = 0; i < in.size(); i++) {
    pts.push_back(in[i]->createStartPoint());
    pts.push_back(in[i]->createEndPoint());
  }
  sort(pts.begin(), pts.end());
  // for (size_t i = 0; i < pts.size(); i++) {
  //   std::cout << pts[i].intervalId();
  //   if (pts[i].isStart())
  //     std::cout << "S";
  //   else
  //     std::cout << "E";
  //   std::cout << ", ";
  //   if ((i + 1) % 20 == 0) std::cout << std::endl;
  // }
  // return;
  int prev_id = -1;
  for (size_t i = 0; i < pts.size(); i++) {
    int curr_id = pts[i].intervalId();
    if (!pts[i].isStart() && prev_id < curr_id) {
      out.push_back(pts[i].getInterval());
      prev_id = curr_id;
    }
  }
}

void DiscordantPairHandler::clusterIntervals(const std::vector<Interval*>& in,
                                             std::vector<IntervalCluster>& clus,
                                             int threshold) {
  std::vector<Point> pts;
  for (size_t i = 0; i < in.size(); i++) {
    pts.push_back(in[i]->createStartPoint());
    pts.push_back(in[i]->createEndPoint());
  }
  sort(pts.begin(), pts.end());
  std::queue<Interval*> q;
  for (size_t i = 0; i < pts.size(); i++) {
    if (pts[i].isStart())
      q.push(pts[i].getInterval());
    else {
      if (q.empty()) continue;
      IntervalCluster clu;
      while (!q.empty()) {
        clu.add(q.front());
        q.pop();
      }
      clus.push_back(clu);
    }
  }

  for (size_t i = 0; i < clus.size(); i++)
    clus[i].removeInvalidIntervals(threshold);
}
