#include "clip-sv.h"
#include "Contig.h"
#include "LeftClipped.h"
#include "LeftClippedCluster.h"
#include "RightClipped.h"
#include "RightClippedCluster.h"
#include "gtest/gtest.h"

class CallDelsTest : public testing::Test {
 protected:
  virtual void SetUp() {
    // donor: AGCATGTTAGATA*GTAGGCAGTCAGCGCCAT*AACTACGCG
    // ref: AGCATGTTAGATA(AGATAGCTGTGCTA)GTAGGCAGTCAGCGCCAT(CTACAGAGC)AACTACGCG
    lclips.push_back(new LeftClip(0, 27, 3, 3, "ATAGTAGGCA"));
    lclips.push_back(new LeftClip(0, 27, 5, 5, "AGATAGTAGG"));
    lclips.push_back(new LeftClip(0, 54, 4, 4, "CCATAACTAC"));
    rclips.push_back(new RightClip(0, 13, 6, 4, "TAGATAGTAG"));
    rclips.push_back(new RightClip(0, 13, 7, 3, "TTAGATAGTA"));
  }

  virtual void TearDown() {
    freeClips(lclips);
    freeClips(rclips);
  }

  std::vector<Clip*> lclips, rclips;
};

// TEST_F(CallDelsTest, countOverlappedReads) {
//   EXPECT_EQ(2, countOverlappedReads(lclips, rclips[0]));      
//   EXPECT_EQ(2, countOverlappedReads(rclips, lclips[0]));
// }

// TEST_F(CallDelsTest, matchClips) {
//   std::map<Clip*, std::vector<Clip*> > matches;
//   matchClips(lclips, rclips, matches);
//   EXPECT_EQ(2, matches[lclips[0]].size());
//   EXPECT_EQ(2, matches[lclips[1]].size());
//   EXPECT_EQ(2, matches[rclips[0]].size());
//   EXPECT_EQ(2, matches[rclips[1]].size());
// }

TEST_F(CallDelsTest, createOverlapGraph) {
  std::vector<Clip*> LCs, RCs;
  extractClipsForDels(lclips, rclips, LCs, RCs, 4, 10);
  std::vector<std::vector<int> > g;
  createOverlapGraph(LCs, RCs, g);
  EXPECT_EQ(2, g.size());
  EXPECT_EQ(2, g[0].size());
  EXPECT_EQ(2, g[1].size());
  EXPECT_EQ(0, g[0][0]);
  EXPECT_EQ(1, g[0][1]);
}

// TEST_F(CallDelsTest, clusterClips) {
//   std::map<int, std::set<Clip*> > clusters;
//   clusterClips(lclips, clusters);
//   clusterClips(rclips, clusters);
//   EXPECT_EQ(1, clusters.count(27));
//   EXPECT_EQ(2, clusters[27].size());
//   EXPECT_EQ(1, clusters.count(13));
//   EXPECT_EQ(2, clusters[13].size());
// }

TEST_F(CallDelsTest, extractClipsForDels) {
  std::vector<Clip*> LCs, RCs;
  extractClipsForDels(lclips, rclips, LCs, RCs, 4, 10);
  EXPECT_EQ(2, LCs.size());
  EXPECT_EQ(2, RCs.size());
}

// TEST_F(CallDelsTest, clusterBreakpoints) {
//   std::vector<Breakpoint> bps;
//   buildBreakpoints(lclips, rclips, bps);
//   std::vector<std::vector<Breakpoint> > clusters;
//   clusterBreakpoints(bps, clusters);
//   EXPECT_EQ(1, clusters.size());
//   EXPECT_EQ(2, clusters[0].size());
// }

// TEST_F(CallDelsTest, overlaps) {
//   Contig c1("AGATAGTAGGCA", Locus("1", 27), 5, 2, false);
//   Contig c2("TTAGATAGTA", Locus("1", 13), 7, 2, true);
//   int offset = 0;
//   EXPECT_TRUE(c1.overlaps(c2, 2, 8, 0.0, offset));
//   offset = 0;
//   EXPECT_TRUE(c2.overlaps(c1, 2, 8, 0.0, offset));
// }

// TEST_F(CallDelsTest, assembleContigs) {
//   Locus a1("1", 27);
//   LeftClippedCluster lcc(a1);
//   SingleClippedCluster& cls1 = lcc;
//   SingleClipped *p1 = new LeftClipped(a1, "ATAGTAGGCA", 20, 0, 3);
//   SingleClipped *p2 = new LeftClipped(a1, "AGATAGTAGG", 20, 0, 5);
//   cls1.add(p1);
//   cls1.add(p2);
//   Contig c1("AGATAGTAGGCA", a1, 5, 2, false);
//   EXPECT_EQ(c1, cls1.contig());
//   delete p1;
//   delete p2;

//   Locus a2("1", 13);
//   RightClippedCluster rcc(a2);
//   SingleClippedCluster& cls2 = rcc;
//   SingleClipped *p3 = new RightClipped(a2, "TAGATAGTAG", 13, 6, 4);
//   SingleClipped *p4 = new RightClipped(a2, "TTAGATAGTA", 13, 7, 3);
//   cls2.add(p3);
//   cls2.add(p4);
//   Contig c2("TTAGATAGTAG", a2, 7, 2, true);
//   EXPECT_EQ("TTAGATAGTAG", c2.sequence());
//   EXPECT_EQ(c2, cls2.contig());
//   delete p3;
//   delete p4;
// }
