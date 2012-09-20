#include "clip-sv.h"
#include "gtest/gtest.h"

class CallDelsTest : public testing::Test {
 protected:
  virtual void SetUp() {
    // donor: AGCATGTTAGATA*GTAGGCAGTCAGCGCCAT*aactacgcg
    // ref: AGCATGTTAGATA-AGATAGCTGTGCTA-GTAGGCAGTCAGCGCCAT-ctacagagc-aactacgcg
    lclips.push_back(new LeftClip(0, 27, 3, 3, "ATAGTAGGCA"));
    lclips.push_back(new LeftClip(0, 27, 5, 5, "AGATAGTAGG"));
    lclips.push_back(new LeftClip(0, 32, 4, 4, "CCATaactac"));
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
  extractClipsForDels(lclips, rclips, LCs, RCs);
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
  extractClipsForDels(lclips, rclips, LCs, RCs);
  EXPECT_EQ(2, LCs.size());
  EXPECT_EQ(2, RCs.size());
}

TEST_F(CallDelsTest, clusterBreakpoints) {
  std::vector<Breakpoint> bps;
  buildBreakpoints(lclips, rclips, bps);
  std::vector<std::vector<Breakpoint> > clusters;
  clusterBreakpoints(bps, clusters);
  EXPECT_EQ(1, clusters.size());
  EXPECT_EQ(2, clusters[0].size());
}
