#include <string>
#include <vector>
#include <algorithm>
#include "gtest/gtest.h"
#include "api/BamReader.h"
#include "clip-sv.h"
#include "LeftClipped.h"
#include "RightClipped.h"
#include "ClippedCreator.h"
#include "LeftClippedCluster.h"
#include "RightClippedCluster.h"
#include "ClusterCreator.h"

// TEST(SingleClipTest, GetSoftClips) {
//   std::string filename = "toy2.bam";
//   BamTools::BamReader reader;
//   if (!reader.Open(filename)) {
//     FAIL() << "Could not open input BAM file.";
//   }
//   BamTools::BamAlignment al;
//   std::vector<int> clipSizes, readPositions, genomePositions;
  
//   reader.GetNextAlignment(al);
//   EXPECT_FALSE(al.GetSoftClips(clipSizes, readPositions, genomePositions));
//   EXPECT_EQ(0, clipSizes.size());
//   EXPECT_EQ(0, readPositions.size());
//   EXPECT_EQ(0, genomePositions.size());
 
//   reader.GetNextAlignment(al);
//   EXPECT_TRUE(al.GetSoftClips(clipSizes, readPositions, genomePositions));
//   EXPECT_EQ("AAAAGATAAGGATA", al.QueryBases);
//   EXPECT_EQ("AGATAA*GGATA", al.AlignedBases);
//   EXPECT_EQ(8, al.Position);
//   EXPECT_EQ(0, al.RefID);
//   EXPECT_EQ(1, clipSizes.size());
//   EXPECT_EQ(1, readPositions.size());
//   EXPECT_EQ(1, genomePositions.size());
//   EXPECT_EQ(3, clipSizes[0]);
//   EXPECT_EQ(3, readPositions[0]);
//   EXPECT_EQ(8, genomePositions[0]);

//   reader.GetNextAlignment(al);
//   clipSizes.clear();
  
//   readPositions.clear();
//   genomePositions.clear();
//   EXPECT_TRUE(al.GetSoftClips(clipSizes, readPositions, genomePositions));
//   EXPECT_EQ(8, al.Position);  
//   EXPECT_EQ(1, clipSizes.size());
//   EXPECT_EQ(1, readPositions.size());
//   EXPECT_EQ(1, genomePositions.size());
//   EXPECT_EQ(3, clipSizes[0]);
//   EXPECT_EQ(11, readPositions[0]);
//   EXPECT_EQ(18, genomePositions[0]);

//   reader.Close();

// }

TEST(SingleClipTest, equals) {
  std::string s1 = "AGGTACT";
  std::string s2 = "ACGTACT";
  EXPECT_FALSE(equals(s1, s2));
  EXPECT_TRUE(equals(s1, s2, 1));
}

TEST(SingleClipTest, toString) {
  LeftClip lcl = LeftClip(0, 27, 3, 3, "ATAGTAGGCA");
  EXPECT_EQ("Clip: 0, 27, 3, 3, ATAGTAGGCA - [L]", lcl.toString());
}

// TEST(SingleClipTest, reverseIterator) {
//   std::vector<int> myvector;
//   for (int i=1; i<=5; i++) myvector.push_back(i);
//   std::vector<int>::reverse_iterator rit;
//   for ( rit=myvector.rbegin() ; rit < myvector.rend();) {
//     if (*rit == 4) {
//       std::vector<int>::iterator tempIt = myvector.erase(--rit.base());
//       rit = std::vector<int>::reverse_iterator(tempIt);
//     }
//     else ++rit;
//   }
//   EXPECT_EQ(4, myvector.size());
//   EXPECT_EQ(3, myvector[2]);
//   EXPECT_EQ(5, myvector[3]);
// }

void createGraph(std::vector<std::vector<int> >& g) {
  std::vector<int> nodes;
  nodes.push_back(0);
  nodes.push_back(3);
  g.push_back(nodes);
  nodes.clear();
  nodes.push_back(2);
  nodes.push_back(4);
  g.push_back(nodes);
  nodes.clear();
  nodes.push_back(0);
  nodes.push_back(3);
  g.push_back(nodes);
  nodes.clear();
  nodes.push_back(1);
  nodes.push_back(2);
  nodes.push_back(4);
  g.push_back(nodes);
}

TEST(SingleClipTest, match) {
  std::vector<std::vector<int> > g;
  int m = 4, n = 5;
  createGraph(g);
  Matching mat(g, m, n);
  EXPECT_EQ(4, mat.match());
  EXPECT_EQ(3, mat.getMateL(0));
  EXPECT_EQ(2, mat.getMateL(1));
  EXPECT_EQ(0, mat.getMateL(2));
  EXPECT_EQ(1, mat.getMateL(3));  
  EXPECT_EQ(2, mat.getMateR(0));
}

TEST(SingleClipTest, findOverlaps) {
  StructVar t1("22", 30, 80);
  StructVar t2("22", 340, 410);
  std::vector<StructVar> cs1;
  StructVar c1("22", 200, 220);
  cs1.push_back(c1);
  StructVar c2("22", 10, 40);
  cs1.push_back(c2);
  StructVar c3("22", 10, 30);
  cs1.push_back(c3);
  StructVar c4("22", 410, 500);
  cs1.push_back(c4);
  StructVar c5("22", 80, 150);
  cs1.push_back(c5);
  StructVar c6("22", 300, 340);
  cs1.push_back(c6);
  std::vector<StructVar> cs2 = cs1;
  sort(cs1.begin(), cs1.end(), cmp1);
  sort(cs2.begin(), cs2.end(), cmp2);
  std::set<StructVar> s1 = findOverlaps(t1, cs1, cs2);
  EXPECT_EQ(3, s1.size());
  EXPECT_EQ(1, s1.count(c2));
  EXPECT_EQ(0, s1.count(c4));
}

TEST(SingleClipTest, groupBreakpoints) {
  std::vector<Breakpoint> bps;
  bps.push_back(Breakpoint(10, 30));
  bps.push_back(Breakpoint(10, 40));
  bps.push_back(Breakpoint(30, 80));
  bps.push_back(Breakpoint(80, 150));
  bps.push_back(Breakpoint(200, 220));
  bps.push_back(Breakpoint(340, 410));
  bps.push_back(Breakpoint(380, 500));
  std::vector<std::vector<Breakpoint> > groups;
  groupBreakpoints(bps, groups);
  EXPECT_EQ(4, groups.size());
  EXPECT_EQ(3, groups[0].size());
  EXPECT_EQ(1, groups[1].size());
  EXPECT_EQ(1, groups[2].size());
  EXPECT_EQ(2, groups[3].size());
}

TEST(SingleClipTest, loadClippeds) {
  std::string filename = "toy2.bam";
  BamTools::BamReader reader;
  reader.Open(filename);
  std::vector<SingleClipped*> lefts, rights;
  loadClippeds(reader, lefts, rights);
  EXPECT_EQ(1, lefts.size());
  EXPECT_EQ(1, rights.size());
}

TEST(SingleClipTest, clusterClippeds) {
  Locus a1("1", 27);
  Locus a2("1", 54);
  StandardClippedCreator<LeftClipped> cliCreator;
  std::vector<SingleClipped*> lefts;
  lefts.push_back(cliCreator.createClipped(a1, "ATAGTAGGCA", 20, 0, 3));
  lefts.push_back(cliCreator.createClipped(a1, "AGATAGTAGG", 20, 0, 5));
  lefts.push_back(cliCreator.createClipped(a2, "CCATAACTAC", 20, 0, 4));
  lefts.push_back(cliCreator.createClipped(a2, "ATAACTACGC", 20, 0, 2));
  StandardClusterCreator<LeftClippedCluster> cluCreator;
  std::vector<SingleClippedCluster*> clus;
  clusterClippeds(lefts, clus, cluCreator, 2);
  EXPECT_EQ(2, clus.size());
  EXPECT_EQ(2, clus[0]->size());
  EXPECT_EQ(2, clus[1]->size());

  Locus a3("1", 13);
  StandardClippedCreator<RightClipped> cliCreator2;
  std::vector<SingleClipped*> rights;
  rights.push_back(cliCreator2.createClipped(a3, "TAGATAGTAG", 20, 6, 4));
  rights.push_back(cliCreator2.createClipped(a3, "TTAGATAGTA", 20, 7, 3));
  StandardClusterCreator<RightClippedCluster> cluCreator2;
  std::vector<SingleClippedCluster*> clus2;
  clusterClippeds(rights, clus2, cluCreator2, 2);
  EXPECT_EQ(1, clus2.size());
  EXPECT_EQ(2, clus2[0]->size());
}

TEST(SingleClipTest, findFirstRegion) {
  std::vector<Contig> cons;
  Locus a1("1", 13);
  Contig c1("TTAGATAGTAG", a1, 7, 2);
  Locus a2("1", 27);
  Contig c2("AGATAGTAGGCA", a2, 5, 2);
  cons.push_back(c2);
  Locus a3("1", 54);
  Contig c3("CCATAACTACGC", a3, 4, 2);
  cons.push_back(c3);

  Region r1;
  EXPECT_TRUE(findFirstRegion(cons.begin(), cons.end(), c1, 0.0, r1));
  Region r2(a1, a2);
  EXPECT_EQ(r2, r1);
}

TEST(SingleClipTest, callDeletions) {
  std::vector<Contig> cons1;
  std::vector<Contig> cons2;
  
  Locus fstA1("1", 13);
  Contig fstC1("TTAGATAGTAG", fstA1, 7, 2);
  cons1.push_back(fstC1);
  Locus fstA2("1", 45);
  Contig fstC2("CAGCGCCATAACTA", fstA2, 9, 2);
  cons1.push_back(fstC2);
  
  Locus sndA1("1", 27);
  Contig sndC1("AGATAGTAGGCA", sndA1, 5, 2);
  cons2.push_back(sndC1);
  Locus sndA2("1", 54);
  Contig sndC2("CCATAACTACGC", sndA2, 4, 2);
  cons2.push_back(sndC2);

  std::vector<Region> calls;
  callDeletions(cons1, cons2, calls, 0.0);
  EXPECT_EQ(2, calls.size());
  EXPECT_EQ(Region(fstA2, sndA2), calls[0]);
  EXPECT_EQ(Region(fstA1, sndA1), calls[1]);
}
