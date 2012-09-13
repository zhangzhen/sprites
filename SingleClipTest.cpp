#include <string>
#include <vector>
#include "gtest/gtest.h"
#include "api/BamReader.h"
#include "clip-sv.h"

TEST(SingleClipTest, GetSoftClips) {
  std::string filename = "toy2.bam";
  BamTools::BamReader reader;
  if (!reader.Open(filename)) {
    FAIL() << "Could not open input BAM file.";
  }
  BamTools::BamAlignment al;
  std::vector<int> clipSizes, readPositions, genomePositions;
  
  reader.GetNextAlignment(al);
  EXPECT_FALSE(al.GetSoftClips(clipSizes, readPositions, genomePositions));
  EXPECT_EQ(0, clipSizes.size());
  EXPECT_EQ(0, readPositions.size());
  EXPECT_EQ(0, genomePositions.size());
 
  reader.GetNextAlignment(al);
  EXPECT_TRUE(al.GetSoftClips(clipSizes, readPositions, genomePositions));
  EXPECT_EQ("AAAAGATAAGGATA", al.QueryBases);
  EXPECT_EQ("AGATAA*GGATA", al.AlignedBases);
  EXPECT_EQ(8, al.Position);
  EXPECT_EQ(0, al.RefID);
  EXPECT_EQ(1, clipSizes.size());
  EXPECT_EQ(1, readPositions.size());
  EXPECT_EQ(1, genomePositions.size());
  EXPECT_EQ(3, clipSizes[0]);
  EXPECT_EQ(3, readPositions[0]);
  EXPECT_EQ(8, genomePositions[0]);

  reader.GetNextAlignment(al);
  clipSizes.clear();
  
  readPositions.clear();
  genomePositions.clear();
  EXPECT_TRUE(al.GetSoftClips(clipSizes, readPositions, genomePositions));
  EXPECT_EQ(8, al.Position);  
  EXPECT_EQ(1, clipSizes.size());
  EXPECT_EQ(1, readPositions.size());
  EXPECT_EQ(1, genomePositions.size());
  EXPECT_EQ(3, clipSizes[0]);
  EXPECT_EQ(11, readPositions[0]);
  EXPECT_EQ(18, genomePositions[0]);

  reader.Close();

}

TEST(SingleClipTest, countMismatches) {
  std::string s1 = "AGGTACT";
  std::string s2 = "ACGTACT";
  EXPECT_EQ(0, countMismatches(s1, s1));  
  EXPECT_EQ(1, countMismatches(s1, s2));
}

TEST(SingleClipTest, countClipsInLengthOneToFive) {
  std::vector<Clip*> lclips;
  lclips.push_back(new LeftClip(0, 27, 3, 3, "ATAGTAGGCA"));
  lclips.push_back(new LeftClip(0, 27, 3, 3, "ATAGTAGGCA"));
  lclips.push_back(new LeftClip(0, 27, 5, 5, "AGATAGTAGG"));
  int stats[5];
  countClipsInLengthOneToFive(lclips, stats);
  EXPECT_EQ(0, stats[0]);
  EXPECT_EQ(0, stats[1]);
  EXPECT_EQ(2, stats[2]);
  EXPECT_EQ(0, stats[3]);
  EXPECT_EQ(1, stats[4]);
  while (!lclips.empty()) delete lclips.back(), lclips.pop_back();  
}

TEST(SingleClipTest, toString) {
  LeftClip lcl = LeftClip(0, 27, 3, 3, "ATAGTAGGCA");
  EXPECT_EQ("Clip: 0, 27, 3, 3, ATAGTAGGCA - [L]", lcl.toString());
}
