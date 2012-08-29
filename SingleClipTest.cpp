#include <string>
#include <vector>
#include "gtest/gtest.h"
#include "api/BamReader.h"

using namespace BamTools;


TEST(SingleClipTest, GetSoftClips) {
  std::string filename = "toy2.bam";
  BamReader reader;
  if (!reader.Open(filename)) {
    FAIL() << "Could not open input BAM file.";
  }
  BamAlignment al;
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
