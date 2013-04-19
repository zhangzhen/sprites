#include "gtest/gtest.h"
#include "SingleClippedCluster2.h"

class SingleClippedClusterTest : public testing::Test {
 protected:
  virtual void SetUp() {
    a1 = new SingleClipped2(0, 16427375, 48, "ACAATTTGCCAGTGGATCTCCGTTTTACCTCTCTATCCTTGCTTTTATAA", "6.55<''9D22::8C)3.BD?-EE?BCC+FD,>CB01AA?BBBB<A8?@@");
    a2 = new SingleClipped2(0, 16427375, 43, "TTGCCAGTGGATATCTAATTTACCTCTCTATGGTTGCTTTTATAATTCCA", "302+2&7-*5B,%,8$9%%%8=9E29;972G+*8C4?#9E?9<E@CEA87");
    a3 = new SingleClipped2(0, 16427375, 32, "TCTTCATTTGACATATCTACCCTTGCTTTTATAATTCCATCAACATCTT", "-%($\"817)+&,%=#52$',.(:B@7F(&=B=5:9*93>?-79@8C75B");
    clu = new SingleClippedCluster2(a1);
  }

  virtual void TearDown() {
    delete a1;
    delete a2;
    delete a3;
    delete clu;
  }
  
  SingleClipped2 *a1, *a2, *a3;
  SingleClippedCluster2 *clu;
};

TEST_F(SingleClippedClusterTest, ConsensusSingleton) {
  ASSERT_EQ("ACAATTTGCCAGTGGATCTCCGTTTTACCTCTCTATCCTTGCTTTTATAA", clu->consensus());
}

TEST_F(SingleClippedClusterTest, ConsensusTwo) {
  clu->add(a2);
  ASSERT_EQ("TTGCCAGTGGATCTCCATTTTACCTCTCTATCCTTGCTTTTATAA", clu->consensus());
}

TEST_F(SingleClippedClusterTest, ConsensusMany) {
  clu->add(a2);
  clu->add(a3);  
  ASSERT_EQ("TTGCCAGTGGATCTCCATTTTACCTCTCTATCCTTGCTTTTATAATTCCA", clu->consensus());
}

TEST_F(SingleClippedClusterTest, ContigTwo) {
  clu->add(a2);
  Contig2* c = clu->contig();
  ASSERT_EQ(0, c->getReferenceId());
  ASSERT_EQ(16427375, c->getClipPosition());
  ASSERT_EQ(43, c->getLocalClipPosition());
  ASSERT_EQ(45, c->length());
}
