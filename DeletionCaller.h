class DeletionCaller {
 public:
  static void callAll(const std::vector<Region2>& in,
               std::vector<Contig>& cons1,
               std::vector<Contig>& cons2,
               std::vector<Region>& calls,
               int minOverlapLen,
               int maxMismatches);
 private:
  static void DeletionCaller::callOne(std::vector<Contig>::iterator first1,
                                      std::vector<Contig>::iterator last1,
                                      std::vector<Contig>::iterator first2,
                                      std::vector<Contig>::iterator last2,
                                      std::vector<Region>& calls,
                                      int minOverlapLen,
                                      int maxMismatches);
};
