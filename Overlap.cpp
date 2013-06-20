// #include <iostream>
#include "Overlap.h"

Overlap::Overlap() {}

Overlap::Overlap(int referenceId, int clipPosition1, int clipPosition2, int numClips1, int numClips2, int length, int numMismatches, int offset) :
    referenceId(referenceId),
    clipPosition1(clipPosition1),
    clipPosition2(clipPosition2),
    numClips1(numClips1),
    numClips2(numClips2),
    length(length),
    numMismatches(numMismatches),
    offset(offset) {
  // std::cout << offset << std::endl;
}

// int Overlap::getLength() const { return length; }

// int Overlap::getNumMismatches() const { return numMismatches; }

int Overlap::deletionLength() const {
  return clipPosition2 + offset - clipPosition1;
}

double Overlap::score() const {
  return 0;
}

Deletion Overlap::getDeletion() const {
  // return Deletion(referenceId, clipPosition1, clipPosition2 + offset, clipPosition1-offset, clipPosition2);
  return Deletion(referenceId, clipPosition1, clipPosition2, offset);
}

bool Overlap::equals(const std::string& s1, const std::string& s2, int maxMismatches, int& numMismatches) {
  int cnt = 0;
  if (s1.size() != s2.size()) return false;
  // assert(s1.size() == s2.size());
  for (size_t i = 0; i < s1.size(); ++i) {
    if (s1[i] != s2[i]) ++cnt;
    if (cnt > maxMismatches) return false;
  }
  numMismatches = cnt;
  return true;
}
