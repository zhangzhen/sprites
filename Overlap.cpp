// #include <iostream>
#include <iomanip>
#include <cassert>
#include "Overlap.h"
#include "SoftClip.h"

Overlap::Overlap() {
}

Overlap::Overlap(const SoftClip *first, const SoftClip *second, int length, int numMismatches, int offset) :
    first(first),
    second(second),
    length(length),
    numMismatches(numMismatches),
    offset(offset) {
  assert(first != NULL && second != NULL && first->referenceId() == second->referenceId());
  // std::cout << offset << std::endl;
}

// int Overlap::getLength() const { return length; }

// int Overlap::getNumMismatches() const { return numMismatches; }

int Overlap::deletionLength() const {
  return second->position() + offset - first->position();
}

double Overlap::score() const {
  return numMismatches / float(length);
}

Deletion Overlap::getDeletion() const {
  // return Deletion(referenceId, clipPosition1, clipPosition2 + offset, clipPosition1-offset, clipPosition2);
  return Deletion(first->referenceId(), first->position(), second->position(), offset);
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

// bug to fix in getBestOverlap
bool Overlap::getHighScoreOverlap(std::vector<Overlap> overlaps, Overlap& ov) {
  if (overlaps.size() == 0) return false;
  ov = overlaps.front();
  double min = ov.score();
  for (auto itr = overlaps.begin(); itr != overlaps.end(); ++itr) {
    if (min > (*itr).score()) {
      min = (*itr).score();
      ov = (*itr);
    }
  }
  return true;
}

std::ostream& operator <<(std::ostream& stream, const Overlap& o) {
  stream << "[CALL] " << o.first->referenceId() << ":" << o.first->position()
         << "-" << o.second->position() << "\t" << o.deletionLength() << "\t"
         << o.score() << std::endl;
  if (o.first->lengthOfLeftPart() >= o.second->lengthOfLeftPart() + o.offset) {
    stream << std::string(o.first->lengthOfLeftPart(), ' ') << '+' << std::endl;
    stream << o.first->sequence() << std::endl;
    stream << std::string(o.first->lengthOfLeftPart() - o.second->lengthOfLeftPart() - o.offset, ' ') << o.second->sequence() << std::endl;
    stream << std::string(o.first->lengthOfLeftPart() - o.offset, ' ') << '^' << std::endl;
  } else {
    stream << std::string(o.second->lengthOfLeftPart() + o.offset, ' ') << '+' << std::endl;
    stream << std::string(o.second->lengthOfLeftPart() - o.first->lengthOfLeftPart() + o.offset, ' ') << o.first->sequence() << std::endl;
    stream << o.second->sequence() << std::endl;
    stream << std::string(o.second->lengthOfLeftPart(), ' ') << '^' << std::endl;
  }
  return stream;
}
