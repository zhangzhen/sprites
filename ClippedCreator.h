#ifndef CLIPPED_CREATOR_H
#define CLIPPED_CREATOR_H

#include "SingleClipped.h"

class ClippedCreator
{
 public:
  virtual SingleClipped* createClipped(const Locus& loc, const std::string& seq, int qual, int start, int clippedLen) = 0;
};

template <class TheClipped>
class StandardClippedCreator : public ClippedCreator {
 public:
  virtual SingleClipped* createClipped(const Locus& loc, const std::string& seq, int qual, int start, int clippedLen);
};

template <class TheClipped>
SingleClipped* StandardClippedCreator<TheClipped>::createClipped(const Locus& loc, const std::string& seq, int qual, int start, int clippedLen) {
  return new TheClipped(loc, seq, qual, start, clippedLen);
}
#endif /* CLIPPED_CREATOR_H */
