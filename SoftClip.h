#ifndef _SOFTCLIP_H_
#define _SOFTCLIP_H_

#include <string>
#include "Overlap.h"

class SoftClip {
 private:
  int refId;
  int pos;
  int clipPos;
  std::string seq;
  std::string quals;

 public:
  SoftClip(int refId, int pos, int clipPos, const std::string& seq, const std::string& quals);
  int referenceId() const;
  int position() const;
  int startPosition() const;
  int endPosition() const;
  const std::string& sequence() const;
  int length() const;
  int lengthOfLeftPart() const;
  int lengthOfRightPart() const;
  char at(int i) const;
  char qual(int i) const;
  int minDeletionLength(const SoftClip& other) const;
  int maxDeletionLength(const SoftClip& other) const;
  bool overlaps(const SoftClip& other, int minOverlapLength, double maxMismatchRate, Overlap& overlap) const;

  static bool compareL(SoftClip* s1, SoftClip* s2);
  static bool compareR(SoftClip* s1, SoftClip* s2);
  static bool compare1(SoftClip* o, const int pos);
  static bool compare2(const int pos, SoftClip* o);

  friend std::ostream& operator <<(std::ostream& stream, const SoftClip& o);
};

#endif /* _SOFTCLIP_H_ */
