#ifndef _SINGLECLIPPED2_H_
#define _SINGLECLIPPED2_H_

#include "Locus.h"

class SingleClipped2 {
 private:
  int refId;
  int pos;
  int clipPos;
  std::string seq;
  std::string quals;
  
 public:
  SingleClipped2(int refId, int pos, int clipPos, const std::string& seq, const std::string& quals);
  virtual ~SingleClipped2();
  int referenceId() const;
  int position() const;
  const std::string& sequence() const;
  int length() const;
  int lengthOfLeftPart() const;
  int lengthOfRightPart() const;
  char at(int i) const;
  char qual(int i) const;

};

#endif /* _SINGLECLIPPED2_H_ */
