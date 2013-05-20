#ifndef _SOFTCLIP_H_
#define _SOFTCLIP_H_

#include <string>

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
  const std::string& sequence() const;
  int length() const;
  int lengthOfLeftPart() const;
  int lengthOfRightPart() const;
  char at(int i) const;
  char qual(int i) const;

};

#endif /* _SOFTCLIP_H_ */
