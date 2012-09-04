#ifndef CLIP_INCLUDED
#define CLIP_INCLUDED

#include <string>

enum ClipType {
  Left, Right
};

class Clip
{
 public:
  Clip() {}
  
  Clip(int refId, int refPosition, int readPosition, int size, const std::string& readSeq);
  virtual ~Clip();

  int getSize() const {
    return size;
  }

  virtual ClipType getType() = 0;
  
  std::string getReadSeq() const {
    return readSeq;
  }
  
 private:
  int refId;
  int refPosition;                      // 0-based position
  int readPosition;
  int size;
  std::string readSeq;
};

class LeftClip : public Clip {
 public:
  LeftClip(int refId, int refPosition, int readPosition, int size, const std::string& readSeq) :
      Clip(refId, refPosition, readPosition, size, readSeq) {}
  virtual ~LeftClip() {}

  virtual ClipType getType();
};

class RightClip : public Clip {
 public:
  RightClip(int refId, int refPosition, int readPosition, int size, const std::string& readSeq) :
      Clip(refId, refPosition, readPosition, size, readSeq) {}
  virtual ~RightClip() {}

  virtual ClipType getType();

};
#endif // CLIP_INCLUDED
