#ifndef CLIP_INCLUDED
#define CLIP_INCLUDED

#include <string>
#include <sstream>

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

  int getReadPosition() const {
    return readPosition;
  }

  int getPosition() const {
    return refPosition;
  }

  virtual ClipType getType() = 0;
  virtual std::string toString() {
    std::stringstream stream;
    stream << "Clip: " << refId << ", " << refPosition << ", " << readPosition << ", "
           << size << ", " << readSeq;
    return stream.str();
  }
  
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
  virtual std::string toString() {
    return Clip::toString() + " - [L]";
  }
};

class RightClip : public Clip {
 public:
  RightClip(int refId, int refPosition, int readPosition, int size, const std::string& readSeq) :
      Clip(refId, refPosition, readPosition, size, readSeq) {}
  virtual ~RightClip() {}

  virtual ClipType getType();
  virtual std::string toString() {
    return Clip::toString() + " - [R]";
  }

};
#endif // CLIP_INCLUDED
