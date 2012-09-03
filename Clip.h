#ifndef CLIP_INCLUDED
#define CLIP_INCLUDED

#include <string>
using namespace std;

class Clip
{
 public:
  Clip() {}
  
  Clip(int refId, int refPosition, int readPosition, int size, const string& readSeq);
  virtual ~Clip();

  int getSize() const {
    return size;
  }
  
  string getReadSeq() const {
    return readSeq;
  }
  
 private:
  int refId;
  int refPosition;                      // 0-based position
  int readPosition;
  int size;
  string readSeq;
};

#endif // CLIP_INCLUDED
