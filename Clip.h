#ifndef CLIP_INCLUDED
#define CLIP_INCLUDED

#include <string>
using namespace std;

class Clip
{
 public:
  Clip(int refId, int refPosition, int readPosition, int size, const string& readSeq);
  virtual ~Clip();

  int getSize();
  
 private:
  int refId;
  int refPosition;                      // 0-based position
  int readPosition;
  int size;
  string readSeq;
};

#endif // CLIP_INCLUDED
