#ifndef CHRREGION_H
#define CHRREGION_H

#include <string>

class EndPoint;

class ChrRegion {
 private:
  int id;
  int referenceId;
  int startPos;
  int endPos;
  int insertSize;
 public:
  ChrRegion(int id, int referenceId, int startPos, int endPos, int insertSize);
  int getId() const;
  int getReferenceId() const;
  int getStartPos() const;
  int getEndPos() const;
  int getInsertSize() const;
  size_t length() const;
  int minDeletionLength(int mean, int std) const;
  int maxDeletionLength(int mean, int std) const;
  // bool overlapsWith(const ChrRegion& other) const;
  // std::string toString() const;

  const EndPoint getStart() const;
  const EndPoint getEnd() const;
  // static bool compare(const ChrRegion& in1, const ChrRegion& in2);
  
  // static bool compareByStart(const ChrRegion& i1, const ChrRegion& i2);
  // static bool compareByEnd(const ChrRegion& i1, const ChrRegion& i2);
  friend std::ostream& operator <<(std::ostream& os, const ChrRegion& self);
  
};

class EndPoint {
 public:
  EndPoint(const ChrRegion *owner, bool start);
  int ownerId() const;
  int position() const;
  bool isStart() const;
  bool operator< (const EndPoint& other) const;
  const ChrRegion* getOwner() const;
    
 private:

  const ChrRegion* owner;
  bool start;
};

#endif /* CHRREGION_H */
