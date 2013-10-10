#ifndef CHRREGION_H
#define CHRREGION_H

#include <string>

class EndPoint;

class ChrRegion {
 private:
  int id;
  std::string name;
  int referenceId;
  int startPos;
  int endPos;
  int insertSize;
  int readLength;
 public:
  ChrRegion(int id, const std::string& name, int referenceId, int startPos, int endPos, int insertSize, int readLength);
  int getId() const;
  std::string getName() const;
  int getReferenceId() const;
  int getStartPos() const;
  int getEndPos() const;
  int getInsertSize() const;
  int getReadLength() const;
  size_t length() const;
  int minDeletionLength(int mean, int std) const;
  int maxDeletionLength(int mean, int std) const;
  // bool overlapsWith(const ChrRegion& other) const;
  std::string toString() const;

  EndPoint getStart();
  EndPoint getEnd();
  // static bool compare(const ChrRegion& in1, const ChrRegion& in2);

  // static bool compareByStart(const ChrRegion& i1, const ChrRegion& i2);
  // static bool compareByEnd(const ChrRegion& i1, const ChrRegion& i2);
  friend std::ostream& operator <<(std::ostream& os, const ChrRegion& cr);
};

class EndPoint {
 public:
  EndPoint(ChrRegion *owner, bool start);
  int ownerId() const;
  int position() const;
  bool isStart() const;
  bool operator< (const EndPoint& other) const;
  ChrRegion* getOwner() const;

 private:

  ChrRegion* owner;
  bool start;
};

#endif /* CHRREGION_H */
