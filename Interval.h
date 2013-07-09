#ifndef INTERVAL_H
#define INTERVAL_H

#include <string>

class EndPoint;

class Interval {
 private:
  int id;
  int referenceId;
  int startPos;
  int endPos;
  int insertSize;
 public:
  Interval(int id, int referenceId, int startPos, int endPos, int insertSize);
  int getId() const;
  int getReferenceId() const;
  int getStartPos() const;
  int getEndPos() const;
  int getInsertSize() const;
  size_t length() const;
  int minDeletionLength(int mean, int std) const;
  int maxDeletionLength(int mean, int std) const;
  // bool overlapsWith(const Interval& other) const;
  // std::string toString() const;

  const EndPoint getStart() const;
  const EndPoint getEnd() const;
  // static bool compare(const Interval& in1, const Interval& in2);
  
  // static bool compareByStart(const Interval& i1, const Interval& i2);
  // static bool compareByEnd(const Interval& i1, const Interval& i2);
  // friend std::ostream& operator <<(std::ostream& os, const Interval& self);
  
};

class EndPoint {
 public:
  EndPoint(const Interval *owner, bool start);
  int ownerId() const;
  int position() const;
  bool isStart() const;
  bool operator< (const EndPoint& other) const;
  const Interval* getOwner() const;
    
 private:

  const Interval* owner;
  bool start;
};

#endif /* INTERVAL_H */
