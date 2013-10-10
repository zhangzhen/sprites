#ifndef _DELETION_H_
#define _DELETION_H_

class Deletion
{
 public:
    Deletion();
  Deletion(int referenceId, int start2, int end1, int offset);
  int getReferenceId() const;
  int getStart1() const;
  int getEnd1() const;
  int getStart2() const;
  int getEnd2() const;
  int length() const;
  bool overlaps(const Deletion& other) const;
  bool operator <(const Deletion& other) const;
  // double score() const;
private:
  int referenceId;
  int start2;
  int end1;
  int offset;
  // double score;
};

#endif /* _DELETION_H_ */
