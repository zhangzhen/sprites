#ifndef _DELETION_H_
#define _DELETION_H_

class Deletion
{
 public:
  Deletion(const Region* region1, const Region* region2, double score);
  virtual ~Deletion();
  int referenceId() const;
  // const std::string& referenceName() const;
  int startPosition1() const;
  int endPosition1() const;
  int startPosition2() const;
  int endPosition2() const;
  double score() const;
 private:
  const Region* region1;
  const Region* region2;
  double score;
};

#endif /* _DELETION_H_ */
