#ifndef _REGIONX_H
#define _REGIONX_H

class RegionX
{
 private:
  int referenceId;
  int startPosition;
  int endPosition;
  
 public:
  RegionX(int referenceId, int startPosition, int endPosition);
  virtual ~RegionX();

  int getReferenceId() const;
  int getStartPosition() const;
  int getEndPosition() const;
  int length() const;

  // bool operator== (const Region& other) const;
  // bool operator!= (const Region& other) const;

 protected:
  bool checkRep();
};

#endif /* _REGIONX_H */
