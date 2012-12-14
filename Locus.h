#ifndef LOCUS_H
#define LOCUS_H

#include <string>

class Locus
{
 private:
  std::string chr;
  int pos;
  
 public:
  Locus();
  Locus(const std::string& chr, int pos);
  virtual ~Locus();
  std::string chrom() const;
  int position() const;
  bool operator== (const Locus& other) const;
  bool operator!= (const Locus& other) const;
  bool operator< (const Locus& other) const;
  friend std::ostream& operator <<(std::ostream& stream, const Locus& self);
};

#endif /* LOCUS_H */
