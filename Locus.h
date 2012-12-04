#ifndef LOCUS_H
#define LOCUS_H

#include <string>

class Locus
{
 private:
  std::string chr;
  int pos;
  
 public:
  Locus(std::string chr, int pos);
  virtual ~Locus();
  std::string chrom() const;
  int position() const;
  bool operator==(const Locus& other) const;
};

#endif /* LOCUS_H */
