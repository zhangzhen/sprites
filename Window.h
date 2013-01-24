#ifndef WINDOW_H
#define WINDOW_H

#include <string>

class Window {
 private:
  std::string name;
  std::string chrom;
  unsigned int start;
  unsigned int end;
  int insertSize;
  
 public:
  Window();
  Window(const std::string& name, const std::string& chrom, unsigned int start, unsigned int end, int insertSize);
  virtual ~Window();
  std::string getName() { return name; }
  std::string getChrom() { return chrom; }
  unsigned int getStart() { return start; }
  unsigned int getEnd() { return end; }
  int getInsertSize() { return insertSize; }
  unsigned int length() { return end - start; }
};

#endif /* WINDOW_H */
