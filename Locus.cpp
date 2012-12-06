#include "Locus.h"

Locus::Locus(const std::string& chr, int pos)
    : chr(chr), pos(pos) {}

Locus::~Locus() {}

std::string Locus::chrom() const { return chr; }

int Locus::position() const { return pos; }

bool Locus::operator== (const Locus& other) const {
  return chr == other.chr &&
      pos == other.pos;
}

bool Locus::operator!= (const Locus& other) const {
  return !(*this == other);
}

bool Locus::operator< (const Locus& other) const {
  if (chr != other.chr) return chr < other.chr;
  return pos < other.pos;
}
