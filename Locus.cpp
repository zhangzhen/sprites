#include "Locus.h"

Locus::Locus(std::string chr, int pos)
    : chr(chr), pos(pos) {}

Locus::~Locus() {}

std::string Locus::chrom() const { return chr; }

int Locus::position() const { return pos; }

bool Locus::operator==(const Locus& other) const {
  return chr == other.chr &&
      pos == other.pos;
}
