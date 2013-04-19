#ifndef _DELETION_H_
#define _DELETION_H_

#include "Deletion.h"

Deletion::Deletion(const Region* region1, const Region* region2, double score) :
    region1(region1), region2(region2), score(score) {}

Deletion::~Deletion() {}

int Deletion::startPosition1() const { return region1->getStartPosition(); }

int Deletion::endPosition1() const { return region1->getEndPosition(); }

int Deletion::startPosition2() const { return region2->getStartPosition(); }

int Deletion::endPosition2() const { return region2->getEndPosition(); }

double Deletion::score() const { return score; }

#endif /* _DELETION_H_ */
