#include "Deletion.h"
#include <cassert>

using namespace std;

Deletion::Deletion()
{
}

Deletion::Deletion(int id,
                   std::string referenceName,
                   int leftBp,
                   int rightBp,
                   int length,
                   std::string alternative,
                   std::string homseq,
                   std::string genotype) :
    id(id),
    referenceName(referenceName),
    leftBp(leftBp),
    rightBp(rightBp),
    length(length),
    alternative(alternative),
    homseq(homseq),
    genotype(genotype) {
    assert(checkRep());
}

Deletion::~Deletion() {
}

int Deletion::getId() const
{
    return id;
}
std::string Deletion::getReferenceName() const
{
    return referenceName;
}

int Deletion::getLeftBp() const
{
    return leftBp;
}

int Deletion::getRightBp() const
{
    return rightBp;
}

int Deletion::getLength() const
{
    return length;
}

std::string Deletion::getAlternative() const
{
    return alternative;
}

std::string Deletion::getHomseq() const
{
    return homseq;
}

std::string Deletion::getGenotype() const
{
    return genotype;
}

std::ostream& operator <<(ostream &stream, const Deletion &del)
{
    stream << del.id << "\t" << del.referenceName << "\t"
           << del.leftBp << "\t" << del.rightBp << "\tDEL\t"
           << del.length << "\t" << del.alternative << "\t"
           << del.homseq << "\t" << del.genotype;
    return stream;
}

bool Deletion::checkRep() const
{
    return (rightBp > leftBp) &&
            (length <= SVLEN_THRESHOLD);
}
