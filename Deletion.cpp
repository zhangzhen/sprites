#include "Deletion.h"
#include "Helper.h"
#include <cassert>

using namespace std;

Deletion::Deletion(std::string referenceName,
                   int leftBp,
                   int rightBp,
                   int length,
                   std::string alternative,
                   std::string homseq,
                   std::string genotype) :
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

bool Deletion::contains(const Deletion &other)
{
    return leftBp <= other.leftBp && rightBp >= other.rightBp;
}

bool Deletion::dovetailsTo(const Deletion &other)
{
    return leftBp < other.leftBp && rightBp < other.rightBp && rightBp > other.leftBp;
}

std::ostream& operator <<(ostream &stream, const Deletion &del)
{
    stream << del.referenceName << "\t"
           << del.leftBp << "\t" << del.rightBp << "\tDEL\t"
           << del.length << "\t" << del.alternative << "\t"
           << del.homseq << "\t" << del.genotype;
    return stream;
}

bool Deletion::checkRep() const
{
    return (rightBp > leftBp) &&
            (length <= Helper::SVLEN_THRESHOLD);
}
