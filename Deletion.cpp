#include "Deletion.h"
#include "Helper.h"
#include <cassert>
#include <sstream>

using namespace std;

Deletion::Deletion(const string &referenceName,
                   int start1,
                   int end1,
                   int start2,
                   int end2,
                   int length,
                   const string& fromTag) :
    referenceName(referenceName),
    start1(start1),
    end1(end1),
    start2(start2),
    end2(end2),
    length(length),
    fromTag(fromTag) {
    assert(checkRep());
}

Deletion::~Deletion() {
}

string Deletion::toBedpe() const {
    stringstream fmt;
    fmt << referenceName << "\t" << start1 - 1 << "\t" << end1 << "\t"
              << referenceName << "\t" << start2 - 1 << "\t" << end2;
    return fmt.str();
}

bool Deletion::overlaps(const Deletion &other) const
{
    if (referenceName != other.referenceName) return false;
    return ((start1-1 >= other.start1-1 && start1-1 <= other.end1) ||
            (other.start1-1 >= start1-1 && other.start1-1 <= end1)) &&
            ((start2-1 >= other.start2-1 && start2-1 <= other.end2) ||
             (other.start2-1 >= start2-1 && other.start2-1 <= end2));
}

bool Deletion::operator<(const Deletion &other) const
{
    if (referenceName != other.referenceName) return referenceName < other.referenceName;
    if (start1 != other.start1) return start1 < other.start1;
    if (start2 != other.start2) return start2 < other.start2;
    if (end1 != other.end1) return end1 < other.end1;
    return end2 < other.end2;
}

bool Deletion::operator==(const Deletion &other) const
{
    return referenceName == other.referenceName &&
            start1 == other.start1 && start2 == other.start2 &&
            end1 == other.end1 && end2 == other.end2;
}

std::ostream& operator <<(ostream &stream, const Deletion &del)
{
    stream << del.toBedpe();
    return stream;
}

bool Deletion::checkRep() const
{
    return (start1 <= end1) &&
            (start2 <= end2) &&
            (length <= Helper::SVLEN_THRESHOLD);
}
