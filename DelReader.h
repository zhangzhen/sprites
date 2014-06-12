#ifndef DELREADER_H
#define DELREADER_H

#include <string>
#include <vector>
#include <set>


struct Del {
    int id;
    std::string referenceName;
    int leftBp;
    int rightBp;
    int length;
    std::string alternative;
    std::string homseq;
    std::string genotype;

    bool hasInsertedSeq() const {
        return alternative.length() > 1;
    }

    bool hasHomseq() const {
        return homseq == "-";
    }

    bool isHomogeneous() const {
        return genotype == "1/1" || genotype == "1|1";
    }
};

class DelReader
{
public:
    static std::vector<Del> readDelsFromFile(const std::string& filename);
};

#endif // DELREADER_H
