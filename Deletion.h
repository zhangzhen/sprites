#ifndef _DELETION_H_
#define _DELETION_H_

#include <string>
#include <iostream>

class Deletion {
public:
    Deletion(int id, std::string referenceName, int leftBp, int rightBp, int length,
             std::string alternative=".", std::string homseq=".", std::string genotype=".");

    virtual ~Deletion();

    bool hasInsertedSeq() const {
        return alternative.length() > 1;
    }

    bool hasHomseq() const {
        return homseq == "-";
    }

    bool isHomogeneous() const {
        return genotype == "1/1" || genotype == "1|1";
    }

    int getId() const;

    std::string getReferenceName() const;

    int getLeftBp() const;

    int getRightBp() const;

    int getLength() const;

    std::string getAlternative() const;

    std::string getHomseq() const;

    std::string getGenotype() const;

    friend std::ostream& operator <<(std::ostream& stream, const Deletion& del);

private:
    int id;
    std::string referenceName;
    int leftBp;
    int rightBp;
    int length;
    std::string alternative;
    std::string homseq;
    std::string genotype;

    bool checkRep() const;

};

#endif /* _DELETION_H_ */
