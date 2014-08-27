#include "DeletionReader.h"
#include <fstream>
#include <iostream>

using namespace std;

std::vector<Deletion> DeletionReader::readDelsFromFile(const std::string &filename)
{
    ifstream input(filename);

    std::vector<Deletion> dels;

    string referenceName;
    int leftBp;
    int rightBp;
    int length;
    string svtype;
    string alternative;
    string homseq;
    string genotype;

    string head;
    getline(input, head);

    while (input >> referenceName >> leftBp >> rightBp >> svtype >> length >> alternative >> homseq >> genotype) {
        Deletion del (referenceName, leftBp, rightBp, length,
                   alternative, homseq, genotype);
        dels.push_back(del);
    }

    return dels;
}
