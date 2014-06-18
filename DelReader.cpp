#include "DelReader.h"
#include <fstream>
#include <iostream>

using namespace std;

std::vector<Del> DelReader::readDelsFromFile(const std::string &filename)
{
    ifstream input(filename);

    std::vector<Del> dels;

    int id;
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

    while (input >> id >> referenceName >> leftBp >> rightBp >> svtype >> length >> alternative >> homseq >> genotype) {
        Del del = {id, referenceName, leftBp, rightBp, length,
                   alternative, homseq == "-" ? "" : homseq, genotype};
        dels.push_back(del);
    }

    return dels;
}
