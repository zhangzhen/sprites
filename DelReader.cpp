#include "DelReader.h"
#include <fstream>

using namespace std;

std::vector<Del> DelReader::readDelsFromFile(const std::string &filename)
{
    ifstream input(filename);

    std::vector<Del> dels;

    string referenceName;
    int leftBp;
    int rightBp;
    int length;

    while (input >> referenceName >> leftBp >> rightBp >> length) {
        Del del = {referenceName, leftBp, rightBp};
        dels.push_back(del);
    }
    return dels;
}
