#include "DeletionReader.h"
#include <fstream>
#include <iostream>

using namespace std;

std::vector<Deletion> DeletionReader::readDelsFromFile(const std::string &filename)
{
    ifstream input(filename);

    std::vector<Deletion> dels;

    string referenceName;
    string referenceName2;
    int start1;
    int end1;
    int start2;
    int end2;
    int length;

    string head;
    getline(input, head);

    while (input >> referenceName >> start1 >> end1 >> referenceName2 >> start2 >> end2 >> length) {
        Deletion del (referenceName, start1, end1, start2,
                   end2, length);
        dels.push_back(del);
    }

    return dels;
}
