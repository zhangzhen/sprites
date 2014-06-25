#ifndef DELREADER_H
#define DELREADER_H

#include "Deletion.h"
#include <string>
#include <vector>
#include <set>


class DeletionReader
{
public:
    static std::vector<Deletion> readDelsFromFile(const std::string& filename);
};

#endif // DELREADER_H
