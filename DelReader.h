#ifndef DELREADER_H
#define DELREADER_H

#include <string>
#include <vector>

struct Del {
    std::string referenceName;
    int leftBp;
    int rightBp;

    int length() {
        return rightBp - leftBp;
    }
};

class DelReader
{
public:
    static std::vector<Del> readDelsFromFile(const std::string& filename);
};

#endif // DELREADER_H
