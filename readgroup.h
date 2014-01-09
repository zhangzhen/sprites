#ifndef READGROUP_H
#define READGROUP_H

#include <string>

#include "library.h"

class ReadGroup
{
public:
    ReadGroup(const std::string& name, Library *lib);

    std::string getName() const;

    int getMeanOfInsertSize() const;
    int getSdOfInsertSize() const;

private:
    std::string name;
    Library *lib;
};

#endif // READGROUP_H
