#include "readgroup.h"

using namespace std;

ReadGroup::ReadGroup(const std::string &name, Library *lib) :
    name(name), lib(lib)
{
}

std::string ReadGroup::getName() const
{
    return name;
}

int ReadGroup::getMeanOfInsertSize() const
{
    return lib->getMeanOfInsertSize();
}

int ReadGroup::getSdOfInsertSize() const
{
    return lib->getSdOfInsertSize();
}
