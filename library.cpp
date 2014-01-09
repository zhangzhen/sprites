#include "library.h"
using namespace std;


Library::Library(const string& name) :
    name(name), meanOfInsertSize(0), sdOfInsertSize(0)
{
}

string Library::getName()
{
    return name;
}

void Library::add(const string &readgroup)
{
    readgroups.insert(readgroup);
}

bool Library::contains(const string &readgroup)
{
    return readgroups.count(readgroup);
}

int Library::getMeanOfInsertSize() const
{
    return meanOfInsertSize;
}

void Library::setMeanOfInsertSize(int value)
{
    meanOfInsertSize = value;
}

int Library::getSdOfInsertSize() const
{
    return sdOfInsertSize;
}

void Library::setSdOfInsertSize(int value)
{
    sdOfInsertSize = value;
}

int Library::getReadLength() const
{
    return readLength;
}

void Library::setReadLength(int value)
{
    readLength = value;
}




