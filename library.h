#ifndef LIBRARY_H
#define LIBRARY_H

#include <unordered_set>
#include <vector>
#include <string>

class Library
{
public:
    Library(const std::string& name);

    std::string getName();

    void add(const std::string& readgroup);

    bool contains(const std::string& readgroup);

    int getMeanOfInsertSize() const;
    void setMeanOfInsertSize(int value);

    int getSdOfInsertSize() const;
    void setSdOfInsertSize(int value);

    int getReadLength() const;
    void setReadLength(int value);

private:
    std::string name;
    int meanOfInsertSize;
    int sdOfInsertSize;
    int readLength;
    std::unordered_set<std::string> readgroups;
};

#endif // LIBRARY_H
