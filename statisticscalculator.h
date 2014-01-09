#ifndef STATISTICSCALCULATOR_H
#define STATISTICSCALCULATOR_H

#include "api/BamReader.h"
#include "library.h"
#include <string>
#include <vector>

class StatisticsCalculator
{
public:
    StatisticsCalculator(std::string filename, Library *lib);
    void calcMean();
    void calcSd();

    ~StatisticsCalculator();
private:
    int mean;
    int sd;
    Library *lib;
    BamTools::BamReader reader;
    std::vector<int> insertsizes;

    void load();
};

#endif // STATISTICSCALCULATOR_H
