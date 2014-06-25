#ifndef BAMSTATCALCULATOR_H
#define BAMSTATCALCULATOR_H

#include "api/BamReader.h"
#include <string>
#include <vector>

class BamStatCalculator
{
public:
    BamStatCalculator(const std::string& filename);
    virtual ~BamStatCalculator();

    int getInsertMean();
    int getInsertSd();

private:
    void loadInserts();
    int mean();
    int sd();

    BamTools::BamReader reader;
    std::vector<int> inserts;
    int insertMean;
    int insertSd;
};

#endif // BAMSTATCALCULATOR_H
