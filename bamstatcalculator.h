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

    int getInsert() const;
    int getStd() const;

private:
    void loadInserts();

    BamTools::BamReader reader;
    std::vector<int> inserts;
    int insert;
    int std;
};

namespace BasicStat {

double mean(const std::vector<int>& v);

double std(const std::vector<int>& v);

}

#endif // BAMSTATCALCULATOR_H
