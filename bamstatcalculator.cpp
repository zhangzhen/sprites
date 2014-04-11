#include "bamstatcalculator.h"
#include "error.h"

using namespace std;
using namespace BamTools;

BamStatCalculator::BamStatCalculator(const string &filename)
{
    if (reader.Open(filename))
        error("Could not open the input BAM file.");
    loadInserts();
    insert = BasicStat::mean(inserts);
    std = BasicStat::std(inserts);
}

BamStatCalculator::~BamStatCalculator()
{
    reader.Close();
}

int BamStatCalculator::getInsert() const
{
    return insert;
}

int BamStatCalculator::getStd() const
{
    return std;
}

void BamStatCalculator::loadInserts()
{
    BamAlignment al;
    int cnt = 0;
    while (reader.GetNextAlignmentCore(al) || cnt <= 10000)
    {
        if (al.IsProperPair() && al.MatePosition > al.Position)
        {
            inserts.push_back(al.InsertSize);
            cnt++;
        }
    }
}


double BasicStat::mean(const std::vector<int> &v)
{
    return accumulate(v.begin(), v.end(), 0.0) / v.size();
}


double BasicStat::std(const std::vector<int> &v)
{
    double m =  mean(v);
    return sqrt(inner_product(v.begin(), v.end(), v.begin(), 0.0) / v.size() - m * m);
}
