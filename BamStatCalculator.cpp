#include "BamStatCalculator.h"
#include "error.h"

#include <numeric>
#include <cmath>
#include <algorithm>

using namespace std;
using namespace BamTools;

BamStatCalculator::BamStatCalculator(const string &filename) :
    insertMean(-1), insertSd(-1)
{
    if (!reader.Open(filename))
        error("Could not open the input BAM file.");
    loadInserts();
}

BamStatCalculator::~BamStatCalculator()
{
    reader.Close();
}

int BamStatCalculator::getInsertMean()
{
    if (insertMean == -1) {
        insertMean = mean();
    }
    return insertMean;
}

int BamStatCalculator::getInsertSd()
{
    if (insertSd == -1) {
        insertSd = sd();
    }
    return insertSd;
}

void BamStatCalculator::loadInserts()
{
    BamAlignment al;
    size_t cnt = 0;
    while (reader.GetNextAlignmentCore(al) && cnt < 10000)
    {
        if (al.IsProperPair() && al.MatePosition > al.Position)
        {
            uint64_t insert = al.MatePosition + al.Length - al.Position;
            if (insert < 10000) {
                inserts.push_back(insert);
                cnt++;
            }
        }
    }
}

int BamStatCalculator::mean()
{
    return accumulate(inserts.begin(), inserts.end(), 0) / inserts.size();
}


int BamStatCalculator::sd()
{
    int m =  getInsertMean();
    vector<int> temp;
    transform(inserts.begin(), inserts.end(), back_inserter(temp), [](int x) { return x*x; });
    uint32_t sum = accumulate(temp.begin(), temp.end(), 0);
    return sqrt( sum / temp.size() - m * m);
}
