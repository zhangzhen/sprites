#include "BamStatCalculator.h"
#include "error.h"

#include <numeric>
#include <cmath>

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
    uint32_t cnt = 0;
    while (reader.GetNextAlignmentCore(al) && cnt <= 10000)
    {
        if (al.IsProperPair() && al.MatePosition > al.Position)
        {
            uint32_t insert = al.MatePosition + al.Length - al.Position;
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
    return sqrt(inner_product(inserts.begin(), inserts.end(), inserts.begin(), 0) / inserts.size() - m * m);
}
