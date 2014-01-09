#include "statisticscalculator.h"
#include "api/BamAlignment.h"
#include "DFinderHelper.h"
#include <algorithm>

using namespace std;
using namespace BamTools;

StatisticsCalculator::StatisticsCalculator(string filename, Library *lib) :
    lib(lib), mean(0), sd(0)
{
    reader.Open(filename);
}

void StatisticsCalculator::calcMean()
{
    if (insertsizes.empty()) load();
    if (!mean) mean = round(::mean(insertsizes));
    lib->setMeanOfInsertSize(mean);
}

void StatisticsCalculator::calcSd()
{
    if (insertsizes.empty()) load();
    if (!sd) sd = round(::sd(insertsizes));
    lib->setSdOfInsertSize(sd);
}

void StatisticsCalculator::load()
{
    vector<int> res;

    BamAlignment al;
    string rg;
    while (reader.GetNextAlignment(al))
    {
        if (!al.GetTag("RG", rg)) continue;
        if (al.InsertSize > 0 && lib->contains(rg))
        {
            res.push_back(al.InsertSize);
        }
    }

    int mn = ::quantile(res, 0.05);
    int mx = ::quantile(res, 0.95);

    copy_if(res.begin(), res.end(), back_inserter(insertsizes), [mn, mx](int i) { return i >= mn && i <= mx; });
}

StatisticsCalculator::~StatisticsCalculator()
{
    reader.Close();
}
