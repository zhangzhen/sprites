#include "FaidxWrapper.h"
#include "error.h"
#include <algorithm>

using namespace std;

FaidxWrapper::FaidxWrapper(const std::string &fasta)
{
    fai = fai_load(fasta.c_str());
    if (fai == NULL) error("Cannot load the indexed fasta.");
}

FaidxWrapper::~FaidxWrapper()
{
    if (fai != NULL) fai_destroy(fai);
}

size_t FaidxWrapper::size()
{
    return faidx_fetch_nseq(fai);
}

string FaidxWrapper::fetch(const string &chrom, int start, int end)
{
    int len;
    char *s = faidx_fetch_seq(fai, (char *)chrom.c_str(), start, end, &len);
    string str(s);
    transform(str.begin(), str.end(),str.begin(), ::toupper);
    return str;
}
