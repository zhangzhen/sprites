#ifndef FAIDXWRAPPER_H
#define FAIDXWRAPPER_H

#include "htslib/faidx.h"
#include <string>

class FaidxWrapper
{
public:
    FaidxWrapper(const std::string& fasta);
    virtual ~FaidxWrapper();
    int size();
    std::string fetch(const std::string& chrom, int start, int end);

private:
    faidx_t *fai;
};

#endif // FAIDXWRAPPER_H
