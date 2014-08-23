#ifndef FAIDXWRAPPER_H
#define FAIDXWRAPPER_H

#include "faidx.h"
#include <string>

class FaidxWrapper
{
public:
    FaidxWrapper(const std::string& fasta);
    virtual ~FaidxWrapper();
    size_t size();
    std::string fetch(const std::string& chrom, int start, int end);

private:
    faidx_t *fai;
};

#endif // FAIDXWRAPPER_H
