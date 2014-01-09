#ifndef _DFINDER_H_
#define _DFINDER_H_

#include <iostream>
#include <iterator>
#include <map>
#include <utility>
#include <unordered_set>

#include "SoftClip.h"
#include "Interval.h"
#include "ChrRegion.h"
#include "ChrRegionCluster.h"
#include "api/BamReader.h"
#include "DFinderHelper.h"
#include "library.h"


template<typename T>
class DeletePtr {
public:
    void operator() (T* ptr) {
        delete ptr;
    }
};

class DFinder
{
public:
    DFinder(const std::string& filename, const std::string& referenceName, int minOverlapLength, double maxMismatchRate);
    ~DFinder();
    void callToFile(const std::string& filename);

private:
    void call(std::vector<Deletion>& calls);
    void getSpanningPairsFor(const SoftClip& softclip, std::vector<SoftClip>& results);

    bool callFor(const SoftClip& softclip, Deletion& del);
    bool callFor(const SoftClip& softclip, const SoftClip& spanner, Deletion& del);
    bool callFor(const SoftClip &softclip, const std::vector<SoftClip> &candidates, Deletion &del);

    bool isValidAlignment(const BamTools::BamAlignment& al);
    bool isCorrectOrientation(const BamTools::BamAlignment& al);
    bool isValidPartnerCandidateFor(const SoftClip& orig, const SoftClip &partner);

    bool findReferenceId(const std::string& name, int& id);

    std::string referenceName;
    int referenceId;

    int minOverlapLength;
    double maxMismatchRate;

    static const int lengthThreshold = 50;
    static const int MapQualityThreshold = 1;

    BamTools::BamReader r1;
    BamTools::BamReader r2;

    std::unordered_set<std::string> readNames;

    std::map<std::string, Library*> libraries;
    std::map<std::string, ReadGroup*> readgroups;

    void findPartnerCandidates(const SoftClip &softclip, int start, int end, std::vector<SoftClip> &softclips, std::vector<SoftClip> &reads);
};

#endif /* _DFINDER_H_ */
