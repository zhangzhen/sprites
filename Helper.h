#ifndef HELPER_H
#define HELPER_H

#include "api/BamReader.h"
#include <string>
#include <set>

//
// Functions
//
std::string stripFilename(const std::string& filename);
std::string stripExtension(const std::string& filename);
std::string stripDirectories(const std::string& filename);
int getOffsetForward(const std::string& s1, const std::string& s2);
int getOffsetReverse(const std::string& s1, const std::string& s2);

template<class T, class Compare>
void cluster(const std::vector<T>& orig, std::vector<std::vector<T> >& clusters, Compare comp) {
    std::vector<T> buffer;

    auto first = orig.begin();
    auto last = orig.end();
    buffer.push_back(*first);
    while (++first != last) {
        if (!comp(*first, buffer[0])) {
            clusters.push_back(buffer);
            buffer.clear();
        }
        buffer.push_back(*first);
    }
    if (!buffer.empty()) clusters.push_back(buffer);
}

namespace Helper {
std::string getReferenceName(BamTools::BamReader& reader, int referenceId);
const int SVLEN_THRESHOLD = -50;
const int CONFLICT_THRESHOLD = 13;

//std::set<std::string> forwardEClipNames;
//std::set<std::string> reverseBClipNames;
}

#endif // HELPER_H
