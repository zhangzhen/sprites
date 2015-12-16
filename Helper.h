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
int numOfTheLongestPrefix(const std::string& s1, const std::string& s2);
int numOfThelongestSuffix(const std::string& s1, const std::string& s2);

int extend(const std::string& read, int offset, int leftOrigin, int rightOrigin);

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

template<class T, class Compare>
void merge(const std::vector<T>& orig, std::vector<T>& results, Compare comp) {
    std::vector<bool> removed(orig.size(), false);

    for (size_t i = 0; i < orig.size() - 1; ++i) {
        for (size_t j = i + 1; j < orig.size(); ++j) {
            if (!removed[j] && comp(orig[i], orig[j])) {
                removed[j] = true;
            }
        }
    }

    for (size_t i = 0; i < removed.size(); ++i) {
        if (!removed[i]) {
            results.push_back(orig[i]);
        }
    }

}

namespace Helper {
std::string getReferenceName(BamTools::BamReader& reader, int referenceId);
const int SVLEN_THRESHOLD = -50;
const int CONFLICT_THRESHOLD = 13;

//std::set<std::string> forwardEClipNames;
//std::set<std::string> reverseBClipNames;
}

#endif // HELPER_H
