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

namespace Helper {
std::string getReferenceName(BamTools::BamReader& reader, int referenceId);
const int SVLEN_THRESHOLD = -50;

//std::set<std::string> forwardEClipNames;
//std::set<std::string> reverseBClipNames;
}

#endif // HELPER_H
