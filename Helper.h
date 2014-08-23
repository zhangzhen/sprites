#ifndef HELPER_H
#define HELPER_H

#include "api/BamReader.h"
#include <string>

//
// Functions
//
std::string stripFilename(const std::string& filename);
std::string stripExtension(const std::string& filename);
std::string stripDirectories(const std::string& filename);

namespace Helper {
std::string getReferenceName(BamTools::BamReader& reader, int referenceId);
}

#endif // HELPER_H
