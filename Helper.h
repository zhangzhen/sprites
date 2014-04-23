#ifndef HELPER_H
#define HELPER_H

#include <string>

//
// Functions
//
std::string stripFilename(const std::string& filename);
std::string stripExtension(const std::string& filename);
std::string stripDirectories(const std::string& filename);

#endif // HELPER_H
