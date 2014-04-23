#include "Helper.h"

// Strip the leading directories and
// the last trailling suffix from a filename
std::string stripFilename(const std::string& filename)
{
    std::string out = stripDirectories(filename);
    return stripExtension(out);
}

// Remove a single file extension from the filename
std::string stripExtension(const std::string& filename)
{
    size_t suffixPos = filename.find_last_of('.');
    if(suffixPos == std::string::npos)
        return filename; // no suffix
    else
        return filename.substr(0, suffixPos);
}

// Strip the leadering directories from a filename
std::string stripDirectories(const std::string& filename)
{
    size_t lastDirPos = filename.find_last_of('/');

    if(lastDirPos == std::string::npos)
        return filename; // no directories
    else
        return filename.substr(lastDirPos + 1);
}

