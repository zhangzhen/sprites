#include "Helper.h"

using namespace std;

// Strip the leading directories and
// the last trailling suffix from a filename
string stripFilename(const string& filename) {
    string out = stripDirectories(filename);
    return stripExtension(out);
}

// Remove a single file extension from the filename
string stripExtension(const string& filename) {
    size_t suffixPos = filename.find_last_of('.');
    if(suffixPos == string::npos)
        return filename; // no suffix
    else
        return filename.substr(0, suffixPos);
}

// Strip the leadering directories from a filename
string stripDirectories(const string& filename) {
    size_t lastDirPos = filename.find_last_of('/');

    if(lastDirPos == string::npos)
        return filename; // no directories
    else
        return filename.substr(lastDirPos + 1);
}



std::string Helper::getReferenceName(BamTools::BamReader &reader, int referenceId) {
    assert(referenceId >= 0 && referenceId < reader.GetReferenceCount());
    return reader.GetReferenceData()[referenceId].RefName;
}


int numOfTheLongestPrefix(const string &s1, const string &s2)
{
    assert(s1.size() == s2.size());
    for (int i = 0; i < s1.size(); i++) {
        if (s1[i] != s2[i]) return i;
    }
    return 0;
}


int numOfThelongestSuffix(const string &s1, const string &s2)
{
    assert(s1.size() == s2.size());
    for (int i = 0; i < s1.size(); i++) {
        if (s1[s1.size() - 1 - i] != s2[s1.size() - 1 - i]) return i;
    }
    return 0;
}
