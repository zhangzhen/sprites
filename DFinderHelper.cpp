#include <sstream>
#include <cassert>
#include <cmath>
#include "DFinderHelper.h"
using namespace std;

string integerToString(int n) {
    ostringstream stream;
    stream << n;
    return stream.str();
}

bool isInconsistent(int a, int b) {
    if (a < b) return double(b) / a > 2;
    return double(a) / b > 2;
}

bool contains(std::pair<int,int> a, std::pair<int,int> b) {
    return (a.first <= b.first && a.second >= b.second) ||
	(b.first <= a.first && b.second >= a.second);
}

int numOfMismatches(const std::string& s1, const std::string& s2) {
    int cnt = 0;
    assert(s1.length() == s2.length());
    for (int i = 0; i < s1.length(); ++i) {
	if (s1[i] != 'N' && s2[i] != 'N' && s1[i] == s2[i])
	    cnt++;
    }
    return cnt;
}

bool overlaps2(const std::string& s1, const std::string& s2, double maxMismatchRate, int& lengthOfOverlap, int& mismatches) {
    int initLength = std::min(s1.length(), s2.length());
    for (int i = initLength; i > 0; --i) {
	int k = ceil(i * maxMismatchRate);
	mismatches = numOfMismatches(s1.substr(s1.length() - i), s2.substr(0, i));
	if (mismatches <= k) {
	    lengthOfOverlap = i;
	    return true;
	}
    }
    return false;
}
