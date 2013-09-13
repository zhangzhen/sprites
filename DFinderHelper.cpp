#include <sstream>
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
