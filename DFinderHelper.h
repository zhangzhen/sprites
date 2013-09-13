#ifndef _DFINDERHELPER_H_
#define _DFINDERHELPER_H_

#include <string>
#include <utility>

/*
 * Function: integerToString
 * Usage: string s = integerToString(n);
 * -------------------------------------
 * Converts an integer into the corresponding string of digits.
 * For example, calling <code>integerToString(123)</code> returns
 * the string <code>"123"</code>.
 */

std::string integerToString(int n);
bool isInconsistent(int a, int b);
bool contains(std::pair<int,int> a, std::pair<int,int> b);

#endif
