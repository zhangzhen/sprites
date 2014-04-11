#include "error.h"

using namespace std;

ErrorException::ErrorException(string msg) : msg(msg) {}

ErrorException::~ErrorException() throw() {}

string ErrorException::getMessage() {
    return msg;
}

const char *ErrorException::what() const throw () {
    return ("Error: " + msg).c_str();
}

void error(string str) {
    throw ErrorException(str);
}
