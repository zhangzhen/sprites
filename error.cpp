#include "error.h"

using namespace std;

ErrorException::ErrorException(std::string msg) : msg(msg) {}

ErrorException::~ErrorException() throw() {}

string ErrorException::getMessage() const {
    return msg;
}

const char *ErrorException::what() const throw () {
    return ("Error: " + msg).c_str();
}

void error(std::string str) {
    throw ErrorException(str);
}
