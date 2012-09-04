#include "error.h"

ErrorException::ErrorException(std::string msg) : msg(msg) {}

ErrorException::~ErrorException() throw() {}

std::string ErrorException::getMessage() {
  return msg;
}

void error(std::string str) {
  throw ErrorException(str);
}
