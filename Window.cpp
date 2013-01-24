#include "Window.h"

Window::Window() {}

Window::Window(const std::string& name,
               const std::string& chrom,
               unsigned int start,
               unsigned int end,
               int insertSize)
    : name(name),
      chrom(chrom),
      start(start),
      end(end),
      insertSize(insertSize) {}

Window::~Window() {}
