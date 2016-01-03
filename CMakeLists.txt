cmake_minimum_required(VERSION 2.8)

project(sprites)

include_directories($ENV{BAMTOOLS_HOME}/include $ENV{HTSLIB_HOME})
link_directories($ENV{BAMTOOLS_HOME}/lib $ENV{HTSLIB_HOME})
add_definitions(-std=c++11)

add_executable(sprites main.cpp error.cpp Helper.cpp
Deletion.cpp Thirdparty/overlapper.cpp BamStatCalculator.cpp ClipReader.cpp clip.cpp FaidxWrapper.cpp range.cpp)
target_link_libraries(sprites libhts.a bamtools pthread)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -gdwarf-3 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-g -O2 -Wall")

message(${CMAKE_CXX_FLAGS_DEBUG})
