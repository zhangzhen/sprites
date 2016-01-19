cmake_minimum_required(VERSION 2.8)

project(sprites)

include_directories($ENV{BAMTOOLS_HOME}/include $ENV{HTSLIB_HOME})
#link_directories($ENV{BAMTOOLS_HOME}/lib $ENV{HTSLIB_HOME})
add_definitions(-std=c++0x)

add_executable(sprites main.cpp error.cpp Helper.cpp
Deletion.cpp Thirdparty/overlapper.cpp BamStatCalculator.cpp ClipReader.cpp clip.cpp FaidxWrapper.cpp range.cpp)
target_link_libraries(sprites $ENV{HTSLIB_HOME}/libhts.a $ENV{BAMTOOLS_HOME}/lib/libbamtools.a pthread z)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O2 -Wall")

if(CMAKE_BUILD_TYPE MATCHES DEBUG)
    message(${CMAKE_CXX_FLAGS_DEBUG})
else(CMAKE_BUILD_TYPE MATCHES DEBUG)
    message(${CMAKE_CXX_FLAGS_RELEASE})
endif(CMAKE_BUILD_TYPE MATCHES DEBUG)

