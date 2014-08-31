cmake_minimum_required(VERSION 2.8)

project(dfinder)

include_directories($ENV{BAMTOOLS_HOME}/include $ENV{SAMTOOLS_HOME})
link_directories($ENV{BAMTOOLS_HOME}/lib $ENV{SAMTOOLS_HOME})
add_definitions(-std=c++0x)

add_executable(dfinder main.cpp error.cpp Helper.cpp
Deletion.cpp Thirdparty/overlapper.cpp Thirdparty/multiple_alignment.cpp
BamStatCalculator.cpp ClipReader.cpp clip.cpp FaidxWrapper.cpp)
target_link_libraries(dfinder bam bamtools z)

add_executable(countReads countReads.cpp Deletion.cpp DeletionReader.cpp SoftClip.cpp SoftClipReader.cpp
SoftClipCounter.cpp error.cpp)
target_link_libraries(countReads bamtools)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")
