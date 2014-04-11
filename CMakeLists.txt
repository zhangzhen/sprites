cmake_minimum_required(VERSION 2.8)

project(dfinder)

include_directories($ENV{BAMTOOLS_HOME}/include)
link_directories($ENV{BAMTOOLS_HOME}/lib)
add_definitions(-std=c++0x)

add_executable(dfinder main.cpp DFinder.cpp SoftClip.cpp error.cpp DFinderHelper.cpp ChrRegion.cpp
ChrRegionCluster.cpp Overlap.cpp Deletion.cpp Thirdparty/overlapper.cpp Thirdparty/multiple_alignment.cpp
bamstatcalculator.cpp caller.cpp softclipreader.cpp)
target_link_libraries(dfinder bamtools "-lz")
