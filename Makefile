# File: Makefile
# Author: Keith Schwarz (htiek@cs.stanford.edu)
#
# A Makefile for building the CS106L GraphViz program.

BAMTOOLS_ROOT = $(HOME)/bamtools
GTEST_ROOT = $(HOME)/gtest-1.6.0

# Disable optimization; turn on debugging.  Feel free to change this to
#
CCFLAGS = -O3 -std=c++0x
#
# If you want to turn on optimization once things get working.
# CCFLAGS = -g -O0

LIBS = -lz -L$(BAMTOOLS_ROOT)/lib
TESTLIBS = -L$(GTEST_ROOT)/lib -lgtest -lpthread $(LIBS)

INCLUDE = -I$(BAMTOOLS_ROOT)/include -I$(GTEST_ROOT)/include

OBJS = DFinder.o SoftClip.o error.o Interval.o IntervalCluster.o Overlap.o Deletion.o $(BAMTOOLS_ROOT)/lib/libbamtools.a
# TESTOBJS = SingleClipTest.o CallDelsTest.o SingleClippedClusterTest.o $(OBJS)

# Builds the main program with the necessary libraries.
all: dfinder #AllTests

dfinder: main.o $(OBJS)
	g++ main.o $(OBJS) -o dfinder $(CCFLAGS) $(LIBS)

# Build object files from sources.
%.o: %.cpp
	g++ $^ -c -o $@ $(CCFLAGS) $(INCLUDE)

# AllTests: AllTests.o $(TESTOBJS)
# 	g++ AllTests.o $(TESTOBJS) -o AllTests $(TESTLIBS)

# Cleans the project by nuking emacs temporary files (*~), object files (*.o),
# and the resulting executable.
clean:
	rm -rf *~ *.o dfinder #AllTests
