# File: Makefile
# Author: Keith Schwarz (htiek@cs.stanford.edu)
#
# A Makefile for building the CS106L GraphViz program.

BAMTOOLS_ROOT = $(HOME)/bamtools
GTEST_ROOT = $(HOME)/gtest-1.6.0

# Disable optimization; turn on debugging.  Feel free to change this to
#
CCFLAGS = -O3
#
# If you want to turn on optimization once things get working.
# CCFLAGS = -g -O0

# Builds the main program with the necessary libraries.
all: clip-sv AllTests

clip-sv: main.o clip-sv.o Clip.o
	g++ main.o clip-sv.o Clip.o -o clip-sv -L$(BAMTOOLS_ROOT)/lib -lbamtools $(CCFLAGS)

# Build object files from sources.
%.o: %.cpp
	g++ $^ -c -o $@ -I$(BAMTOOLS_ROOT)/include -I$(GTEST_ROOT)/include $(CCFLAGS)

AllTests: AllTests.o SingleClipTest.o clip-sv.o Clip.o
	g++ AllTests.o SingleClipTest.o clip-sv.o Clip.o -o AllTests -L$(GTEST_ROOT)/lib -L$(BAMTOOLS_ROOT)/lib -lgtest -lbamtools -lpthread

# Cleans the project by nuking emacs temporary files (*~), object files (*.o),
# and the resulting executable.
clean:
	rm -rf *~ *.o clip-sv AllTests
