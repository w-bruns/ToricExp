#
# Makefile for ToricExp
#

CXX = gcc
CXXFLAGS = -O3


all: RandTriang Support ChPoly Triangulate

%: %.c
	$(CXX) $(CXXFLAGS) $< -o $@

clean: 
	rm -f *.o RandTriang Support ChPoly Triangulate
