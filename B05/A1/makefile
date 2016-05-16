CC=g++
CXX=g++
CXXFLAGS=-std=c++11 -fopenmp -g
LDLIBS=-fopenmp

all: a1show

initializeSquareGrid.o: initializeSquareGrid.cpp

initializeSquareGrid: initializeSquareGrid.o

a1.dat: initializeSquareGrid
	./initializeSquareGrid

a1show: a1.dat
	#python a1.py

clean:
	rm -f a1 *.o *.dat *.pdf *.png out.txt
