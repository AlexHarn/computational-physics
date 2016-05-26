CC=g++
CXX=g++
CXXFLAGS=-std=c++11 -fopenmp -g -Wno-ignored-attributes -Wno-deprecated-declarations
LDLIBS=-fopenmp

all: a1show

initializeSquareGrid.o: initializeSquareGrid.cpp initializeSquareGrid.h

mdsimulation.o: mdsimulation.cpp mdsimulation.h initializeSquareGrid.h ljforces.h

ljforces.o: ljforces.cpp ljforces.h periodRB.h

periodRB.o: periodRB.cpp periodRB.h

mdsimulation: mdsimulation.o initializeSquareGrid.o ljforces.o periodRB.o

mdsimulation.dat: mdsimulation
	./mdsimulation

a1show: mdsimulation.dat
	#python a1.py

clean:
	rm -f mdsimulation *.o *.dat *.pdf *.png out.txt
