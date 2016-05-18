CC=g++
CXX=g++
CXXFLAGS=-std=c++11 -fopenmp -g
LDLIBS=-fopenmp

all: a1show

initializeSquareGrid.o: initializeSquareGrid.cpp initializeSquareGrid.h

mdsimulation.o: mdsimulation.cpp mdsimulation.h initializeSquareGrid.h

#Dpendulum.o: Dpendulum.cpp Dpendulum.h rungkutt.h

mdsimulation: mdsimulation.o initializeSquareGrid.o

mdsimulation.dat: mdsimulation
	./mdsimulation

a1show: mdsimulation.dat
	#python a1.py

clean:
	rm -f a1 *.o *.dat *.pdf *.png out.txt
