CC=g++
CXX=g++
CXXFLAGS=-std=c++11 -fopenmp -g -Wno-ignored-attributes -Wno-deprecated-declarations
LDLIBS=-fopenmp

all: a1show

main.o: main.cpp 

main: main.o

mdsimulation.dat: main 
	./main

a1show: mdsimulation.dat
	python a1.py

clean:
	rm -f main *.o *.dat *.pdf *.png out.txt
