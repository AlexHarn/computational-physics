CC=g++
CXX=g++
CXXFLAGS=-std=c++11 -fopenmp -g -Wno-ignored-attributes -Wno-deprecated-declarations
LDLIBS=-fopenmp

all: run

poisson.o: poisson.cpp poisson.h

main.o: main.cpp poisson.h

main: main.o poisson.o

run: main
	./main
	python plot.py

clean:
	rm -f  main *.o *.dat *.pdf out.txt
