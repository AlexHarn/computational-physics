CC=g++
CXX=g++
CXXFLAGS=-std=c++11 -fopenmp -g -Wno-ignored-attributes -Wno-deprecated-declarations
LDLIBS=-fopenmp

all: run

main.o: main.cpp

main: main.o

run: main
	./main
	python plot.py

clean:
	rm -f  main *.o *.dat *.pdf out.txt
