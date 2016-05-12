CC=g++
CXX=g++
CXXFLAGS=-std=c++11 -fopenmp -g
LDLIBS=-fopenmp

all: a1show

a1.o: a1.cpp

a1: a1.o

a1.dat: a1
	./a1

a1show: a1.dat
	python a1.py

clean:
	rm -f a1 *.o *.dat *.pdf *.png out.txt
