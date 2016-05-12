CC=g++
CXX=g++
CXXFLAGS=-std=c++11 -fopenmp -g
LDLIBS=-fopenmp

all: a2show

a2.o: a2.cpp

a2: a2.o

a2.dat: a2
	./a2

a2show: a2.dat
	python a2.py

clean:
	rm -f a2 *.o *.dat *.pdf *.png out.txt
