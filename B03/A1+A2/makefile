CXXFLAGS=-std=c++11 -fopenmp -g
CC=g++
CXX=g++

all: a2show

rungkutt.o: rungkutt.cpp rungkutt.h

a2.o: a2.cpp rungkutt.h

a2: rungkutt.o a2.o

a2.dat: a2
	./a2

a2show: a2.dat
	python a2.py

clean:
	rm -f a2 *.o *.dat *.pdf
