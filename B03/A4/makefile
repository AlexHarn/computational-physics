CXXFLAGS=-std=c++11 -fopenmp -g
LDLIBS=-fopenmp
CC=g++
CXX=g++

all: a4show

rungkutt.o: rungkutt.cpp rungkutt.h

kepler.o: kepler.cpp kepler.h rungkutt.h

a4.o: a4.cpp kepler.h

a4: a4.o rungkutt.o kepler.o

a4.dat: a4
	./a4

a4show: a4.dat
	python a4.py

clean:
	rm -f a4 *.o *.dat *.dat1 *.pdf *.png out.txt
