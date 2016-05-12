CC=g++
CXX=g++
CXXFLAGS=-std=c++11 -fopenmp -g
LDLIBS=-fopenmp

all: a2abgabeshow

a2abgabe.o: a2abgabe.cpp

a2abgabe: a2abgabe.o

a2abgabe.dat: a2abgabe
	./a2abgabe

a2abgabeshow: a2abgabe.dat
	python a2.py

clean:
	rm -f a2 *.o *.dat *.pdf *.png out.txt
