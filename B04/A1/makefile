CC=g++
CXX=g++
CXXFLAGS=-std=c++11 -fopenmp -g
LDLIBS=-fopenmp

all: a1show

a1abgabe.o: a1abgabe.cpp

a1abgabe: a1abgabe.o

a1abgabe.dat: a1abgabe
	./a1abgabe

a1show: a1abgabe.dat
	python a1.py

clean:
	rm -f a1 *.o *.dat *.pdf *.png out.txt
