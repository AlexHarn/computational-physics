CC=g++
CXX=g++
CXXFLAGS=-std=c++11 -fopenmp -g -Wno-ignored-attributes -Wno-deprecated-declarations $(shell pkg-config eigen3 --cflags) -I "C:\Program Files\cpplibraries\eigen3"
LDLIBS=-fopenmp

all: plot 

main.o: main.cpp

main: main.o

run: main
	./main

plot: run
	python plot.py

clean:
	rm -f  main *.o *.dat *.pdf out.txt
