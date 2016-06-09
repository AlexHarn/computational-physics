CC=g++
CXX=g++
CXXFLAGS=-std=c++11 -fopenmp -g -Wno-ignored-attributes -Wno-deprecated-declarations #-I "C:\Program Files\cpplibraries" # -I Teil nur für Sebastians Windows einkommentieren
LDLIBS=-fopenmp

all: run

lingen.o: lingen.cpp lingen.h

main.o: main.cpp lingen.h

main: main.o lingen.o

run: main
	./main
	python plot.py

clean:
	rm -f  main *.o *.dat *.pdf out.txt
