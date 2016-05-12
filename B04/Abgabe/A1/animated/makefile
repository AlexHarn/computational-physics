CC=g++
CXX=g++
CXXFLAGS=-std=c++11 -fopenmp -g
LDLIBS=-fopenmp -lsfml-graphics -lsfml-window -lsfml-system

all: run

rungkutt.o: rungkutt.cpp rungkutt.h

main.o: main.cpp Dpendulum.h

Dpendulum.o: Dpendulum.cpp Dpendulum.h rungkutt.h

main: main.o Dpendulum.o rungkutt.o

run: main
	./main

clean:
	rm -f main *.o
