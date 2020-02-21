CC = g++
.DEFAULT_GOAL := all

tsp.o: tsp.cpp matrix.hpp
	$(CC) -c tsp.cpp

matrix.o: matrix.cpp matrix.hpp
	$(CC) -c matrix.cpp

ga_frame.o: ga_frame.cpp ga_frame.hpp
	$(CC) -c ga_frame.cpp

tsp: tsp.o matrix.o ga_frame.o
	$(CC) -g -o tsp tsp.o matrix.o ga_frame.o

all: tsp

clean:
	rm tsp tsp.o matrix.o ga_frame.o solution.csv