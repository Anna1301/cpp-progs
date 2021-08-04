CC = g++ --std=c++17

all: prog 

prog: main.o test.o matrix.o
	$(CC) main.o test.o matrix.o -o prog

test.o: test.cpp test.h
	$(CC) -c test.cpp

matrix.o: matrix.cpp matrix.h
	$(CC) -c matrix.cpp

main.o: main.cpp
	$(CC) -c main.cpp

test: prog
	./prog

clean:
	rm *.o prog
