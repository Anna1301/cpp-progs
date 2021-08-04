CC = g++ --std=c++17

all: prog 

prog: test.o alctr.o
	$(CC) test.o alctr.o -o prog

alctr.o: alctr.cpp
	$(CC) -c alctr.cpp

test.o: test.cpp
	$(CC) -c test.cpp

test: prog
	./prog

clean:
	rm *.o prog
