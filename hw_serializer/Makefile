CC = g++ --std=c++17

all: prog 

prog: main.o serializer.o
	$(CC) main.o serializer.o -o prog

serializer.o: serializer.cpp serializer.h
	$(CC) -c serializer.cpp

main.o: main.cpp
	$(CC) -c main.cpp

test: prog
	./prog

clean:
	rm *.o prog
