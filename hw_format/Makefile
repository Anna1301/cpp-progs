CC = --std=c++17

prog: main.o
	g++ main.o -o prog $(CC)
main.o: main.cpp
	g++ -c main.cpp $(CC)
format.o: format.h
	g++ -c format.h $(CC)
clean:
	rm *.o
test:
	./prog