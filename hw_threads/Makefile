prog: main.o
	g++ main.o -o prog -pthread --std=c++17
main.o: main.cpp tpool.h
	g++ -c main.cpp --std=c++17
clean:
	rm *.o
test:
	./prog