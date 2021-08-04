#include <iostream>
#include "serializer.h"
#include <cassert>
#include <sstream>

struct Data
{
	uint64_t a;
	bool b;
	uint64_t c;

	template <class Serializer>
	Error serialize(Serializer &A)
	{
		return A(a, b, c);
	}

	template <class Deserializer>
	Error deserialize(Deserializer &A)
	{
		return A(a, b, c);
	}
};

void RightTest()
{
	Data x { 1, true, 2 };

	std::stringstream str;

	Serializer serializer(str);
	serializer.save(x);

	Data y { 1, false, 0 };

	Deserializer deserializer(str);
	const Error err = deserializer.load(y);

	assert(err == Error::NoError);

	assert(x.a == y.a);
	assert(x.b == y.b);
	assert(x.c == y.c);
}

void WrongTest1()
{
	Data y { 1, false, 0 };

	std::stringstream str("100000000000000000000000000000000000000000000000000000000000000000000000 true 1");

	Deserializer deserializer(str);
	const Error err = deserializer.load(y);

	assert(err == Error::CorruptedArchive);
}

void WrongTest2()
{
	Data y { 1, false, 0 };

	std::stringstream str("-1 true 1"); 

	Deserializer deserializer(str);
	const Error err = deserializer.load(y);

	assert(err == Error::CorruptedArchive);
}

void WrongTest3()
{
        Data y { 1, false, 0 };

	std::stringstream str("1err true 1");

        Deserializer deserializer(str);
        const Error err = deserializer.load(y);

	assert(err == Error::CorruptedArchive);
}

int main()
{
	RightTest();
	WrongTest1();
	WrongTest2();
	WrongTest3();

	std::cout << "success" << std::endl;

	return 0;
}
