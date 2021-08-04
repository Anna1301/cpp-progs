#include "vector.h"
#include <iostream>
#include <cassert>

void Test()
{
	Vector<int> vec;

	assert(vec.capacity() == 4);

	for(int i = 5; i < 15; i++)
	{
		vec.push_back(i);
	}

	assert(vec.size() == 10);

	vec.pop_back();

	assert(vec.size() == 9);

	Vector<int> vec2(vec);

	assert(vec2 == vec);

	vec.emplace_back(7);

	assert(vec[9] == 7);

	vec.resize(50);

	assert(vec.size() == 50);

	vec.clear();

	assert(vec.size() == 0);
}

int main()
{
	Test();

	std::cout << "success" << std::endl;

	return 0;
}
