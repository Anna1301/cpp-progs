#include <iostream>
#include "matrix.h"
#include <cassert>
#include <cmath>
#include <algorithm>
#include "test.h"

int main()
{
	CreationTest();
	EqualityTest();
	AssignmentTest();
	MultiplicationTest();

	std::cout << "success" << std::endl;

	return 0;

}
