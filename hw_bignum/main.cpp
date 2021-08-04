#include "bignum.h"
#include <iostream>
#include <cassert>

void Test()
{
	std::string s1 = "11111111111111111111111111111111111111111111111111111111111111111111111111111111111";
	BigNum num1(s1);

	std::string s2 = "22222222222222222222222222222222222222222222222222222222222222222222222222222222222";
	BigNum num2(s2);

	std::string s3 = "33333333333333333333333333333333333333333333333333333333333333333333333333333333333";

	assert(num1 < num2);

	BigNum num3(s3); 

	num1 = num1 + num2;

	assert(num3 == num1);

	num3 = num1 - num2;

	BigNum num4(s1);

	assert(num3 == num4);

	s1 = "1000000000000000000000000";
	s2 = "1111111111111111111111111";
	BigNum num5(s1);
	BigNum num6(s2);
	s3 = "1111111111111111111111111000000000000000000000000";
	BigNum num7(s3);

	assert(num7 == num5 * num6);

	s1 = "-100";
	s2 = "100";

	BigNum num8(s1);
	BigNum num9(s2);
	BigNum num10("0");

	assert(num10 == num8 + num9);
}

int main()
{
	Test();

	std::cout << "success" << std::endl;

	return 0;
}
