#include <iostream>
#include "format.h"
#include <cassert>
#include <stdexcept>

void BasicTest()
{
	try
	{
		std::string s = format("{0} + {0} = {1}", "one", 2);
		assert(s == "one + one = 2");
	}
	catch(...)
	{
		std::cout << "Error" << std::endl;
	}
}

void FirstTest()
{
	bool IsError = false;
	try
	{
		std::string s = format("{0} + {0} = {2}", "one", 2);   
	}
	catch(const std::exception& err)
	{                        	
		IsError = true;
		std::cerr << "First test error: " << err.what() << std::endl;
	}
	catch(...)  
	{
		std::cout << "Unknown error" << std::endl;
	}
	assert(IsError == true);
}


void SecondTest()
{
	bool IsError = false;  
	try
	{
		std::string s = format("{0} + {0} = {}", "one", 2);
	}
	catch(const std::exception& err)
	{
		IsError = true;   
		std::cerr << "Second test error: " << err.what() << std::endl;
	}
	catch(...)
	{
		std::cout << "Unknown error" << std::endl;
	}
	assert(IsError == true);
}

void ThirdTest()
{
	bool IsError = false;
	try
	{
		std::string s = format("{0} { + {0} = {1}", "one", 2);
	}
	catch(const std::exception& err)
	{
		IsError = true;
		std::cerr << "Third test error: " << err.what() << std::endl;
	}
	catch(...)
	{
		std::cout << "Unknown error" << std::endl;
	}
	assert(IsError == true);
}


int main()
{
	BasicTest();
	FirstTest();
	SecondTest();
	ThirdTest();

	std::cout << "success" << std::endl;

	return 0;
}

