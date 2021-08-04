#include <vector>
#include <iostream>

int main()
{
	std::vector<int> vec;
	vec.push_back(2);
	vec.push_back(2);
	std::cout << vec.size() << "  " << vec.capacity() << std::endl;
	for(int i = 0; i < vec.size(); ++i)
	{
		std::cout << vec[i] << " "; 
	}
	std::cout << std::endl;
	vec.emplace_back(5);
	for(int i = 0; i < vec.size(); ++i)
	{
		std::cout << vec[i] << " ";
	}
	std::cout << std::endl;
}
