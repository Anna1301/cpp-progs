#include "tpool.h"
#include <numeric>
#include <iostream>
#include <unistd.h>
#include <iostream>
#include <cassert>
#include <cmath>


int64_t from_binary_to_decimal(const std::vector<int> &vec, size_t i, size_t j) 
{
	int res = 0;
	for (size_t l = i; l < j; ++l) 
	{
		res += vec[l]*pow(2, l);
	}
	return res;
}


int main()
{
	vector<int> vec(20, 1);
	int64_t multiple_num = 0;
	int64_t single_num = 0;

	size_t thr_num = 4;
	ThreadPool pool(thr_num);

	vector<future<int64_t>> futures;

	for (size_t i = 0; i != vec.size(); i += vec.size() / thr_num) 
	{
		futures.push_back(pool.exec(from_binary_to_decimal, vec, i, i + vec.size() / thr_num));
	}

	for (auto& f : futures) 
	{
		multiple_num += f.get();
	}

	single_num = from_binary_to_decimal(vec, 0, vec.size());

	assert(multiple_num == single_num);

	std::cout << "success" << std::endl;

	return 0;
}	
