#include <iostream>
#include <sstream>
#include <vector>
#include <cctype>
#include "exc.h"

using std::to_string;

	template<typename T>
std::string to_string(const T &t)
{
	std::stringstream ss;
	ss << t;
	return ss.str();
}


std::string format(const std::string &s)
{
	return s;
}

	template<class... Args>
std::string format(const std::string &s, const Args &...args)
{
	std::vector<std::string> strs = {to_string(args)...};

	static const char FORMAT_SYMBOL = '{';
	std::string res;
	std::string buf;
	bool arg = false;
	for (size_t i = 0; i <= s.size(); ++i)
	{
		bool last = (i == s.size());
		if(arg && last)
		{
			throw WrongUsingBrackets("Wrong using brackets");
		}
		if(arg)
		{
			if(isdigit(s[i]))
			{
				buf += s[i];
			}
			else
			{
				if(s[i] == '}')
				{	
					size_t num = 0;
					if (!buf.empty())
					{
						try
						{
							num = std::stoi(buf);
						}
						catch(std::out_of_range)
						{
							std::cout << "index out of range" << std::endl;
						}
					}
					else
						throw CorruptedIndex("Empty index");

					if (num >= 0 && num < strs.size())
						res += strs[num];
					else
						throw CorruptedIndex("Wrong index or arguments");
					buf.clear();
					arg = false;
				}
				else
				{
					throw CorruptedIndex("Wrong index");
				}
			}
		}
		else
		{
			if(s[i] == FORMAT_SYMBOL)
			{
				arg = true;
			}
			else
			{
				if(s[i] == '}')
					throw WrongUsingBrackets("Wrong using brackets");
				if(!last)
					res += s[i];
			}
		}
	}
	return res;
}
