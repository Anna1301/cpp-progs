#include <iostream>
#include "serializer.h"

bool Long_is_positive_digit(const std::string &str)
{
	bool res = true;
	int i = 0;
	while(res == true && i < str.size())
	{
		if(!std::isdigit(str[i]))
		{
			res = false;
		}
		i++;
	}
	return res;
}



template <>
Error Serializer::save(uint64_t &object)
{
	out_ << object << Separator;

	if (out_.fail())
	{
		return Error::CorruptedArchive;
	}
	else
	{
		return Error::NoError;
	}
}

template <>
Error Serializer::save<bool>(bool &object)
{
	if (!object)
	{
		out_ << "false" << Separator;
	}
	else
	{
		out_ << "true" << Separator;
	}

	if (out_.fail())
	{
		return Error::CorruptedArchive;
	}
	else
	{
		return Error::NoError;
	}
}

template <>
Error Deserializer::load(uint64_t &object)
{
	std::string str;
	in_ >> str;
	if(!Long_is_positive_digit(str))
	{
		return Error::CorruptedArchive;
	}
	else
	{
		try
		{
			object = std::stoull(str);
		}
		catch(std::invalid_argument)
		{
			return Error::CorruptedArchive;
		}
		catch(std::out_of_range)
		{
			return Error::CorruptedArchive;
		}

	}

	if (in_.fail())
	{
		return Error::CorruptedArchive;
	}
	else
	{
		return Error::NoError;
	}
}

template <>
Error Deserializer::load<bool>(bool &object)
{
	std::string str;
	in_ >> str;
	if (str == "false")
	{
		object = false;
	}
	else if (str == "true")
	{
		object = true;
	}
	else
	{
		return Error::CorruptedArchive;
	}
	return Error::NoError;
}


