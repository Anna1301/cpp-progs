#pragma once
#include <iostream>


enum class Error
{
	NoError,
	CorruptedArchive
};

class Serializer
{
	private:  
	static constexpr char Separator = ' ';
	std::ostream &out_;

	public:
	explicit Serializer(std::ostream &out)
		: out_(out)
	{
	}

	template <class T>
	Error save(T &object)
	{
		return object.serialize(*this);
	}

	template <class... ArgsT>
	Error operator()(ArgsT &...args)
	{
		return process(args...);
	}

	private:
	template <typename T>
	Error process(T &object)
	{
		return save(object);
	}

	template <typename T, typename... ArgsT>
	Error process(T &object, ArgsT &...args)
	{
		if(save(object) == Error::NoError)
		{
			return process(args...);
		}
		else 
		{
			return Error::CorruptedArchive;
		}
	}	
	// process использует variadic templates
};

template <>
Error Serializer::save<uint64_t>(uint64_t &object);

template <>
Error Serializer::save<bool>(bool &object);

class Deserializer
{
	private:
	static constexpr char Separator = ' ';
	std::istream &in_;

	public:
	Deserializer(std::istream &in)
		: in_(in)
	{
	}

	template <class T>
	Error load(T &object)
	{
		return object.deserialize(*this);
	}

	template <class... ArgsT>
	Error operator()(ArgsT &...args)
	{
		return process(args...);
	}

	private:
	template <typename T>
	Error process(T &object)
	{
		return load(object);
	}

	template <typename T, typename... ArgsT>
	Error process(T &object, ArgsT &...args)
	{
		if(load(object) == Error::NoError)
		{
			return process(args...);
		}
		else
		{
			return Error::CorruptedArchive;
		}
	}
};

template <>
Error Deserializer::load<uint64_t>(uint64_t &object);

template <>
Error Deserializer::load<bool>(bool &object);


