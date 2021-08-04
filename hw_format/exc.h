#pragma once

#include <stdexcept>
#include <string>

class WrongUsingBrackets : public std::logic_error {
	public:
		WrongUsingBrackets(const std::string& text) : std::logic_error(text) {}
};

class CorruptedArguments : public std::logic_error {
	public:
		CorruptedArguments(const std::string& text) : std::logic_error(text) {}
};

class CorruptedIndex : public std::logic_error {
	public:
		CorruptedIndex(const std::string& text) : std::logic_error(text) {}
};

