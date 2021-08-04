#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <functional>

using callBackStart = std::function<void(std::ofstream &)>;
using callBackNum = std::function<void(std::string &, std:: ofstream &)>;
using callBackStr = std::function<void(std::string &, std:: ofstream &)>; 
using callBackFinish = std::function<void(std::ofstream &)>;

class Parser
{
	private:
		std::string line;

	public:
		Parser()
		{
			line = "";
		}

		void MakeParser(std::string &line_)
		{
			line = line_;
		}

		void ParseLine(callBackStart call_back_start, callBackNum call_back_num, callBackStr call_back_str,
				callBackFinish call_back_finish, std::ofstream &out);

};
