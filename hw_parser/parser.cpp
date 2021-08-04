#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "parser.h"
#include <functional>

void Parser:: ParseLine(callBackStart call_back_start, callBackNum call_back_num, callBackStr call_back_str,
		callBackFinish call_back_finish, std::ofstream &out)
{
	int len = line.length();
	int current_pos = 0;
	int i;
	std::string token = "";

	call_back_start(out);

	while(current_pos < len)
	{
		if(!isspace(line[current_pos]) && current_pos != len-1)
		{
			token += line[current_pos];
			++current_pos;			
		}
		else
		{
			if(isdigit(line[current_pos]))
			{
				call_back_num(token, out);
			}
			else
			{
				call_back_str(token, out);
			}

			token = "";
			++current_pos;
		}

	}

	call_back_finish(out);
}

