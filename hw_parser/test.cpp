#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "parser.h"
#include <cassert>
#include <functional>


void ParseNum(std::string &token, std::ofstream &out)
{
	out << token << "\n";
}

void  ParseStr(std::string &token, std::ofstream &out)
{
	out << token << "\n";
}

void Start_fun(std::ofstream &out)
{
	out << "Start parcing line\n" ;
}

void Finish_fun(std::ofstream &out)
{
	out << "Finish parcing line\n";
}



void DefaultWorkTest()
{
	callBackNum num_parcer = &ParseNum;
	callBackStart start_parcer = &Start_fun;
	callBackFinish finish_parcer = &Finish_fun;
	callBackStr str_parcer = &ParseStr;

	Parser p;
	std::string line;
	std::ifstream in("textin.txt");
	std::ofstream fout;
	fout.open("textout.txt");

	if(in.is_open())
	{
		while (getline(in,line))
		{
			p.MakeParser(line);
			p.ParseLine(start_parcer, num_parcer, str_parcer, finish_parcer, fout);
		}
	}

	in.close();
	fout.close();

	std::ifstream checkin("textout.txt");
	if(checkin.is_open())
	{
		getline(checkin, line);
		assert(line == "Start parcing line");

		getline(checkin, line);
		assert(line == "pbpbpbpbpbp");

		getline(checkin, line); 
		assert(line == "45");

		getline(checkin, line);
		assert(line == "7");

		getline(checkin, line);   
		assert(line == "Finish parcing line");

		getline(checkin, line);  
		assert(line == "Start parcing line");

		getline(checkin, line); 
		assert(line == "3gufi");

		getline(checkin, line);          
		assert(line == "Finish parcing line");

		assert(!getline(checkin,line));  
	}
	checkin.close();
}



int main()
{
	DefaultWorkTest();

	std::cout << "success" << std::endl;

	return 0;
}
