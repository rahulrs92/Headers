#pragma once

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>

using namespace std;

class bash
{
public:

	static stringstream to_ss(string ic,string in = "")
	{
		system((ic + " > cpp_sys_output" + in).data());
		fstream fi("cpp_sys_output" + in,fstream::in);
		string s;
		stringstream ss;
		while(1)
		{
			getline(fi,s);
			if(fi.eof()) break; else ss<<endl;
			ss<<s;
		}
		fi.close();
		system(("rm cpp_sys_output" + in).data());
		return ss;
	}
	static string to_string(string ic,string in = "") { return to_ss(ic,in).str(); }
	static int to_int(string ic,string in = "") { int oi; to_ss(ic,in)>>oi; return oi; }
};
