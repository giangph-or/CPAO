#include <iostream>
#include <vector>
#include <fstream>
#include "ProblemReader.h"

using namespace std;
class Report
{
public:
	int nProducts;
	string sol_file;
	string report_file;
	void create_report(vector<string> instance, vector<string> noPay, vector<string> budget, string nProduct);
};
#pragma once

