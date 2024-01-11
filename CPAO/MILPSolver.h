#include <iostream>
#include <math.h>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#define IL_STD
#include <ilcplex/ilocplex.h>
#include "ProblemReader.h"
using namespace std;

class MILPSolver {
public:
    Data data;
    string output;
    string log_file;
    string out_res_csv;
    double time_limit;
    double obj_val;
    double gap;
    double time_for_solve;
    char var_name[1000];
    MILPSolver();
    MILPSolver(Data data, double time_limit, string outfile);
    void solve(Data data);
};
