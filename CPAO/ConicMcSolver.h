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

class ConicMcSolver {
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
    ConicMcSolver();
    ConicMcSolver(Data data, double time_limit, string outfile);
    vector<int> find_bound_y(Data data, int i, int budget);
    double calculate_sum_utility(Data data, int budget, int i, int j);
    double calculate_master_obj(Data data, vector<int> x);
    void solve(Data data, int budget);
};
