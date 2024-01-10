#include <iostream>
#include <cmath>
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

class CuttingPlaneSolver {
public:
    Data data;
    const char* output;
    string log;
    string out_res_csv;
    double time_limit;
    double obj_val;
    double gap;
    double time_for_solve;
    char var_name[1000];
    CuttingPlaneSolver();
    CuttingPlaneSolver(Data data, const char* output, string log, string out_res_csv, double time_limit);
    vector<vector<double>> create_optimal_sub_intervals(Data data, vector<int> initial_x);
    bool solve(Data data, vector<int> initial_x);
};
