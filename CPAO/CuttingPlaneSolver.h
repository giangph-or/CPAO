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

class CuttingPlaneSolver {
public:
    Data data;
    string output;
    string log_file;
    string out_res_csv;
    double time_limit;
    double master_obj_val;
    double gap;
    double time_for_solve;
    char var_name[1000];
    CuttingPlaneSolver();
    CuttingPlaneSolver(Data data, double time_limit, string outfile);
    vector<vector<double>> create_optimal_sub_intervals(Data data);
    vector<double> calculate_y(Data data, vector<int> x, double alpha);
    vector<double> calculate_z(Data data, vector<int> x, double alpha);
    double calculate_original_obj(Data data, vector<int> x, double alpha);
    double calculate_master_obj(Data data, vector<int> x);
    void solve(Data data);
};
