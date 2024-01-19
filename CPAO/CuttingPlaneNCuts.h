#include <iostream>
#include <math.h>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <C:/gurobi1100/win64/include/gurobi_c++.h>
#include "ProblemReader.h"
using namespace std;

class CuttingPlaneNCuts {
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
    CuttingPlaneNCuts();
    CuttingPlaneNCuts(Data data, double time_limit, string outfile);
    vector<double> calculate_y(Data data, vector<int> x, vector<double> alpha);
    vector<double> calculate_z(Data data, vector<int> x);
    double calculate_original_obj(Data data, vector<int> x, vector<double> alpha);
    double calculate_master_obj(Data data, vector<int> x);
    double calculate_original_obj_tmp(Data data, vector<int> x, vector<double> alpha, int i);
    double calculate_bound_y(Data data, int budget, int i, double alpha);
    double calculate_bound_z(Data data, int budget, int i);
    double calculate_optimal_bound_y(Data data, int budget, int i, double alpha);
    vector<int> greedy(Data data, int budget, vector<double> alpha);
    void solve_build_in(Data data, int budget, int nCuts);
};
