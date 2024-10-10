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

class Conic {
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
    Conic();
    Conic(Data data, double time_limit, string outfile);
    double calculate_master_obj(Data data, vector<int> x);
    double calculate_optimal_bound_denominator(Data data, int i);
    vector<vector<double>> calculate_bound_y_in(Data data);
    vector<vector<double>> calculate_bound_y_notin(Data data);
    vector<vector<double>> subset_bound_y_in(Data data);
    vector<vector<double>> subset_bound_y_notin(Data data);
    void solve(Data data);
};