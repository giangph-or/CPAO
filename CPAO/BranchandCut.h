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

class CBBranchandCut : public GRBCallback {
public:
    int products;
    int customers;
    vector<double> noPay;
    vector<vector<double>> util;
    vector<vector<double>> ren;
    GRBVar* x;
    GRBVar* y;
    CBBranchandCut();
    CBBranchandCut(GRBVar* x, GRBVar* y, int products, int customers, vector<double> noPay, vector<vector<double>> util, vector<vector<double>> ren);
protected:
    void callback();
};

class BranchandCut {
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
    BranchandCut();
    BranchandCut(Data data, double time_limit, string outfile);
    double calculate_optimal_bound_denominator(Data data, int i);
    double calculate_master_obj(Data data, vector<int> x);
    vector<vector<double>> calculate_bound_y_in(Data data);
    vector<vector<double>> calculate_bound_y_notin(Data data);
    vector<vector<double>> subset_bound_y_in(Data data);
    vector<vector<double>> subset_bound_y_notin(Data data);
    void solve_multicut_bi(Data data);
    void solve_multicut_milp(Data data);
};

