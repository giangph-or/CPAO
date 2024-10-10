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

class CBMulticut : public GRBCallback {
public:
    int products;
    int customers;
    int cuts;
    vector<double> noPay;
    vector<double> al;
    vector<double> frac;
    vector<vector<double>> util;
    vector<vector<double>> ren;
    GRBVar* x;
    GRBVar* y;
    GRBVar* cut_customers;
    GRBVar* theta;
    CBMulticut();
    CBMulticut(GRBVar* x, GRBVar* y, GRBVar* cut_customers, GRBVar* theta, int products, int customers, vector<double> noPay, vector<vector<double>> util, vector<vector<double>> ren, vector<double> al, vector<double> frac, int cuts);
protected:
    void callback();
};

class BCMulticut {
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
    BCMulticut();
    BCMulticut(Data data, double time_limit, string outfile);
    double calculate_z(Data data, vector<double> alpha, vector<double> denominator);
    double calculate_optimal_bound_denominator(Data data, int i);
    double calculate_master_obj(Data data, vector<int> x);
    vector<vector<double>> calculate_bound_y_in(Data data);
    vector<vector<double>> calculate_bound_y_notin(Data data);
    vector<vector<double>> subset_bound_y_in(Data data);
    vector<vector<double>> subset_bound_y_notin(Data data);
    void solve_multicut_milp(Data data, int number_cuts);
    void solve_multicut_bi(Data data, int number_cuts);
};
