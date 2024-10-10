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

class CuttingPlane {
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
    CuttingPlane();
    CuttingPlane(Data data, double time_limit, string outfile);
    vector<double> calculate_denominator(Data data, vector<int> x);
    double calculate_z(Data data, vector<double> alpha, vector<double> denominator);
    double calculate_optimal_bound_denominator(Data data, int i);
    double calculate_original_obj(Data data, vector<int> x, vector<double> alpha);
    double calculate_master_obj(Data data, vector<int> x);
    double calculate_original_obj_tmp(Data data, vector<int> x, vector<double> alpha, int i);
    vector<int> greedy(Data data, vector<double> alpha);
    vector<int> greedy_general(Data data, vector<double> alpha);
    vector<vector<double>> calculate_bound_y_in(Data data);
    vector<vector<double>> calculate_bound_y_notin(Data data);
    vector<vector<double>> subset_bound_y_in(Data data);
    vector<vector<double>> subset_bound_y_notin(Data data);
    void solve_milp(Data data);
    void solve_multicut_milp(Data data, int number_cut);
    void solve_bi(Data data);
    void solve_multicut_bi(Data data, int number_cuts);
    void solve_multi_multicut(Data data, int number_cut, int number_cut_theta);
};
