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

class CBPolymatroid : public GRBCallback {
public:
    Data data;
    vector<int> index_w;
    vector<double> alpha;
    GRBVar* s;
    GRBVar** w;
    CBPolymatroid();
    CBPolymatroid(GRBVar* s, GRBVar** w, Data data, vector<int> index_w, vector<double> alpha);
protected:
    void callback();
};

class CEFlog {
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
    CEFlog();
    CEFlog(Data data, double time_limit, string outfile);
    double calculate_master_obj(Data data, vector<int> x);
    double calculate_optimal_bound_denominator(Data data, int i);
    void solve(Data data);
};
