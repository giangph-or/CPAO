#include<iostream>
#include "ProblemReader.h"
#include "CuttingPlaneSolver.h"
#include "MILPSolver.h"
#include "ConicMcSolver.h"
#include <string>
#include <iomanip>
using namespace std;

int main(int argc, char* argv[]) {
    string instance_name = argv[1];
    string no_pay = argv[2];
    string model = argv[3];
    string budget = argv[4];
    string instance_file = "AO_data//" + instance_name + ".dat";
    double time_limit = 3600;
    double noPay = stod(no_pay);
    Data data;
    data.read_data(instance_file, noPay);
    //data.print_data();
    if (model == "CP") {
        string out_file = "AO_result_cp//" + instance_name + "_" + budget + ".txt";
        CuttingPlaneSolver cpoa(data, time_limit, out_file);
        cpoa.solve(data, stoi(budget));
    }
    if (model == "MILP") {
        string out_file = "AO_result_milp//" + instance_name + "_" + budget + ".txt";
        MILPSolver milpoa(data, time_limit, out_file);
        milpoa.solve(data, stoi(budget));
    }
    if (model == "Conic") {
        string out_file = "AO_result_conic//" + instance_name + "_" + budget + ".txt";
        ConicMcSolver conicmcoa(data, time_limit, out_file);
        conicmcoa.solve(data, stoi(budget));
    }
}