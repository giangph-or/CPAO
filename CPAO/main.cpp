#include<iostream>
#include "ProblemReader.h"
#include "CuttingPlaneSolver.h"
#include "CuttingPlaneGurobi.h"
#include "MILPSolver.h"
#include "MILPGurobi.h"
#include "ConicMcSolver.h"
#include "ConicMcGurobi.h"
#include <string>
#include <iomanip>
using namespace std;

int main(int argc, char* argv[]) {
    string instance_name = argv[1];
    string no_pay = argv[2];
    string model = argv[4];
    string budget = argv[3];
    string instance_file = "AO_data//" + instance_name + ".dat";
    double time_limit = 3600;
    double noPay = stod(no_pay);
    Data data;
    data.read_data(instance_file, noPay);
    //data.print_data();
    if (model == "CP") {
        /*string out_file = "AO_result_cp//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlaneSolver cpoa(data, time_limit, out_file);
        cpoa.solve(data, stoi(budget));*/
        string out_file = "AO_result_cp_gurobi//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlaneGurobi cpoa(data, time_limit, out_file);
        cpoa.solve(data, stoi(budget));
    }
    if (model == "MILP") {
        //string out_file = "AO_result_milp//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        //MILPSolver milpoa(data, time_limit, out_file);
        //milpoa.solve(data, stoi(budget));
        string out_file = "AO_result_milp_gurobi//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        MILPGurobi milpoa(data, time_limit, out_file);
        milpoa.solve(data, stoi(budget));
    }
    if (model == "Conic") {
        /*string out_file = "AO_result_conic//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        ConicMcSolver conicmcoa(data, time_limit, out_file);
        conicmcoa.solve(data, stoi(budget));*/
        string out_file = "AO_result_conic_gurobi//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        ConicMcGurobi conicmcoa(data, time_limit, out_file);
        conicmcoa.solve(data, stoi(budget));
    }
}