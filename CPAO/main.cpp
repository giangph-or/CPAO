#include<iostream>
#include "ProblemReader.h"
#include "CuttingPlaneSolver.h"
#include "CuttingPlaneGurobi.h"
#include "CuttingPlaneNCuts.h"
#include "CuttingPlaneGeneral.h"
#include "BCGurobi.h"
#include "BCNCuts.h"
#include "BCGeneral.h"
#include "MILPSolver.h"
#include "MILPGurobi.h"
#include "ConicMcSolver.h"
#include "ConicMcGurobi.h"
#include "ConicMcGeneral.h"
#include <string>
#include <iomanip>
using namespace std;

int main(int argc, char* argv[]) {
    string instance_name = argv[1];
    string no_pay = argv[2];
    string model = argv[4];
    string budget = argv[3];
    //string instance_file = "AO_data//" + instance_name + ".dat";
    string instance_file = "AO_data//" + instance_name + ".dat";
    double time_limit = 3600;
    double noPay = stod(no_pay);
    Data data;
    //data.read_data(instance_file, noPay);
    data.read_general_data(instance_file, noPay);
    //data.print_data_general();
    if (model == "CP") {
        /*string out_file = "AO_result_cp//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlaneSolver cpoa(data, time_limit, out_file);
        cpoa.solve(data, stoi(budget));*/
        //string out_file = "AO_result_cp_gurobi//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        string out_file = "AO_result_cp_gurobi_build_in//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlaneGurobi cpoa(data, time_limit, out_file);
        //cpoa.solve(data, stoi(budget));
        cpoa.solve_build_in(data, stoi(budget));
    }
    if (model == "CPNCuts") {
        /*string out_file = "AO_result_cp//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlaneSolver cpoa(data, time_limit, out_file);
        cpoa.solve(data, stoi(budget));*/
        //string out_file = "AO_result_cp_gurobi//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        string out_file = "AO_result_cp_gurobi_ncuts//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlaneNCuts cpoa(data, time_limit, out_file);
        //cpoa.solve(data, stoi(budget));
        cpoa.solve_build_in(data, stoi(budget), 5);
    }
    if (model == "CPGeneral") {
        /*string out_file = "AO_result_cp//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlaneSolver cpoa(data, time_limit, out_file);
        cpoa.solve(data, stoi(budget));*/
        //string out_file = "AO_result_cp_gurobi//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        string out_file = "AO_result_cp_gurobi_general//" + instance_name + "_" + no_pay + ".txt";
        CuttingPlaneGeneral cpoa(data, time_limit, out_file);
        //cpoa.solve(data, stoi(budget));
        cpoa.solve_build_in(data, 5);
    }
    if (model == "BC") {
        /*string out_file = "AO_result_cp//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlaneSolver cpoa(data, time_limit, out_file);
        cpoa.solve(data, stoi(budget));*/
        //string out_file = "AO_result_cp_gurobi//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        string out_file = "AO_result_bc_gurobi_build_in//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        BCGurobi bcoa(data, time_limit, out_file);
        //cpoa.solve(data, stoi(budget));
        bcoa.solve_build_in(data, stoi(budget));
    }
    if (model == "BCNCuts") {
        /*string out_file = "AO_result_cp//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlaneSolver cpoa(data, time_limit, out_file);
        cpoa.solve(data, stoi(budget));*/
        //string out_file = "AO_result_cp_gurobi//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        string out_file = "AO_result_bc_gurobi_ncuts//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        BCNCuts bcoa(data, time_limit, out_file);
        //cpoa.solve(data, stoi(budget));
        bcoa.solve_build_in(data, stoi(budget), 5);
    }
    if (model == "BCGeneral") {
        /*string out_file = "AO_result_cp//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlaneSolver cpoa(data, time_limit, out_file);
        cpoa.solve(data, stoi(budget));*/
        //string out_file = "AO_result_cp_gurobi//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        string out_file = "AO_result_bc_gurobi_general//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        BCGeneral bcoa(data, time_limit, out_file);
        //cpoa.solve(data, stoi(budget));
        bcoa.solve_build_in(data, 5);
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
    if (model == "ConicGeneral") {
        /*string out_file = "AO_result_conic//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        ConicMcSolver conicmcoa(data, time_limit, out_file);
        conicmcoa.solve(data, stoi(budget));*/
        string out_file = "AO_result_conic_general//" + instance_name + "_" + no_pay + ".txt";
        ConicMcGeneral conicmcoa(data, time_limit, out_file);
        conicmcoa.solve(data);
    }
}