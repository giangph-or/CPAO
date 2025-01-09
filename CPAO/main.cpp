#include<iostream>
#include "ProblemReader.h"
#include "Conic.h"
#include "CEFlog.h"
#include "CuttingPlane.h"
#include "BranchandCut.h"
#include "BCMulticut.h"
#include <string>
#include <iomanip>
#include <chrono>
using namespace std;

int main(int argc, char* argv[]) {
    string instance_name = argv[1];
    string no_pay = argv[2];
    string model = argv[3];
    string instance_file = "AO_data//" + instance_name + ".dat";
    double time_limit = 3600;
    double noPay = stod(no_pay);
    string budget = argv[4];
    Data data_Sen;
    data_Sen.read_data_Sen(instance_file, noPay, stod(budget));
    //data_Sen.read_data_Sen_general(instance_file, noPay, stod(budget), stod(subBudget));
    //data_Sen.print_data();

    //CP and Conic models on Sen et al (2018) dataset
    if (model == "Conic") {
        string out_file = "result_conic//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        Conic conicmcoa(data_Sen, time_limit, out_file);
        conicmcoa.solve(data_Sen);

        //Run from .lp file

        //vector<double> alpha(data_Sen.number_customers, -1);
        //for (int i = 0; i < data_Sen.number_customers; ++i)
        //    for (int j = 0; j < data_Sen.number_products; ++j)
        //        if (data_Sen.revenue[i][j] > alpha[i])
        //            alpha[i] = data_Sen.revenue[i][j];

        //double total_alpha = 0;
        //for (int i = 0; i < data_Sen.number_customers; ++i)
        //    total_alpha += data_Sen.fraction[i] * alpha[i];

        //instance_file = "lpmodel//" + instance_name + "_" + no_pay + "_cap" + budget + ".lp";
        //GRBEnv env = GRBEnv(true);
        //env.start();
        //GRBModel model = GRBModel(env, instance_file);
        //model.set(GRB_DoubleParam_TimeLimit, 600);
        //model.set(GRB_IntParam_Threads, 1);
        //auto start = chrono::steady_clock::now();        
        //model.optimize();
        //auto end = chrono::steady_clock::now();
        //std::chrono::duration<double> solvingTime = end - start;
        //string out_file = "result_conic//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        //ofstream report_results(out_file, ofstream::out);
        //report_results.precision(10);
        //report_results << model.get(GRB_DoubleAttr_ObjVal) << " " << total_alpha - model.get(GRB_DoubleAttr_ObjVal) << " " << solvingTime.count();
        //report_results.close();
    }

    if (model == "CEFlog") {
        string out_file = "result_ceflog//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CEFlog ceflog(data_Sen, time_limit, out_file);
        ceflog.solve(data_Sen);
    }

    if (model == "MILP") {
        //Run from .lp file
        string budget = argv[4];
        instance_file = "Sen_data//" + instance_name + ".dat";
        Data data_Sen;
        data_Sen.read_data_Sen(instance_file, noPay, stod(budget));

        vector<double> alpha(data_Sen.number_customers, -1);
        for (int i = 0; i < data_Sen.number_customers; ++i)
            for (int j = 0; j < data_Sen.number_products; ++j)
                if (data_Sen.revenue[i][j] > alpha[i])
                    alpha[i] = data_Sen.revenue[i][j];

        double total_alpha = 0;
        for (int i = 0; i < data_Sen.number_customers; ++i)
            total_alpha += alpha[i];

        instance_file = "milpmodel//" + instance_name + "_" + no_pay + "_cap" + budget + ".lp";
        GRBEnv env = GRBEnv(true);
        env.start();
        GRBModel model = GRBModel(env, instance_file);
        model.set(GRB_DoubleParam_TimeLimit, 600);
        auto start = chrono::steady_clock::now();
        model.optimize();
        auto end = chrono::steady_clock::now();
        std::chrono::duration<double> solvingTime = end - start;
        string out_file = "result_milp//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        ofstream report_results(out_file, ofstream::out);
        report_results.precision(10);
        report_results << model.get(GRB_DoubleAttr_ObjVal) << " " << total_alpha - model.get(GRB_DoubleAttr_ObjVal) << " " << solvingTime.count();
        report_results.close();
    }

    if (model == "CPLi") {
        //string subBudget = argv[5];
        string out_file = "result_cpli//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlane suboa(data_Sen, time_limit, out_file);
        suboa.solve_milp(data_Sen);
    }

    if (model == "CPBi") {
        //string subBudget = argv[5];
        string out_file = "result_cpbi//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlane suboa(data_Sen, time_limit, out_file);
        suboa.solve_bi(data_Sen);
    }

    if (model == "CPLiSB") {
        //string subBudget = argv[5];
        string out_file = "result_cplisb//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlane suboa(data_Sen, time_limit, out_file);
        suboa.solve_multicut_milp(data_Sen, 20);
    }

    if (model == "CPBiSB") {
        //string subBudget = argv[5];
        string out_file = "result_cpbisb//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlane suboa(data_Sen, time_limit, out_file);
        suboa.solve_multicut_bi(data_Sen, 20);
    }

    if (model == "CPLiMulti") {
        //string subBudget = argv[5];
        string out_file = "result_cplimulti//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlane suboa(data_Sen, time_limit, out_file);
        suboa.solve_multi_multicut_milp(data_Sen, 20, 20);
    }

    if (model == "CPBiMulti") {
        //string subBudget = argv[5];
        string out_file = "result_cpbimulti//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlane suboa(data_Sen, time_limit, out_file);
        suboa.solve_multi_multicut_bi(data_Sen, 20, 20);
    }

    if (model == "BCLi") {
        //string subBudget = argv[5];
        string out_file = "result_bcli//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        BranchandCut bcsuboa(data_Sen, time_limit, out_file);
        bcsuboa.solve_multicut_milp(data_Sen);
    }

    if (model == "BCBi") {
        //string subBudget = argv[5];
        string out_file = "result_bcbi//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        BranchandCut bcconic(data_Sen, time_limit, out_file);
        bcconic.solve_multicut_bi(data_Sen);
    }

    if (model == "BCLiSB") {
        //string subBudget = argv[5];
        string out_file = "result_bclisb//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        BCMulticut bcmulti(data_Sen, time_limit, out_file);
        bcmulti.solve_multicut_milp(data_Sen, 20);
    }

    if (model == "BCBiSB") {
        //string subBudget = argv[5];
        string out_file = "result_bcbisb//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        BCMulticut bcmulti(data_Sen, time_limit, out_file);
        bcmulti.solve_multicut_bi(data_Sen, 20);
    }
}