#include<iostream>
#include "ProblemReader.h"
#include "CuttingPlaneSolver.h"
#include <string>
#include <iomanip>
using namespace std;

int main(int argc, char* argv[]) {
    string instance_name = argv[1];
    string no_pay = argv[2];
    string instance_file = "AO_data//" + instance_name + ".dat";
    double noPay = stod(no_pay);
    Data data;
    data.read_data(instance_file, noPay);
    //data.print_data();
    CuttingPlaneSolver cpoa;
    bool solved = cpoa.solve(data);
}