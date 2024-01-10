#include<iostream>
#include "ProblemReader.h"
#include <string>
#include <iomanip>
using namespace std;

int main(int argc, char* argv[]) {
    string instance_name = argv[1];
    string instance_file = "AO_data//" + instance_name + ".dat";
    Data data;
    data.read_data(instance_file);
    data.print_data();
}