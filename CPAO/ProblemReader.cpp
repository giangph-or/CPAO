#include "ProblemReader.h"
#include <math.h>

void Data::read_data(string facilities_data) {
    fstream input;
    input.open(facilities_data, ios::in);
    input >> number_products;
    input >> number_customers;
    utilities.resize(number_customers);
    double utility_data;
    for(int j = 0; j < number_customers; ++j)
        for (int i = 0; i < number_products; ++i) {
            input >> utility_data;
            utilities[j].push_back(utility_data);
        }
    double revenue_data;
    for (int i = 0; i < number_products; ++i) {
        input >> revenue_data;
        revenue.push_back(revenue_data);
    }
    input.close();
}

void Data::print_data() {
    cout << "nProducts = " << number_products << endl;
    cout << "nCustomers = " << number_customers << endl;
    cout << "utility = " << endl;
    for (int j = 0; j < number_customers; ++j) {
        cout << "[ ";
        for (int i = 0; i < number_products; ++i)
            cout << utilities[j][i] << " ";
        cout << "]" << endl;
    }
    cout << "revenue = [ ";
    for (int i = 0; i < number_products; ++i) {
        cout << revenue[i] << " ";
    }
    cout << "]" << endl;
}