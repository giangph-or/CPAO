#include "ProblemReader.h"
#include <math.h>

void Data::read_data(string facilities_data, double noPay) {
    fstream input;
    input.open(facilities_data, ios::in);
    input >> number_products;
    input >> number_customers;

    utilities.resize(number_customers);
    double utility_data;
    for(int i = 0; i < number_customers; ++i)
        for (int j = 0; j < number_products; ++j) {
            input >> utility_data;
            utilities[i].push_back(utility_data);
        }

    revenue.resize(number_customers);
    double revenue_data;
    for (int i = 0; i < number_customers; ++i)
        for (int j = 0; j < number_products; ++j) {
            input >> revenue_data;
            revenue[i].push_back(revenue_data);
        }

    input.close();
    no_purchase.resize(number_customers, noPay);
}

void Data::print_data() {
    cout << "nProducts = " << number_products << endl;
    cout << "nCustomers = " << number_customers << endl;
    cout << "utility = " << endl;
    for (int i = 0; i < number_customers; ++i) {
        cout << "[ ";
        for (int j = 0; j < number_products; ++j)
            cout << utilities[i][j] << " ";
        cout << "]" << endl;
    }

    cout << "revenue = [ ";
    for (int i = 0; i < number_customers; ++i) {
        cout << "[ ";
        for (int j = 0; j < number_products; ++j)
            cout << revenue[i][j] << " ";
        cout << "]" << endl;
    }

    cout << "noPurchase = [ ";
    for (int i = 0; i < number_customers; ++i) {
        cout << no_purchase[i] << " ";
    }
    cout << "]" << endl;
}