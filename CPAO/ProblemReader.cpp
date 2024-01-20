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

void Data::read_general_data(string facilities_data, double noPay) {
    fstream input;
    input.open(facilities_data, ios::in);
    input >> number_products;
    input >> number_customers;
    input >> total_capacity;
    input >> capacity_each_set;
    input >> number_sets;

    utilities.resize(number_customers);
    double utility_data;
    for (int i = 0; i < number_customers; ++i)
        for (int j = 0; j < number_products; ++j) {
            input >> utility_data;
            utilities[i].push_back(utility_data);
        }

    vector<double> tmp_revenue(number_products);
    for (int j = 0; j < number_products; ++j)
        input >> tmp_revenue[j];

    revenue.resize(number_customers);
    for (int i = 0; i < number_customers; ++i)
         revenue[i]= tmp_revenue;

    in_set.resize(number_products);
    for (int j = 0; j < number_products; ++j)
        input >> in_set[j];
    for (int j = 0; j < number_products; ++j)
        in_set[j] -= 1;

    cost.resize(number_products);
    for (int j = 0; j < number_products; ++j)
        input >> cost[j];

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

    cout << "\nrevenue = [ ";
    for (int i = 0; i < number_customers; ++i) {
        cout << "[ ";
        for (int j = 0; j < number_products; ++j)
            cout << revenue[i][j] << " ";
        cout << "]" << endl;
    }

    cout << "\nnoPurchase = [ ";
    for (int i = 0; i < number_customers; ++i) {
        cout << no_purchase[i] << " ";
    }
    cout << "]" << endl;
}

void Data::print_data_general() {
    cout << "nProducts = " << number_products << endl;
    cout << "nCustomers = " << number_customers << endl;
    cout << "utility = " << endl;
    for (int i = 0; i < number_customers; ++i) {
        cout << "[ ";
        for (int j = 0; j < number_products; ++j)
            cout << utilities[i][j] << " ";
        cout << "]" << endl;
    }

    cout << "\nrevenue = [ ";
    for (int i = 0; i < number_customers; ++i) {
        cout << "[ ";
        for (int j = 0; j < number_products; ++j)
            cout << revenue[i][j] << " ";
        cout << "]" << endl;
    }

    cout << "\ncost = [ ";
    for (int j = 0; j < number_products; ++j)
        cout << cost[j] << " ";
    cout << "]" << endl;

    cout << "\nin_set = [ ";
    for (int j = 0; j < number_products; ++j)
        cout << in_set[j] << " ";
    cout << "]" << endl;

    cout << "\nnoPurchase = [ ";
    for (int i = 0; i < number_customers; ++i) {
        cout << no_purchase[i] << " ";
    }
    cout << "]" << endl;
}