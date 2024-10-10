#include "ProblemReader.h"
#include <math.h>

void Data::read_data(string data, double noPay) {
    fstream input;
    input.open(data, ios::in);
    input >> number_products;
    input >> number_customers;
    input >> number_sets;
    input >> capacity_each_set;

    utilities.resize(number_customers);
    double utility_data;
    for (int i = 0; i < number_customers; ++i)
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

    cost.resize(number_products);
    for (int j = 0; j < number_products; ++j)
        input >> cost[j];

    in_set.resize(number_products);
    int set_id;
    for (int j = 0; j < number_products; ++j)
        for (int s = 0; s < number_sets; ++s) {
            input >> set_id;
            in_set[j].push_back(set_id);
        }

    input.close();
    no_purchase.resize(number_customers, noPay);
}

void Data::read_data_Sen(string data, double noPay, double budget) {
    fstream input;
    input.open(data, ios::in);
    input >> number_products;
    input >> number_customers;
    number_sets = 1;
    capacity_each_set = budget;

    utilities.resize(number_customers);
    double utility_data;
    for (int i = 0; i < number_customers; ++i)
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

    cost.resize(number_products);
    for (int j = 0; j < number_products; ++j)
        cost[j] = 1;

    //cost.resize(number_products);
    //for (int j = 0; j < number_products; ++j)
    //    input >> cost[j];

    //200_20 and 500_50
    //fraction.resize(number_customers, 1);

    //100_100 and 1000_100 and larger number of customers
    fraction.resize(number_customers);
    for (int i = 0; i < number_customers; ++i)
        input >> fraction[i];

    in_set.resize(number_products);
    int set_id;
    for (int j = 0; j < number_products; ++j)
        for (int s = 0; s < number_sets; ++s)
            in_set[j].push_back(1);

    input.close();
    no_purchase.resize(number_customers, noPay);
}

void Data::read_data_Sen_general(string data, double noPay, double budget, double subBudget) {
    fstream input;
    input.open(data, ios::in);
    input >> number_products;
    input >> number_customers;
    number_sets = 1;
    capacity_each_set = budget;
    sub_capacity_each_set = subBudget;

    utilities.resize(number_customers);
    double utility_data;
    for (int i = 0; i < number_customers; ++i)
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

    cost.resize(number_products);
    for (int j = 0; j < number_products; ++j)
        cost[j] = 1;

    //cost.resize(number_products);
    //for (int j = 0; j < number_products; ++j)
    //    input >> cost[j];

    fraction.resize(number_customers, 1);

    //100_100 and 1000_100
    //fraction.resize(number_customers);
    //for (int i = 0; i < number_customers; ++i)
    //    input >> fraction[i];

    fraction2.resize(number_products);
    for (int j = 0; j < number_products; ++j)
        input >> fraction2[j];

    in_set.resize(number_products);
    for (int j = 0; j < number_products; j++)
        in_set[j].resize(5, 0);
    int count = 0;
    for (int j = 0; j < number_products; ++j)
        for (int s = 0; s < 5; ++s)
            if (s == j / 40)
                in_set[j][s] = 1;
            else
                in_set[j][s] = 0;

    input.close();
    no_purchase.resize(number_customers, noPay);
}

vector<vector<double>> Data::compute_utilities(double alpha, double beta, vector<double> incumbent_utilities, vector<vector<double>> util) {
    //cout << "Utility: " << endl;
    vector<vector<double>> re;
    for (int i = 0; i < incumbent_utilities.size(); ++i)
        incumbent_utilities[i] = alpha * beta * incumbent_utilities[i];

    re.resize(incumbent_utilities.size());
    for (int i = 0; i < incumbent_utilities.size(); ++i)
        for (int j = 0; j < util[i].size(); ++j) {
            util[i][j] = beta * util[i][j];
            util[i][j] = exp(util[i][j] - incumbent_utilities[i]);
            re[i].push_back(util[i][j]);
        }
    return re;
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

    cout << "\ngroup = [ ";
    for (int i = 0; i < number_sets; ++i) {
        cout << "[ ";
        for (int j = 0; j < number_products; ++j) {
            if (in_set[j][i] == 1)
                cout << j << " ";
        }
        cout << "]" << endl;
    }
    cout << "]" << endl;
}