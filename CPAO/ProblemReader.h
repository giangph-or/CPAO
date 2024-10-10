#include <iostream>
#include <vector>
#include <fstream>

using namespace std;
class Data
{
public:
	int number_products;
	int number_customers;
	int capacity_each_set;
	int sub_capacity_each_set;
	int number_sets;
	vector<vector<double>> utilities;
	vector<vector<double>> revenue;
	vector<double> no_purchase;
	vector<vector<int>> in_set;
	vector<double> cost;
	vector<double> fraction;
	vector<double> fraction2;
	void read_data(string data, double no_purchase);
	void read_data_Sen(string data, double no_purchase, double budget);
	void read_data_Sen_general(string data, double no_purchase, double budget, double subBudget);
	vector<vector<double>> compute_utilities(double alpha, double beta, vector<double> incum, vector<vector<double>> util);
	void print_data();
};
#pragma once