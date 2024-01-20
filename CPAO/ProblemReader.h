#include <iostream>
#include <vector>
#include <fstream>

using namespace std;
class Data
{
public:
	int number_products;
	int number_customers;
	vector<vector<double>> utilities;
	vector<vector<double>> revenue;
	vector<double> no_purchase;
	vector<int> in_set;
	vector<double> cost;
	int capacity_each_set;
	int number_sets;
	int total_capacity;
	void read_data(string facilities_data, double noPay);
	void read_general_data(string facilities_data, double noPay);
	void print_data();
	void print_data_general();
};
#pragma once
