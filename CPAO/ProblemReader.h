#include<iostream>
#include<vector>
#include<fstream>

using namespace std;
class Data
{
public:
	int number_products;
	int number_customers;
	vector<vector<double>> utilities;
	vector<double> revenue;
	void read_data(string facilities_data);
	void print_data();
};
#pragma once
