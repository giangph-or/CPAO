#include "BCMulticutSen.h"
#include <chrono>
#include <algorithm>
#include <sstream>
#include <cstdlib>
#include <cassert>

BCMulticutSen::BCMulticutSen() {

}

BCMulticutSen::BCMulticutSen(Data data, double time_limit, string outfile) {
	this->data = data;
	this->time_limit = time_limit;
	this->out_res_csv = outfile;
}

CBMulticutSen::CBMulticutSen() {

}

CBMulticutSen::CBMulticutSen(GRBVar* varx, GRBVar* vary, GRBVar* varcut, int p, int c, vector<double> n, vector<vector<double>> u, vector<vector<double>> r, vector<double> a, vector<double> f, int cut) {
	x = varx;
	y = vary;
	cut_customers = varcut;
	products = p;
	customers = c;
	noPay = n;
	util = u;
	ren = r;
	al = a;
	cuts = cut;
	frac = f;
}

double BCMulticutSen::calculate_z(Data data, vector<double> alpha, vector<double> denominator) {
	double z = 0;
	for (int i = 0; i < data.number_customers; ++i)
		z += alpha[i] * data.no_purchase[i] / denominator[i];

	return z;
}

double BCMulticutSen::calculate_optimal_bound_denominator(Data data, int i) {
	GRBEnv env = GRBEnv(true);
	env.start();

	GRBModel model = GRBModel(env);

	//Decison variables: x_j = 1 if product j is chosen, 0 otherwise
	GRBVar* x = 0;
	x = model.addVars(data.number_products, GRB_BINARY);

	//Set capacity constraints
	for (int s = 0; s < data.number_sets; ++s) {
		GRBLinExpr sum;
		for (int j = 0; j < data.number_products; ++j)
			if (data.in_set[j][s] == 1)
				sum += data.cost[j] * x[j];
		model.addConstr(sum <= data.capacity_each_set, "ct_set_cap" + to_string(s));
	}

	//GRBLinExpr cost;
	//for (int j = 0; j < data.number_products; ++j)
	//	cost += data.fraction2[j] * x[j];
	//model.addConstr(cost <= data.capacity_each_set, "ct_set_cap");

	GRBLinExpr obj;
	for (int j = 0; j < data.number_products; ++j)
		obj += data.utilities[i][j] * x[j];
	model.setObjective(obj, GRB_MAXIMIZE);

	model.set(GRB_IntParam_OutputFlag, 0);

	model.optimize();

	return model.get(GRB_DoubleAttr_ObjVal) + data.no_purchase[i];
}

double BCMulticutSen::calculate_master_obj(Data data, vector<int> x) {
	double obj = 0;
	for (int i = 0; i < data.number_customers; ++i) {
		double ts = 0, ms = data.no_purchase[i];
		for (int j = 0; j < data.number_products; ++j) {
			ts += data.revenue[i][j] * x[j] * data.utilities[i][j];
			ms += x[j] * data.utilities[i][j];
		}
		obj += data.fraction[i] * ts / ms;
	}

	return obj;
}

vector<vector<double>> BCMulticutSen::calculate_bound_y_in(Data data) {
	vector<vector<double>> lb_in(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i)
		lb_in[i].resize(data.number_products);

	for (int i = 0; i < data.number_customers; ++i) {
		for (int k = 0; k < data.number_products; ++k) {
			vector<pair<double, int>> u(data.number_products);
			for (int j = 0; j < data.number_products; ++j)
				u[j] = make_pair(data.utilities[i][j], j);

			sort(u.begin(), u.end(), greater<pair<double, int>>());
			for (int j = 0; j < data.number_products; ++j)
				if (u[j].second == k) {
					u.erase(u.begin() + j);
					break;
				}

			lb_in[i][k] = data.utilities[i][k];

			int count = 0;
			while (count < data.capacity_each_set) {
				lb_in[i][k] += u[count].first;
				count++;
			}

			////General
			//double cost = 0;
			//for (int a = 0; a < u.size(); ++a)
			//	if (cost + data.fraction2[u[a].second] <= data.capacity_each_set) {
			//		lb_in[i][k] += u[a].first;
			//		cost += data.fraction2[u[a].second];
			//	}
		}
	}
	return lb_in;
}

vector<vector<double>> BCMulticutSen::calculate_bound_y_notin(Data data) {
	vector<vector<double>> lb_notin(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i)
		lb_notin[i].resize(data.number_products);

	for (int i = 0; i < data.number_customers; ++i) {
		for (int k = 0; k < data.number_products; ++k) {
			vector<pair<double, int>> u(data.number_products);
			for (int j = 0; j < data.number_products; ++j)
				u[j] = make_pair(data.utilities[i][j], j);

			sort(u.begin(), u.end(), greater<pair<double, int>>());
			for (int j = 0; j < data.number_products; ++j)
				if (u[j].second == k) {
					u.erase(u.begin() + j);
					break;
				}

			lb_notin[i][k] = 0;

			int count = 0;
			while (count < data.capacity_each_set) {
				lb_notin[i][k] += u[count].first;
				count++;
			}

			////General
			//double cost = 0;
			//for (int a = 0; a < u.size(); ++a)
			//	if (cost + data.fraction2[u[a].second] <= data.capacity_each_set) {
			//		lb_notin[i][k] += u[a].first;
			//		cost += data.fraction2[u[a].second];
			//	}
		}
	}
	return lb_notin;
}

vector<vector<double>> BCMulticutSen::subset_bound_y_in(Data data) {
	vector<vector<double>> lb_in(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i)
		lb_in[i].resize(data.number_products);

	for (int i = 0; i < data.number_customers; ++i) {
		for (int k = 0; k < data.number_products; ++k) {
			vector<pair<double, int>> u(data.number_products);
			for (int j = 0; j < data.number_products; ++j)
				u[j] = make_pair(data.utilities[i][j], j);

			sort(u.begin(), u.end(), greater<pair<double, int>>());
			for (int j = 0; j < data.number_products; ++j)
				if (u[j].second == k) {
					u.erase(u.begin() + j);
					break;
				}

			lb_in[i][k] = data.utilities[i][k];
			vector<int> in(5, data.sub_capacity_each_set);
			for (int a = 0; a < u.size(); ++a)
				for (int s = 0; s < 5; ++s)
					if (in[s] > 0 && data.in_set[u[a].second][s] == 1) {
						lb_in[i][k] += u[a].first;
						in[s]--;
					}
		}
	}
	return lb_in;
}

vector<vector<double>> BCMulticutSen::subset_bound_y_notin(Data data) {
	vector<vector<double>> lb_notin(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i)
		lb_notin[i].resize(data.number_products);

	for (int i = 0; i < data.number_customers; ++i) {
		for (int k = 0; k < data.number_products; ++k) {
			vector<pair<double, int>> u(data.number_products);
			for (int j = 0; j < data.number_products; ++j)
				u[j] = make_pair(data.utilities[i][j], j);

			sort(u.begin(), u.end(), greater<pair<double, int>>());
			for (int j = 0; j < data.number_products; ++j)
				if (u[j].second == k) {
					u.erase(u.begin() + j);
					break;
				}

			lb_notin[i][k] = 0;
			vector<int> in(5, data.sub_capacity_each_set);
			for (int a = 0; a < u.size(); ++a)
				for (int s = 0; s < 5; ++s)
					if (in[s] > 0 && data.in_set[u[a].second][s] == 1) {
						lb_notin[i][k] += u[a].first;
						in[s]--;
					}
		}
	}
	return lb_notin;
}

void BCMulticutSen::solve_multicut_milp(Data data, int number_cuts) {
	vector<double> alpha(data.number_customers, -1);
	for (int i = 0; i < data.number_customers; ++i)
		for (int j = 0; j < data.number_products; ++j)
			if (data.revenue[i][j] > alpha[i])
				alpha[i] = data.revenue[i][j];

	vector<vector<double>> bound_in = calculate_bound_y_in(data);
	vector<vector<double>> bound_notin = calculate_bound_y_notin(data);
	////General
	//vector<vector<double>> subset_bound_in = subset_bound_y_in(data);
	//vector<vector<double>> subset_bound_notin = subset_bound_y_notin(data);

	vector<double> upper_bound_denominator(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i)
		upper_bound_denominator[i] = calculate_optimal_bound_denominator(data, i);

	//auto start = chrono::steady_clock::now();

	GRBEnv env = GRBEnv(true);
	env.start();

	GRBModel model = GRBModel(env);

	//cout << "Decison variables : x\n" << endl;
	GRBVar* x;
	x = new GRBVar[data.number_products];
	for (int j = 0; j < data.number_products; ++j)
		x[j] = model.addVar(0, 1, 0, GRB_BINARY, "x_" + to_string(j));

	//cout << "Slack variables : y_i\n" << endl;
	vector<double> upper_bound_y(data.number_customers, 0);
	vector<double> lower_bound_y(data.number_customers, 0);
	for (int i = 0; i < data.number_customers; ++i) {
		lower_bound_y[i] = 1 / upper_bound_denominator[i];
		upper_bound_y[i] = 1 / data.no_purchase[i];
	}

	GRBVar* y;
	y = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		y[i] = model.addVar(lower_bound_y[i], upper_bound_y[i], 0, GRB_CONTINUOUS, "y_" + to_string(i));

	GRBVar** z;
	z = new GRBVar * [data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		z[i] = new GRBVar[data.number_products];
	for (int i = 0; i < data.number_customers; ++i)
		for (int j = 0; j < data.number_products; ++j)
			z[i][j] = model.addVar(0, upper_bound_y[i], 0, GRB_CONTINUOUS, "z_" + to_string(i) + "_" + to_string(j));

	vector<double> upper_bound_cut_customers(number_cuts, 0);
	vector<double> lower_bound_cut_customers(number_cuts, 0);
	for (int c = 0; c < number_cuts; ++c)
		for (int i = 0; i < data.number_customers; ++i)
			if (i % number_cuts == c) {
				lower_bound_cut_customers[c] += data.fraction[i] * data.no_purchase[i] * alpha[i] * 1 / upper_bound_denominator[i];
				upper_bound_cut_customers[c] += data.fraction[i] * data.no_purchase[i] * alpha[i] * 1 / data.no_purchase[i];
			}

	GRBVar* cut_customers;
	cut_customers = new GRBVar[number_cuts];
	for (int c = 0; c < number_cuts; ++c)
		cut_customers[c] = model.addVar(lower_bound_cut_customers[c], upper_bound_cut_customers[c], 0, GRB_CONTINUOUS, "cut_customers_" + to_string(c));

	//Set capacity constraints
	for (int s = 0; s < data.number_sets; ++s) {
		GRBLinExpr sum;
		for (int j = 0; j < data.number_products; ++j)
			if (data.in_set[j][s] == 1)
				sum += data.cost[j] * x[j];
		model.addConstr(sum <= data.capacity_each_set, "ct_set_cap" + to_string(s));
	}

	////General
	//for (int s = 0; s < 5; ++s) {
	//	GRBLinExpr sum;
	//	for (int j = 0; j < data.number_products; ++j)
	//		if (data.in_set[j][s] == 1)
	//			sum += data.cost[j] * x[j];
	//	model.addConstr(sum <= data.sub_capacity_each_set, "ct_set_cap" + to_string(s));
	//}

	//GRBLinExpr cost;
	//for (int j = 0; j < data.number_products; ++j)
	//	cost += data.fraction2[j] * x[j];
	//model.addConstr(cost <= data.capacity_each_set, "ct_set_cap");

	//cout << "Constraints related to y and z\n" << endl;
	for (int i = 0; i < data.number_customers; ++i) {
		GRBLinExpr sum_z;
		for (int j = 0; j < data.number_products; ++j)
			sum_z += data.utilities[i][j] * z[i][j];

		model.addConstr(data.no_purchase[i] * y[i] + sum_z >= 1);
	}

	//cout << "z_ij = y_i * x_j\n" << endl;
	for (int i = 0; i < data.number_customers; ++i)
		for (int j = 0; j < data.number_products; ++j) {
			//model.addConstr(y[i] - z[i][j] <= upper_bound_y[i] - upper_bound_y[i] * x[j]);
			//model.addConstr(z[i][j] <= upper_bound_y[i] * x[j]);
			model.addConstr(z[i][j] <= y[i]);
		}

	//McCornick constraints
	for (int i = 0; i < data.number_customers; ++i)
		for (int j = 0; j < data.number_products; ++j) {
			model.addConstr(z[i][j] <= x[j] * (1 / (data.no_purchase[i] + data.utilities[i][j])));

			//if (bound_in[i][j] >= subset_bound_in[i][j])
			model.addConstr(z[i][j] >= x[j] * (1 / (data.no_purchase[i] + bound_in[i][j])));
			//else
				//model.addConstr(z[i][j] >= x[j] * (1 / (data.no_purchase[i] + subset_bound_in[i][j])));

			//if (bound_notin[i][j] >= subset_bound_notin[i][j])
			model.addConstr(z[i][j] <= y[i] - (1 - x[j]) * (1 / (data.no_purchase[i] + bound_notin[i][j])));
			//else
				//model.addConstr(z[i][j] <= y[i] - (1 - x[j]) * (1 / (data.no_purchase[i] + subset_bound_notin[i][j])));

			model.addConstr(z[i][j] >= y[i] - (1 - x[j]) * (1 / data.no_purchase[i]));
		}

	for (int c = 0; c < number_cuts; ++c) {
		GRBLinExpr sum_y;
		for (int i = 0; i < data.number_customers; ++i)
			if (i % number_cuts == c)
				sum_y += data.fraction[i] * data.no_purchase[i] * alpha[i] * y[i];

		model.addConstr(cut_customers[c] >= sum_y, "ct_sub_gradient_y_" + to_string(c));
	}

	//Objective
	GRBLinExpr obj;
	for (int c = 0; c < number_cuts; ++c)
		obj += cut_customers[c];
	for (int i = 0; i < data.number_customers; ++i)
		for (int j = 0; j < data.number_products; ++j)
			obj += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * z[i][j];

	model.setObjective(obj, GRB_MINIMIZE);

	//auto time_before_cut = chrono::steady_clock::now();
	//chrono::duration<double> before_cut = time_before_cut - start;

	//double run_time = time_limit - before_cut.count();

	double obj_val_cplex = 0.0;

	model.write("bc.lp");
	model.set(GRB_IntParam_LazyConstraints, 1);
	model.set(GRB_IntParam_PreCrush, 1);
	//model.set(GRB_DoubleParam_TimeLimit, run_time);
	model.set(GRB_DoubleParam_TimeLimit, time_limit);
	model.set(GRB_IntParam_Threads, 8);
	//model.set(GRB_IntParam_OutputFlag, 0);

	CBMulticutSen cb = CBMulticutSen(x, y, cut_customers, data.number_products, data.number_customers, data.no_purchase, data.utilities, data.revenue, alpha, data.fraction, number_cuts);
	model.setCallback(&cb);

	auto start = chrono::steady_clock::now();

	model.optimize();

	auto end = chrono::steady_clock::now();
	chrono::duration<double> total_time = end - start;
	time_for_solve = total_time.count();

	vector<int> sol_x(data.number_products);

	if (model.get(GRB_IntAttr_SolCount) > 0) {
		//update obj, variables
		obj_val_cplex = model.get(GRB_DoubleAttr_ObjVal);
		cout << "\nResult product list: " << endl;
		for (int j = 0; j < data.number_products; ++j)
			if (x[j].get(GRB_DoubleAttr_X) > 0.5) {
				sol_x[j] = 1;
				cout << j << " ";
			}
			else sol_x[j] = 0;
		cout << endl;

		cout << "Gurobi obj = " << std::setprecision(7) << fixed << obj_val_cplex << endl;
		master_obj_val = calculate_master_obj(data, sol_x);
		cout << "Master obj = " << std::setprecision(7) << fixed << master_obj_val << endl;

		//check time
		auto time_now = std::chrono::steady_clock::now(); //get now time
		std::chrono::duration<double> after_cut = time_now - start;
		cout << "Time now: " << after_cut.count() << endl;
		cout << "--- --- --- --- --- --- ---" << endl;
	}
	else {
		cout << "No solution found..." << endl;
		auto end = chrono::steady_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;
		time_for_solve = elapsed_seconds.count();
	}

	cout << "\nObjective value: " << setprecision(5) << master_obj_val << endl;
	cout << "Solution: ";
	for (int j = 0; j < data.number_products; ++j)
		if (sol_x[j] == 1)
			cout << j << " ";
	cout << "\nTotal time: " << time_for_solve << " seconds" << endl;

	ofstream report_results(out_res_csv, ofstream::out);
	report_results.precision(10);
	report_results << obj_val_cplex << " " << master_obj_val << " " << time_for_solve << endl;
	for (int j = 0; j < data.number_products; ++j)
		if (sol_x[j] == 1)
			report_results << j << " ";
	report_results.close();
}

void BCMulticutSen::solve_multicut_bi(Data data, int number_cuts) {
	vector<double> alpha(data.number_customers, -1);
	for (int i = 0; i < data.number_customers; ++i)
		for (int j = 0; j < data.number_products; ++j)
			if (data.revenue[i][j] > alpha[i])
				alpha[i] = data.revenue[i][j];

	vector<vector<double>> bound_in = calculate_bound_y_in(data);
	vector<vector<double>> bound_notin = calculate_bound_y_notin(data);
	////General
	//vector<vector<double>> subset_bound_in = subset_bound_y_in(data);
	//vector<vector<double>> subset_bound_notin = subset_bound_y_notin(data);

	vector<double> upper_bound_denominator(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i)
		upper_bound_denominator[i] = calculate_optimal_bound_denominator(data, i);

	//auto start = chrono::steady_clock::now();

	GRBEnv env = GRBEnv(true);
	env.start();

	GRBModel model = GRBModel(env);

	//cout << "Decison variables : x\n" << endl;
	GRBVar* x;
	x = new GRBVar[data.number_products];
	for (int j = 0; j < data.number_products; ++j)
		x[j] = model.addVar(0, 1, 0, GRB_BINARY, "x_" + to_string(j));

	//cout << "Slack variables : y_i\n" << endl;
	vector<double> upper_bound_y(data.number_customers, 0);
	vector<double> lower_bound_y(data.number_customers, 0);
	for (int i = 0; i < data.number_customers; ++i) {
		lower_bound_y[i] = 1 / upper_bound_denominator[i];
		upper_bound_y[i] = 1 / data.no_purchase[i];
	}

	GRBVar* y;
	y = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		y[i] = model.addVar(lower_bound_y[i], upper_bound_y[i], 0, GRB_CONTINUOUS, "y_" + to_string(i));

	vector<double> upper_bound_cut_customers(number_cuts, 0);
	vector<double> lower_bound_cut_customers(number_cuts, 0);
	for (int c = 0; c < number_cuts; ++c)
		for (int i = 0; i < data.number_customers; ++i)
			if (i % number_cuts == c) {
				lower_bound_cut_customers[c] += data.fraction[i] * data.no_purchase[i] * alpha[i] * 1 / upper_bound_denominator[i];
				upper_bound_cut_customers[c] += data.fraction[i] * data.no_purchase[i] * alpha[i] * 1 / data.no_purchase[i];
			}

	GRBVar* cut_customers;
	cut_customers = new GRBVar[number_cuts];
	for (int c = 0; c < number_cuts; ++c)
		cut_customers[c] = model.addVar(lower_bound_cut_customers[c], upper_bound_cut_customers[c], 0, GRB_CONTINUOUS, "cut_customers_" + to_string(c));

	//Set capacity constraints
	for (int s = 0; s < data.number_sets; ++s) {
		GRBLinExpr sum;
		for (int j = 0; j < data.number_products; ++j)
			if (data.in_set[j][s] == 1)
				sum += data.cost[j] * x[j];
		model.addConstr(sum <= data.capacity_each_set, "ct_set_cap" + to_string(s));
	}

	////General
	//for (int s = 0; s < 5; ++s) {
	//	GRBLinExpr sum;
	//	for (int j = 0; j < data.number_products; ++j)
	//		if (data.in_set[j][s] == 1)
	//			sum += data.cost[j] * x[j];
	//	model.addConstr(sum <= data.sub_capacity_each_set, "ct_set_cap" + to_string(s));
	//}

	//GRBLinExpr cost;
	//for (int j = 0; j < data.number_products; ++j)
	//	cost += data.fraction2[j] * x[j];
	//model.addConstr(cost <= data.capacity_each_set, "ct_set_cap");

	//cout << "Constraints related to y and z\n" << endl;
	for (int i = 0; i < data.number_customers; ++i) {
		GRBQuadExpr sum_z;
		for (int j = 0; j < data.number_products; ++j)
			sum_z += data.utilities[i][j] * y[i] * x[j];

		model.addQConstr(data.no_purchase[i] * y[i] + sum_z >= 1);
	}

	//McCornick constraints
	for (int i = 0; i < data.number_customers; ++i)
		for (int j = 0; j < data.number_products; ++j) {
			model.addQConstr(y[i] * x[j] <= x[j] * (1 / (data.no_purchase[i] + data.utilities[i][j])));

			//if (bound_in[i][j] >= subset_bound_in[i][j])
			model.addQConstr(y[i] * x[j] >= x[j] * (1 / (data.no_purchase[i] + bound_in[i][j])));
			//else
				//model.addQConstr(y[i] * x[j] >= x[j] * (1 / (data.no_purchase[i] + subset_bound_in[i][j])));

			//if (bound_notin[i][j] >= subset_bound_notin[i][j])
			model.addQConstr(y[i] * x[j] <= y[i] - (1 - x[j]) * (1 / (data.no_purchase[i] + bound_notin[i][j])));
			//else
				//model.addQConstr(y[i] * x[j] <= y[i] - (1 - x[j]) * (1 / (data.no_purchase[i] + subset_bound_notin[i][j])));

			model.addQConstr(y[i] * x[j] >= y[i] - (1 - x[j]) * (1 / data.no_purchase[i]));
		}

	for (int c = 0; c < number_cuts; ++c) {
		GRBLinExpr sum_y;
		for (int i = 0; i < data.number_customers; ++i)
			if (i % number_cuts == c)
				sum_y += data.fraction[i] * data.no_purchase[i] * alpha[i] * y[i];

		model.addConstr(cut_customers[c] >= sum_y, "ct_sub_gradient_y_" + to_string(c));
	}

	//Objective
	GRBQuadExpr obj;
	for (int c = 0; c < number_cuts; ++c)
		obj += cut_customers[c];
	for (int i = 0; i < data.number_customers; ++i)
		for (int j = 0; j < data.number_products; ++j)
			obj += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * y[i] * x[j];

	model.setObjective(obj, GRB_MINIMIZE);

	//auto time_before_cut = chrono::steady_clock::now();
	//chrono::duration<double> before_cut = time_before_cut - start;

	//double run_time = time_limit - before_cut.count();

	double obj_val_cplex = 0.0;

	model.write("bc.lp");
	model.set(GRB_IntParam_LazyConstraints, 1);
	model.set(GRB_IntParam_PreCrush, 1);
	model.set(GRB_IntParam_MIQCPMethod, 1);
	//model.set(GRB_DoubleParam_TimeLimit, run_time);
	model.set(GRB_DoubleParam_TimeLimit, time_limit);
	model.set(GRB_IntParam_Threads, 8);
	//model.set(GRB_IntParam_OutputFlag, 0);

	CBMulticutSen cb = CBMulticutSen(x, y, cut_customers, data.number_products, data.number_customers, data.no_purchase, data.utilities, data.revenue, alpha, data.fraction, number_cuts);
	model.setCallback(&cb);

	auto start = chrono::steady_clock::now();

	model.optimize();

	auto end = chrono::steady_clock::now();
	chrono::duration<double> total_time = end - start;
	time_for_solve = total_time.count();

	vector<int> sol_x(data.number_products);

	if (model.get(GRB_IntAttr_SolCount) > 0) {
		//update obj, variables
		obj_val_cplex = model.get(GRB_DoubleAttr_ObjVal);
		cout << "\nResult product list: " << endl;
		for (int j = 0; j < data.number_products; ++j)
			if (x[j].get(GRB_DoubleAttr_X) > 0.5) {
				sol_x[j] = 1;
				cout << j << " ";
			}
			else sol_x[j] = 0;
		cout << endl;

		cout << "Gurobi obj = " << std::setprecision(7) << fixed << obj_val_cplex << endl;
		master_obj_val = calculate_master_obj(data, sol_x);
		cout << "Master obj = " << std::setprecision(7) << fixed << master_obj_val << endl;

		//check time
		auto time_now = std::chrono::steady_clock::now(); //get now time
		std::chrono::duration<double> after_cut = time_now - start;
		cout << "Time now: " << after_cut.count() << endl;
		cout << "--- --- --- --- --- --- ---" << endl;
	}
	else {
		cout << "No solution found..." << endl;
		auto end = chrono::steady_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;
		time_for_solve = elapsed_seconds.count();
	}

	cout << "\nObjective value: " << setprecision(5) << master_obj_val << endl;
	cout << "Solution: ";
	for (int j = 0; j < data.number_products; ++j)
		if (sol_x[j] == 1)
			cout << j << " ";
	cout << "\nTotal time: " << time_for_solve << " seconds" << endl;

	ofstream report_results(out_res_csv, ofstream::out);
	report_results.precision(10);
	report_results << obj_val_cplex << " " << master_obj_val << " " << time_for_solve << endl;
	for (int j = 0; j < data.number_products; ++j)
		if (sol_x[j] == 1)
			report_results << j << " ";
	report_results.close();
}

void CBMulticutSen::callback() {
	try {
		if (where == GRB_CB_MIPSOL) {
			double* initial_x = new double[products];
			double* initial_y = new double[customers];
			initial_x = getSolution(x, products);
			initial_y = getSolution(y, customers);

			vector<double> initial_denominator(customers, 0);
			for (int i = 0; i < customers; ++i) {
				initial_denominator[i] += noPay[i];
				for (int j = 0; j < products; ++j)
					initial_denominator[i] += initial_x[j] * util[i][j];
			}

			//cout << "Outer-cuts y\n" << endl;
			vector<double> partitial_y(customers, 0);
			for (int i = 0; i < customers; ++i)
				partitial_y[i] = 1 / initial_denominator[i];

			vector<vector<double>> subgradient_y(customers);
			for (int i = 0; i < customers; ++i)
				subgradient_y[i].resize(products, 0);

			for (int i = 0; i < customers; ++i)
				for (int j = 0; j < products; ++j)
					subgradient_y[i][j] -= util[i][j] / (initial_denominator[i] * initial_denominator[i]);

			for (int c = 0; c < cuts; ++c) {
				GRBLinExpr grad_y;
				double part_y = 0;
				for (int i = 0; i < customers; ++i)
					if (i % cuts == c) {
						part_y += frac[i] * noPay[i] * al[i] * partitial_y[i];
						for (int j = 0; j < products; ++j)
							grad_y += frac[i] * noPay[i] * al[i] * subgradient_y[i][j] * (x[j] - initial_x[j]);
					}

				addLazy(cut_customers[c] >= part_y + grad_y);
			}

			////cout << "Calculate total utility\n" << endl;
			//vector<double> sum_uti_customer(customers, 0);
			//for (int i = 0; i < customers; ++i) {
			//	sum_uti_customer[i] += noPay[i];
			//	for (int j = 0; j < products; ++j)
			//		sum_uti_customer[i] += util[i][j];
			//}

			////cout << "Submodular-cuts\n" << endl;
			//for (int i = 0; i < customers; ++i) {
			//	if (initial_y[i] < partitial_y[i]) {
			//		GRBLinExpr submodular_cut_a_z, submodular_cut_b_z;
			//		for (int j = 0; j < products; ++j)
			//			if (initial_x[j] == 1) {
			//				submodular_cut_a_z += (1 - x[j]) * util[i][j] /
			//					(sum_uti_customer[i] * (sum_uti_customer[i] - util[i][j]));
			//				submodular_cut_b_z += (1 - x[j]) * util[i][j] /
			//					(initial_denominator[i] * (initial_denominator[i] - util[i][j]));
			//			}
			//			else {
			//				submodular_cut_a_z -= x[j] * util[i][j] /
			//					(initial_denominator[i] * (initial_denominator[i] + util[i][j]));
			//				submodular_cut_b_z -= x[j] * util[i][j] /
			//					(noPay[i] * (noPay[i] + util[i][j]));
			//			}

			//		submodular_cut_a_z += partitial_y[i];
			//		addLazy(y[i] >= submodular_cut_a_z);
			//		submodular_cut_b_z += partitial_y[i];
			//		addLazy(y[i] >= submodular_cut_b_z);
			//		//cout << "Submodular-cuts " << i << "\n" << endl;
			//	}
			//}
		}
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during callback" << endl;
	}
}