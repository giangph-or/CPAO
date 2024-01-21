#include "CuttingPlaneGeneral.h"
#include <chrono>
#include <algorithm>
#include <sstream>
#include <cstdlib>
#include <cassert>

CuttingPlaneGeneral::CuttingPlaneGeneral() {

}

CuttingPlaneGeneral::CuttingPlaneGeneral(Data data, double time_limit, string outfile) {
	this->data = data;
	this->time_limit = time_limit;
	this->out_res_csv = outfile;
}

vector<double> CuttingPlaneGeneral::calculate_y(Data data, vector<int> x, vector<double> alpha) {
	vector<double> y(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i) {
		double tmp_y = alpha[i] * data.no_purchase[i];
		for (int j = 0; j < data.number_products; ++j)
			tmp_y += (alpha[i] - data.revenue[i][j]) * x[j] * data.utilities[i][j];
		y[i] = log(tmp_y);
	}

	return y;
}

vector<double> CuttingPlaneGeneral::calculate_z(Data data, vector<int> x) {
	vector<double> z(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i) {
		double tmp_z = data.no_purchase[i];
		for (int j = 0; j < data.number_products; ++j)
			tmp_z += x[j] * data.utilities[i][j];
		z[i] = -log(tmp_z);
	}

	return z;
}

double  CuttingPlaneGeneral::calculate_original_obj(Data data, vector<int> x, vector<double> alpha) {
	double obj = 0;
	for (int i = 0; i < data.number_customers; ++i) {
		double ts = alpha[i] * data.no_purchase[i], ms = data.no_purchase[i];
		for (int j = 0; j < data.number_products; ++j) {
			ts += (alpha[i] - data.revenue[i][j]) * x[j] * data.utilities[i][j];
			ms += x[j] * data.utilities[i][j];
		}
		obj += ts / ms;
	}

	return obj;
}

double CuttingPlaneGeneral::calculate_master_obj(Data data, vector<int> x) {
	double obj = 0;
	for (int i = 0; i < data.number_customers; ++i) {
		double ts = 0, ms = data.no_purchase[i];
		for (int j = 0; j < data.number_products; ++j) {
			ts += data.revenue[i][j] * x[j] * data.utilities[i][j];
			ms += x[j] * data.utilities[i][j];
		}
		obj += ts / ms;
	}

	return obj;
}

double CuttingPlaneGeneral::calculate_original_obj_tmp(Data data, vector<int> x, vector<double> alpha, int candidate) {
	double obj = 0;
	for (int i = 0; i < data.number_customers; ++i) {
		double ts = alpha[i] * data.no_purchase[i] + (alpha[i] - data.revenue[i][candidate]) * data.utilities[i][candidate],
			ms = data.no_purchase[i] + data.utilities[i][candidate];
		for (int j = 0; j < data.number_products; ++j) {
			ts += data.revenue[i][j] * x[j] * data.utilities[i][j];
			ms += x[j] * data.utilities[i][j];
		}

		obj += ts / ms;
	}

	return obj;
}

vector<int> CuttingPlaneGeneral::greedy(Data data, vector<double> alpha) {
	auto start = chrono::steady_clock::now();
	vector<int> chosen(data.number_products, 0);
	double obj = calculate_original_obj(data, chosen, alpha);
	//cout << "Obj: " << obj << endl;
	double current_capacity = 0;
	vector<int> set_capacity(data.number_sets, 0);

	while (current_capacity < data.total_capacity) {
		double extra_ratio = 0;
		int inserted_product = -1;
		double min = 999999;
		for (int j = 0; j < data.number_products; ++j) {
			if (chosen[j] == 0 && set_capacity[data.in_set[j]] + 1 <= data.capacity_each_set && current_capacity + data.cost[j] <= data.total_capacity) {
				extra_ratio = calculate_original_obj_tmp(data, chosen, alpha, j) - obj;
				if (min > extra_ratio) {
					min = extra_ratio;
					inserted_product = j;
				}
			}
		}
		if (inserted_product != -1) {
			chosen[inserted_product] = 1;
			cout << "Best-Insertion: " << inserted_product << endl;
			obj = calculate_original_obj(data, chosen, alpha);
			current_capacity += data.cost[inserted_product];
			cout << "Current capacity: " << current_capacity << endl;
			set_capacity[data.in_set[inserted_product]]++;
		}
		else break;
	}

	auto end = chrono::steady_clock::now();
	std::chrono::duration<double> greedy_time = end - start;

	cout << "Greedy solution: " << endl;
	for (int j = 0; j < data.number_products; j++)
		if (chosen[j] == 1)
			cout << j << " ";
	cout << endl;
	cout << "Sub obj = " << obj << endl;
	double master_obj = 0;
	for (int i = 0; i < data.number_customers; i++)
		master_obj += alpha[i];
	cout << "Master obj = " << master_obj - obj << endl;
	cout << "Time: " << greedy_time.count() << endl;

	return chosen;
}

double CuttingPlaneGeneral::calculate_bound_y_total(Data data, int i, double alpha) {
	vector<pair<double, int>> u(data.number_products);
	for (int j = 0; j < data.number_products; ++j)
		u[j] = make_pair(data.utilities[i][j], j);

	sort(u.begin(), u.end(), greater<pair<double, int>>());

	double total_cost = 0;
	double current_capacity = 0;
	double total_cap_obj = 0;
	int count = 0;
	while (count < data.number_products) {
		if (current_capacity + data.cost[u[count].second] <= data.total_capacity) {
			current_capacity += data.cost[u[count].second];
			total_cap_obj += u[count].first * (alpha - data.revenue[i][count]);
		}
		count++;
	}

	return total_cap_obj + alpha * data.no_purchase[i];
}

double CuttingPlaneGeneral::calculate_bound_y_set(Data data, int i, double alpha) {
	vector<pair<double, int>> u(data.number_products);
	for (int j = 0; j < data.number_products; ++j)
		u[j] = make_pair(data.utilities[i][j], j);

	sort(u.begin(), u.end(), greater<pair<double, int>>());

	double set_cap_obj = 0;
	vector<int> set_capacity(data.number_sets, 0);
	int count = 0;
	while (count < data.number_products) {
		if (set_capacity[data.in_set[count]] + 1 <= data.capacity_each_set) {
			set_capacity[data.in_set[count]]++;
			set_cap_obj += u[count].first * (alpha - data.revenue[i][count]);
		}
		count++;
	}

	return set_cap_obj + alpha * data.no_purchase[i];
}

double CuttingPlaneGeneral::calculate_bound_z_total(Data data, int i) {
	vector<pair<double, int>> u(data.number_products);
	for (int j = 0; j < data.number_products; ++j)
		u[j] = make_pair(data.utilities[i][j], j);

	sort(u.begin(), u.end(), greater<pair<double, int>>());

	double current_capacity = 0;
	double total_cap_obj = 0;
	int count = 0;
	while (count < data.number_products) {
		if (current_capacity + data.cost[u[count].second] <= data.total_capacity) {
			current_capacity += data.cost[u[count].second];
			total_cap_obj += u[count].first;
		}
		count++;
	}

	return total_cap_obj + data.no_purchase[i];
}

double CuttingPlaneGeneral::calculate_bound_z_set(Data data, int i) {
	vector<pair<double, int>> u(data.number_products);
	for (int j = 0; j < data.number_products; ++j)
		u[j] = make_pair(data.utilities[i][j], j);

	sort(u.begin(), u.end(), greater<pair<double, int>>());

	double set_cap_obj = 0;
	vector<int> set_capacity(data.number_sets, 0);
	int count = 0;
	while (count < data.number_products) {
		if (set_capacity[data.in_set[count]] + 1 <= data.capacity_each_set) {
			set_capacity[data.in_set[count]]++;
			set_cap_obj += u[count].first;
		}
		count++;
	}

	return set_cap_obj + data.no_purchase[i];
}

double CuttingPlaneGeneral::calculate_optimal_bound_y(Data data, int i, double alpha) {
	GRBEnv env = GRBEnv(true);
	env.start();

	GRBModel model = GRBModel(env);

	//Decison variables: x_j = 1 if product j is chosen, 0 otherwise
	GRBVar* x = 0;
	x = model.addVars(data.number_products, GRB_BINARY);

	//Budget constraint
	GRBLinExpr capacity = 0;
	for (int j = 0; j < data.number_products; ++j) {
		capacity += data.cost[j] * x[j];
	}
	model.addConstr(capacity <= data.total_capacity);

	GRBLinExpr obj;
	for (int j = 0; j < data.number_products; ++j)
		obj += data.utilities[i][j] * x[j] * (alpha - data.revenue[i][j]);
	model.setObjective(obj, GRB_MAXIMIZE);

	model.optimize();

	return model.get(GRB_DoubleAttr_ObjVal) + alpha * data.no_purchase[i];
}

double CuttingPlaneGeneral::calculate_optimal_bound_z(Data data, int i) {
	GRBEnv env = GRBEnv(true);
	env.start();

	GRBModel model = GRBModel(env);

	//Decison variables: x_j = 1 if product j is chosen, 0 otherwise
	GRBVar* x = 0;
	x = model.addVars(data.number_products, GRB_BINARY);

	//Budget constraint
	GRBLinExpr capacity = 0;
	for (int j = 0; j < data.number_products; ++j) {
		capacity += data.cost[j] * x[j];
	}
	model.addConstr(capacity <= data.total_capacity);

	GRBLinExpr obj;
	for (int j = 0; j < data.number_products; ++j)
		obj += data.utilities[i][j] * x[j];
	model.setObjective(obj, GRB_MAXIMIZE);

	model.optimize();

	return model.get(GRB_DoubleAttr_ObjVal) + data.no_purchase[i];
}

void CuttingPlaneGeneral::solve_build_in(Data data, int nCuts) {
	auto start = chrono::steady_clock::now(); //get start time

	vector<int> index;
	for (int i = 0; i < data.number_customers; ++i) {
		index.push_back(i);
	}
	vector<vector<int>> group(nCuts);

	int size = data.number_customers;
	srand(time(NULL));
	for (int t = 0; t < nCuts; ++t) {
		for (int i = 0; i < data.number_customers / nCuts; ++i) {
			int res = rand() % size;
			group[t].push_back(index[res]);
			index.erase(index.begin() + res);
			size--;
			index.resize(size);
		}
	}

	vector<double> alpha(data.number_customers, -1);
	for (int i = 0; i < data.number_customers; ++i)
		for (int j = 0; j < data.number_products; ++j)
			if (data.revenue[i][j] > alpha[i])
				alpha[i] = data.revenue[i][j];

	//Calculate initial_x, initial_y, initial_z
	//vector<int> initial_x(data.number_products, 0);
	vector<int> initial_x = greedy(data, alpha);

	vector<int> best_x = initial_x;
	double best_obj = calculate_master_obj(data, best_x);

	vector<double> initial_y = calculate_y(data, initial_x, alpha);
	vector<double> initial_z = calculate_z(data, initial_x);
	vector<double> initial_theta(nCuts, 0);

	GRBEnv env = GRBEnv(true);
	env.start();

	GRBModel model = GRBModel(env);

	//cout << "Decison variables : x_j\n" << endl;
	GRBVar* x;
	x = new GRBVar[data.number_products];
	for (int j = 0; j < data.number_products; ++j)
		x[j] = model.addVar(0, 1, 0, GRB_BINARY, "x_" + to_string(j));

	//cout << "Slack variables : y_i\n" << endl;
	GRBVar* y;
	GRBVar* u;
	y = new GRBVar[data.number_customers];
	u = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i) {
		double bound;
		//double bound_total = calculate_optimal_bound_y(data, i, alpha[i]);
		double bound_total = calculate_bound_y_total(data, i, alpha[i]);
		double bound_set = calculate_bound_y_set(data, i, alpha[i]);
		if (bound_set > bound_total) bound = bound_set;
		else bound = bound_total;
		y[i] = model.addVar(log(alpha[i] * data.no_purchase[i]), log(bound), 0, GRB_CONTINUOUS, "y_" + to_string(i));
		u[i] = model.addVar(alpha[i] * data.no_purchase[i], bound, 0, GRB_CONTINUOUS, "u_" + to_string(i));
	}

	//cout << "Slack variables : z_i\n" << endl;
	GRBVar* z;
	z = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i) {
		double bound;
		double bound_total = calculate_bound_z_total(data, i);
		double bound_set = calculate_bound_z_set(data, i);
		if (bound_set > bound_total) bound = bound_set;
		else bound = bound_total;
		z[i] = model.addVar(-log(bound), -log(data.no_purchase[i]), 0, GRB_CONTINUOUS, "z_" + to_string(i));
	}

	//cout << "Slack variables : theta_i\n" << endl;
	GRBVar* theta;
	theta = new GRBVar[nCuts];
	for (int c = 0; c < nCuts; ++c)
		theta[c] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "theta_" + to_string(c));

	for (int i = 0; i < data.number_customers; ++i)
		model.addGenConstrExp(y[i], u[i]);

	for (int i = 0; i < data.number_customers; ++i) {
		GRBLinExpr sum_x;
		for (int j = 0; j < data.number_products; ++j)
			sum_x += (alpha[i] - data.revenue[i][j]) * x[j] * data.utilities[i][j];
		sum_x += alpha[i] * data.no_purchase[i];
		model.addConstr(u[i] >= sum_x);
	}

	//Total capacity constraint
	GRBLinExpr capacity;
	for (int j = 0; j < data.number_products; ++j) {
		capacity += data.cost[j] * x[j];
	}
	model.addConstr(capacity <= data.total_capacity, "ct_budget");

	//Set capacity constraints
	for (int s = 0; s < data.number_sets; ++s) {
		GRBLinExpr sum;
		for (int j = 0; j < data.number_products; ++j)
			if (data.in_set[j] == s)
				sum += x[j];
		model.addConstr(sum <= data.capacity_each_set, "ct_set_cap" + to_string(s));
	}

	//Objective
	GRBLinExpr obj;
	for (int c = 0; c < nCuts; ++c)
		obj += theta[c];
	model.setObjective(obj, GRB_MINIMIZE);

	auto time_before_cut = chrono::steady_clock::now();
	chrono::duration<double> before_cut = time_before_cut - start;

	double run_time = time_limit - before_cut.count();

	int num_iterative = 0;
	double stop_param = 1e-4;
	double sub_obj = 1.0;
	double obj_val_cplex = 0.0;
	double best_sub_obj = 0;

	//cout << "Add cut iteratively" << endl;
	while (sub_obj > obj_val_cplex + stop_param) {

		//compute gradient e^{y+z} at initial_x, initial_y, initial_z and set up constraints related to theta
		for (int c = 0; c < nCuts; ++c) {
			GRBLinExpr sum;
			for (int i = 0; i < group[c].size(); ++i) {
				sum += exp(initial_y[group[c][i]] + initial_z[group[c][i]]) * (1 + y[group[c][i]] - initial_y[group[c][i]] + z[group[c][i]] - initial_z[group[c][i]]);
			}
			model.addConstr(theta[c] >= sum, "ct_sub_gradient_e_" + to_string(c));
		}

		//compute gradient e^{-z} at initial_x, initial_y, initial_z and set up constraints related to e^{-z}
		for (int i = 0; i < data.number_customers; ++i) {
			GRBLinExpr sum = 0;
			for (int j = 0; j < data.number_products; ++j)
				sum += x[j] * data.utilities[i][j];
			sum += data.no_purchase[i];
			model.addConstr(exp(-initial_z[i]) * (1 - z[i] + initial_z[i]) <= sum, "ct_sub_gradient_z_" + to_string(i));
		}

		//solve
		num_iterative++;
		cout << "Remaining time: " << run_time << endl;

		model.write("cutting_plane.lp");
		model.set(GRB_DoubleParam_TimeLimit, run_time);
		//model.set(GRB_IntParam_MIPFocus, 3);
		model.set(GRB_IntParam_FuncPieces, 1);
		model.set(GRB_DoubleParam_FuncPieceLength, 1e-2);
		model.set(GRB_DoubleParam_FuncPieceError, 1e-3);
		//model.set(GRB_DoubleParam_MIPGap, 1e-3);
		//model.set(GRB_IntParam_OutputFlag, 0);
		//model.set(GRB_IntParam_Threads, 1);

		model.optimize();

		if (model.get(GRB_IntAttr_SolCount) > 0) {
			cout << "\nIteration " << num_iterative << endl;
			//update obj, variables
			obj_val_cplex = model.get(GRB_DoubleAttr_ObjVal);
			cout << "\nResult product list: " << endl;
			for (int j = 0; j < data.number_products; ++j)
				if (x[j].get(GRB_DoubleAttr_X) > 0.5) {
					initial_x[j] = 1;
					cout << j << " ";
				}
				else initial_x[j] = 0;
			cout << endl;

			for (int i = 0; i < data.number_customers; ++i) {
				initial_y[i] = y[i].get(GRB_DoubleAttr_X);
				initial_z[i] = z[i].get(GRB_DoubleAttr_X);
			}

			for (int c = 0; c < nCuts; ++c)
				initial_theta[c] = theta[c].get(GRB_DoubleAttr_X);

			//check the in equation related to theta_i and e^{y_i + z_i} for next iteration
			sub_obj = 0;

			for (int i = 0; i < data.number_customers; ++i)
				sub_obj += exp(initial_y[i] + initial_z[i]);

			cout << "Sub obj = " << std::setprecision(7) << fixed << sub_obj << endl;
			cout << "Gurobi obj = " << std::setprecision(7) << fixed << obj_val_cplex << endl;
			master_obj_val = calculate_master_obj(data, initial_x);
			cout << "Master obj = " << std::setprecision(7) << fixed << master_obj_val << endl;

			if (master_obj_val >= best_obj) {
				best_obj = master_obj_val;
				best_x = initial_x;
				best_sub_obj = obj_val_cplex;
			}

			//check time
			auto time_now = std::chrono::steady_clock::now(); //get now time
			std::chrono::duration<double> after_cut = time_now - start;
			cout << "Time now: " << after_cut.count() << endl;
			cout << "--- --- --- --- --- --- ---" << endl;

			if (after_cut.count() > time_limit) break;
			run_time = time_limit - after_cut.count();
		}
		else {
			cout << "Iteration " << num_iterative << ". No solution found..." << endl;
			auto end = chrono::steady_clock::now();
			chrono::duration<double> elapsed_seconds = end - start;
			time_for_solve = elapsed_seconds.count();
			break;
		}
	}
	auto end = chrono::steady_clock::now();
	chrono::duration<double> total_time = end - start;
	time_for_solve = total_time.count();

	cout << "\nObjective value: " << setprecision(5) << best_obj << endl;
	cout << "Solution: ";
	for (int j = 0; j < data.number_products; ++j)
		if (best_x[j] == 1)
			cout << j << " ";
	cout << "\nTotal time: " << time_for_solve << " seconds" << endl;

	ofstream report_results(out_res_csv, ofstream::out);
	report_results.precision(10);
	report_results << best_sub_obj << " " << best_obj << " " << time_for_solve << endl;
	for (int j = 0; j < data.number_products; ++j)
		if (best_x[j] == 1)
			report_results << j << " ";
	report_results.close();
}