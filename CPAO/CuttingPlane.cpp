#include "CuttingPlane.h"
#include <chrono>
#include <algorithm>
#include <sstream>
#include <cstdlib>
#include <cassert>

CuttingPlane::CuttingPlane() {

}

CuttingPlane::CuttingPlane(Data data, double time_limit, string outfile) {
	this->data = data;
	this->time_limit = time_limit;
	this->out_res_csv = outfile;
}

vector<double> CuttingPlane::calculate_denominator(Data data, vector<int> x) {
	vector<double> denominator(data.number_customers, 0);
	for (int i = 0; i < data.number_customers; ++i) {
		denominator[i] += data.no_purchase[i];
		for (int j = 0; j < data.number_products; ++j)
			denominator[i] += x[j] * data.utilities[i][j];
	}

	return denominator;
}

double CuttingPlane::calculate_z(Data data, vector<double> alpha, vector<double> denominator) {
	double z = 0;
	for (int i = 0; i < data.number_customers; ++i)
		z += alpha[i] * data.no_purchase[i] / denominator[i];

	return z;
}

double CuttingPlane::calculate_optimal_bound_denominator(Data data, int i) {
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

	GRBLinExpr obj;
	for (int j = 0; j < data.number_products; ++j)
		obj += data.utilities[i][j] * x[j];
	model.setObjective(obj, GRB_MAXIMIZE);

	model.set(GRB_IntParam_OutputFlag, 0);

	model.optimize();

	return model.get(GRB_DoubleAttr_ObjVal) + data.no_purchase[i];
}

double  CuttingPlane::calculate_original_obj(Data data, vector<int> x, vector<double> alpha) {
	double obj = 0;
	for (int i = 0; i < data.number_customers; ++i) {
		double ts = alpha[i] * data.no_purchase[i], ms = data.no_purchase[i];
		for (int j = 0; j < data.number_products; ++j) {
			ts += (alpha[i] - data.revenue[i][j]) * x[j] * data.utilities[i][j];
			ms += x[j] * data.utilities[i][j];
		}
		obj += data.fraction[i] * ts / ms;
	}

	return obj;
}

double CuttingPlane::calculate_master_obj(Data data, vector<int> x) {
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

double CuttingPlane::calculate_original_obj_tmp(Data data, vector<int> x, vector<double> alpha, int candidate) {
	double obj = 0;
	for (int i = 0; i < data.number_customers; ++i) {
		double ts = alpha[i] * data.no_purchase[i] + (alpha[i] - data.revenue[i][candidate]) * data.utilities[i][candidate],
			ms = data.no_purchase[i] + data.utilities[i][candidate];
		for (int j = 0; j < data.number_products; ++j) {
			ts += (alpha[i] - data.revenue[i][j]) * x[j] * data.utilities[i][j];
			ms += x[j] * data.utilities[i][j];
		}

		obj += data.fraction[i] * ts / ms;
	}

	return obj;
}

vector<int> CuttingPlane::greedy(Data data, vector<double> alpha) {
	auto start = chrono::steady_clock::now();
	vector<int> chosen(data.number_products, 0);
	double obj = calculate_original_obj(data, chosen, alpha);
	//cout << "Obj: " << obj << endl;
	vector<double> set_capacity(data.number_sets, 0);

	while (true) {
		double obj_tmp = 0;
		int inserted_product = -1;
		double min = 999999;
		for (int j = 0; j < data.number_products; ++j)
			if (chosen[j] == 0) {
				int check_in_set = 0;
				int check_valid = 0;
				for (int s = 0; s < data.number_sets; ++s) {
					if (data.in_set[j][s] == 1) {
						check_in_set++;
						if (set_capacity[s] + data.cost[j] <= data.capacity_each_set)
							check_valid++;
						else break;
					}
				}
				if (check_in_set != 0 && check_in_set == check_valid) {
					obj_tmp = calculate_original_obj_tmp(data, chosen, alpha, j);
					if (min > obj_tmp) {
						min = obj_tmp;
						inserted_product = j;
					}
				}
				if (check_in_set == 0) {
					obj_tmp = calculate_original_obj_tmp(data, chosen, alpha, j);
					if (min > obj_tmp) {
						min = obj_tmp;
						inserted_product = j;
					}
				}
			}

		if (inserted_product != -1) {
			chosen[inserted_product] = 1;
			obj = calculate_original_obj(data, chosen, alpha);
			//cout << "Best-Insertion: " << inserted_product << " " << obj << endl;
			for (int s = 0; s < data.number_sets; ++s)
				if (data.in_set[inserted_product][s] == 1)
					set_capacity[s] += data.cost[inserted_product];
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
	cout << "Master obj = " << calculate_master_obj(data, chosen) << endl;
	cout << "Time: " << greedy_time.count() << endl;

	return chosen;
}

vector<int> CuttingPlane::greedy_general(Data data, vector<double> alpha) {
	auto start = chrono::steady_clock::now();
	vector<int> chosen(data.number_products, 0);
	double obj = calculate_original_obj(data, chosen, alpha);
	//cout << "Obj: " << obj << endl;
	vector<double> set_capacity(5, 0);

	double cost = 0;
	while (true) {
		double obj_tmp = 0;
		int inserted_product = -1;
		double min = 999999;
		for (int j = 0; j < data.number_products; ++j)
			if (chosen[j] == 0 && cost + data.fraction2[j] <= data.capacity_each_set) {
				int check_valid = 0;
				for (int s = 0; s < 5; ++s) {
					if (data.in_set[j][s] == 1) {
						if (set_capacity[s] + 1 <= data.sub_capacity_each_set)
							check_valid++;
						else break;
					}
				}
				if (check_valid == 1) {
					obj_tmp = calculate_original_obj_tmp(data, chosen, alpha, j);
					if (min > obj_tmp) {
						min = obj_tmp;
						inserted_product = j;
					}
				}
			}

		if (inserted_product != -1) {
			chosen[inserted_product] = 1;
			obj = calculate_original_obj(data, chosen, alpha);
			//cout << "Best-Insertion: " << inserted_product << " " << obj << endl;
			for (int s = 0; s < 5; ++s)
				if (data.in_set[inserted_product][s] == 1)
					set_capacity[s]++;
			cost += data.fraction2[inserted_product];
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
	cout << "Master obj = " << calculate_master_obj(data, chosen) << endl;
	cout << "Time: " << greedy_time.count() << endl;

	return chosen;
}

vector<vector<double>> CuttingPlane::calculate_bound_y_in(Data data) {
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
			//while (count < calculate_bound_y(data)) {
			while (count < data.capacity_each_set - 1) {
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

vector<vector<double>> CuttingPlane::calculate_bound_y_notin(Data data) {
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

vector<vector<double>> CuttingPlane::subset_bound_y_in(Data data) {
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

vector<vector<double>> CuttingPlane::subset_bound_y_notin(Data data) {
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

void CuttingPlane::solve_milp(Data data) {
	//auto start = chrono::steady_clock::now();

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

	auto start = chrono::steady_clock::now();

	//Calculate initial_x
	//vector<int> initial_x(data.number_products, 0);
	vector<int> initial_x = greedy(data, alpha);
	//vector<int> initial_x = greedy_general(data, alpha);
	vector<double> initial_y(data.number_customers, 0);
	vector<double> initial_theta(data.number_products, 0);

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

	vector<double> upper_bound_theta(data.number_products, 0);
	vector<double> lower_bound_theta(data.number_products, 0);
	for (int j = 0; j < data.number_products; ++j)
		for (int i = 0; i < data.number_customers; ++i) {
			lower_bound_theta[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / upper_bound_denominator[i];
			upper_bound_theta[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / data.no_purchase[i];
		}

	GRBVar* theta;
	theta = new GRBVar[data.number_products];
	for (int j = 0; j < data.number_products; ++j)
		theta[j] = model.addVar(lower_bound_theta[j], upper_bound_theta[j], 0, GRB_CONTINUOUS, "theta_" + to_string(j));

	GRBVar* z;
	z = new GRBVar[data.number_products];
	for (int j = 0; j < data.number_products; ++j)
		z[j] = model.addVar(0, upper_bound_theta[j], 0, GRB_CONTINUOUS, "z_" + to_string(j));

	//Set capacity constraints
	for (int s = 0; s < data.number_sets; ++s) {
		GRBLinExpr sum;
		for (int j = 0; j < data.number_products; ++j)
			if (data.in_set[j][s] == 1)
				sum += data.cost[j] * x[j];
		model.addConstr(sum <= data.capacity_each_set, "ct_set_cap" + to_string(s));
	}

	////General
	//for (int s = 0; s < 10; ++s) {
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

	//cout << "Constraints related to y and theta\n" << endl;
	for (int j = 0; j < data.number_products; ++j) {
		GRBLinExpr sum;
		for (int i = 0; i < data.number_customers; ++i)
			sum += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * y[i];
		model.addConstr(theta[j] >= sum);
	}

	//cout << "z_j = theta_j * x_j\n" << endl;
	for (int j = 0; j < data.number_products; ++j) {
		//model.addConstr(theta[j] - z[j] <= upper_bound_theta[j] - upper_bound_theta[j] * x[j]);
		//model.addConstr(z[j] <= upper_bound_theta[j] * x[j]);
		model.addConstr(z[j] <= theta[j]);
	}

	vector<double> boundx1(data.number_products, 0);
	vector<double> boundx0(data.number_products, 0);
	vector<double> boundin(data.number_products, 0);
	vector<double> boundnotin(data.number_products, 0);
	for (int j = 0; j < data.number_products; ++j)
		for (int i = 0; i < data.number_customers; ++i) {
			boundx1[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + data.utilities[i][j]);
			boundx0[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / data.no_purchase[i];
			boundin[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + bound_in[i][j]);
			boundnotin[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + bound_notin[i][j]);
		}

	//McCornick constraints
	for (int j = 0; j < data.number_products; ++j) {
		model.addConstr(z[j] <= x[j] * boundx1[j]);

		//if (bound_in[i][j] >= subset_bound_in[i][j])
		model.addConstr(z[j] >= x[j] * boundin[j]);
		//else
			//model.addConstr(z[j] >= x[j] * (1 / (data.no_purchase[i] + subset_bound_in[i][j])));

		//if (bound_notin[i][j] >= subset_bound_notin[i][j])
		model.addConstr(z[j] <= theta[j] - (1 - x[j]) * boundnotin[j]);
		//else
			//model.addConstr(z[j] <= y[i] - (1 - x[j]) * (1 / (data.no_purchase[i] + subset_bound_notin[i][j])));

		model.addConstr(z[j] >= theta[j] - (1 - x[j]) * boundx0[j]);
	}

	//Objective
	GRBLinExpr obj;
	for (int i = 0; i < data.number_customers; ++i)
		obj += data.fraction[i] * data.no_purchase[i] * alpha[i] * y[i];
	for (int j = 0; j < data.number_products; ++j)
		obj += z[j];

	model.setObjective(obj, GRB_MINIMIZE);

	auto time_before_cut = chrono::steady_clock::now();
	chrono::duration<double> before_cut = time_before_cut - start;

	double run_time = time_limit - before_cut.count();

	int num_iterative = 0;
	double stop_param = 1e-4;
	double sub_obj = 1.0;
	double obj_val_cplex = 0.0;
	double best_sub_obj = 0;
	vector<int> best_x = initial_x;
	double best_obj = calculate_master_obj(data, best_x);

	while (sub_obj > obj_val_cplex + stop_param) {
		vector<double> initial_denominator = calculate_denominator(data, initial_x);

		//cout << "Outer-cuts y\n" << endl;
		vector<double> partitial_y(data.number_customers, 0);
		for (int i = 0; i < data.number_customers; ++i)
			partitial_y[i] = 1 / initial_denominator[i];

		vector<vector<double>> subgradient_y(data.number_customers);
		for (int i = 0; i < data.number_customers; ++i)
			subgradient_y[i].resize(data.number_products, 0);

		for (int i = 0; i < data.number_customers; ++i)
			for (int j = 0; j < data.number_products; ++j)
				subgradient_y[i][j] -= data.utilities[i][j] / (initial_denominator[i] * initial_denominator[i]);

		GRBConstr* ct_y = new GRBConstr[data.number_customers];

		for (int i = 0; i < data.number_customers; ++i)
			/*if (initial_y[i] < partitial_y[i])*/ {
			GRBLinExpr grad_y;
			for (int j = 0; j < data.number_products; ++j)
				grad_y += subgradient_y[i][j] * (x[j] - initial_x[j]);

			ct_y[i] = model.addConstr(y[i] >= partitial_y[i] + grad_y, "ct_sub_gradient_y_" + to_string(i));
		}

		//cout << "Outer-cuts theta\n" << endl;
		vector<double> partitial_theta(data.number_products, 0);
		for (int j = 0; j < data.number_products; ++j)
			for (int i = 0; i < data.number_customers; ++i)
				partitial_theta[j] +=
				data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j])
				/ initial_denominator[i];

		vector<vector<double>> subgradient_theta(data.number_products);
		for (int j = 0; j < data.number_products; ++j)
			subgradient_theta[j].resize(data.number_products, 0);

		for (int j = 0; j < data.number_products; ++j)
			for (int j1 = 0; j1 < data.number_products; ++j1)
				for (int i = 0; i < data.number_customers; ++i)
					subgradient_theta[j][j1] -=
					data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1]
					/ (initial_denominator[i] * initial_denominator[i]);

		GRBConstr* ct_theta = new GRBConstr[data.number_products];

		for (int j = 0; j < data.number_products; ++j)
			/*if(initial_theta[j] < partitial_theta[j])*/ {
			GRBLinExpr grad_theta;
			for (int j1 = 0; j1 < data.number_products; ++j1)
				grad_theta += subgradient_theta[j][j1] * (x[j1] - initial_x[j1]);

			ct_theta[j] = model.addConstr(theta[j] >= partitial_theta[j] + grad_theta, "ct_sub_gradient_theta_" + to_string(j));
		}

		////cout << "Calculate total utility\n" << endl;
		//vector<double> sum_uti_customer(data.number_customers, 0);
		//for (int i = 0; i < data.number_customers; ++i) {
		//	sum_uti_customer[i] += data.no_purchase[i];
		//	for (int j = 0; j < data.number_products; ++j)
		//		sum_uti_customer[i] += data.utilities[i][j];
		//}

		////cout << "Submodular-cuts y\n" << endl;
		//for (int i = 0; i < data.number_customers; ++i) {
		//	//if (initial_y[i] < partitial_y[i]) {
		//		GRBLinExpr submodular_cut_a_z, submodular_cut_b_z;
		//		for (int j = 0; j < data.number_products; ++j)
		//			if (initial_x[j] == 1) {
		//				submodular_cut_a_z += (1 - x[j]) * data.utilities[i][j] /
		//					(sum_uti_customer[i] * (sum_uti_customer[i] - data.utilities[i][j]));
		//				submodular_cut_b_z += (1 - x[j]) * data.utilities[i][j] /
		//					(initial_denominator[i] * (initial_denominator[i] - data.utilities[i][j]));
		//			}
		//			else {
		//				submodular_cut_a_z -= x[j] * data.utilities[i][j] /
		//					(initial_denominator[i] * (initial_denominator[i] + data.utilities[i][j]));
		//				submodular_cut_b_z -= x[j] * data.utilities[i][j] /
		//					(data.no_purchase[i] * (data.no_purchase[i] + data.utilities[i][j]));
		//			}

		//		submodular_cut_a_z += partitial_y[i];
		//		model.addConstr(y[i] >= submodular_cut_a_z, "ct_sub_modular_a_z_" + to_string(i));
		//		submodular_cut_b_z += partitial_y[i];
		//		model.addConstr(y[i] >= submodular_cut_b_z, "ct_sub_modular_b_z_" + to_string(i));
		//	//}
		//}

		////cout << "Submodular-cuts theta\n" << endl;
		//for (int j = 0; j < data.number_products; ++j) {
		//	//if (initial_theta[j] < partitial_theta[j]) {
		//		GRBLinExpr submodular_cut_a_z, submodular_cut_b_z;
		//		for (int j1 = 0; j1 < data.number_products; ++j1)
		//			if (initial_x[j1] == 1) {
		//				for (int i = 0; i < data.number_customers; ++i) {
		//					submodular_cut_a_z +=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(sum_uti_customer[i] * (sum_uti_customer[i] - data.utilities[i][j1]));
		//					submodular_cut_b_z +=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(initial_denominator[i] * (initial_denominator[i] - data.utilities[i][j1]));
		//				}
		//			}
		//			else {
		//				for (int i = 0; i < data.number_customers; ++i) {
		//					submodular_cut_a_z -=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(initial_denominator[i] * (initial_denominator[i] + data.utilities[i][j1]));
		//					submodular_cut_b_z -=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(data.no_purchase[i] * (data.no_purchase[i] + data.utilities[i][j1]));
		//				}
		//			}

		//		submodular_cut_a_z += partitial_theta[j];
		//		model.addConstr(theta[j] >= submodular_cut_a_z, "ct_sub_modular_a_z_" + to_string(j));
		//		submodular_cut_b_z += partitial_theta[j];
		//		model.addConstr(theta[j] >= submodular_cut_b_z, "ct_sub_modular_b_z_" + to_string(j));
		//	//}
		//}

		//solve
		num_iterative++;
		cout << "Remaining time: " << run_time << endl;

		model.write("submodular.lp");
		model.set(GRB_DoubleParam_TimeLimit, run_time);
		model.set(GRB_IntParam_Threads, 8);
		//model.set(GRB_DoubleParam_MIPGap, 1e-3);
		//model.set(GRB_IntParam_MIQCPMethod, 1);
		//model.set(GRB_IntParam_OutputFlag, 0);

		model.optimize();
		for (int i = 0; i < data.number_customers; ++i)
			if (initial_y[i] >= partitial_y[i])
				model.remove(ct_y[i]);
		for (int j = 0; j < data.number_products; ++j)
			if (initial_theta[j] >= partitial_theta[j])
				model.remove(ct_theta[j]);

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

			for (int i = 0; i < data.number_customers; ++i)
				initial_y[i] = y[i].get(GRB_DoubleAttr_X);
			for (int j = 0; j < data.number_products; ++j)
				initial_theta[j] = theta[j].get(GRB_DoubleAttr_X);

			initial_denominator = calculate_denominator(data, initial_x);

			//check the in equation related to theta_j, x_j and y_j for next iteration
			sub_obj = 0;
			for (int i = 0; i < data.number_customers; ++i)
				sub_obj += data.fraction[i] * data.no_purchase[i] * alpha[i] / initial_denominator[i];
			for (int j = 0; j < data.number_products; ++j)
				for (int i = 0; i < data.number_customers; ++i)
					sub_obj += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * initial_x[j]
					/ initial_denominator[i];

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

void CuttingPlane::solve_multicut_milp(Data data, int number_cuts) {
	//auto start = chrono::steady_clock::now();

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

	auto start = chrono::steady_clock::now();

	//Calculate initial_x
	//vector<int> initial_x(data.number_products, 0);
	vector<int> initial_x = greedy(data, alpha);
	//vector<int> initial_x = greedy_general(data, alpha);
	double greedy_obj = calculate_master_obj(data, initial_x);
	vector<double> initial_y(data.number_customers, 0);
	vector<double> initial_theta(data.number_products, 0);
	vector<double> initial_cut(number_cuts, 0);

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

	vector<double> upper_bound_theta(data.number_products, 0);
	vector<double> lower_bound_theta(data.number_products, 0);
	for (int j = 0; j < data.number_products; ++j)
		for (int i = 0; i < data.number_customers; ++i) {
			lower_bound_theta[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / upper_bound_denominator[i];
			upper_bound_theta[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / data.no_purchase[i];
		}

	GRBVar* theta;
	theta = new GRBVar[data.number_products];
	for (int j = 0; j < data.number_products; ++j)
		theta[j] = model.addVar(lower_bound_theta[j], upper_bound_theta[j], 0, GRB_CONTINUOUS, "theta_" + to_string(j));

	GRBVar* z;
	z = new GRBVar[data.number_products];
	for (int j = 0; j < data.number_products; ++j)
		z[j] = model.addVar(0, upper_bound_theta[j], 0, GRB_CONTINUOUS, "z_" + to_string(j));

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

	//cout << "Constraints related to y and theta\n" << endl;
	for (int j = 0; j < data.number_products; ++j) {
		GRBLinExpr sum;
		for (int i = 0; i < data.number_customers; ++i)
			sum += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * y[i];
		model.addConstr(theta[j] >= sum);
	}

	//cout << "z_j = theta_j * x_j\n" << endl;
	for (int j = 0; j < data.number_products; ++j) {
		//model.addConstr(theta[j] - z[j] <= upper_bound_theta[j] - upper_bound_theta[j] * x[j]);
		//model.addConstr(z[j] <= upper_bound_theta[j] * x[j]);
		model.addConstr(z[j] <= theta[j]);
	}

	vector<double> boundx1(data.number_products, 0);
	vector<double> boundx0(data.number_products, 0);
	vector<double> boundin(data.number_products, 0);
	vector<double> boundnotin(data.number_products, 0);
	for (int j = 0; j < data.number_products; ++j)
		for (int i = 0; i < data.number_customers; ++i) {
			boundx1[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + data.utilities[i][j]);
			boundx0[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / data.no_purchase[i];
			boundin[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + bound_in[i][j]);
			boundnotin[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + bound_notin[i][j]);
		}

	//vector<double> subset_boundin(data.number_products, 0);
	//vector<double> subset_boundnotin(data.number_products, 0);
	//for (int j = 0; j < data.number_products; ++j)
	//	for (int i = 0; i < data.number_customers; ++i) {
	//		subset_boundin[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + subset_bound_in[i][j]);
	//		subset_boundnotin[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + subset_bound_notin[i][j]);
	//	}

	//McCornick constraints
	for (int j = 0; j < data.number_products; ++j) {
		model.addConstr(z[j] <= x[j] * boundx1[j]);

		//if (boundin[j] >= subset_boundin[j])
		model.addConstr(z[j] >= x[j] * boundin[j]);
		//else
			//model.addConstr(z[j] >= x[j] * subset_boundin[j]);

		//if (boundnotin[j] >= subset_boundnotin[j])
		model.addConstr(z[j] <= theta[j] - (1 - x[j]) * boundnotin[j]);
		//else
			//model.addConstr(z[j] <= theta[j] - (1 - x[j]) * subset_boundnotin[j]);

		model.addConstr(z[j] >= theta[j] - (1 - x[j]) * boundx0[j]);
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
	for (int j = 0; j < data.number_products; ++j)
		obj += z[j];

	model.setObjective(obj, GRB_MINIMIZE);

	auto time_before_cut = chrono::steady_clock::now();
	chrono::duration<double> before_cut = time_before_cut - start;

	double run_time = time_limit - before_cut.count();

	int num_iterative = 0;
	double stop_param = 1e-4;
	double sub_obj = 1.0;
	double obj_val_cplex = 0.0;
	double best_sub_obj = 0;
	vector<int> best_x = initial_x;
	double best_obj = calculate_master_obj(data, best_x);

	while (sub_obj > obj_val_cplex + stop_param) {
		//cout << "Calculate total utility\n" << endl;
		vector<double> sum_uti_customer(data.number_customers, 0);
		for (int i = 0; i < data.number_customers; ++i) {
			sum_uti_customer[i] += data.no_purchase[i];
			for (int j = 0; j < data.number_products; ++j)
				sum_uti_customer[i] += data.utilities[i][j];
		}

		vector<double> initial_denominator = calculate_denominator(data, initial_x);

		//cout << "Outer-cuts y\n" << endl;
		vector<double> partitial_y(data.number_customers, 0);
		for (int i = 0; i < data.number_customers; ++i)
			partitial_y[i] = 1 / initial_denominator[i];

		vector<vector<double>> subgradient_y(data.number_customers);
		for (int i = 0; i < data.number_customers; ++i)
			subgradient_y[i].resize(data.number_products, 0);

		for (int i = 0; i < data.number_customers; ++i)
			for (int j = 0; j < data.number_products; ++j)
				subgradient_y[i][j] -= data.utilities[i][j] / (initial_denominator[i] * initial_denominator[i]);

		for (int c = 0; c < number_cuts; ++c) {
			GRBLinExpr grad_y;
			double part_y = 0;
			for (int i = 0; i < data.number_customers; ++i)
				if (i % number_cuts == c) {
					part_y += data.fraction[i] * data.no_purchase[i] * alpha[i] * partitial_y[i];
					for (int j = 0; j < data.number_products; ++j)
						grad_y += data.fraction[i] * data.no_purchase[i] * alpha[i] * subgradient_y[i][j] * (x[j] - initial_x[j]);
				}
			//if(initial_cut[c] < part_y)
			model.addConstr(cut_customers[c] >= part_y + grad_y, "ct_sub_gradient_cut_" + to_string(c));
		}

		//cout << "Outer-cuts theta\n" << endl;
		vector<double> partitial_theta(data.number_products, 0);
		for (int j = 0; j < data.number_products; ++j)
			for (int i = 0; i < data.number_customers; ++i)
				partitial_theta[j] +=
				data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j])
				/ initial_denominator[i];

		vector<vector<double>> subgradient_theta(data.number_products);
		for (int j = 0; j < data.number_products; ++j)
			subgradient_theta[j].resize(data.number_products, 0);

		for (int j = 0; j < data.number_products; ++j)
			for (int j1 = 0; j1 < data.number_products; ++j1)
				for (int i = 0; i < data.number_customers; ++i)
					subgradient_theta[j][j1] -=
					data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1]
					/ (initial_denominator[i] * initial_denominator[i]);


		for (int j = 0; j < data.number_products; ++j)
			if (initial_theta[j] < partitial_theta[j]) {
				GRBLinExpr grad_theta;
				for (int j1 = 0; j1 < data.number_products; ++j1)
					grad_theta += subgradient_theta[j][j1] * (x[j1] - initial_x[j1]);

				model.addConstr(theta[j] >= partitial_theta[j] + grad_theta, "ct_sub_gradient_theta_" + to_string(j));
			}

		////cout << "Submodular-cuts y\n" << endl;
		//for (int c = 0; c < number_cuts; ++c) {
		//	GRBLinExpr submodular_cut_a_z, submodular_cut_b_z;
		//	double part_y = 0;
		//	for (int i = 0; i < data.number_customers; ++i)
		//		if (i % number_cuts == c) {
		//			part_y += partitial_y[i];
		//			for (int j = 0; j < data.number_products; ++j)
		//				if (initial_x[j] == 1) {
		//					submodular_cut_a_z += (1 - x[j]) * data.utilities[i][j] /
		//						(sum_uti_customer[i] * (sum_uti_customer[i] - data.utilities[i][j]));
		//					submodular_cut_b_z += (1 - x[j]) * data.utilities[i][j] /
		//						(initial_denominator[i] * (initial_denominator[i] - data.utilities[i][j]));
		//				}
		//				else {
		//					submodular_cut_a_z -= x[j] * data.utilities[i][j] /
		//						(initial_denominator[i] * (initial_denominator[i] + data.utilities[i][j]));
		//					submodular_cut_b_z -= x[j] * data.utilities[i][j] /
		//						(data.no_purchase[i] * (data.no_purchase[i] + data.utilities[i][j]));
		//				}
		//	}
		//	if (initial_cut[c] < part_y) {
		//		submodular_cut_a_z += part_y;
		//		model.addConstr(cut_customers[c] >= submodular_cut_a_z, "ct_sub_modular_a_z_" + to_string(c));
		//		submodular_cut_b_z += part_y;
		//		model.addConstr(cut_customers[c] >= submodular_cut_b_z, "ct_sub_modular_b_z_" + to_string(c));
		//	}
		//}

		////cout << "Submodular-cuts theta\n" << endl;
		//for (int j = 0; j < data.number_products; ++j) {
		//	if (initial_theta[j] < partitial_theta[j]) {
		//		GRBLinExpr submodular_cut_a_z, submodular_cut_b_z;
		//		for (int j1 = 0; j1 < data.number_products; ++j1)
		//			if (initial_x[j1] == 1) {
		//				for (int i = 0; i < data.number_customers; ++i) {
		//					submodular_cut_a_z +=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(sum_uti_customer[i] * (sum_uti_customer[i] - data.utilities[i][j1]));
		//					submodular_cut_b_z +=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(initial_denominator[i] * (initial_denominator[i] - data.utilities[i][j1]));
		//				}
		//			}
		//			else {
		//				for (int i = 0; i < data.number_customers; ++i) {
		//					submodular_cut_a_z -=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(initial_denominator[i] * (initial_denominator[i] + data.utilities[i][j1]));
		//					submodular_cut_b_z -=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(data.no_purchase[i] * (data.no_purchase[i] + data.utilities[i][j1]));
		//				}
		//			}

		//		submodular_cut_a_z += partitial_theta[j];
		//		model.addConstr(theta[j] >= submodular_cut_a_z, "ct_sub_modular_a_z_" + to_string(j));
		//		submodular_cut_b_z += partitial_theta[j];
		//		model.addConstr(theta[j] >= submodular_cut_b_z, "ct_sub_modular_b_z_" + to_string(j));
		//	}
		//}

		//solve
		num_iterative++;
		cout << "Remaining time: " << run_time << endl;

		model.write("submodular.lp");
		model.set(GRB_DoubleParam_TimeLimit, run_time);
		model.set(GRB_IntParam_Threads, 8);
		//model.set(GRB_DoubleParam_MIPGap, 1e-3);
		//model.set(GRB_IntParam_MIQCPMethod, 1);
		//model.set(GRB_IntParam_OutputFlag, 0);

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

			for (int i = 0; i < data.number_customers; ++i)
				initial_y[i] = y[i].get(GRB_DoubleAttr_X);
			for (int j = 0; j < data.number_products; ++j)
				initial_theta[j] = theta[j].get(GRB_DoubleAttr_X);
			for (int c = 0; c < number_cuts; ++c)
				initial_cut[c] = cut_customers[c].get(GRB_DoubleAttr_X);

			initial_denominator = calculate_denominator(data, initial_x);

			//check the in equation related to theta_j, x_j and y_j for next iteration
			sub_obj = 0;
			for (int i = 0; i < data.number_customers; ++i)
				sub_obj += data.fraction[i] * data.no_purchase[i] * alpha[i] / initial_denominator[i];
			for (int j = 0; j < data.number_products; ++j)
				for (int i = 0; i < data.number_customers; ++i)
					sub_obj += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * initial_x[j]
					/ initial_denominator[i];

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
	report_results << best_sub_obj << " " << best_obj << " " << greedy_obj << " " << time_for_solve << endl;
	for (int j = 0; j < data.number_products; ++j)
		if (best_x[j] == 1)
			report_results << j << " ";
	report_results.close();
}

void CuttingPlane::solve_bi(Data data) {
	//auto start = chrono::steady_clock::now();

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

	auto start = chrono::steady_clock::now();

	//Calculate initial_x
	//vector<int> initial_x(data.number_products, 0);
	vector<int> initial_x = greedy(data, alpha);
	//vector<int> initial_x = greedy_general(data, alpha);
	vector<double> initial_y(data.number_customers, 0);
	vector<double> initial_theta(data.number_products, 0);

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

	vector<double> upper_bound_theta(data.number_products, 0);
	vector<double> lower_bound_theta(data.number_products, 0);
	for (int j = 0; j < data.number_products; ++j)
		for (int i = 0; i < data.number_customers; ++i) {
			lower_bound_theta[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / upper_bound_denominator[i];
			upper_bound_theta[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / data.no_purchase[i];
		}

	GRBVar* theta;
	theta = new GRBVar[data.number_products];
	for (int j = 0; j < data.number_products; ++j)
		theta[j] = model.addVar(lower_bound_theta[j], upper_bound_theta[j], 0, GRB_CONTINUOUS, "theta_" + to_string(j));

	//Set capacity constraints
	for (int s = 0; s < data.number_sets; ++s) {
		GRBLinExpr sum;
		for (int j = 0; j < data.number_products; ++j)
			if (data.in_set[j][s] == 1)
				sum += data.cost[j] * x[j];
		model.addConstr(sum <= data.capacity_each_set, "ct_set_cap" + to_string(s));
	}

	////General
	//for (int s = 0; s < 10; ++s) {
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

	//cout << "Constraints related to y and theta\n" << endl;
	for (int j = 0; j < data.number_products; ++j) {
		GRBLinExpr sum;
		for (int i = 0; i < data.number_customers; ++i)
			sum += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * y[i];
		model.addConstr(theta[j] >= sum);
	}

	vector<double> boundx1(data.number_products, 0);
	vector<double> boundx0(data.number_products, 0);
	vector<double> boundin(data.number_products, 0);
	vector<double> boundnotin(data.number_products, 0);
	for (int j = 0; j < data.number_products; ++j)
		for (int i = 0; i < data.number_customers; ++i) {
			boundx1[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + data.utilities[i][j]);
			boundx0[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / data.no_purchase[i];
			boundin[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + bound_in[i][j]);
			boundnotin[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + bound_notin[i][j]);
		}

	//McCornick constraints
	for (int j = 0; j < data.number_products; ++j) {
		model.addQConstr(theta[j] * x[j] <= x[j] * boundx1[j]);

		//if (bound_in[i][j] >= subset_bound_in[i][j])
		model.addQConstr(theta[j] * x[j] >= x[j] * boundin[j]);
		//else
			//model.addConstr(z[j] >= x[j] * (1 / (data.no_purchase[i] + subset_bound_in[i][j])));

		//if (bound_notin[i][j] >= subset_bound_notin[i][j])
		model.addQConstr(theta[j] * x[j] <= theta[j] - (1 - x[j]) * boundnotin[j]);
		//else
			//model.addConstr(z[j] <= y[i] - (1 - x[j]) * (1 / (data.no_purchase[i] + subset_bound_notin[i][j])));

		model.addQConstr(theta[j] * x[j] >= theta[j] - (1 - x[j]) * boundx0[j]);
	}

	//Objective
	GRBQuadExpr obj;
	for (int i = 0; i < data.number_customers; ++i)
		obj += data.fraction[i] * data.no_purchase[i] * alpha[i] * y[i];
	for (int j = 0; j < data.number_products; ++j)
		obj += theta[j] * x[j];

	model.setObjective(obj, GRB_MINIMIZE);

	auto time_before_cut = chrono::steady_clock::now();
	chrono::duration<double> before_cut = time_before_cut - start;

	double run_time = time_limit - before_cut.count();

	int num_iterative = 0;
	double stop_param = 1e-4;
	double sub_obj = 1.0;
	double obj_val_cplex = 0.0;
	double best_sub_obj = 0;
	vector<int> best_x = initial_x;
	double best_obj = calculate_master_obj(data, best_x);

	while (sub_obj > obj_val_cplex + stop_param) {
		vector<double> initial_denominator = calculate_denominator(data, initial_x);

		//cout << "Outer-cuts y\n" << endl;
		vector<double> partitial_y(data.number_customers, 0);
		for (int i = 0; i < data.number_customers; ++i)
			partitial_y[i] = 1 / initial_denominator[i];

		vector<vector<double>> subgradient_y(data.number_customers);
		for (int i = 0; i < data.number_customers; ++i)
			subgradient_y[i].resize(data.number_products, 0);

		for (int i = 0; i < data.number_customers; ++i)
			for (int j = 0; j < data.number_products; ++j)
				subgradient_y[i][j] -= data.utilities[i][j] / (initial_denominator[i] * initial_denominator[i]);

		GRBConstr* ct_y = new GRBConstr[data.number_customers];

		for (int i = 0; i < data.number_customers; ++i)
			/*if (initial_y[i] < partitial_y[i])*/ {
			GRBLinExpr grad;
			for (int j = 0; j < data.number_products; ++j)
				grad += subgradient_y[i][j] * (x[j] - initial_x[j]);

			ct_y[i] = model.addConstr(y[i] >= partitial_y[i] + grad, "ct_sub_gradient_y_" + to_string(i));
		}

		//cout << "Outer-cuts theta\n" << endl;
		vector<double> partitial_theta(data.number_products, 0);
		for (int j = 0; j < data.number_products; ++j)
			for (int i = 0; i < data.number_customers; ++i)
				partitial_theta[j] +=
				data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j])
				/ initial_denominator[i];

		vector<vector<double>> subgradient_theta(data.number_products);
		for (int j = 0; j < data.number_products; ++j)
			subgradient_theta[j].resize(data.number_products, 0);

		for (int j = 0; j < data.number_products; ++j)
			for (int j1 = 0; j1 < data.number_products; ++j1)
				for (int i = 0; i < data.number_customers; ++i)
					subgradient_theta[j][j1] -=
					data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1]
					/ (initial_denominator[i] * initial_denominator[i]);

		GRBConstr* ct_theta = new GRBConstr[data.number_products];

		for (int j = 0; j < data.number_products; ++j)
			/*if (initial_theta[j] < partitial_theta[j])*/ {
			GRBLinExpr grad_theta;
			for (int j1 = 0; j1 < data.number_products; ++j1)
				grad_theta += subgradient_theta[j][j1] * (x[j1] - initial_x[j1]);

			ct_theta[j] = model.addConstr(theta[j] >= partitial_theta[j] + grad_theta, "ct_sub_gradient_theta_" + to_string(j));
		}

		////cout << "Calculate total utility\n" << endl;
		//vector<double> sum_uti_customer(data.number_customers, 0);
		//for (int i = 0; i < data.number_customers; ++i) {
		//	sum_uti_customer[i] += data.no_purchase[i];
		//	for (int j = 0; j < data.number_products; ++j)
		//		sum_uti_customer[i] += data.utilities[i][j];
		//}

		////cout << "Submodular-cuts y\n" << endl;
		//for (int i = 0; i < data.number_customers; ++i) {
		//	//if (initial_y[i] < partitial_y[i]) {
		//		GRBLinExpr submodular_cut_a_z, submodular_cut_b_z;
		//		for (int j = 0; j < data.number_products; ++j)
		//			if (initial_x[j] == 1) {
		//				submodular_cut_a_z += (1 - x[j]) * data.utilities[i][j] /
		//					(sum_uti_customer[i] * (sum_uti_customer[i] - data.utilities[i][j]));
		//				submodular_cut_b_z += (1 - x[j]) * data.utilities[i][j] /
		//					(initial_denominator[i] * (initial_denominator[i] - data.utilities[i][j]));
		//			}
		//			else {
		//				submodular_cut_a_z -= x[j] * data.utilities[i][j] /
		//					(initial_denominator[i] * (initial_denominator[i] + data.utilities[i][j]));
		//				submodular_cut_b_z -= x[j] * data.utilities[i][j] /
		//					(data.no_purchase[i] * (data.no_purchase[i] + data.utilities[i][j]));
		//			}

		//		submodular_cut_a_z += partitial_y[i];
		//		model.addConstr(y[i] >= submodular_cut_a_z, "ct_sub_modular_a_z_" + to_string(i));
		//		submodular_cut_b_z += partitial_y[i];
		//		model.addConstr(y[i] >= submodular_cut_b_z, "ct_sub_modular_b_z_" + to_string(i));
		//	//}
		//}

		////cout << "Submodular-cuts theta\n" << endl;
		//for (int j = 0; j < data.number_products; ++j) {
		//	//if (initial_theta[j] < partitial_theta[j]) {
		//		GRBLinExpr submodular_cut_a_z, submodular_cut_b_z;
		//		for (int j1 = 0; j1 < data.number_products; ++j1)
		//			if (initial_x[j1] == 1) {
		//				for (int i = 0; i < data.number_customers; ++i) {
		//					submodular_cut_a_z +=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(sum_uti_customer[i] * (sum_uti_customer[i] - data.utilities[i][j1]));
		//					submodular_cut_b_z +=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(initial_denominator[i] * (initial_denominator[i] - data.utilities[i][j1]));
		//				}
		//			}
		//			else {
		//				for (int i = 0; i < data.number_customers; ++i) {
		//					submodular_cut_a_z -=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(initial_denominator[i] * (initial_denominator[i] + data.utilities[i][j1]));
		//					submodular_cut_b_z -=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(data.no_purchase[i] * (data.no_purchase[i] + data.utilities[i][j1]));
		//				}
		//			}

		//		submodular_cut_a_z += partitial_theta[j];
		//		model.addConstr(theta[j] >= submodular_cut_a_z, "ct_sub_modular_a_z_" + to_string(j));
		//		submodular_cut_b_z += partitial_theta[j];
		//		model.addConstr(theta[j] >= submodular_cut_b_z, "ct_sub_modular_b_z_" + to_string(j));
		//	//}
		//}

		//solve
		num_iterative++;
		cout << "Remaining time: " << run_time << endl;

		model.write("submodular.lp");
		model.set(GRB_DoubleParam_TimeLimit, run_time);
		model.set(GRB_IntParam_Threads, 8);
		//model.set(GRB_DoubleParam_MIPGap, 1e-3);
		model.set(GRB_IntParam_MIQCPMethod, 1);
		//model.set(GRB_IntParam_MIQCPMethod, 1);
		//model.set(GRB_IntParam_OutputFlag, 0);

		model.optimize();
		for (int i = 0; i < data.number_customers; ++i)
			if (initial_y[i] >= partitial_y[i])
				model.remove(ct_y[i]);
		for (int j = 0; j < data.number_products; ++j)
			if (initial_theta[j] >= partitial_theta[j])
				model.remove(ct_theta[j]);

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

			for (int i = 0; i < data.number_customers; ++i)
				initial_y[i] = y[i].get(GRB_DoubleAttr_X);
			for (int j = 0; j < data.number_products; ++j)
				initial_theta[j] = theta[j].get(GRB_DoubleAttr_X);

			initial_denominator = calculate_denominator(data, initial_x);

			//check the in equation related to theta_j, x_j and y_j for next iteration
			sub_obj = 0;
			for (int i = 0; i < data.number_customers; ++i)
				sub_obj += data.fraction[i] * data.no_purchase[i] * alpha[i] / initial_denominator[i];
			for (int j = 0; j < data.number_products; ++j)
				for (int i = 0; i < data.number_customers; ++i)
					sub_obj += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * initial_x[j]
					/ initial_denominator[i];

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

void CuttingPlane::solve_multicut_bi(Data data, int number_cuts) {
	//auto start = chrono::steady_clock::now();

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

	auto start = chrono::steady_clock::now();

	//Calculate initial_x
	//vector<int> initial_x(data.number_products, 0);
	vector<int> initial_x = greedy(data, alpha);
	//vector<int> initial_x = greedy_general(data, alpha);
	double greedy_obj = calculate_master_obj(data, initial_x);
	vector<double> initial_y(data.number_customers, 0);
	vector<double> initial_theta(data.number_products, 0);
	vector<double> initial_cut(number_cuts, 0);

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

	vector<double> upper_bound_theta(data.number_products, 0);
	vector<double> lower_bound_theta(data.number_products, 0);
	for (int j = 0; j < data.number_products; ++j)
		for (int i = 0; i < data.number_customers; ++i) {
			lower_bound_theta[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / upper_bound_denominator[i];
			upper_bound_theta[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / data.no_purchase[i];
		}

	GRBVar* theta;
	theta = new GRBVar[data.number_products];
	for (int j = 0; j < data.number_products; ++j)
		theta[j] = model.addVar(lower_bound_theta[j], upper_bound_theta[j], 0, GRB_CONTINUOUS, "theta_" + to_string(j));

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

	//cout << "Constraints related to y and theta\n" << endl;
	for (int j = 0; j < data.number_products; ++j) {
		GRBLinExpr sum;
		for (int i = 0; i < data.number_customers; ++i)
			sum += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * y[i];
		model.addConstr(theta[j] >= sum);
	}

	vector<double> boundx1(data.number_products, 0);
	vector<double> boundx0(data.number_products, 0);
	vector<double> boundin(data.number_products, 0);
	vector<double> boundnotin(data.number_products, 0);
	for (int j = 0; j < data.number_products; ++j)
		for (int i = 0; i < data.number_customers; ++i) {
			boundx1[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + data.utilities[i][j]);
			boundx0[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / data.no_purchase[i];
			boundin[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + bound_in[i][j]);
			boundnotin[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + bound_notin[i][j]);
		}

	//vector<double> subset_boundin(data.number_products, 0);
	//vector<double> subset_boundnotin(data.number_products, 0);
	//for (int j = 0; j < data.number_products; ++j)
	//	for (int i = 0; i < data.number_customers; ++i) {
	//		subset_boundin[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + subset_bound_in[i][j]);
	//		subset_boundnotin[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + subset_bound_notin[i][j]);
	//	}

	//McCornick constraints
	for (int j = 0; j < data.number_products; ++j) {
		model.addQConstr(theta[j] * x[j] <= x[j] * boundx1[j]);

		//if (boundin[j] >= subset_boundin[j])
		model.addQConstr(theta[j] * x[j] >= x[j] * boundin[j]);
		//else
			//model.addQConstr(theta[j] * x[j] >= x[j] * subset_boundin[j]);

		//if (boundnotin[j] >= subset_boundnotin[j])
		model.addQConstr(theta[j] * x[j] <= theta[j] - (1 - x[j]) * boundnotin[j]);
		//else
			//model.addQConstr(theta[j] * x[j] <= theta[j] - (1 - x[j]) * subset_boundnotin[j]);

		model.addQConstr(theta[j] * x[j] >= theta[j] - (1 - x[j]) * boundx0[j]);
	}

	for (int c = 0; c < number_cuts; ++c) {
		GRBLinExpr sum_y;
		for (int i = 0; i < data.number_customers; ++i)
			if (i % number_cuts == c)
				sum_y += data.fraction[i] * data.no_purchase[i] * alpha[i] * y[i];

		model.addConstr(cut_customers[c] >= sum_y, "ct_sub_gradient_y_" + to_string(c));
	}

	//for (int i = 0; i < data.number_customers; ++i) {
	//	GRBQuadExpr sum_z;
	//	for (int j = 0; j < data.number_products; ++j)
	//		sum_z += data.utilities[i][j] * y[i] * x[j];

	//	model.addQConstr(data.no_purchase[i] * y[i] + sum_z >= 1);
	//}

	//Objective
	GRBQuadExpr obj;
	for (int c = 0; c < number_cuts; ++c)
		obj += cut_customers[c];
	for (int j = 0; j < data.number_products; ++j)
		obj += theta[j] * x[j];

	model.setObjective(obj, GRB_MINIMIZE);

	auto time_before_cut = chrono::steady_clock::now();
	chrono::duration<double> before_cut = time_before_cut - start;

	double run_time = time_limit - before_cut.count();

	int num_iterative = 0;
	double stop_param = 1e-4;
	double sub_obj = 1.0;
	double obj_val_cplex = 0.0;
	double best_sub_obj = 0;
	vector<int> best_x = initial_x;
	double best_obj = calculate_master_obj(data, best_x);

	while (sub_obj > obj_val_cplex + stop_param) {
		vector<double> initial_denominator = calculate_denominator(data, initial_x);

		//cout << "Outer-cuts y\n" << endl;
		vector<double> partitial_y(data.number_customers, 0);
		for (int i = 0; i < data.number_customers; ++i)
			partitial_y[i] = 1 / initial_denominator[i];

		vector<vector<double>> subgradient_y(data.number_customers);
		for (int i = 0; i < data.number_customers; ++i)
			subgradient_y[i].resize(data.number_products, 0);

		for (int i = 0; i < data.number_customers; ++i)
			for (int j = 0; j < data.number_products; ++j)
				subgradient_y[i][j] -= data.utilities[i][j] / (initial_denominator[i] * initial_denominator[i]);

		for (int c = 0; c < number_cuts; ++c) {
			GRBLinExpr grad_y;
			double part_y = 0;
			for (int i = 0; i < data.number_customers; ++i)
				if (i % number_cuts == c) {
					part_y += data.fraction[i] * data.no_purchase[i] * alpha[i] * partitial_y[i];
					for (int j = 0; j < data.number_products; ++j)
						grad_y += data.fraction[i] * data.no_purchase[i] * alpha[i] * subgradient_y[i][j] * (x[j] - initial_x[j]);
				}

			//if(initial_cut[c] < part_y)
			model.addConstr(cut_customers[c] >= part_y + grad_y, "ct_sub_gradient_cut_" + to_string(c));
		}

		//cout << "Outer-cuts theta\n" << endl;
		vector<double> partitial_theta(data.number_products, 0);
		for (int j = 0; j < data.number_products; ++j)
			for (int i = 0; i < data.number_customers; ++i)
				partitial_theta[j] +=
				data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j])
				/ initial_denominator[i];

		vector<vector<double>> subgradient_theta(data.number_products);
		for (int j = 0; j < data.number_products; ++j)
			subgradient_theta[j].resize(data.number_products, 0);

		for (int j = 0; j < data.number_products; ++j)
			for (int j1 = 0; j1 < data.number_products; ++j1)
				for (int i = 0; i < data.number_customers; ++i)
					subgradient_theta[j][j1] -=
					data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1]
					/ (initial_denominator[i] * initial_denominator[i]);

		for (int j = 0; j < data.number_products; ++j)
			if (initial_theta[j] < partitial_theta[j]) {
				GRBLinExpr grad_theta;
				for (int j1 = 0; j1 < data.number_products; ++j1)
					grad_theta += subgradient_theta[j][j1] * (x[j1] - initial_x[j1]);

				model.addConstr(theta[j] >= partitial_theta[j] + grad_theta, "ct_sub_gradient_theta_" + to_string(j));
			}

		////cout << "Calculate total utility\n" << endl;
		//vector<double> sum_uti_customer(data.number_customers, 0);
		//for (int i = 0; i < data.number_customers; ++i) {
		//	sum_uti_customer[i] += data.no_purchase[i];
		//	for (int j = 0; j < data.number_products; ++j)
		//		sum_uti_customer[i] += data.utilities[i][j];
		//}

		////cout << "Submodular-cuts y\n" << endl;
		//for (int c = 0; c < number_cuts; ++c) {
		//	GRBLinExpr submodular_cut_a_z, submodular_cut_b_z;
		//	double part_y = 0;
		//	for (int i = 0; i < data.number_customers; ++i)
		//		if (i % number_cuts == c) {
		//			part_y += partitial_y[i];
		//			for (int j = 0; j < data.number_products; ++j)
		//				if (initial_x[j] == 1) {
		//					submodular_cut_a_z += (1 - x[j]) * data.utilities[i][j] /
		//						(sum_uti_customer[i] * (sum_uti_customer[i] - data.utilities[i][j]));
		//					submodular_cut_b_z += (1 - x[j]) * data.utilities[i][j] /
		//						(initial_denominator[i] * (initial_denominator[i] - data.utilities[i][j]));
		//				}
		//				else {
		//					submodular_cut_a_z -= x[j] * data.utilities[i][j] /
		//						(initial_denominator[i] * (initial_denominator[i] + data.utilities[i][j]));
		//					submodular_cut_b_z -= x[j] * data.utilities[i][j] /
		//						(data.no_purchase[i] * (data.no_purchase[i] + data.utilities[i][j]));
		//				}
		//		}

		//	submodular_cut_a_z += part_y;
		//	model.addConstr(cut_customers[c] >= submodular_cut_a_z, "ct_sub_modular_a_z_" + to_string(c));
		//	submodular_cut_b_z += part_y;
		//	model.addConstr(cut_customers[c] >= submodular_cut_b_z, "ct_sub_modular_b_z_" + to_string(c));
		//}

		////cout << "Submodular-cuts theta\n" << endl;
		//for (int j = 0; j < data.number_products; ++j) {
		//	if (initial_theta[j] < partitial_theta[j]) {
		//		GRBLinExpr submodular_cut_a_z, submodular_cut_b_z;
		//		for (int j1 = 0; j1 < data.number_products; ++j1)
		//			if (initial_x[j1] == 1) {
		//				for (int i = 0; i < data.number_customers; ++i) {
		//					submodular_cut_a_z +=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(sum_uti_customer[i] * (sum_uti_customer[i] - data.utilities[i][j1]));
		//					submodular_cut_b_z +=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(initial_denominator[i] * (initial_denominator[i] - data.utilities[i][j1]));
		//				}
		//			}
		//			else {
		//				for (int i = 0; i < data.number_customers; ++i) {
		//					submodular_cut_a_z -=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(initial_denominator[i] * (initial_denominator[i] + data.utilities[i][j1]));
		//					submodular_cut_b_z -=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(data.no_purchase[i] * (data.no_purchase[i] + data.utilities[i][j1]));
		//				}
		//			}

		//		submodular_cut_a_z += partitial_theta[j];
		//		model.addConstr(theta[j] >= submodular_cut_a_z, "ct_sub_modular_a_z_" + to_string(j));
		//		submodular_cut_b_z += partitial_theta[j];
		//		model.addConstr(theta[j] >= submodular_cut_b_z, "ct_sub_modular_b_z_" + to_string(j));
		//	}
		//}

		//solve
		num_iterative++;
		cout << "Remaining time: " << run_time << endl;

		model.write("submodular.lp");
		model.set(GRB_DoubleParam_TimeLimit, run_time);
		model.set(GRB_IntParam_Threads, 8);
		//model.set(GRB_DoubleParam_MIPGap, 1e-3);
		model.set(GRB_IntParam_MIQCPMethod, 1);
		//model.set(GRB_IntParam_MIQCPMethod, 1);
		//model.set(GRB_IntParam_OutputFlag, 0);

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

			for (int i = 0; i < data.number_customers; ++i)
				initial_y[i] = y[i].get(GRB_DoubleAttr_X);
			for (int j = 0; j < data.number_products; ++j)
				initial_theta[j] = theta[j].get(GRB_DoubleAttr_X);
			for (int c = 0; c < number_cuts; ++c)
				initial_cut[c] = cut_customers[c].get(GRB_DoubleAttr_X);

			initial_denominator = calculate_denominator(data, initial_x);

			//check the in equation related to theta_j, x_j and y_j for next iteration
			sub_obj = 0;
			for (int i = 0; i < data.number_customers; ++i)
				sub_obj += data.fraction[i] * data.no_purchase[i] * alpha[i] / initial_denominator[i];
			for (int j = 0; j < data.number_products; ++j)
				for (int i = 0; i < data.number_customers; ++i)
					sub_obj += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * initial_x[j]
					/ initial_denominator[i];

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
	report_results << best_sub_obj << " " << best_obj << " " << greedy_obj << " " << time_for_solve << endl;
	for (int j = 0; j < data.number_products; ++j)
		if (best_x[j] == 1)
			report_results << j << " ";
	report_results.close();
}

void CuttingPlane::solve_multi_multicut_milp(Data data, int number_cuts, int number_cuts_theta) {
	//auto start = chrono::steady_clock::now();

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

	auto start = chrono::steady_clock::now();

	//Calculate initial_x
	//vector<int> initial_x(data.number_products, 0);
	vector<int> initial_x = greedy(data, alpha);
	//vector<int> initial_x = greedy_general(data, alpha);
	double greedy_obj = calculate_master_obj(data, initial_x);
	vector<double> initial_y(data.number_customers, 0);
	vector<vector<double>> initial_theta(data.number_products);
	for (int j = 0; j < data.number_products; ++j)
		initial_theta[j].resize(number_cuts_theta, 0);
	vector<double> initial_cut(number_cuts, 0);

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

	vector<vector<double>> upper_bound_theta(data.number_products);
	vector<vector<double>> lower_bound_theta(data.number_products);
	for (int j = 0; j < data.number_products; ++j) {
		lower_bound_theta[j].resize(number_cuts_theta, 0);
		upper_bound_theta[j].resize(number_cuts_theta, 0);
	}

	for (int j = 0; j < data.number_products; ++j)
		for (int c = 0; c < number_cuts_theta; ++c)
			for (int i = 0; i < data.number_customers; ++i)
				if (i % number_cuts_theta == c) {
					lower_bound_theta[j][c] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / upper_bound_denominator[i];
					upper_bound_theta[j][c] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / data.no_purchase[i];
				}

	GRBVar** theta_multicut;
	theta_multicut = new GRBVar * [data.number_products];
	for (int j = 0; j < data.number_products; ++j)
		theta_multicut[j] = new GRBVar[number_cuts_theta];
	for (int j = 0; j < data.number_products; ++j)
		for (int c = 0; c < number_cuts_theta; ++c)
			theta_multicut[j][c] = model.addVar(lower_bound_theta[j][c], upper_bound_theta[j][c], 0, GRB_CONTINUOUS, "theta_" + to_string(j));

	vector<double> upper_bound_cut_customers(number_cuts, 0);
	vector<double> lower_bound_cut_customers(number_cuts, 0);
	for (int c = 0; c < number_cuts; ++c)
		for (int i = 0; i < data.number_customers; ++i)
			if (i % number_cuts == c) {
				lower_bound_cut_customers[c] += data.fraction[i] * data.no_purchase[i] * alpha[i] * 1 / upper_bound_denominator[i];
				upper_bound_cut_customers[c] += data.fraction[i] * data.no_purchase[i] * alpha[i] * 1 / data.no_purchase[i];
			}

	GRBVar** z;
	z = new GRBVar * [data.number_products];
	for (int j = 0; j < data.number_products; ++j)
		z[j] = new GRBVar[number_cuts_theta];
	for (int j = 0; j < data.number_products; ++j)
		for (int c = 0; c < number_cuts_theta; ++c)
			z[j][c] = model.addVar(0, upper_bound_theta[j][c], 0, GRB_CONTINUOUS, "z_" + to_string(j));

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

	//cout << "Constraints related to y and theta\n" << endl;
	for (int j = 0; j < data.number_products; ++j) {
		for (int c = 0; c < number_cuts_theta; ++c) {
			GRBLinExpr sum;
			for (int i = 0; i < data.number_customers; ++i)
				if (i % number_cuts_theta == c)
					sum += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * y[i];
			model.addConstr(theta_multicut[j][c] >= sum);
		}
	}

	//cout << "Constraints related to z and theta\n" << endl;
	for (int j = 0; j < data.number_products; ++j)
		for (int c = 0; c < number_cuts_theta; ++c)
			model.addConstr(z[j][c] <= theta_multicut[j][c]);

	vector<vector<double>> boundx1(data.number_products);
	vector<vector<double>> boundx0(data.number_products);
	vector<vector<double>> boundin(data.number_products);
	vector<vector<double>> boundnotin(data.number_products);
	for (int j = 0; j < data.number_products; ++j) {
		boundx1[j].resize(number_cuts_theta, 0);
		boundx0[j].resize(number_cuts_theta, 0);
		boundin[j].resize(number_cuts_theta, 0);
		boundnotin[j].resize(number_cuts_theta, 0);
	}
	for (int j = 0; j < data.number_products; ++j)
		for (int c = 0; c < number_cuts; ++c)
			for (int i = 0; i < data.number_customers; ++i)
				if (i % number_cuts_theta == c) {
					boundx1[j][c] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + data.utilities[i][j]);
					boundx0[j][c] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / data.no_purchase[i];
					boundin[j][c] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + bound_in[i][j]);
					boundnotin[j][c] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + bound_notin[i][j]);
				}

	//vector<double> subset_boundin(data.number_products, 0);
	//vector<double> subset_boundnotin(data.number_products, 0);
	//for (int j = 0; j < data.number_products; ++j)
	//	for (int i = 0; i < data.number_customers; ++i) {
	//		subset_boundin[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + subset_bound_in[i][j]);
	//		subset_boundnotin[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + subset_bound_notin[i][j]);
	//	}

	//McCornick constraints
	for (int j = 0; j < data.number_products; ++j)
		for (int c = 0; c < number_cuts_theta; ++c) {
			model.addQConstr(z[j][c] <= x[j] * boundx1[j][c]);

			//if (boundin[j] >= subset_boundin[j])
			model.addQConstr(z[j][c] >= x[j] * boundin[j][c]);
			//else
				//model.addQConstr(z[j][c] >= x[j] * subset_boundin[j]);

			//if (boundnotin[j] >= subset_boundnotin[j])
			model.addQConstr(z[j][c] * x[j] <= theta_multicut[j][c] - (1 - x[j]) * boundnotin[j][c]);
			//else
				//model.addQConstr(z[j][c] <= theta_multicut[j][c] - (1 - x[j]) * subset_boundnotin[j]);

			model.addQConstr(z[j][c] * x[j] >= theta_multicut[j][c] - (1 - x[j]) * boundx0[j][c]);
		}

	for (int c = 0; c < number_cuts; ++c) {
		GRBLinExpr sum_y;
		for (int i = 0; i < data.number_customers; ++i)
			if (i % number_cuts == c)
				sum_y += data.fraction[i] * data.no_purchase[i] * alpha[i] * y[i];

		model.addConstr(cut_customers[c] >= sum_y, "ct_sub_gradient_y_" + to_string(c));
	}

	//for (int i = 0; i < data.number_customers; ++i) {
	//	GRBQuadExpr sum_z;
	//	for (int j = 0; j < data.number_products; ++j)
	//		sum_z += data.utilities[i][j] * y[i] * x[j];

	//	model.addQConstr(data.no_purchase[i] * y[i] + sum_z >= 1);
	//}

	//Objective
	GRBQuadExpr obj;
	for (int c = 0; c < number_cuts; ++c)
		obj += cut_customers[c];
	for (int j = 0; j < data.number_products; ++j)
		for (int c = 0; c < number_cuts_theta; ++c)
			obj += z[j][c];

	model.setObjective(obj, GRB_MINIMIZE);

	auto time_before_cut = chrono::steady_clock::now();
	chrono::duration<double> before_cut = time_before_cut - start;

	double run_time = time_limit - before_cut.count();

	int num_iterative = 0;
	double stop_param = 1e-4;
	double sub_obj = 1.0;
	double obj_val_cplex = 0.0;
	double best_sub_obj = 0;
	vector<int> best_x = initial_x;
	double best_obj = calculate_master_obj(data, best_x);

	while (sub_obj > obj_val_cplex + stop_param) {
		vector<double> initial_denominator = calculate_denominator(data, initial_x);

		//cout << "Outer-cuts y\n" << endl;
		vector<double> partitial_y(data.number_customers, 0);
		for (int i = 0; i < data.number_customers; ++i)
			partitial_y[i] = 1 / initial_denominator[i];

		vector<vector<double>> subgradient_y(data.number_customers);
		for (int i = 0; i < data.number_customers; ++i)
			subgradient_y[i].resize(data.number_products, 0);

		for (int i = 0; i < data.number_customers; ++i)
			for (int j = 0; j < data.number_products; ++j)
				subgradient_y[i][j] -= data.utilities[i][j] / (initial_denominator[i] * initial_denominator[i]);

		for (int c = 0; c < number_cuts; ++c) {
			GRBLinExpr grad_y;
			double part_y = 0;
			for (int i = 0; i < data.number_customers; ++i)
				if (i % number_cuts == c) {
					part_y += data.fraction[i] * data.no_purchase[i] * alpha[i] * partitial_y[i];
					for (int j = 0; j < data.number_products; ++j)
						grad_y += data.fraction[i] * data.no_purchase[i] * alpha[i] * subgradient_y[i][j] * (x[j] - initial_x[j]);
				}

			//if(initial_cut[c] < part_y)
			model.addConstr(cut_customers[c] >= part_y + grad_y, "ct_sub_gradient_cut_" + to_string(c));
		}

		//cout << "Outer-cuts theta\n" << endl;
		vector<vector<double>> partitial_theta_multicut(data.number_products);
		for (int j = 0; j < data.number_products; ++j)
			partitial_theta_multicut[j].resize(number_cuts_theta, 0);

		for (int j = 0; j < data.number_products; ++j)
			for (int c = 0; c < number_cuts_theta; ++c)
				for (int i = 0; i < data.number_customers; ++i)
					if (i % number_cuts_theta == c)
						partitial_theta_multicut[j][c] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / initial_denominator[i];

		vector<vector<vector<double>>> subgradient_theta(data.number_products);
		for (int j = 0; j < data.number_products; ++j)
			subgradient_theta[j].resize(number_cuts_theta);

		for (int j = 0; j < data.number_products; ++j)
			for (int c = 0; c < number_cuts_theta; ++c)
				subgradient_theta[j][c].resize(data.number_products, 0);

		for (int j = 0; j < data.number_products; ++j)
			for (int c = 0; c < number_cuts_theta; ++c)
				for (int j1 = 0; j1 < data.number_products; ++j1)
					for (int i = 0; i < data.number_customers; ++i)
						if (i % number_cuts_theta == c)
							subgradient_theta[j][c][j1] -= data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1]
							/ (initial_denominator[i] * initial_denominator[i]);

		for (int j = 0; j < data.number_products; ++j)
			for (int c = 0; c < number_cuts_theta; ++c)
				if (initial_theta[j][c] < partitial_theta_multicut[j][c]) {
					GRBLinExpr grad_theta;
					for (int j1 = 0; j1 < data.number_products; ++j1)
						grad_theta += subgradient_theta[j][c][j1] * (x[j1] - initial_x[j1]);

					model.addConstr(theta_multicut[j][c] >= partitial_theta_multicut[j][c] + grad_theta, "ct_sub_gradient_theta_" + to_string(j));
				}

		////cout << "Calculate total utility\n" << endl;
		//vector<double> sum_uti_customer(data.number_customers, 0);
		//for (int i = 0; i < data.number_customers; ++i) {
		//	sum_uti_customer[i] += data.no_purchase[i];
		//	for (int j = 0; j < data.number_products; ++j)
		//		sum_uti_customer[i] += data.utilities[i][j];
		//}

		////cout << "Submodular-cuts y\n" << endl;
		//for (int c = 0; c < number_cuts; ++c) {
		//	GRBLinExpr submodular_cut_a_z, submodular_cut_b_z;
		//	double part_y = 0;
		//	for (int i = 0; i < data.number_customers; ++i)
		//		if (i % number_cuts == c) {
		//			part_y += partitial_y[i];
		//			for (int j = 0; j < data.number_products; ++j)
		//				if (initial_x[j] == 1) {
		//					submodular_cut_a_z += (1 - x[j]) * data.utilities[i][j] /
		//						(sum_uti_customer[i] * (sum_uti_customer[i] - data.utilities[i][j]));
		//					submodular_cut_b_z += (1 - x[j]) * data.utilities[i][j] /
		//						(initial_denominator[i] * (initial_denominator[i] - data.utilities[i][j]));
		//				}
		//				else {
		//					submodular_cut_a_z -= x[j] * data.utilities[i][j] /
		//						(initial_denominator[i] * (initial_denominator[i] + data.utilities[i][j]));
		//					submodular_cut_b_z -= x[j] * data.utilities[i][j] /
		//						(data.no_purchase[i] * (data.no_purchase[i] + data.utilities[i][j]));
		//				}
		//		}

		//	submodular_cut_a_z += part_y;
		//	model.addConstr(cut_customers[c] >= submodular_cut_a_z, "ct_sub_modular_a_z_" + to_string(c));
		//	submodular_cut_b_z += part_y;
		//	model.addConstr(cut_customers[c] >= submodular_cut_b_z, "ct_sub_modular_b_z_" + to_string(c));
		//}

		////cout << "Submodular-cuts theta\n" << endl;
		//for (int j = 0; j < data.number_products; ++j) {
		//	if (initial_theta[j] < partitial_theta[j]) {
		//		GRBLinExpr submodular_cut_a_z, submodular_cut_b_z;
		//		for (int j1 = 0; j1 < data.number_products; ++j1)
		//			if (initial_x[j1] == 1) {
		//				for (int i = 0; i < data.number_customers; ++i) {
		//					submodular_cut_a_z +=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(sum_uti_customer[i] * (sum_uti_customer[i] - data.utilities[i][j1]));
		//					submodular_cut_b_z +=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(initial_denominator[i] * (initial_denominator[i] - data.utilities[i][j1]));
		//				}
		//			}
		//			else {
		//				for (int i = 0; i < data.number_customers; ++i) {
		//					submodular_cut_a_z -=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(initial_denominator[i] * (initial_denominator[i] + data.utilities[i][j1]));
		//					submodular_cut_b_z -=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(data.no_purchase[i] * (data.no_purchase[i] + data.utilities[i][j1]));
		//				}
		//			}

		//		submodular_cut_a_z += partitial_theta[j];
		//		model.addConstr(theta[j] >= submodular_cut_a_z, "ct_sub_modular_a_z_" + to_string(j));
		//		submodular_cut_b_z += partitial_theta[j];
		//		model.addConstr(theta[j] >= submodular_cut_b_z, "ct_sub_modular_b_z_" + to_string(j));
		//	}
		//}

		//solve
		num_iterative++;
		cout << "Remaining time: " << run_time << endl;

		model.write("submodular.lp");
		model.set(GRB_DoubleParam_TimeLimit, run_time);
		model.set(GRB_IntParam_Threads, 8);
		//model.set(GRB_DoubleParam_MIPGap, 1e-3);
		model.set(GRB_IntParam_MIQCPMethod, 1);
		//model.set(GRB_IntParam_MIQCPMethod, 1);
		//model.set(GRB_IntParam_OutputFlag, 0);

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

			for (int i = 0; i < data.number_customers; ++i)
				initial_y[i] = y[i].get(GRB_DoubleAttr_X);
			for (int j = 0; j < data.number_products; ++j)
				for (int c = 0; c < number_cuts_theta; ++c)
					initial_theta[j][c] = theta_multicut[j][c].get(GRB_DoubleAttr_X);
			for (int c = 0; c < number_cuts; ++c)
				initial_cut[c] = cut_customers[c].get(GRB_DoubleAttr_X);

			initial_denominator = calculate_denominator(data, initial_x);

			//check the in equation related to theta_j, x_j and y_j for next iteration
			sub_obj = 0;
			for (int i = 0; i < data.number_customers; ++i)
				sub_obj += data.fraction[i] * data.no_purchase[i] * alpha[i] / initial_denominator[i];
			for (int j = 0; j < data.number_products; ++j)
				for (int i = 0; i < data.number_customers; ++i)
					sub_obj += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * initial_x[j]
					/ initial_denominator[i];

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
	report_results << best_sub_obj << " " << best_obj << " " << greedy_obj << " " << time_for_solve << endl;
	for (int j = 0; j < data.number_products; ++j)
		if (best_x[j] == 1)
			report_results << j << " ";
	report_results.close();
}

void CuttingPlane::solve_multi_multicut_bi(Data data, int number_cuts, int number_cuts_theta) {
	//auto start = chrono::steady_clock::now();

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

	auto start = chrono::steady_clock::now();

	//Calculate initial_x
	//vector<int> initial_x(data.number_products, 0);
	vector<int> initial_x = greedy(data, alpha);
	//vector<int> initial_x = greedy_general(data, alpha);
	double greedy_obj = calculate_master_obj(data, initial_x);
	vector<double> initial_y(data.number_customers, 0);
	vector<vector<double>> initial_theta(data.number_products);
	for (int j = 0; j < data.number_products; ++j)
		initial_theta[j].resize(number_cuts_theta, 0);
	vector<double> initial_cut(number_cuts, 0);

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

	vector<vector<double>> upper_bound_theta(data.number_products);
	vector<vector<double>> lower_bound_theta(data.number_products);
	for (int j = 0; j < data.number_products; ++j) {
		lower_bound_theta[j].resize(number_cuts_theta, 0);
		upper_bound_theta[j].resize(number_cuts_theta, 0);
	}

	for (int j = 0; j < data.number_products; ++j)
		for (int c = 0; c < number_cuts_theta; ++c)
			for (int i = 0; i < data.number_customers; ++i)
				if (i % number_cuts_theta == c) {
					lower_bound_theta[j][c] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / upper_bound_denominator[i];
					upper_bound_theta[j][c] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / data.no_purchase[i];
				}

	GRBVar** theta_multicut;
	theta_multicut = new GRBVar * [data.number_products];
	for (int j = 0; j < data.number_products; ++j)
		theta_multicut[j] = new GRBVar[number_cuts_theta];
	for (int j = 0; j < data.number_products; ++j)
		for (int c = 0; c < number_cuts_theta; ++c)
			theta_multicut[j][c] = model.addVar(lower_bound_theta[j][c], upper_bound_theta[j][c], 0, GRB_CONTINUOUS, "theta_" + to_string(j));

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

	//cout << "Constraints related to y and theta\n" << endl;
	for (int j = 0; j < data.number_products; ++j) {
		for (int c = 0; c < number_cuts_theta; ++c) {
			GRBLinExpr sum;
			for (int i = 0; i < data.number_customers; ++i)
				if (i % number_cuts_theta == c)
					sum += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * y[i];
			model.addConstr(theta_multicut[j][c] >= sum);
		}
	}

	vector<vector<double>> boundx1(data.number_products);
	vector<vector<double>> boundx0(data.number_products);
	vector<vector<double>> boundin(data.number_products);
	vector<vector<double>> boundnotin(data.number_products);
	for (int j = 0; j < data.number_products; ++j) {
		boundx1[j].resize(number_cuts_theta, 0);
		boundx0[j].resize(number_cuts_theta, 0);
		boundin[j].resize(number_cuts_theta, 0);
		boundnotin[j].resize(number_cuts_theta, 0);
	}
	for (int j = 0; j < data.number_products; ++j)
		for (int c = 0; c < number_cuts; ++c)
			for (int i = 0; i < data.number_customers; ++i)
				if (i % number_cuts_theta == c) {
					boundx1[j][c] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + data.utilities[i][j]);
					boundx0[j][c] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / data.no_purchase[i];
					boundin[j][c] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + bound_in[i][j]);
					boundnotin[j][c] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + bound_notin[i][j]);
				}

	//vector<double> subset_boundin(data.number_products, 0);
	//vector<double> subset_boundnotin(data.number_products, 0);
	//for (int j = 0; j < data.number_products; ++j)
	//	for (int i = 0; i < data.number_customers; ++i) {
	//		subset_boundin[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + subset_bound_in[i][j]);
	//		subset_boundnotin[j] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / (data.no_purchase[i] + subset_bound_notin[i][j]);
	//	}

	//McCornick constraints
	for (int j = 0; j < data.number_products; ++j)
		for (int c = 0; c < number_cuts_theta; ++c) {
			model.addQConstr(theta_multicut[j][c] * x[j] <= x[j] * boundx1[j][c]);

			//if (boundin[j] >= subset_boundin[j])
			model.addQConstr(theta_multicut[j][c] * x[j] >= x[j] * boundin[j][c]);
			//else
				//model.addQConstr(theta_multicut[j][c] * x[j] >= x[j] * subset_boundin[j]);

			//if (boundnotin[j] >= subset_boundnotin[j])
			model.addQConstr(theta_multicut[j][c] * x[j] <= theta_multicut[j][c] - (1 - x[j]) * boundnotin[j][c]);
			//else
				//model.addQConstr(theta_multicut[j][c] * x[j] <= theta_multicut[j][c] - (1 - x[j]) * subset_boundnotin[j]);

			model.addQConstr(theta_multicut[j][c] * x[j] >= theta_multicut[j][c] - (1 - x[j]) * boundx0[j][c]);
		}

	for (int c = 0; c < number_cuts; ++c) {
		GRBLinExpr sum_y;
		for (int i = 0; i < data.number_customers; ++i)
			if (i % number_cuts == c)
				sum_y += data.fraction[i] * data.no_purchase[i] * alpha[i] * y[i];

		model.addConstr(cut_customers[c] >= sum_y, "ct_sub_gradient_y_" + to_string(c));
	}

	//for (int i = 0; i < data.number_customers; ++i) {
	//	GRBQuadExpr sum_z;
	//	for (int j = 0; j < data.number_products; ++j)
	//		sum_z += data.utilities[i][j] * y[i] * x[j];

	//	model.addQConstr(data.no_purchase[i] * y[i] + sum_z >= 1);
	//}

	//Objective
	GRBQuadExpr obj;
	for (int c = 0; c < number_cuts; ++c)
		obj += cut_customers[c];
	for (int j = 0; j < data.number_products; ++j)
		for (int c = 0; c < number_cuts_theta; ++c)
			obj += theta_multicut[j][c] * x[j];

	model.setObjective(obj, GRB_MINIMIZE);

	auto time_before_cut = chrono::steady_clock::now();
	chrono::duration<double> before_cut = time_before_cut - start;

	double run_time = time_limit - before_cut.count();

	int num_iterative = 0;
	double stop_param = 1e-4;
	double sub_obj = 1.0;
	double obj_val_cplex = 0.0;
	double best_sub_obj = 0;
	vector<int> best_x = initial_x;
	double best_obj = calculate_master_obj(data, best_x);

	while (sub_obj > obj_val_cplex + stop_param) {
		vector<double> initial_denominator = calculate_denominator(data, initial_x);

		//cout << "Outer-cuts y\n" << endl;
		vector<double> partitial_y(data.number_customers, 0);
		for (int i = 0; i < data.number_customers; ++i)
			partitial_y[i] = 1 / initial_denominator[i];

		vector<vector<double>> subgradient_y(data.number_customers);
		for (int i = 0; i < data.number_customers; ++i)
			subgradient_y[i].resize(data.number_products, 0);

		for (int i = 0; i < data.number_customers; ++i)
			for (int j = 0; j < data.number_products; ++j)
				subgradient_y[i][j] -= data.utilities[i][j] / (initial_denominator[i] * initial_denominator[i]);

		for (int c = 0; c < number_cuts; ++c) {
			GRBLinExpr grad_y;
			double part_y = 0;
			for (int i = 0; i < data.number_customers; ++i)
				if (i % number_cuts == c) {
					part_y += data.fraction[i] * data.no_purchase[i] * alpha[i] * partitial_y[i];
					for (int j = 0; j < data.number_products; ++j)
						grad_y += data.fraction[i] * data.no_purchase[i] * alpha[i] * subgradient_y[i][j] * (x[j] - initial_x[j]);
				}

			//if(initial_cut[c] < part_y)
			model.addConstr(cut_customers[c] >= part_y + grad_y, "ct_sub_gradient_cut_" + to_string(c));
		}

		//cout << "Outer-cuts theta\n" << endl;
		vector<vector<double>> partitial_theta_multicut(data.number_products);
		for (int j = 0; j < data.number_products; ++j)
			partitial_theta_multicut[j].resize(number_cuts_theta, 0);

		for (int j = 0; j < data.number_products; ++j)
			for (int c = 0; c < number_cuts_theta; ++c)
				for (int i = 0; i < data.number_customers; ++i)
					if (i % number_cuts_theta == c)
						partitial_theta_multicut[j][c] += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) / initial_denominator[i];

		vector<vector<vector<double>>> subgradient_theta(data.number_products);
		for (int j = 0; j < data.number_products; ++j)
			subgradient_theta[j].resize(number_cuts_theta);

		for (int j = 0; j < data.number_products; ++j)
			for (int c = 0; c < number_cuts_theta; ++c)
				subgradient_theta[j][c].resize(data.number_products, 0);

		for (int j = 0; j < data.number_products; ++j)
			for (int c = 0; c < number_cuts_theta; ++c)
				for (int j1 = 0; j1 < data.number_products; ++j1)
					for (int i = 0; i < data.number_customers; ++i)
						if (i % number_cuts_theta == c)
							subgradient_theta[j][c][j1] -= data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1]
							/ (initial_denominator[i] * initial_denominator[i]);

		for (int j = 0; j < data.number_products; ++j)
			for (int c = 0; c < number_cuts_theta; ++c)
				if (initial_theta[j][c] < partitial_theta_multicut[j][c]) {
					GRBLinExpr grad_theta;
					for (int j1 = 0; j1 < data.number_products; ++j1)
						grad_theta += subgradient_theta[j][c][j1] * (x[j1] - initial_x[j1]);

					model.addConstr(theta_multicut[j][c] >= partitial_theta_multicut[j][c] + grad_theta, "ct_sub_gradient_theta_" + to_string(j));
				}

		////cout << "Calculate total utility\n" << endl;
		//vector<double> sum_uti_customer(data.number_customers, 0);
		//for (int i = 0; i < data.number_customers; ++i) {
		//	sum_uti_customer[i] += data.no_purchase[i];
		//	for (int j = 0; j < data.number_products; ++j)
		//		sum_uti_customer[i] += data.utilities[i][j];
		//}

		////cout << "Submodular-cuts y\n" << endl;
		//for (int c = 0; c < number_cuts; ++c) {
		//	GRBLinExpr submodular_cut_a_z, submodular_cut_b_z;
		//	double part_y = 0;
		//	for (int i = 0; i < data.number_customers; ++i)
		//		if (i % number_cuts == c) {
		//			part_y += partitial_y[i];
		//			for (int j = 0; j < data.number_products; ++j)
		//				if (initial_x[j] == 1) {
		//					submodular_cut_a_z += (1 - x[j]) * data.utilities[i][j] /
		//						(sum_uti_customer[i] * (sum_uti_customer[i] - data.utilities[i][j]));
		//					submodular_cut_b_z += (1 - x[j]) * data.utilities[i][j] /
		//						(initial_denominator[i] * (initial_denominator[i] - data.utilities[i][j]));
		//				}
		//				else {
		//					submodular_cut_a_z -= x[j] * data.utilities[i][j] /
		//						(initial_denominator[i] * (initial_denominator[i] + data.utilities[i][j]));
		//					submodular_cut_b_z -= x[j] * data.utilities[i][j] /
		//						(data.no_purchase[i] * (data.no_purchase[i] + data.utilities[i][j]));
		//				}
		//		}

		//	submodular_cut_a_z += part_y;
		//	model.addConstr(cut_customers[c] >= submodular_cut_a_z, "ct_sub_modular_a_z_" + to_string(c));
		//	submodular_cut_b_z += part_y;
		//	model.addConstr(cut_customers[c] >= submodular_cut_b_z, "ct_sub_modular_b_z_" + to_string(c));
		//}

		////cout << "Submodular-cuts theta\n" << endl;
		//for (int j = 0; j < data.number_products; ++j) {
		//	if (initial_theta[j] < partitial_theta[j]) {
		//		GRBLinExpr submodular_cut_a_z, submodular_cut_b_z;
		//		for (int j1 = 0; j1 < data.number_products; ++j1)
		//			if (initial_x[j1] == 1) {
		//				for (int i = 0; i < data.number_customers; ++i) {
		//					submodular_cut_a_z +=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(sum_uti_customer[i] * (sum_uti_customer[i] - data.utilities[i][j1]));
		//					submodular_cut_b_z +=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(initial_denominator[i] * (initial_denominator[i] - data.utilities[i][j1]));
		//				}
		//			}
		//			else {
		//				for (int i = 0; i < data.number_customers; ++i) {
		//					submodular_cut_a_z -=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(initial_denominator[i] * (initial_denominator[i] + data.utilities[i][j1]));
		//					submodular_cut_b_z -=
		//						(1 - x[j1]) * data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * data.utilities[i][j1] /
		//						(data.no_purchase[i] * (data.no_purchase[i] + data.utilities[i][j1]));
		//				}
		//			}

		//		submodular_cut_a_z += partitial_theta[j];
		//		model.addConstr(theta[j] >= submodular_cut_a_z, "ct_sub_modular_a_z_" + to_string(j));
		//		submodular_cut_b_z += partitial_theta[j];
		//		model.addConstr(theta[j] >= submodular_cut_b_z, "ct_sub_modular_b_z_" + to_string(j));
		//	}
		//}

		//solve
		num_iterative++;
		cout << "Remaining time: " << run_time << endl;

		model.write("submodular.lp");
		model.set(GRB_DoubleParam_TimeLimit, run_time);
		model.set(GRB_IntParam_Threads, 8);
		//model.set(GRB_DoubleParam_MIPGap, 1e-3);
		model.set(GRB_IntParam_MIQCPMethod, 1);
		//model.set(GRB_IntParam_MIQCPMethod, 1);
		//model.set(GRB_IntParam_OutputFlag, 0);

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

			for (int i = 0; i < data.number_customers; ++i)
				initial_y[i] = y[i].get(GRB_DoubleAttr_X);
			for (int j = 0; j < data.number_products; ++j)
				for (int c = 0; c < number_cuts_theta; ++c)
					initial_theta[j][c] = theta_multicut[j][c].get(GRB_DoubleAttr_X);
			for (int c = 0; c < number_cuts; ++c)
				initial_cut[c] = cut_customers[c].get(GRB_DoubleAttr_X);

			initial_denominator = calculate_denominator(data, initial_x);

			//check the in equation related to theta_j, x_j and y_j for next iteration
			sub_obj = 0;
			for (int i = 0; i < data.number_customers; ++i)
				sub_obj += data.fraction[i] * data.no_purchase[i] * alpha[i] / initial_denominator[i];
			for (int j = 0; j < data.number_products; ++j)
				for (int i = 0; i < data.number_customers; ++i)
					sub_obj += data.fraction[i] * data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * initial_x[j]
					/ initial_denominator[i];

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
	report_results << best_sub_obj << " " << best_obj << " " << greedy_obj << " " << time_for_solve << endl;
	for (int j = 0; j < data.number_products; ++j)
		if (best_x[j] == 1)
			report_results << j << " ";
	report_results.close();
}