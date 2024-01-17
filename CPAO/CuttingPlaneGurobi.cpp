#include "CuttingPlaneGurobi.h"
#include <chrono>
#include <algorithm>
#include <sstream>
#include <cstdlib>
#include <cassert>

CuttingPlaneGurobi::CuttingPlaneGurobi() {

}

CuttingPlaneGurobi::CuttingPlaneGurobi(Data data, double time_limit, string outfile) {
	this->data = data;
	this->time_limit = time_limit;
	this->out_res_csv = outfile;
}

vector<double> CuttingPlaneGurobi::calculate_y(Data data, vector<int> x, vector<double> alpha) {
	vector<double> y(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i) {
		double tmp_y = alpha[i] * data.no_purchase[i];
		for (int j = 0; j < data.number_products; ++j)
			tmp_y += (alpha[i] - data.revenue[i][j]) * x[j] * data.utilities[i][j];
		y[i] = log(tmp_y);
	}

	return y;
}

vector<double> CuttingPlaneGurobi::calculate_z(Data data, vector<int> x) {
	vector<double> z(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i) {
		double tmp_z = data.no_purchase[i];
		for (int j = 0; j < data.number_products; ++j)
			tmp_z += x[j] * data.utilities[i][j];
		z[i] = -log(tmp_z);
	}

	return z;
}

double  CuttingPlaneGurobi::calculate_original_obj(Data data, vector<int> x, vector<double> alpha) {
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

double CuttingPlaneGurobi::calculate_master_obj(Data data, vector<int> x) {
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

double CuttingPlaneGurobi::calculate_original_obj_tmp(Data data, vector<int> x, vector<double> alpha, int candidate) {
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

vector<int> CuttingPlaneGurobi::greedy(Data data, int budget, vector<double> alpha) {
	auto start = chrono::steady_clock::now();
	vector<int> chosen(data.number_products, 0);
	double obj = calculate_original_obj(data, chosen, alpha);
	int number_chosen = 0;

	while (number_chosen < budget) {
		double extra_ratio = 0;
		int inserted_product = -1;
		double min = 999999;
		for (int j = 0; j < data.number_products; ++j) {
			if (chosen[j] == 0) {
				extra_ratio = calculate_original_obj_tmp(data, chosen, alpha, j) - obj;
				if (min > extra_ratio) {
					min = extra_ratio;
					inserted_product = j;
				}
			}
		}
		if (inserted_product != -1) {
			chosen[inserted_product] = 1;
			//cout << "Best-Insertion: " << inserted_product << endl;
			obj = calculate_original_obj(data, chosen, alpha);
			number_chosen++;
		}
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

double CuttingPlaneGurobi::calculate_bound_y(Data data, int budget, int i, double alpha) {
	double sum = 0;
	vector<double> u(data.number_products);
	for (int j = 0; j < data.number_products; ++j)
		u[j] = (alpha - data.revenue[i][j]) * data.utilities[i][j];

	sort(u.begin(), u.end(), greater<double>());

	for (int j = 0; j < budget; ++j)
		sum += u[j];

	return sum;
}

double CuttingPlaneGurobi::calculate_optimal_bound_y(Data data, int budget, int i, double alpha) {
	GRBEnv env = GRBEnv(true);
	env.start();

	GRBModel model = GRBModel(env);

	//Decison variables: x_j = 1 if product j is chosen, 0 otherwise
	GRBVar* x = 0;
	x = model.addVars(data.number_products, GRB_BINARY);

	//Budget constraint
	GRBLinExpr capacity = 0;
	for (int j = 0; j < data.number_products; ++j) {
		capacity += x[j];
	}
	model.addConstr(capacity <= budget);

	GRBLinExpr obj;
	for (int j = 0; j < data.number_products; ++j)
		obj += data.utilities[i][j] * x[j] * (alpha - data.revenue[i][j]);
	model.setObjective(obj, GRB_MAXIMIZE);

	model.optimize();

	return model.get(GRB_DoubleAttr_ObjVal);
}

vector<vector<double>> CuttingPlaneGurobi::create_sub_intervals(Data data, int budget, vector<double> alpha, int number_itervals) {
	vector<vector<double>> c(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i)
		c[i].resize(number_itervals + 1, 0);
	for (int i = 0; i < data.number_customers; ++i) {
		c[i][0] = log(alpha[i] * data.no_purchase[i]);
		double upper = alpha[i] * data.no_purchase[i] + calculate_bound_y(data, budget, i, alpha[i]);
		c[i][c[i].size() - 1] = log(upper);
		for (int k = 1; k < c[i].size() - 1; ++k)
			c[i][k] = c[i][0] + k * (c[i][c[i].size() - 1] - c[i][0]) / number_itervals;
	}

	return c;
}

double CuttingPlaneGurobi::next_approximate_point(double b, double epsilon) {
	double x = b + 0.0001;
	int count = 0;
	while (exp(b) + (exp(x) - exp(b)) / (x - b) * (log((exp(x) - exp(b)) / (x - b)) - b) - (exp(x) - exp(b)) / (x - b) <= epsilon) {
		x += 0.0001;
		count++;
	}
	return x;
}

vector<vector<double>> CuttingPlaneGurobi::optimal_sub_intervals(Data data, int budget, vector<double> alpha, double epsilon) {
	vector<vector<double>> c(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i) {
		c[i].push_back(log(alpha[i] * data.no_purchase[i]));
		double upper = log(alpha[i] * data.no_purchase[i] + calculate_bound_y(data, budget, i, alpha[i]));
		while (true) {
			double x = next_approximate_point(c[i][c[i].size() - 1], epsilon);
			if (x >= upper) break;
			else c[i].push_back(x);
		}
		c[i].push_back(upper);
		cout << c[i].size() << " ";
	}
	cout << endl;
	return c;
}

void CuttingPlaneGurobi::solve(Data data, int budget) {
	auto start = chrono::steady_clock::now(); //get start time

	vector<double> alpha(data.number_customers, -1);
	for (int i = 0; i < data.number_customers; ++i)
		for (int j = 0; j < data.number_products; ++j)
			if (data.revenue[i][j] > alpha[i])
				alpha[i] = data.revenue[i][j];

	//Calculate initial_x, initial_y, initial_z
	vector<int> initial_x = greedy(data, budget, alpha);

	vector<int> best_x = initial_x;
	double best_obj = calculate_master_obj(data, best_x);

	vector<double> initial_y = calculate_y(data, initial_x, alpha);
	vector<double> initial_z = calculate_z(data, initial_x);

	//create bounds c^i_k for e^{y_i}
	//vector<vector<double>> c = create_sub_intervals(data, budget, alpha, 200);
	vector<vector<double>> c = optimal_sub_intervals(data, budget, alpha, 0.00005);
	vector<int> number_sub_intervals(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i)
		number_sub_intervals[i] = c[i].size() - 1;

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
	y = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		y[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "y_" + to_string(i));

	//cout << "Slack variables : z_i\n" << endl;
	GRBVar* z;
	z = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		z[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "z_" + to_string(i));

	//cout << "Decision variables: r_{ik}\n" << endl;
	GRBVar** r;
	r = new GRBVar* [data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		r[i] = new GRBVar[number_sub_intervals[i]];
	for (int i = 0; i < data.number_customers; ++i)
		for (int k = 0; k < number_sub_intervals[i]; ++k)
			r[i][k] = model.addVar(0, 1, 0, GRB_BINARY, "r_" + to_string(i) + "_" + to_string(k));

	//cout << "Decision variables: s_{ik}\n" << endl;
	GRBVar** s;
	s = new GRBVar* [data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		s[i] = new GRBVar[number_sub_intervals[i]];
	for (int i = 0; i < data.number_customers; ++i)
		for (int k = 0; k < number_sub_intervals[i]; ++k)
			s[i][k] = model.addVar(0, 1, 0, GRB_CONTINUOUS,"s_" + to_string(i) + "_" + to_string(k));

	//cout << "Slack variables : theta_i\n" << endl;
	GRBVar* theta = 0;
	theta = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		theta[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "theta_" + to_string(i));

	//Constraints related to e^{y_i}
	//exp(c[i][0]) = alpha * no_purchase[i] => remove form both sides
	for (int i = 0; i < data.number_customers; ++i) {
		GRBLinExpr sum_r, sum_x;
		for (int k = 0; k < number_sub_intervals[i]; ++k) {
			sum_r += (exp(c[i][k + 1]) - exp(c[i][k])) * r[i][k];
		}
		for (int j = 0; j < data.number_products; ++j)
			sum_x += (alpha[i] - data.revenue[i][j]) * x[j] * data.utilities[i][j];

		model.addConstr(sum_r >= sum_x, "ct_e_y_" + to_string(i));
	}

	//Constraints related to y_i
	for (int i = 0; i < data.number_customers; ++i) {
		GRBLinExpr sum_r;
		for (int k = 0; k < number_sub_intervals[i]; ++k) {
			sum_r += (c[i][k + 1] - c[i][k]) * r[i][k];
		}

		model.addConstr(y[i] == c[i][0] + sum_r, "ct_y_" + to_string(i));
	}

	//Constraints related to s_{ik} and r_{ik}
	for (int i = 0; i < data.number_customers; ++i) {
		for (int k = 0; k < number_sub_intervals[i]; ++k) {
			model.addConstr(r[i][k] >= s[i][k], "ct_rs_" + to_string(i) + "_" + to_string(k));

			if (k < number_sub_intervals[i] - 1) {
				model.addConstr(s[i][k] >= s[i][k + 1], "ct_ss_" + to_string(i) + "_" + to_string(k));
				model.addConstr(s[i][k] >= r[i][k + 1], "ct_sr_" + to_string(i) + "_" + to_string(k));
			}
		}
	}

	//Budget constraint
	GRBLinExpr capacity;
	for (int j = 0; j < data.number_products; ++j) {
		capacity += x[j];
	}
	model.addConstr(capacity <= budget, "ct_budget");

	//Objective
	GRBLinExpr obj;
	for (int i = 0; i < data.number_customers; ++i)
		obj += theta[i];
	model.setObjective(obj, GRB_MINIMIZE);

	auto time_before_cut = chrono::steady_clock::now();
	chrono::duration<double> before_cut = time_before_cut - start;

	double run_time = time_limit - before_cut.count();

	int num_iterative = 0;
	double stop_param = 1e-4;
	double sub_obj = 1.0;
	double obj_val_cplex = 0.0;

	//Add cut iteratively
	while (sub_obj > stop_param + obj_val_cplex) {

		//compute gradient e^{y+z} at initial_x, initial_y, initial_z and set up constraints related to theta
		for (int i = 0; i < data.number_customers; ++i)
			model.addConstr(theta[i] >= exp(initial_y[i] + initial_z[i]) * (1 + y[i] - initial_y[i] + z[i] - initial_z[i]), "ct_sub_gradient_y+z_" + to_string(i));

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
		model.set(GRB_IntParam_MIPFocus, 3);
		//model.set(GRB_IntParam_MIPFocus, 1);
		model.set(GRB_IntParam_OutputFlag, 0);
		
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
				//cout << initial_y[i] << " " << initial_z[i] << endl;
			}

			//sub_obj = calculate_original_obj(data, initial_x, alpha);
			sub_obj = 0;
			for (int i = 0; i < data.number_customers; ++i) {
				sub_obj += exp(initial_y[i] + initial_z[i]);
			}

			cout << "Sub obj = " << std::setprecision(7) << fixed << sub_obj << endl;
			cout << "Gurobi obj = " << std::setprecision(7) << fixed << obj_val_cplex << endl;
			master_obj_val = calculate_master_obj(data, initial_x);
			cout << "Master obj = " << std::setprecision(7) << fixed << master_obj_val << endl;

			if (master_obj_val > best_obj) {
				best_obj = master_obj_val;
				best_x = initial_x;
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
	report_results << best_obj << " " << time_for_solve << endl;
	for (int j = 0; j < data.number_products; ++j)
		if (best_x[j] == 1)
			report_results << j << " ";
	report_results.close();
}

void CuttingPlaneGurobi::solve_build_in(Data data, int budget) {
	auto start = chrono::steady_clock::now(); //get start time

	vector<double> alpha(data.number_customers, -1);
	for (int i = 0; i < data.number_customers; ++i)
		for (int j = 0; j < data.number_products; ++j)
			if (data.revenue[i][j] > alpha[i])
				alpha[i] = data.revenue[i][j];

	//Calculate initial_x, initial_y, initial_z
	vector<int> initial_x = greedy(data, budget, alpha);

	vector<int> best_x = initial_x;
	double best_obj = calculate_master_obj(data, best_x);

	vector<double> initial_y = calculate_y(data, initial_x, alpha);
	vector<double> initial_z = calculate_z(data, initial_x);

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
	y = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		y[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "y_" + to_string(i));

	//cout << "Slack variables : u_i\n" << endl;
	GRBVar* u;
	u = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		u[i] = model.addVar(alpha[i] * data.no_purchase[i], alpha[i] * data.no_purchase[i] + calculate_bound_y(data, budget, i, alpha[i]), 0, GRB_CONTINUOUS, "u_" + to_string(i));

	//cout << "Slack variables : z_i\n" << endl;
	GRBVar* z;
	z = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		z[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "z_" + to_string(i));

	//cout << "Slack variables : theta_i\n" << endl;
	GRBVar* theta = 0;
	theta = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		theta[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "theta_" + to_string(i));

	for (int i = 0; i < data.number_customers; ++i)
		model.addGenConstrExp(y[i], u[i]);
	
	for (int i = 0; i < data.number_customers; ++i) {
		GRBLinExpr sum_x;
		for (int j = 0; j < data.number_products; ++j)
			sum_x += (alpha[i] - data.revenue[i][j]) * x[j] * data.utilities[i][j];
		sum_x += alpha[i] * data.no_purchase[i];
		model.addConstr(u[i] >= sum_x);
	}

	//Budget constraint
	GRBLinExpr capacity;
	for (int j = 0; j < data.number_products; ++j) {
		capacity += x[j];
	}
	model.addConstr(capacity <= budget, "ct_budget");

	//Objective
	GRBLinExpr obj;
	for (int i = 0; i < data.number_customers; ++i)
		obj += theta[i];
	model.setObjective(obj, GRB_MINIMIZE);

	auto time_before_cut = chrono::steady_clock::now();
	chrono::duration<double> before_cut = time_before_cut - start;

	double run_time = time_limit - before_cut.count();

	int num_iterative = 0;
	double stop_param = 1e-4;
	double sub_obj = 1.0;
	double obj_val_cplex = 0.0;

	//Add cut iteratively
	while (sub_obj > stop_param + obj_val_cplex) {

		//compute gradient e^{y+z} at initial_x, initial_y, initial_z and set up constraints related to theta
		for (int i = 0; i < data.number_customers; ++i)
			model.addConstr(theta[i] >= exp(initial_y[i] + initial_z[i]) * (1 + y[i] - initial_y[i] + z[i] - initial_z[i]), "ct_sub_gradient_y+z_" + to_string(i));

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
		model.set(GRB_IntParam_MIPFocus, 3);
		model.set(GRB_IntParam_FuncPieces, 1);
		model.set(GRB_DoubleParam_FuncPieceLength, 1e-2);
		model.set(GRB_IntParam_OutputFlag, 0);

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
				//cout << initial_y[i] << " " << initial_z[i] << endl;
			}

			//sub_obj = calculate_original_obj(data, initial_x, alpha);
			sub_obj = 0;
			for (int i = 0; i < data.number_customers; ++i) {
				sub_obj += exp(initial_y[i] + initial_z[i]);
			}

			cout << "Sub obj = " << std::setprecision(7) << fixed << sub_obj << endl;
			cout << "Gurobi obj = " << std::setprecision(7) << fixed << obj_val_cplex << endl;
			master_obj_val = calculate_master_obj(data, initial_x);
			cout << "Master obj = " << std::setprecision(7) << fixed << master_obj_val << endl;

			if (master_obj_val > best_obj) {
				best_obj = master_obj_val;
				best_x = initial_x;
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
	report_results << best_obj << " " << time_for_solve << endl;
	for (int j = 0; j < data.number_products; ++j)
		if (best_x[j] == 1)
			report_results << j << " ";
	report_results.close();
}