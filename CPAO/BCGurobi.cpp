#include "BCGurobi.h"
#include <chrono>
#include <algorithm>
#include <sstream>
#include <cstdlib>
#include <cassert>

BCGurobi::BCGurobi() {

}

BCGurobi::BCGurobi(Data data, double time_limit, string outfile) {
	this->data = data;
	this->time_limit = time_limit;
	this->out_res_csv = outfile;
}

CB::CB() {

}

CB::CB(GRBVar* varx, GRBVar* vary, GRBVar* varz, GRBVar* vartheta, int p, int c, vector<double> n, vector<vector<double>> u, vector<vector<double>> r) {
	cout << "Set up object CB..." << endl;
	x = varx;
	y = vary;
	z = varz;
	theta = vartheta;
	products = p;
	customers = c;
	noPay = n;
	util = u;
	ren = r;
}

vector<double> BCGurobi::calculate_y(Data data, vector<int> x, vector<double> alpha) {
	vector<double> y(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i) {
		double tmp_y = alpha[i] * data.no_purchase[i];
		for (int j = 0; j < data.number_products; ++j)
			tmp_y += (alpha[i] - data.revenue[i][j]) * x[j] * data.utilities[i][j];
		y[i] = log(tmp_y);
	}

	return y;
}

vector<double> BCGurobi::calculate_z(Data data, vector<int> x) {
	vector<double> z(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i) {
		double tmp_z = data.no_purchase[i];
		for (int j = 0; j < data.number_products; ++j)
			tmp_z += x[j] * data.utilities[i][j];
		z[i] = -log(tmp_z);
	}

	return z;
}

double  BCGurobi::calculate_original_obj(Data data, vector<int> x, vector<double> alpha) {
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

double BCGurobi::calculate_master_obj(Data data, vector<int> x) {
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

double BCGurobi::calculate_original_obj_tmp(Data data, vector<int> x, vector<double> alpha, int candidate) {
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

vector<int> BCGurobi::greedy(Data data, int budget, vector<double> alpha) {
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

	return chosen;
}

double BCGurobi::calculate_bound_y(Data data, int budget, int i, double alpha) {
	double sum = 0;
	vector<double> u(data.number_products);
	for (int j = 0; j < data.number_products; ++j)
		u[j] = (alpha - data.revenue[i][j]) * data.utilities[i][j];

	sort(u.begin(), u.end(), greater<double>());

	for (int j = 0; j < budget; ++j)
		sum += u[j];
	sum += alpha * data.no_purchase[i];

	return sum;
}

double BCGurobi::calculate_bound_z(Data data, int budget, int i) {
	double sum = 0;
	vector<double> u(data.number_products);
	for (int j = 0; j < data.number_products; ++j)
		u[j] = data.utilities[i][j];

	sort(u.begin(), u.end(), greater<double>());

	for (int j = 0; j < budget; ++j)
		sum += u[j];
	sum += data.no_purchase[i];

	return sum;
}

double BCGurobi::calculate_optimal_bound_y(Data data, int budget, int i, double alpha) {
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

void BCGurobi::solve_build_in(Data data, int budget) {
	auto start = chrono::steady_clock::now(); //get start time

	vector<double> alpha(data.number_customers, -1);
	for (int i = 0; i < data.number_customers; ++i)
		for (int j = 0; j < data.number_products; ++j)
			if (data.revenue[i][j] > alpha[i])
				alpha[i] = data.revenue[i][j];

	vector<int> initial_x = greedy(data, budget, alpha);
	vector<double> initial_y = calculate_y(data, initial_x, alpha);
	vector<double> initial_z = calculate_z(data, initial_x);

	GRBEnv env = GRBEnv(true);
	env.start();

	GRBModel model = GRBModel(env);

	model.set(GRB_IntParam_LazyConstraints, 1);

	//cout << "Decison variables : x_j\n" << endl;
	GRBVar* x;
	x = new GRBVar[data.number_products];
	for (int j = 0; j < data.number_products; ++j)
		x[j] = model.addVar(0, 1, 0, GRB_BINARY, "x_" + to_string(j));

	//cout << "Slack variables : y_i\n" << endl;
	GRBVar* y;
	y = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		y[i] = model.addVar(log(alpha[i] * data.no_purchase[i]), log(calculate_bound_y(data, budget, i, alpha[i])), 0, GRB_CONTINUOUS, "y_" + to_string(i));

	//cout << "Slack variables : u_i\n" << endl;
	GRBVar* u;
	u = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		u[i] = model.addVar(alpha[i] * data.no_purchase[i], calculate_bound_y(data, budget, i, alpha[i]), 0, GRB_CONTINUOUS, "u_" + to_string(i));

	//cout << "Slack variables : z_i\n" << endl;
	GRBVar* z;
	z = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		z[i] = model.addVar(-log(data.no_purchase[i]), -log(calculate_bound_z(data, budget, i)), 0, GRB_CONTINUOUS, "z_" + to_string(i));

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

	//Objective
	GRBLinExpr obj;
	for (int i = 0; i < data.number_customers; ++i)
		obj += theta[i];
	model.setObjective(obj, GRB_MINIMIZE);

	auto time_before_cut = chrono::steady_clock::now();
	chrono::duration<double> before_cut = time_before_cut - start;

	double run_time = time_limit - before_cut.count();

	double obj_val_cplex = 0.0;

	model.write("branch_and_cut.lp");
	model.set(GRB_DoubleParam_TimeLimit, run_time);
	model.set(GRB_IntParam_MIPFocus, 2);
	model.set(GRB_IntParam_FuncPieces, 1);
	model.set(GRB_DoubleParam_FuncPieceLength, 1e-2);
	//model.set(GRB_IntParam_OutputFlag, 0);

	CB cb = CB(x, y, z, theta, data.number_products, data.number_customers, data.no_purchase, data.utilities, data.revenue);
	model.setCallback(&cb);

	model.optimize();

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

	auto end = chrono::steady_clock::now();
	chrono::duration<double> total_time = end - start;
	time_for_solve = total_time.count();

	cout << "\nObjective value: " << setprecision(5) << master_obj_val << endl;
	cout << "Solution: ";
	for (int j = 0; j < data.number_products; ++j)
		if (sol_x[j] == 1)
			cout << j << " ";
	cout << "\nTotal time: " << time_for_solve << " seconds" << endl;

	ofstream report_results(out_res_csv, ofstream::out);
	report_results.precision(10);
	report_results << master_obj_val << " " << time_for_solve << endl;
	for (int j = 0; j < data.number_products; ++j)
		if (sol_x[j] == 1)
			report_results << j << " ";
	report_results.close();
}

void CB::callback() {
	try {
		if (where == GRB_CB_MIPSOL) {
			double* initial_x = new double[products];
			double* initial_y = new double[customers];
			double* initial_z = new double[customers];
			double* initial_theta = new double[customers];
			initial_x = getSolution(x, products);
			initial_y = getSolution(y, customers);
			initial_z = getSolution(z, customers);
			initial_theta = getSolution(theta, customers);

			int check_theta = 0;
			for (int i = 0; i < customers; ++i)
				if (initial_theta[i] < exp(initial_y[i] + initial_z[i])) {
					check_theta = 1;
					break;
				}
			
			int check_z = 0;
			for (int i = 0; i < customers; ++i) {
				double sum_x = 0;
				for (int j = 0; j < products; ++j)
					sum_x += initial_x[j] * util[i][j];
				sum_x += noPay[i];
				if (exp(-initial_z[i]) < sum_x) {
					check_z = 1;
					break;
				}
			}

			if (check_z == 1 || check_theta == 1) {
				//compute gradient e^{y+z} at initial_x, initial_y, initial_z and set up constraints related to theta
				for (int i = 0; i < customers; ++i)
					addLazy(theta[i] >= exp(initial_y[i] + initial_z[i]) * (1 + y[i] - initial_y[i] + z[i] - initial_z[i]));

				//compute gradient e^{-z} at initial_x, initial_y, initial_z and set up constraints related to e^{-z}
				for (int i = 0; i < customers; ++i) {
					GRBLinExpr sum = 0;
					for (int j = 0; j < products; ++j)
						sum += x[j] *util[i][j];
					sum += noPay[i];
					addLazy(exp(-initial_z[i]) * (1 - z[i] + initial_z[i]) <= sum);
				}
			}
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

double CB::calculate_master_obj(double* x) {
	double obj = 0;
	for (int i = 0; i < customers; ++i) {
		double ts = 0, ms = noPay[i];
		for (int j = 0; j < products; ++j) {
			ts += ren[i][j] * x[j] * util[i][j];
			ms += x[j] * util[i][j];
		}
		obj += ts / ms;
	}
	return obj;
}