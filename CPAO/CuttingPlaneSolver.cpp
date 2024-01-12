#include "CuttingPlaneSolver.h"
#include <chrono>
#include <algorithm>

CuttingPlaneSolver::CuttingPlaneSolver() {

}

CuttingPlaneSolver::CuttingPlaneSolver(Data data, double time_limit, string outfile) {
	this->data = data;
	this->time_limit = time_limit;
	this->out_res_csv = outfile;
}

vector<int> CuttingPlaneSolver::greedy(Data data, int budget) {
	vector<int> chosen(data.number_products, 0);
	double obj = 0;
	int number_chosen = 0;

	while (number_chosen < budget) {
		double extra_ratio = 0;
		int inserted_product = -1;
		double max = 0;
		for (int j = 0; j < data.number_products; ++j) {
			if (chosen[j] == 0) {
				extra_ratio = calculate_master_obj_tmp(data, chosen, j) - obj;
				if (max < extra_ratio) {
					max = extra_ratio;
					inserted_product = j;
					break;
				}
			}
		}
		if (inserted_product != -1) {
			double distance;
			chosen[inserted_product] = 1;
			//cout << "Best-Insertion: " << inserted_product << endl;
			obj += max;
			number_chosen++;
		}
	}

	cout << "greedy solution: " << endl;
	for (int j = 0; j < data.number_products; j++)
		cout << chosen[j] << " ";
	cout << endl;
	cout << "master obj = " << obj << endl;
	cout << endl;

	return chosen;
}

double CuttingPlaneSolver::calculate_sum_utility(Data data, int budget, int i, double alpha) {
	double sum = 0;
	vector<double> u(data.number_products);
	for (int j = 0; j < data.number_products; ++j)
		u[j] = (alpha - data.revenue[i][j]) * data.utilities[i][j];

	sort(u.begin(), u.end(), greater<double>());

	for (int j = 0; j < budget; ++j)
		sum += u[j];

	return sum;
}

vector<vector<double>> CuttingPlaneSolver::create_optimal_sub_intervals(Data data, int budget, vector<double> alpha, int number_itervals) {
	vector<vector<double>> c(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i)
		c[i].resize(number_itervals + 1, 0);
	for (int i = 0; i < data.number_customers; ++i) {
		c[i][0] = log(alpha[i] * data.no_purchase[i]);
		double upper = alpha[i] * data.no_purchase[i] + calculate_sum_utility(data, budget, i, alpha[i]);
		c[i][c[i].size() - 1] = log(upper);
		for (int k = 1; k < c[i].size() - 1; ++k)
			c[i][k] = c[i][0] + k * (c[i][c[i].size() - 1] - c[i][0]) / number_itervals;
	}

	return c;
}

vector<double> CuttingPlaneSolver::calculate_y(Data data, vector<int> x, vector<double> alpha) {
	vector<double> y(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i) {
		double tmp_y = alpha[i] * data.no_purchase[i];
		for (int j = 0; j < data.number_products; ++j)
			tmp_y += (alpha[i] - data.revenue[i][j]) * x[j] * data.utilities[i][j];
		y[i] = log(tmp_y);
	}

	return y;
}

vector<double> CuttingPlaneSolver::calculate_z(Data data, vector<int> x, vector<double> alpha) {
	vector<double> z(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i) {
		double tmp_z = data.no_purchase[i];
		for (int j = 0; j < data.number_products; ++j)
			tmp_z += x[j] * data.utilities[i][j];
		z[i] = -log(tmp_z);
	}

	return z;
}

double  CuttingPlaneSolver::calculate_original_obj(Data data, vector<int> x, vector<double> alpha) {
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

double CuttingPlaneSolver::calculate_master_obj(Data data, vector<int> x) {
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

double CuttingPlaneSolver::calculate_master_obj_tmp(Data data, vector<int> x, int candidate) {
	double obj = 0;
	for (int i = 0; i < data.number_customers; ++i) {
		double ts = data.revenue[i][candidate] * data.utilities[i][candidate], 
			ms = data.no_purchase[i] + data.utilities[i][candidate];
		for (int j = 0; j < data.number_products; ++j)
			if (j != candidate) {
				ts += data.revenue[i][j] * x[j] * data.utilities[i][j];
				ms += x[j] * data.utilities[i][j];
			}

		obj += ts / ms;
	}

	return obj;
}

void CuttingPlaneSolver::solve(Data data, int budget) {
	auto start = chrono::steady_clock::now(); //get start time

	double best_obj = 0;
	vector<int> best_x;

	vector<double> alpha(data.number_customers, -1);
	for (int i = 0; i < data.number_customers; ++i)
		for (int j = 0; j < data.number_products; ++j)
			if (data.revenue[i][j] > alpha[i])
				alpha[i] = data.revenue[i][j];

	//Calculate initial_x, initial_y, initial_z
	vector<int> initial_x = greedy(data, budget);
	vector<double> initial_y = calculate_y(data, initial_x, alpha);
	vector<double> initial_z = calculate_z(data, initial_x, alpha);

	//create bounds c^i_k for e^{y_i}
	vector<vector<double>> c = create_optimal_sub_intervals(data, budget, alpha, 200);
	vector<int> number_sub_intervals(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i)
		number_sub_intervals[i] = c[i].size() - 1;

	//ILOSTLBEGIN
	IloEnv env;
	IloModel model(env);

	//Decison variables: x_j = 1 if product j is chosen, 0 otherwise
	IloIntVarArray x(env, data.number_products);
	for (int j = 0; j < data.number_products; ++j) {
		sprintf_s(var_name, "x(%d)", j);
		x[j] = IloIntVar(env, 0, 1, var_name);
	}

	//Slack variables: y_i
	IloNumVarArray y(env, data.number_customers);
	for (int i = 0; i < data.number_customers; ++i) {
		sprintf_s(var_name, "y(%d)", i);
		y[i] = IloNumVar(env, -INFINITY, INFINITY, var_name);
	}

	//Slack variables: z_i
	IloNumVarArray z(env, data.number_customers);
	for (int i = 0; i < data.number_customers; ++i) {
		sprintf_s(var_name, "z(%d)", i);
		z[i] = IloNumVar(env, -INFINITY, INFINITY, var_name);
	}

	//Decision variables: r_{ik}
	IloArray<IloIntVarArray> r(env, data.number_customers);
	for (int i = 0; i < data.number_customers; ++i) {
		r[i] = IloIntVarArray(env, number_sub_intervals[i]);
		for (int k = 0; k < number_sub_intervals[i]; ++k) {
			sprintf_s(var_name, "r(%d,%d)", i, k);
			r[i][k] = IloIntVar(env, 0, 1, var_name);
		}
	}

	//Decision variables: s_{ik}
	IloArray<IloNumVarArray> s(env, data.number_customers);
	for (int i = 0; i < data.number_customers; ++i) {
		s[i] = IloNumVarArray(env, number_sub_intervals[i]);
		for (int k = 0; k < number_sub_intervals[i]; ++k) {
			sprintf_s(var_name, "s(%d,%d)", i, k);
			s[i][k] = IloNumVar(env, 0, 1, var_name);
		}
	}

	//Slack variable: theta_i
	IloNumVarArray theta(env, data.number_customers);
	for (int i = 0; i < data.number_customers; ++i) {
		sprintf_s(var_name, "theta(%d)", i);
		theta[i] = IloNumVar(env, -INFINITY, INFINITY, var_name);
	}

	//Constraints related to e^{y_i}
	//exp(c[i][0]) = alpha * no_purchase[i] => remove form both sides
	for (int i = 0; i < data.number_customers; ++i) {
		IloExpr sum_r(env), sum_x(env);
		for (int k = 0; k < number_sub_intervals[i]; ++k) {
			//sum_r += sigma[i][k] * (c[i][k + 1] - c[i][k]) * r[i][k];
			sum_r += (exp(c[i][k + 1]) - exp(c[i][k])) * r[i][k];
		}
		for (int j = 0; j < data.number_products; ++j)
			sum_x += (alpha[i] - data.revenue[i][j]) * x[j] * data.utilities[i][j];

		IloConstraint constraint;
		constraint = IloConstraint(sum_r >= sum_x);
		sprintf_s(var_name, "ct_e(%d)", i);
		constraint.setName(var_name);
		model.add(constraint);
	}

	//Constraints related to y_i
	for (int i = 0; i < data.number_customers; ++i) {
		IloExpr sum_r(env);
		for (int k = 0; k < number_sub_intervals[i]; ++k) {
			sum_r += (c[i][k + 1] - c[i][k]) * r[i][k];
		}

		IloConstraint constraint;
		constraint = IloConstraint(y[i] == c[i][0] + sum_r);
		sprintf_s(var_name, "ct_y(%d)", i);
		constraint.setName(var_name);
		model.add(constraint);
	}

	//Constraints related to s_{ik} and r_{ik}
	for (int i = 0; i < data.number_customers; ++i) {
		for (int k = 0; k < number_sub_intervals[i]; ++k) {
			IloConstraint constraint_rs;
			constraint_rs = IloConstraint(r[i][k] >= s[i][k]);
			sprintf_s(var_name, "ct_rs(%d,%d)", i, k);
			constraint_rs.setName(var_name);
			model.add(constraint_rs);

			if (k < number_sub_intervals[i] - 1) {
				IloConstraint constraint_s;
				constraint_s = IloConstraint(s[i][k] >= s[i][k + 1]);
				sprintf_s(var_name, "ct_s(%d,%d)", i, k);
				constraint_s.setName(var_name);
				model.add(constraint_s);

				IloConstraint constraint_sr;
				constraint_sr = IloConstraint(s[i][k] >= r[i][k + 1]);
				sprintf_s(var_name, "ct_sr(%d,%d)", i, k);
				constraint_sr.setName(var_name);
				model.add(constraint_sr);
			}
		}
	}

	//Budget constraint
	IloExpr capacity(env);
	for (int j = 0; j < data.number_products; ++j) {
		capacity += x[j];
	}
	IloConstraint constraint;
	constraint = IloConstraint(capacity <= budget);
	sprintf_s(var_name, "ct_budget");
	constraint.setName(var_name);
	model.add(constraint);

	//Objective
	IloExpr obj(env);
	for (int i = 0; i < data.number_customers; ++i)
		obj += theta[i];
	model.add(IloMinimize(env, obj));

	auto time_now = chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = time_now - start;

	IloCplex cplex(model);
	IloNum tol = cplex.getParam(IloCplex::EpInt);
	cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 1e-4);
	IloNum run_time = time_limit - elapsed_seconds.count();
	//cplex.setParam(IloCplex::TiLim, run_time);
	//cplex.setParam(IloCplex::Threads, 8);
	cplex.exportModel("cpao.lp");
	string log_file;
	ofstream logfile(log_file);
	cplex.setOut(logfile);

	int num_iterative = 0;
	double stop_param = 0.0001;
	double sub_obj = 1.0;
	double obj_val_cplex = 0.0;

	//Add cut iteratively
	while (sub_obj > stop_param + obj_val_cplex) {

		//compute gradient e^{y+z} at initial_x, initial_y, initial_z and set up constraints related to theta
		for (int i = 0; i < data.number_customers; ++i) {
			IloConstraint constraint;
			constraint = IloConstraint(theta[i] >= exp(initial_y[i] + initial_z[i]) * (1 + y[i] - initial_y[i] + z[i] - initial_z[i]));
			sprintf_s(var_name, "ct_theta(%d)", i);
			constraint.setName(var_name);
			model.add(constraint);
		}

		//compute gradient e^{-z} at initial_x, initial_y, initial_z and set up constraints related to e^{-z}
		for (int i = 0; i < data.number_customers; ++i) {
			IloExpr sum(env);
			for (int j = 0; j < data.number_products; ++j)
				sum += x[j] * data.utilities[i][j];
			sum += data.no_purchase[i];
			IloConstraint constraint;
			constraint = IloConstraint(exp(-initial_z[i]) * (1 - z[i] + initial_z[i]) <= sum);
			sprintf_s(var_name, "ct_ez(%d)", i);
			constraint.setName(var_name);
			model.add(constraint);
		}

		//solve
		num_iterative++;
		cplex.setParam(IloCplex::Param::TimeLimit, run_time);
		//cplex.setParam(IloCplex::Threads, 8);
		//cplex.exportModel("cpao_.lp");

		if (cplex.solve()) {
			cout << "\nIteration " << num_iterative << ": Solved!" << endl;
			//update obj, variables
			obj_val_cplex = cplex.getObjValue();
			cout << "\nResult product list: " << endl;
			for (int j = 0; j < data.number_products; ++j)
				if (cplex.getValue(x[j]) > 1 - tol) {
					initial_x[j] = 1;
					cout << j << " ";
				}
				else initial_x[j] = 0;
			cout << endl;

			for (int i = 0; i < data.number_customers; ++i) {
				initial_y[i] = cplex.getValue(y[i]);
				initial_z[i] = cplex.getValue(z[i]);
			}

			sub_obj = 0;
			for (int i = 0; i < data.number_customers; ++i) {
				sub_obj += exp(initial_y[i] + initial_z[i]);
			}

			cout << "\nsub obj = " << std::setprecision(5) << fixed << sub_obj << endl;
			cout << "cplex obj = " << std::setprecision(5) << fixed << obj_val_cplex << endl;
			master_obj_val = calculate_master_obj(data, initial_x);
			cout << "master obj = " << std::setprecision(5) << fixed << master_obj_val << endl;

			if (master_obj_val > best_obj) {
				best_obj = master_obj_val;
				best_x = initial_x;
			}

			//check time
			auto time_now = std::chrono::steady_clock::now(); //get now time
			std::chrono::duration<double> total_time = time_now - start;
			cout << "time now: " << total_time.count() << endl;
			cout << "--- --- --- --- --- --- ---" << endl;

			if (total_time.count() > time_limit) break;
			run_time = time_limit - total_time.count();
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
	chrono::duration<double> elapsed_seconds = end - start;
	time_for_solve = elapsed_seconds.count();

	cout << "\nObjective value: " << setprecision(5) << best_obj << endl;
	cout << "Total time: " << time_for_solve << " seconds" << endl;

	ofstream report_results(out_res_csv, ofstream::out);
	report_results.precision(10);
	report_results << best_obj << " " << time_for_solve << endl;
	for (int j = 0; j < data.number_products; ++j)
		if (best_x[j] == 1)
			report_results << j << " ";
	report_results.close();
	env.end();
}