#include "CEFlog.h"
#include <chrono>
#include <algorithm>

int multiply = 10;

CEFlog::CEFlog() {

}

CEFlog::CEFlog(Data data, double time_limit, string outfile) {
	this->data = data;
	this->time_limit = time_limit;
	this->out_res_csv = outfile;
}

double CEFlog::calculate_master_obj(Data data, vector<int> x) {
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

	//double obj = 0;
	//for (int i = 0; i < data.number_customers; ++i) {
	//	double ts = 0, ms = data.no_purchase[i];
	//	for (int j = 0; j < data.number_products; ++j) {
	//		ts += data.revenue[i][j] * x[j] * data.utilities[i][j];
	//		ms += x[j] * data.utilities[i][j];
	//	}
	//	obj += ts / ms;
	//}
	//return obj;
}

double CEFlog::calculate_optimal_bound_denominator(Data data, int i) {
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
	//    cost += data.fraction2[j] * x[j];
	//model.addConstr(cost <= data.capacity_each_set, "ct_set_cap");

	GRBLinExpr obj;
	for (int j = 0; j < data.number_products; ++j)
		obj += round(data.utilities[i][j] * multiply) * x[j];
	model.setObjective(obj, GRB_MAXIMIZE);

	model.set(GRB_IntParam_OutputFlag, 0);

	model.optimize();

	return model.get(GRB_DoubleAttr_ObjVal) + round(data.no_purchase[i] * multiply);
}

void CEFlog::solve(Data data) {
	vector<double> alpha(data.number_customers, -1);
	for (int i = 0; i < data.number_customers; ++i)
		for (int j = 0; j < data.number_products; ++j)
			if (data.revenue[i][j] > alpha[i])
				alpha[i] = data.revenue[i][j];

	vector<double> upper_bound_denominator(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i)
		upper_bound_denominator[i] = calculate_optimal_bound_denominator(data, i);

	auto start = chrono::steady_clock::now(); //get start time

	GRBEnv env = GRBEnv(true);
	env.start();

	GRBModel model = GRBModel(env);

	//cout << "Decison variables : x_j\n" << endl;
	GRBVar* x;
	x = new GRBVar[data.number_products];
	for (int j = 0; j < data.number_products; ++j)
		x[j] = model.addVar(0, 1, 0, GRB_BINARY, "x_" + to_string(j));

	vector<double> upper_bound_y(data.number_customers, 0);
	vector<double> lower_bound_y(data.number_customers, 0);
	for (int i = 0; i < data.number_customers; ++i) {
		lower_bound_y[i] = 1 / upper_bound_denominator[i];
		upper_bound_y[i] = 1 / round(data.no_purchase[i] * multiply);
	}

	//cout << "Slack variables : y_i\n" << endl;
	GRBVar* y;
	y = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		y[i] = model.addVar(lower_bound_y[i], upper_bound_y[i], 0, GRB_CONTINUOUS, "y_" + to_string(i));
	//y[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "y_" + to_string(i));

	vector<int> index_w(data.number_customers, 0);
	for (int i = 0; i < data.number_customers; ++i) {
		//double frac = 0;
		//if (round(data.fraction[i] * multiply) == 0) frac = 1;
		//else frac = round(data.fraction[i] * multiply);
		double sum_uti = 0;
		for (int j = 0; j < data.number_products; ++j)
			sum_uti += round(data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * multiply);

		index_w[i] = log2(sum_uti) + 1;
		//cout << sum_uti << " " << index_w[i] << endl;
	}

	//cout << "Slack variables: w_{ik}\n" << endl;
	GRBVar** w;
	w = new GRBVar * [data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		w[i] = new GRBVar[index_w[i] + 1];
	for (int i = 0; i < data.number_customers; ++i)
		for (int k = 0; k <= index_w[i]; ++k)
			w[i][k] = model.addVar(0, 1, 0, GRB_BINARY, "w_" + to_string(i) + "_" + to_string(k));

	//cout << "Slack variables: z_{ik}\n" << endl;
	GRBVar** z;
	z = new GRBVar * [data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		z[i] = new GRBVar[index_w[i] + 1];
	for (int i = 0; i < data.number_customers; ++i)
		for (int k = 0; k <= index_w[i]; ++k)
			z[i][k] = model.addVar(0, upper_bound_y[i], 0, GRB_CONTINUOUS, "z_" + to_string(i) + "_" + to_string(k));
	//z[i][k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "z_" + to_string(i) + "_" + to_string(k));

	for (int i = 0; i < data.number_customers; ++i) {
		model.addConstr(z[i][0] == 0);
		model.addConstr(w[i][0] == 0);
	}

	//cout << "Slack variables : r_i\n" << endl;
	GRBVar* r;
	r = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		r[i] = model.addVar(1 / upper_bound_y[i], 1 / lower_bound_y[i], 0, GRB_CONTINUOUS, "r_" + to_string(i));
	//r[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "r_" + to_string(i));

	//cout << "Slack variables : s_i\n" << endl;
	GRBVar* s;
	s = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		s[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "r_" + to_string(i));

	GRBVar* t;
	t = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i)
		t[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "t_" + to_string(i));

	//cout << "Constraints related to t and x\n" << endl;
	for (int i = 0; i < data.number_customers; ++i) {
		//double frac = 0;
		//if (round(data.fraction[i] * multiply) == 0) frac = 1;
		//else frac = round(data.fraction[i] * multiply);
		GRBLinExpr sum = 0;
		for (int k = 1; k <= index_w[i]; ++k)
			sum += pow(2, k - 1) * z[i][k];
		model.addConstr(t[i] >= sum + round(data.no_purchase[i] * alpha[i] * multiply) * y[i]);
	}

	//cout << "Constraints related to r and x\n" << endl;
	for (int i = 0; i < data.number_customers; ++i) {
		GRBLinExpr sum_x = 0;
		for (int j = 0; j < data.number_products; ++j)
			sum_x += round(data.utilities[i][j] * multiply) * x[j];

		model.addConstr(r[i] == sum_x + round(data.no_purchase[i] * multiply));
	}

	//cout << "Constraints related to r, w and z\n" << endl;
	for (int i = 0; i < data.number_customers; ++i)
		for (int k = 1; k <= index_w[i]; ++k)
			model.addQConstr(z[i][k] * r[i] >= w[i][k] * w[i][k]);

	//cout << "Constraints related to y and r\n" << endl;
	for (int i = 0; i < data.number_customers; ++i)
		model.addQConstr(y[i] * r[i] >= 1);

	//cout << "Constraints related to x and w\n" << endl;
	for (int i = 0; i < data.number_customers; ++i) {
		GRBLinExpr sum_x = 0, sum_w = 0;
		for (int k = 1; k <= index_w[i]; ++k)
			sum_w += pow(2, k - 1) * w[i][k];
		//double frac = 0;
		//if (round(data.fraction[i] * multiply) == 0) frac = 1;
		//else frac = round(data.fraction[i] * multiply);
		for (int j = 0; j < data.number_products; ++j)
			sum_x += round(data.utilities[i][j] * (alpha[i] - data.revenue[i][j]) * multiply) * x[j];

		model.addConstr(sum_x == sum_w);
	}

	//cout << "Constraints related to z, y and w\n" << endl;
	for (int i = 0; i < data.number_customers; ++i)
		for (int k = 1; k <= index_w[i]; ++k) {
			model.addConstr(z[i][k] <= upper_bound_y[i] * w[i][k]);
			model.addConstr(z[i][k] <= y[i]);
			model.addConstr(upper_bound_y[i] * (y[i] - z[i][k]) <= 1 - w[i][k]);
			model.addConstr(z[i][k] >= lower_bound_y[i] * w[i][k]);
			model.addConstr(z[i][k] >= y[i] + upper_bound_y[i] * (w[i][k] - 1));
		}

	////cout << "Polymatroid constraints" << endl;
	//for (int i = 0; i < data.number_customers; ++i) {
	//	GRBQuadExpr sum = round(data.no_purchase[i] * alpha[i] * multiply);
	//	for (int k = 1; k <= index_w[i]; ++k)
	//		sum += pow(2, k - 1) * w[i][k] * w[i][k];
	//	model.addQConstr(t[i] * r[i] >= sum);
	//}

	for (int i = 0; i < data.number_customers; ++i)
		model.addQConstr(t[i] * r[i] >= s[i] * s[i]);

	//for (int i = 0; i < data.number_customers; ++i) {
	//	vector<double> gamma(index_w[i] + 1);
	//	vector<double> lambda(index_w[i] + 1);
	//	//double frac = 0;
	//	//if (round(data.fraction[i] * multiply) == 0) frac = 1;
	//	//else frac = round(data.fraction[i] * multiply);
	//	gamma[0] = round(data.no_purchase[i] * alpha[i] * multiply);
	//	for (int k = 1; k <= index_w[i]; ++k)
	//		gamma[k] = pow(2, k - 1) + gamma[k - 1];
	//	for (int k = 1; k <= index_w[i]; ++k)
	//		lambda[k] = sqrt(gamma[k]) - sqrt(gamma[k - 1]);

	//	GRBLinExpr sum = sqrt(gamma[0]);
	//	for (int k = 1; k <= index_w[i]; ++k)
	//		sum += lambda[k] * w[i][k];

	//	model.addConstr(s[i] >= sum);
	//}

	for (int i = 0; i < data.number_customers; ++i) {
		int num_polymatroid_cuts = 0;
		vector<int> permutation;
		for (int k = 0; k <= index_w[i]; ++k)
			permutation.push_back(k);

		while (num_polymatroid_cuts < data.number_products) {
			random_shuffle(permutation.begin() + 1, permutation.end());
			vector<double> gamma(index_w[i] + 1);
			vector<double> lambda(index_w[i] + 1);
			//double frac = 0;
			//if (round(data.fraction[i] * multiply) == 0) frac = 1;
			//else frac = round(data.fraction[i] * multiply);
			gamma[0] = round(data.no_purchase[i] * alpha[i] * multiply);
			for (int k = 1; k <= index_w[i]; ++k)
				gamma[permutation[k]] = pow(2, permutation[k] - 1) + gamma[permutation[k - 1]];
			for (int k = 1; k <= index_w[i]; ++k)
				lambda[permutation[k]] = sqrt(gamma[permutation[k]]) - sqrt(gamma[permutation[k - 1]]);

			GRBLinExpr sum = sqrt(gamma[0]);
			for (int k = 1; k <= index_w[i]; ++k)
				sum += lambda[permutation[k]] * w[i][permutation[k]];

			model.addConstr(s[i] >= sum);
			//cout << "Add " << num_polymatroid_cuts + 1 << " polymatroid cuts for " << i << "\n";

			num_polymatroid_cuts++;
		}
	}

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
	//    GRBLinExpr sum;
	//    for (int j = 0; j < data.number_products; ++j)
	//        if (data.in_set[j][s] == 1)
	//            sum += data.cost[j] * x[j];
	//    model.addConstr(sum <= data.sub_capacity_each_set, "ct_set_cap" + to_string(s));
	//}

	//GRBLinExpr cost;
	//for (int j = 0; j < data.number_products; ++j)
	//    cost += data.fraction2[j] * x[j];
	//model.addConstr(cost <= data.capacity_each_set, "ct_set_cap");

	//cout << "Objective\n" << endl;
	GRBLinExpr obj;
	for (int i = 0; i < data.number_customers; ++i)
		obj += t[i];

	model.setObjective(obj, GRB_MINIMIZE);

	auto time_now = chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = time_now - start;

	model.set(GRB_DoubleParam_TimeLimit, time_limit - elapsed_seconds.count());
	model.set(GRB_IntParam_MIQCPMethod, 1);
	model.set(GRB_IntParam_Threads, 8);
	//model.set(GRB_DoubleParam_MIPGap, 1e-6);
	//model.set(GRB_IntParam_LazyConstraints, 1);
	//model.set(GRB_IntParam_PreCrush, 1);
	model.write("ceflog.lp");
	//model.set(GRB_IntParam_OutputFlag, 0);

	//CBPolymatroid cb = CBPolymatroid(s, w, data, index_w, alpha);
	//model.setCallback(&cb);

	model.optimize();

	vector<int> x_sol(data.number_products);

	double master_obj = 0;

	if (model.get(GRB_IntAttr_SolCount) != 0) {
		cout << "\nResult product list: " << endl;
		for (int j = 0; j < data.number_products; ++j)
			if (x[j].get(GRB_DoubleAttr_X) > 0.5) {
				x_sol[j] = 1;
				cout << j << " ";
			}
			else x_sol[j] = 0;
		cout << endl;

		cout << "Conic obj = " << std::setprecision(5) << fixed << model.get(GRB_DoubleAttr_ObjVal) << endl;
		master_obj = calculate_master_obj(data, x_sol);
		cout << "Master obj = " << std::setprecision(5) << fixed << master_obj << endl;

		//check time
		auto time_now = std::chrono::steady_clock::now(); //get now time
		std::chrono::duration<double> total_time = time_now - start;
		cout << "Time now: " << total_time.count() << endl;
		cout << "--- --- --- --- --- --- ---" << endl;
	}
	else {
		cout << "No solution found..." << endl;
		auto end = chrono::steady_clock::now();
		elapsed_seconds = end - start;
		time_for_solve = elapsed_seconds.count();
	}

	auto end = chrono::steady_clock::now();
	elapsed_seconds = end - start;
	time_for_solve = elapsed_seconds.count();

	ofstream report_results(out_res_csv, ofstream::out);
	report_results.precision(10);
	report_results << model.get(GRB_DoubleAttr_ObjVal) << " " << master_obj << " " << model.get(GRB_DoubleAttr_MIPGap) << " " << time_for_solve << endl;
	for (int j = 0; j < data.number_products; ++j)
		if (x_sol[j] == 1)
			report_results << j << " ";
	report_results.close();
}

CBPolymatroid::CBPolymatroid() {

}

CBPolymatroid::CBPolymatroid(GRBVar* vars, GRBVar** varw, Data d, vector<int> i, vector<double> a) {
	s = vars;
	w = varw;
	data = d;
	index_w = i;
	alpha = a;
}

void CBPolymatroid::callback() {
	try {
		if (where == GRB_CB_MIPSOL) {
			for (int i = 0; i < data.number_customers; ++i) {
				int num_polymatroid_cuts = 0;
				vector<int> permutation;
				for (int k = 0; k <= index_w[i]; ++k)
					permutation.push_back(k);

				while (num_polymatroid_cuts < 10) {
					random_shuffle(permutation.begin() + 1, permutation.end());
					vector<double> gamma(index_w[i] + 1);
					vector<double> lambda(index_w[i] + 1);
					//double frac = 0;
					//if (round(data.fraction[i] * multiply) == 0) frac = 1;
					//else frac = round(data.fraction[i] * multiply);
					gamma[0] = round(data.no_purchase[i] * alpha[i] * multiply);
					for (int k = 1; k <= index_w[i]; ++k)
						gamma[permutation[k]] = pow(2, permutation[k] - 1) + gamma[permutation[k - 1]];
					for (int k = 1; k <= index_w[i]; ++k)
						lambda[permutation[k]] = sqrt(gamma[permutation[k]]) - sqrt(gamma[permutation[k - 1]]);

					GRBLinExpr sum = sqrt(gamma[0]);
					for (int k = 1; k <= index_w[i]; ++k)
						sum += lambda[permutation[k]] * w[i][permutation[k]];

					addLazy(s[i] >= sum);
					cout << "Add polymatroid cut for " << i << "\n";

					num_polymatroid_cuts++;
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
