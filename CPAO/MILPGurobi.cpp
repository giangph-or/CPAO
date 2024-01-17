#include "MILPGurobi.h"
#include <chrono>
#include <algorithm>

MILPGurobi::MILPGurobi() {

}

MILPGurobi::MILPGurobi(Data data, double time_limit, string outfile) {
    this->data = data;
    this->time_limit = time_limit;
    this->out_res_csv = outfile;
}

vector<int> MILPGurobi::find_bound_y(Data data, int i, int budget) {
    GRBEnv env = GRBEnv(true);
    env.start();

    GRBModel model = GRBModel(env);

    //Decison variables: x_j = 1 if product j is chosen, 0 otherwise
    GRBVar* x;
    x = new GRBVar[data.number_products];
    for (int j = 0; j < data.number_products; ++j)
        x[j] = model.addVar(0, 1, 0, GRB_BINARY, "x_" + to_string(j));

    //Budget constraint
    GRBLinExpr capacity;
    for (int j = 0; j < data.number_products; ++j) {
        capacity += x[j];
    }
    model.addConstr(capacity <= budget);

    GRBLinExpr obj;
    for (int j = 0; j < data.number_products; ++j)
        obj += data.utilities[i][j] * x[j];
    model.setObjective(obj, GRB_MAXIMIZE);

    model.optimize();

    vector<int> x_sol(data.number_products);

    if (model.get(GRB_IntAttr_SolCount) != 0) {
        for (int j = 0; j < data.number_products; ++j)
            if (x[j].get(GRB_DoubleAttr_X) > 0.5)
                x_sol[j] = 1;
            else
                x_sol[j] = 0;
    }
    return x_sol;
}

double MILPGurobi::calculate_sum_utility(Data data, int budget, int i, int w) {
    vector<pair<double, int>> u(data.number_products);
    for (int j = 0; j < data.number_products; ++j) {
        u[j].first = data.utilities[i][j];
        u[j].second = j;
    }
    sort(u.begin(), u.end(), greater<pair<double, int>>());
    double sum = 0;
    int k = 0, count = 0;
    while (count < budget) {
        if (u[k].second != w) {
            sum += u[k].first;
            count++;
        }
        k++;
    }
    return sum;
}

double MILPGurobi::calculate_master_obj(Data data, vector<int> x) {
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

void MILPGurobi::solve(Data data, int budget) {
    auto start = chrono::steady_clock::now(); //get start time

    vector<vector<int>> bound(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        bound[i] = find_bound_y(data, i, budget);

    vector<vector<double>> y_l(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        y_l[i].resize(data.number_products);
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            if (bound[i][j] == 1) {
                y_l[i][j] = 1.0 / (data.no_purchase[i] + data.utilities[i][j] + calculate_sum_utility(data, budget - 1, i, j));
            }
            else {
                y_l[i][j] = 1.0 / (data.no_purchase[i] + calculate_sum_utility(data, budget, i, j));
            }

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
        y[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "y_" + to_string(i));

    //cout << "Slack variables: z_{ij}\n" << endl;
    GRBVar** z;
    z = new GRBVar * [data.number_customers];
    for (int i = 0; i < data.number_customers; ++i)
        z[i] = new GRBVar[data.number_products];
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            z[i][j] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "z_" + to_string(i) + "_" + to_string(j));

    //Constraints related to y and z
    for (int i = 0; i < data.number_customers; ++i) {
        GRBLinExpr sum_z;
        for (int j = 0; j < data.number_products; ++j)
            sum_z += z[i][j] * data.utilities[i][j];

        model.addConstr(sum_z + y[i] * data.no_purchase[i] == 1);
    }

    //Constraints related to x, y and z
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            model.addConstr(data.no_purchase[i] * (y[i] - z[i][j]) <= 1 - x[j]);

    //Bound z
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            model.addConstr(z[i][j] <= y[i]);

    //Constraints related to x and z
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            model.addConstr(data.no_purchase[i] * z[i][j] <= x[j]);  

    //McCornick constraints
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            if (bound[i][j] == 1) {
                model.addConstr(z[i][j] * (data.no_purchase[i] + data.utilities[i][j]) <= x[j]);
                model.addConstr(z[i][j] >= y_l[i][j] * x[j]);
            }
            else model.addConstr(z[i][j] <= y[i] - y_l[i][j] * (1 - x[j]));

    //Budget constraint
    GRBLinExpr capacity;
    for (int j = 0; j < data.number_products; ++j) {
        capacity += x[j];
    }
    model.addConstr(capacity <= budget);

    //Objective
    GRBLinExpr obj;
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            obj += data.revenue[i][j] * z[i][j] * data.utilities[i][j];
    model.setObjective(obj, GRB_MAXIMIZE);

    auto time_now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = time_now - start;

    double run_time = time_limit - elapsed_seconds.count();
    
    model.write("milp.lp");
    model.set(GRB_DoubleParam_TimeLimit, time_limit - elapsed_seconds.count());
    model.set(GRB_IntParam_OutputFlag, 0);

    model.optimize();

    vector<int> x_sol(data.number_products);
    //vector<double> y_sol(data.number_customers);
    //vector<vector<double>> z_sol(data.number_customers);
    //for (int i = 0; i < data.number_customers; ++i)
    //    z_sol[i].resize(data.number_products);

    if (model.get(GRB_IntAttr_SolCount) != 0) {
        obj_val = model.get(GRB_DoubleAttr_ObjVal);
        cout << "\nResult product list: " << endl;
        for (int j = 0; j < data.number_products; ++j)
            if (x[j].get(GRB_DoubleAttr_X) > 0.5) {
                x_sol[j] = 1;
                cout << j << " ";
            }
            else x_sol[j] = 0;
        cout << endl;

        //for (int i = 0; i < data.number_customers; ++i) {
        //    y_sol[i] = y[i].get(GRB_DoubleAttr_X);
        //    for (int j = 0; j < data.number_products; ++j)
        //        z_sol[i][j] = z[i][j].get(GRB_DoubleAttr_X);
        //}

        cout << "MILP obj = " << std::setprecision(5) << fixed << obj_val << endl;

        //check time
        auto time_now = std::chrono::steady_clock::now(); //get now time
        std::chrono::duration<double> total_time = time_now - start;
        cout << "time now: " << total_time.count() << endl;
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
    report_results << obj_val << " " << time_for_solve << endl;
    for (int j = 0; j < data.number_products; ++j)
        if (x_sol[j] == 1)
            report_results << j << " ";
    report_results.close();
}