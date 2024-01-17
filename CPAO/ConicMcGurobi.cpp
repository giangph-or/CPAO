#include "ConicMcGurobi.h"
#include <chrono>
#include <algorithm>

ConicMcGurobi::ConicMcGurobi() {

}

ConicMcGurobi::ConicMcGurobi(Data data, double time_limit, string outfile) {
    this->data = data;
    this->time_limit = time_limit;
    this->out_res_csv = outfile;
}

vector<int> ConicMcGurobi::find_bound_y(Data data, int i, int budget) {
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

    model.set(GRB_IntParam_OutputFlag, 0);

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

double ConicMcGurobi::calculate_sum_utility(Data data, int budget, int i) {
    vector<pair<double, int>> u(data.number_products);
    for (int j = 0; j < data.number_products; ++j) {
        u[j].first = data.utilities[i][j];
        u[j].second = j;
    }
    sort(u.begin(), u.end(), greater<pair<double, int>>());
    double sum = 0;
    int k = 0, count = 0;
    while (count < budget) {
        sum += u[count].first;
        count++;
    }
    return sum;
}

double  ConicMcGurobi::calculate_master_obj(Data data, vector<int> x) {
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

void ConicMcGurobi::solve(Data data, int budget) {
    auto start = chrono::steady_clock::now(); //get start time

    vector<double> alpha(data.number_customers, -1);
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            if (data.revenue[i][j] > alpha[i])
                alpha[i] = data.revenue[i][j];

    vector<vector<int>> bound(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        bound[i] = find_bound_y(data, i, budget);

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
    z = new GRBVar* [data.number_customers];
    for (int i = 0; i < data.number_customers; ++i)
        z[i] = new GRBVar[data.number_products];
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            z[i][j] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "z_" + to_string(i) + "_" + to_string(j));

    //cout << "Slack variables: w_i\n" << endl;
    GRBVar* w = 0;
    w = new GRBVar[data.number_customers];
    for (int i = 0; i < data.number_customers; ++i)
        w[i] = model.addVar(data.no_purchase[i], GRB_INFINITY, 0, GRB_CONTINUOUS, "w_" + to_string(i));

    //cout << "Constraints related to w and x\n" << endl;
    for (int i = 0; i < data.number_customers; ++i) {
        GRBLinExpr sum_x;
        for (int j = 0; j < data.number_products; ++j)
            sum_x += x[j] * data.utilities[i][j];
        model.addConstr(sum_x + data.no_purchase[i] == w[i]);
    }

    //cout << "Constraints related to x, w and z\n" << endl;
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            model.addQConstr(z[i][j] * w[i] >= x[j] * x[j]);

    //cout << "Constraints related to y and w\n" << endl;
    for (int i = 0; i < data.number_customers; ++i)
        model.addQConstr(y[i] * w[i] >= 1);

    //cout << "Constraints related to y and z\n" << endl;
    for (int i = 0; i < data.number_customers; ++i) {
        GRBLinExpr sum_z;
        for (int j = 0; j < data.number_products; ++j)
            sum_z += data.utilities[i][j] * z[i][j];

        model.addConstr(data.no_purchase[i] * y[i] + sum_z >= 1);
    }

    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            model.addConstr(z[i][j] <= x[j] * (1.0 / (data.no_purchase[i] + data.utilities[i][j])));

    //McCornick constraints
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            if (bound[i][j] == 1) {
                model.addConstr(z[i][j] >= x[j] * 1.0 / (data.no_purchase[i] + calculate_sum_utility(data, budget, i)));
                if(budget < data.number_products)
                    model.addConstr(z[i][j] <= y[i] - (1 - x[j]) * 1.0 / (data.no_purchase[i] + calculate_sum_utility(data, budget + 1, i) - data.utilities[i][j]));
                else
                    model.addConstr(z[i][j] <= y[i] - (1 - x[j]) * 1.0 / (data.no_purchase[i] + calculate_sum_utility(data, budget, i) - data.utilities[i][j]));
            }
            else {
                model.addConstr(z[i][j] <= y[i] - (1 - x[j]) * 1.0 / (data.no_purchase[i] + calculate_sum_utility(data, budget, i)));
                model.addConstr(z[i][j] >= x[j] * 1.0 / (data.no_purchase[i] + data.utilities[i][j] + calculate_sum_utility(data, budget - 1, i)));
            }

    //cout << "Budget constraint\n" << endl;
    GRBLinExpr capacity;
    for (int j = 0; j < data.number_products; ++j) {
        capacity += x[j];
    }
    model.addConstr(capacity <= budget);

    //cout << "Objective\n" << endl;
    GRBLinExpr obj;
    for (int i = 0; i < data.number_customers; ++i) {
        obj += alpha[i] * data.no_purchase[i] * y[i];
        for (int j = 0; j < data.number_products; ++j)
            obj += (alpha[i] - data.revenue[i][j]) * z[i][j] * data.utilities[i][j];
    }
    model.setObjective(obj, GRB_MINIMIZE);

    auto time_now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = time_now - start;

    model.set(GRB_DoubleParam_TimeLimit, time_limit - elapsed_seconds.count());
    model.set(GRB_IntParam_MIQCPMethod, 1);
    model.write("conic.lp");
    //model.set(GRB_IntParam_OutputFlag, 0);

    model.optimize();

    double obj_value;
    vector<int> x_sol(data.number_products);
    vector<double> y_sol(data.number_customers);
    vector<vector<double>> z_sol(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        z_sol[i].resize(data.number_products);

    double master_obj = 0;

    if (model.get(GRB_IntAttr_SolCount) != 0) {
        obj_value = model.get(GRB_DoubleAttr_ObjVal);
        cout << "\nResult product list: " << endl;
        for (int j = 0; j < data.number_products; ++j)
            if (x[j].get(GRB_DoubleAttr_X) > 0.5) {
                x_sol[j] = 1;
                cout << j << " ";
            }
            else x_sol[j] = 0;
        cout << endl;

        for (int i = 0; i < data.number_customers; ++i) {
            y_sol[i] = y[i].get(GRB_DoubleAttr_X);
            for (int j = 0; j < data.number_products; ++j)
                z_sol[i][j] = z[i][j].get(GRB_DoubleAttr_X);
        }

        cout << "Conic obj = " << std::setprecision(5) << fixed << obj_value << endl;
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
    report_results << master_obj << " " << time_for_solve << endl;
    for (int j = 0; j < data.number_products; ++j)
        if (x_sol[j] == 1)
            report_results << j << " ";
    report_results.close();
}