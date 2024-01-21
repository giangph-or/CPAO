#include "ConicMcGeneral.h"
#include <chrono>
#include <algorithm>

ConicMcGeneral::ConicMcGeneral() {

}

ConicMcGeneral::ConicMcGeneral(Data data, double time_limit, string outfile) {
    this->data = data;
    this->time_limit = time_limit;
    this->out_res_csv = outfile;
}

double ConicMcGeneral::lb_in_total_cap(Data data, int i, int k) {
    double current_capacity = data.cost[k];

    vector<pair<double, int>> u(data.number_products);
    for (int j = 0; j < data.number_products; ++j)
        u[j] = make_pair(data.utilities[i][j], j);

    sort(u.begin(), u.end(), greater<pair<double, int>>());
    for (int j = 0; j < data.number_products; ++j)
        if (u[j].second == k) {
            u.erase(u.begin() + j);
            break;
        }

    double total_cap_obj = data.utilities[i][k];
    int count = 0;
    while (count < data.number_products) {
        if (current_capacity + data.cost[u[count].second] <= data.total_capacity) {
            current_capacity += data.cost[u[count].second];
            total_cap_obj += u[count].first;
        }
        count++;
    }

    return total_cap_obj;
}

double ConicMcGeneral::lb_in_set_cap(Data data, int i, int k) {
    double set_cap_obj = data.utilities[i][k];

    vector<pair<double, int>> u(data.number_products);
    for (int j = 0; j < data.number_products; ++j)
        u[j] = make_pair(data.utilities[i][j], j);

    sort(u.begin(), u.end(), greater<pair<double, int>>());
    for (int j = 0; j < data.number_products; ++j)
        if (u[j].second == k) {
            u.erase(u.begin() + j);
            break;
        }

    vector<int> set_capacity(data.number_sets, 0);
    set_capacity[data.in_set[k]] = 1;
    int count = 0;
    while (count < data.number_products) {
        if (set_capacity[data.in_set[count]] + 1 <= data.capacity_each_set) {
            set_capacity[data.in_set[count]]++;
            set_cap_obj += u[count].first;
        }
        count++;
    }

    return set_cap_obj;
}

double ConicMcGeneral::lb_notin_total_cap(Data data, int i, int k) {
    double current_capacity = 0;

    vector<pair<double, int>> u(data.number_products);
    for (int j = 0; j < data.number_products; ++j)
        u[j] = make_pair(data.utilities[i][j], j);

    sort(u.begin(), u.end(), greater<pair<double, int>>());
    for (int j = 0; j < data.number_products; ++j)
        if (u[j].second == k) {
            u.erase(u.begin() + j);
            break;
        }

    double total_cap_obj = 0;
    int count = 0;
    while (count < data.number_products) {
        if (current_capacity + data.cost[u[count].second] <= data.total_capacity) {
            current_capacity += data.cost[u[count].second];
            total_cap_obj += u[count].first; 
        }
        count++;
    }

    return total_cap_obj;
}

double ConicMcGeneral::lb_notin_set_cap(Data data, int i, int k) {
    double set_cap_obj = 0;

    vector<pair<double, int>> u(data.number_products);
    for (int j = 0; j < data.number_products; ++j)
        u[j] = make_pair(data.utilities[i][j], j);

    sort(u.begin(), u.end(), greater<pair<double, int>>());
    for (int j = 0; j < data.number_products; ++j)
        if (u[j].second == k) {
            u.erase(u.begin() + j);
            break;
        }

    vector<int> set_capacity(data.number_sets, 0);
    int count = 0;
    while (count < data.number_products) {
        if (set_capacity[data.in_set[count]] + 1 <= data.capacity_each_set) {
            set_capacity[data.in_set[count]]++;
            set_cap_obj += u[count].first;
        }
        count++;
    }

    return set_cap_obj;
}

vector<vector<double>> ConicMcGeneral::lb_in(Data data) {
    vector<vector<double>> lb(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        lb[i].resize(data.number_products, 0);
    
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j) {
            double lb_total = lb_in_total_cap(data, i, j);
            double lb_cap = lb_in_set_cap(data, i, j);
            if (lb_total <= lb_cap)
                lb[i][j] = lb_total;
            else
                lb[i][j] = lb_cap;
        }

    return lb;
}

vector<vector<double>> ConicMcGeneral::lb_notin(Data data) {
    vector<vector<double>> lb(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        lb[i].resize(data.number_products, 0);

    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j) {
            double lb_total = lb_notin_total_cap(data, i, j);
            double lb_cap = lb_notin_set_cap(data, i, j);
            if (lb_total <= lb_cap)
                lb[i][j] = lb_total;
            else
                lb[i][j] = lb_cap;
        }

    return lb;
}

double ConicMcGeneral::calculate_master_obj(Data data, vector<int> x) {
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

void ConicMcGeneral::solve(Data data) {
    auto start = chrono::steady_clock::now(); //get start time

    vector<double> alpha(data.number_customers, -1);
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            if (data.revenue[i][j] > alpha[i])
                alpha[i] = data.revenue[i][j];

    vector<vector<double>> bound(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        bound[i].resize(data.number_products);

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

    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            model.addConstr(z[i][j] >= y[i] - (1 - x[j]) * (1.0 / data.no_purchase[i]));

    vector<vector<double>> bound_in = lb_in(data);
    vector<vector<double>> bound_notin = lb_notin(data);

    //McCornick constraints
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j){
            model.addConstr(z[i][j] >= x[j] * 1.0 / (data.no_purchase[i] + bound_in[i][j]));
            model.addConstr(z[i][j] <= y[i] - (1 - x[j]) * 1.0 / (data.no_purchase[i] + bound_notin[i][j]));
        }

    //cout << "Budget constraint\n" << endl;
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
    //model.set(GRB_DoubleParam_MIPGap, 1e-3);
    model.write("conic.lp");
    //model.set(GRB_IntParam_OutputFlag, 0);
    //model.set(GRB_IntParam_Threads, 1);

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
    report_results << obj_value << " " << master_obj << " " << time_for_solve << endl;
    for (int j = 0; j < data.number_products; ++j)
        if (x_sol[j] == 1)
            report_results << j << " ";
    report_results.close();
}