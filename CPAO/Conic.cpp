#include "Conic.h"
#include <chrono>
#include <algorithm>

Conic::Conic() {

}

Conic::Conic(Data data, double time_limit, string outfile) {
    this->data = data;
    this->time_limit = time_limit;
    this->out_res_csv = outfile;
}

double Conic::calculate_master_obj(Data data, vector<int> x) {
    //double obj = 0;
    //for (int i = 0; i < data.number_customers; ++i) {
    //    double ts = 0, ms = data.no_purchase[i];
    //    for (int j = 0; j < data.number_products; ++j) {
    //        ts += data.revenue[i][j] * x[j] * data.utilities[i][j];
    //        ms += x[j] * data.utilities[i][j];
    //    }
    //    obj += data.fraction[i] * ts / ms;
    //}
    //return obj;

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

double Conic::calculate_optimal_bound_denominator(Data data, int i) {
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
        obj += data.utilities[i][j] * x[j];
    model.setObjective(obj, GRB_MAXIMIZE);

    model.set(GRB_IntParam_OutputFlag, 0);

    model.optimize();

    return model.get(GRB_DoubleAttr_ObjVal) + data.no_purchase[i];
}

vector<vector<double>> Conic::calculate_bound_y_in(Data data) {
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
            while (count < data.capacity_each_set - 1) {
                lb_in[i][k] += u[count].first;
                count++;
            }

            ////General
            //double cost = 0;
            //for (int a = 0; a < u.size(); ++a)
            //    if (cost + data.fraction2[u[a].second] <= data.capacity_each_set) {
            //        lb_in[i][k] += u[a].first;
            //        cost += data.fraction2[u[a].second];
            //    }
        }
    }
    return lb_in;
}

vector<vector<double>> Conic::calculate_bound_y_notin(Data data) {
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
            //    if (cost + data.fraction2[u[a].second] <= data.capacity_each_set) {
            //        lb_notin[i][k] += u[a].first;
            //        cost += data.fraction2[u[a].second];
            //    }
        }
    }
    return lb_notin;
}

vector<vector<double>> Conic::subset_bound_y_in(Data data) {
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
            vector<int> in(10, data.sub_capacity_each_set);

            for (int a = 0; a < u.size(); ++a)
                for (int s = 0; s < 10; ++s)
                    if (in[s] > 0 && data.in_set[u[a].second][s] == 1) {
                        lb_in[i][k] += u[a].first;
                        in[s]--;
                    }
        }
    }
    return lb_in;
}

vector<vector<double>> Conic::subset_bound_y_notin(Data data) {
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
            vector<int> in(10, data.sub_capacity_each_set);
            for (int a = 0; a < u.size(); ++a)
                for (int s = 0; s < 10; ++s)
                    if (in[s] > 0 && data.in_set[u[a].second][s] == 1) {
                        lb_notin[i][k] += u[a].first;
                        in[s]--;
                    }
        }
    }
    return lb_notin;
}

void Conic::solve(Data data) {
    vector<double> alpha(data.number_customers, -1);
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            if (data.revenue[i][j] > alpha[i])
                alpha[i] = data.revenue[i][j];

    vector<double> upper_bound_denominator(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        upper_bound_denominator[i] = calculate_optimal_bound_denominator(data, i);

    vector<vector<double>> bound_in = calculate_bound_y_in(data);
    vector<vector<double>> bound_notin = calculate_bound_y_notin(data);

    //General
    //vector<vector<double>> subset_bound_in = subset_bound_y_in(data);
    //vector<vector<double>> subset_bound_notin = subset_bound_y_notin(data);

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
        upper_bound_y[i] = 1 / data.no_purchase[i];
    }

    //cout << "Slack variables : y_i\n" << endl;
    GRBVar* y;
    y = new GRBVar[data.number_customers];
    for (int i = 0; i < data.number_customers; ++i)
        y[i] = model.addVar(lower_bound_y[i], upper_bound_y[i], 0, GRB_CONTINUOUS, "y_" + to_string(i));
        //y[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "y_" + to_string(i));

    //cout << "Slack variables: z_{ij}\n" << endl;
    GRBVar** z;
    z = new GRBVar * [data.number_customers];
    for (int i = 0; i < data.number_customers; ++i)
        z[i] = new GRBVar[data.number_products];
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            //z[i][j] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "z_" + to_string(i) + "_" + to_string(j));
            z[i][j] = model.addVar(0, upper_bound_y[i], 0, GRB_CONTINUOUS, "z_" + to_string(i) + "_" + to_string(j));

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

    //cout << "McCornick constraints" << endl;
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

    cout << "Objective\n" << endl;
    GRBLinExpr obj;
    for (int i = 0; i < data.number_customers; ++i) {
        obj += data.fraction[i] * alpha[i] * data.no_purchase[i] * y[i];
        for (int j = 0; j < data.number_products; ++j)
            obj += data.fraction[i] * (alpha[i] - data.revenue[i][j]) * z[i][j] * data.utilities[i][j];
    }

    model.setObjective(obj, GRB_MINIMIZE);

    auto time_now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = time_now - start;

    model.set(GRB_DoubleParam_TimeLimit, time_limit - elapsed_seconds.count());
    model.set(GRB_IntParam_MIQCPMethod, 1);
    model.set(GRB_IntParam_Threads, 8);
    //model.set(GRB_IntParam_LazyConstraints, 1);
    //model.set(GRB_IntParam_PreCrush, 1);
    //model.set(GRB_DoubleParam_MIPGap, 1e-8);
    model.write("conic.lp");
    //model.set(GRB_IntParam_OutputFlag, 0);

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
