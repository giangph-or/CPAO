#include "MILPSolver.h"
#include <chrono>
#include <algorithm>

MILPSolver::MILPSolver() {

}

MILPSolver::MILPSolver(Data data, double time_limit, string outfile) {
    this->data = data;
    this->time_limit = time_limit;
    this->out_res_csv = outfile;
}

vector<int> MILPSolver::find_bound_y(Data data, int i, int budget) {
    IloEnv env;
    IloModel model(env);

    //Decison variables: x_j = 1 if product j is chosen, 0 otherwise
    IloIntVarArray x(env, data.number_products);
    for (int j = 0; j < data.number_products; ++j) {
        sprintf_s(var_name, "x(%d)", j);
        x[j] = IloIntVar(env, 0, 1, var_name);
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

    IloExpr obj(env);
    for (int j = 0; j < data.number_products; ++j)
        obj += data.utilities[i][j] * x[j];
    model.add(IloMaximize(env, obj));

    IloCplex cplex(model);
    IloNum tol = cplex.getParam(IloCplex::EpInt);
    cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 1e-8);
    IloNum run_time = time_limit;
    cplex.setParam(IloCplex::TiLim, run_time);
    string log_file;
    ofstream logfile(log_file);
    cplex.setOut(logfile);

    vector<int> x_sol(data.number_products);

    if (cplex.solve()) {
        for (int j = 0; j < data.number_products; ++j)
            if (cplex.getValue(x[j]) > 1 - tol)
                x_sol[j] = 1;
            else
                x_sol[j] = 0;
    }
    return x_sol;
}

double MILPSolver::calculate_sum_utility(Data data, int budget, int i, int w) {
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

double MILPSolver::calculate_master_obj(Data data, vector<int> x) {
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

void MILPSolver::solve(Data data, int budget) {
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
        y[i] = IloNumVar(env, 0, INFINITY, var_name);
    }

    //Slack variables: z_{ij}
    IloArray<IloNumVarArray> z(env, data.number_customers);
    for (int i = 0; i < data.number_customers; ++i) {
        z[i] = IloNumVarArray(env, data.number_products);
        for (int j = 0; j < data.number_products; ++j) {
            sprintf_s(var_name, "z(%d,%d)", i, j);
            z[i][j] = IloNumVar(env, 0, INFINITY, var_name);
        }
    }

    //Constraints related to y and z
    for (int i = 0; i < data.number_customers; ++i) {
        IloExpr sum_z(env);
        for (int j = 0; j < data.number_products; ++j)
            sum_z += z[i][j] * data.utilities[i][j];

        IloConstraint constraint;
        //constraint = IloConstraint(exp(c[i][0]) + sum_r >= alpha * data.no_purchase[i] + sum_x);
        constraint = IloConstraint(sum_z + y[i] * data.no_purchase[i] == 1);
        sprintf_s(var_name, "ct_yz(%d)", i);
        constraint.setName(var_name);
        model.add(constraint);
    }

    //Constraints related to x, y and z
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j) {
            IloConstraint constraint;
            constraint = IloConstraint(data.no_purchase[i] * (y[i] - z[i][j]) <= 1 - x[j]);
            sprintf_s(var_name, "ct_xyz(%d,%d)", i, j);
            constraint.setName(var_name);
            model.add(constraint);
        }

    //Bound z
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j) {
            IloConstraint constraint;
            constraint = IloConstraint(z[i][j] <= y[i]);
            sprintf_s(var_name, "ct_z(%d,%d)", i, j);
            constraint.setName(var_name);
            model.add(constraint);
        }

    //Constraints related to x and z
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j) {
            IloConstraint constraint;
            constraint = IloConstraint(data.no_purchase[i] * z[i][j] <= x[j]);
            sprintf_s(var_name, "ct_xz(%d,%d)", i, j);
            constraint.setName(var_name);
            model.add(constraint);
        }

    //McCornick constraints
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            if (bound[i][j] == 1) {
                IloConstraint constraint;
                constraint = IloConstraint(z[i][j] * (data.no_purchase[i] + data.utilities[i][j]) <= x[j]);
                sprintf_s(var_name, "ct_zy_u1(%d)", i);
                constraint.setName(var_name);
                model.add(constraint);

                IloConstraint constraint_;
                constraint_ = IloConstraint(z[i][j] >= y_l[i][j] * x[j]);
                sprintf_s(var_name, "ct_zy_l1(%d)", i);
                constraint_.setName(var_name);
                model.add(constraint_);
            }
            else {
                IloConstraint constraint;
                constraint = IloConstraint(z[i][j] <= y[i] - y_l[i][j] * (1 - x[j]));
                sprintf_s(var_name, "ct_zy_l2(%d)", i);
                constraint.setName(var_name);
                model.add(constraint);
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
        for(int j = 0; j < data.number_products; ++j)
            obj += data.revenue[i][j] * z[i][j] * data.utilities[i][j];
    model.add(IloMaximize(env, obj));

    auto time_now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = time_now - start;

    IloCplex cplex(model);
    IloNum tol = cplex.getParam(IloCplex::EpInt);
    cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 1e-8);
    IloNum run_time = time_limit - elapsed_seconds.count();
    cplex.setParam(IloCplex::TiLim, run_time);
    cplex.setParam(IloCplex::Threads, 4);
    //cplex.exportModel("milp_ao.lp");
    string log_file;
    ofstream logfile(log_file);
    cplex.setOut(logfile);

    vector<int> x_sol(data.number_products);
    vector<double> y_sol(data.number_customers);
    vector<vector<double>> z_sol(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        z_sol[i].resize(data.number_products);

    if (cplex.solve()) {
        obj_val = cplex.getObjValue();
        cout << "\nResult product list: " << endl;
        for (int j = 0; j < data.number_products; ++j)
            if (cplex.getValue(x[j]) > 0.5) {
                x_sol[j] = 1;
                cout << j << " ";
            }
            else x_sol[j] = 0;
        cout << endl;

        for (int i = 0; i < data.number_customers; ++i) {
            y_sol[i] = cplex.getValue(y[i]);
            for(int j = 0; j < data.number_products; ++j)
                z_sol[i][j] = cplex.getValue(z[i][j]);
        }

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
        if(x_sol[j] == 1)
            report_results << j << " ";
    report_results.close();
    env.end();
}