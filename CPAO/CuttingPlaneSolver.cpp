#include "CuttingPlaneSolver.h"
#include<chrono>

CuttingPlaneSolver::CuttingPlaneSolver() {

}

CuttingPlaneSolver::CuttingPlaneSolver(Data data, double time_limit) {
    this->data = data;
    this->time_limit = time_limit;
}

vector<vector<double>> CuttingPlaneSolver::create_optimal_sub_intervals(Data data) {
    double alpha = -1;
    for (int j = 0; j < data.number_products; ++j)
        if (data.revenue[j] > alpha)
            alpha = data.revenue[j];

    vector<vector<double>> c(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        c[i].resize(10, 0);
    for (int i = 0; i < data.number_customers; ++i) {
        c[i][0] = log(alpha * data.no_purchase[i]);
        double upper = 0;
        for (int j = 0; j < data.number_products; ++j)
            upper += (alpha - data.revenue[j]) * data.utilities[i][j];
        c[i][c[i].size() - 1] = log(upper);
        for (int k = 1; k < c[i].size() - 1; ++k)
            c[i][k] = c[i][0] + k * (c[i][c[i].size() - 1] - c[i][0]) / 9;
    }

    //cout << "c^i_k = " << endl;
    //for (int i = 0; i < data.number_customers; ++i) {
    //    for (int k = 0; k < c[i].size(); ++k)
    //        cout << c[i][k] << " ";
    //    cout << endl;
    //}

    return c;
}

vector<double> CuttingPlaneSolver::calculate_y(Data data, vector<int> x, double alpha) {
    cout << "y = " << endl;
    vector<double> y(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i) {
        double tmp_y = alpha * data.no_purchase[i];
        for (int j = 0; j < data.number_products; ++j)
            tmp_y += (alpha - data.revenue[j]) * x[j] * data.utilities[i][j];
        y[i] = log(tmp_y);
        cout << y[i] << " ";
    }
    cout << endl;
    return y;
}

vector<double> CuttingPlaneSolver::calculate_z(Data data, vector<int> x, double alpha) {
    cout << "z = " << endl;
    vector<double> z(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i) {
        double tmp_z = data.no_purchase[i];
        for (int j = 0; j < data.number_products; ++j)
            tmp_z += x[j] * data.utilities[i][j];
        z[i] = -log(tmp_z);
        cout << z[i] << " ";
    }
    cout << endl;
    return z;
}

double  CuttingPlaneSolver::calculate_original_obj(Data data, vector<int> x, double alpha) {
    double obj = 0;
    for (int i = 0; i < data.number_customers; ++i) {
        double ts = alpha * data.no_purchase[i], ms = data.no_purchase[i];
        for (int j = 0; j < data.number_products; ++j) {
            ts += (alpha - data.revenue[j]) * data.utilities[i][j];
            ms += x[i] * data.utilities[i][j];
        }
        obj += ts / ms;
    }
    return obj;
}

bool CuttingPlaneSolver::solve(Data data) {
    double alpha = -1;
    for (int j = 0; j < data.number_products; ++j)
        if (data.revenue[j] > alpha)
            alpha = data.revenue[j];

    //create bounds c^i_k for e^{y_i}
    vector<vector<double>> c = create_optimal_sub_intervals(data);
    vector<int> number_sub_intervals(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        number_sub_intervals[i] = c[i].size() - 1;
    //create sigma^i_k
    vector<vector<double>> sigma (data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        for (int k = 0; k < number_sub_intervals[i]; ++k)
            sigma[i].push_back((exp(c[i][k + 1]) - exp(c[i][k])) / (c[i][k + 1] - c[i][k]));

    //cout << "sigma^i_k = " << endl;
    //for (int i = 0; i < data.number_customers; ++i) {
    //    for (int k = 0; k < number_sub_intervals[i]; ++k)
    //        cout << sigma[i][k] << " ";
    //    cout << endl;
    //}

    //ILOSTLBEGIN
    IloEnv env;
    IloModel model(env);

    //Decison variables: x_j = 1 if product j is chosen, 0 otherwise
    IloIntVarArray x(env, data.number_products);
    for (int j = 0; j < data.number_products; ++j) {
        sprintf_s(var_name, "x(%d)", j);
        x[j] = IloIntVar(env, 0, 1, var_name);
        //cout << "x_" << j << endl;
    }

    //Slack variables: y_i
    IloNumVarArray y(env, data.number_customers);
    for (int i = 0; i < data.number_customers; ++i) {
        sprintf_s(var_name, "y(%d)", i);
        y[i] = IloNumVar(env, -INFINITY, INFINITY, var_name);
        //cout << "y_" << i << endl;
    }

    //Slack variables: z_i
    IloNumVarArray z(env, data.number_customers);
    for (int i = 0; i < data.number_customers; ++i) {
        sprintf_s(var_name, "z(%d)", i);
        z[i] = IloNumVar(env, -INFINITY, INFINITY, var_name);
        //cout << "z_" << i << endl;
    }

    //Decision variables: r_{ik}
    IloArray<IloIntVarArray> r(env, data.number_customers);
    for (int i = 0; i < data.number_customers; ++i) {
        r[i] = IloIntVarArray(env, number_sub_intervals[i]);
        for (int k = 0; k < number_sub_intervals[i]; ++k) {
            sprintf_s(var_name, "r(%d,%d)", i, k);
            r[i][k] = IloIntVar(env, 0, 1, var_name);
            //cout << "r_" << i << "_" << k << endl;
        }
    }

    //Decision variables: s_{ik}
    IloArray<IloNumVarArray> s(env, data.number_customers);
    for (int i = 0; i < data.number_customers; ++i) {
        s[i] = IloNumVarArray(env, number_sub_intervals[i]);
        for (int k = 0; k < number_sub_intervals[i]; ++k) {
            sprintf_s(var_name, "s(%d,%d)", i, k);
            s[i][k] = IloNumVar(env, 0, 1, var_name);
            //cout << "s_" << i << "_" << k << endl;
        }
    }

    //Slack variable: theta_i
    IloNumVarArray theta(env, data.number_customers);
    for (int i = 0; i < data.number_customers; ++i) {
        sprintf_s(var_name, "theta(%d)", i);
        theta[i] = IloNumVar(env, -INFINITY, INFINITY, var_name);
        //cout << "theta_" << i << endl;
    }

    //Constraints related to e^{y_i}
    for (int i = 0; i < data.number_customers; ++i) {
        IloExpr sum_r(env), sum_x(env);
        for (int k = 0; k < number_sub_intervals[i]; ++k) {
            sum_r += sigma[i][k] * (c[i][k + 1] - c[i][k]) * r[i][k];
        }
        for (int j = 0; j < data.number_products; ++j)
            sum_x += (alpha - data.revenue[j]) * x[j] * data.utilities[i][j];

        IloConstraint constraint;
        constraint = IloConstraint(exp(c[i][0]) + sum_r >= alpha * data.no_purchase[i] + sum_x);
        sprintf_s(var_name, "ct_e(%d)", i);
        constraint.setName(var_name);
        model.add(constraint);
        //cout << "ct_e_" << i << endl;
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
        //cout << "ct_y_" << i << endl;
    }

    //Constraints related to s_{ik} and r_{ik}
    for (int i = 0; i < data.number_customers; ++i) {
        for (int k = 0; k < number_sub_intervals[i]; ++k) {
            IloConstraint constraint_rs;
            constraint_rs = IloConstraint(r[i][k] >= s[i][k]);
            sprintf_s(var_name, "ct_rs(%d,%d)", i, k);
            constraint_rs.setName(var_name);
            model.add(constraint_rs);
            //cout << "ct_rs_" << i << "_" << k << endl;

            if (k < number_sub_intervals[i] - 1) {
                IloConstraint constraint_s;
                constraint_s = IloConstraint(s[i][k] >= s[i][k + 1]);
                sprintf_s(var_name, "ct_s(%d,%d)", i, k);
                constraint_s.setName(var_name);
                model.add(constraint_s);
                //cout << "ct_s_" << i << "_" << k << endl;

                IloConstraint constraint_sr;
                constraint_sr = IloConstraint(s[i][k] >= r[i][k + 1]);
                sprintf_s(var_name, "ct_sr(%d,%d)", i, k);
                constraint_sr.setName(var_name);
                model.add(constraint_sr);
                //cout << "ct_sr_" << i << "_" << k << endl;
            }
        }
    }

    //Objective
    IloExpr obj(env);
    for (int i = 0; i < data.number_customers; ++i)
        obj += theta[i];
    model.add(IloMinimize(env, obj));

    IloCplex cplex(model);
    IloNum tol = cplex.getParam(IloCplex::EpInt);
    cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 1e-8);
    IloNum run_time = time_limit;
    cplex.setParam(IloCplex::TiLim, run_time);
    cplex.setParam(IloCplex::Threads, 4);
    cplex.exportModel("cpao.lp");
    //string log_file;
    //ofstream logfile(log_file);
    //cplex.setOut(logfile);

    auto start = chrono::steady_clock::now(); //get start time
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;

    bool check_solve = false;
    int num_iterative = 0;
    double stop_param = 0.0001;
    double original_obj = 1.0;
    double obj_val_cplex = 0.0;

    //Calculate initial_x, initial_y, initial_z
    vector<int> initial_x(data.number_products, 1);
    vector<double> initial_y = calculate_y(data, initial_x, alpha);
    vector<double> initial_z = calculate_z(data, initial_x, alpha);
    
    cout << endl;
    cout << "Add cut..." << endl;

    //Add cut iteratively
    while (original_obj > stop_param + obj_val_cplex) {

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
            constraint = IloConstraint(exp(-initial_z[i]) - exp(-initial_z[i]) * (z[i] - initial_z[i]) <= sum);
            sprintf_s(var_name, "ct_ez(%d)", i);
            constraint.setName(var_name);
            model.add(constraint);
        }

        //solve
        num_iterative++;
        cplex.setParam(IloCplex::Param::TimeLimit, run_time);

        if (cplex.solve()) {
            cout << "Solved!" << endl;
            //update obj, variables
            obj_val_cplex = cplex.getObjValue();
            initial_x.resize(data.number_products);
            for (int j = 0; j < data.number_products; ++j)
                if (cplex.getValue(x[j]) > 0.5) {
                    initial_x[j] = 1;
                    cout << j << " ";
                }
                else initial_x[j] = 0;
            cout << endl;

            initial_y.resize(data.number_customers);
            initial_z.resize(data.number_customers);
            for (int i = 0; i < data.number_customers; ++i) {
                initial_y[i] = cplex.getValue(y[i]);
                initial_z[i] = cplex.getValue(z[i]);
            }
            
            //update original obj
            original_obj = calculate_original_obj(data, initial_x, alpha);

            cout << "G is: " << std::setprecision(5) << fixed << original_obj << endl;
            cout << "obj is: " << std::setprecision(5) << fixed << obj_val_cplex << endl;

            //check time
            auto time_now = std::chrono::steady_clock::now(); //get now time
            std::chrono::duration<double> total_time = time_now - start;
            cout << "time now: " << total_time.count() << endl;
            cout << "--- --- --- --- --- --- ---" << endl;

            if (total_time.count() > time_limit) break;
            run_time = time_limit - total_time.count();
        }
        else {
            cout << "No solution found..." << endl;
            end = chrono::steady_clock::now();
            elapsed_seconds = end - start;
            time_for_solve = elapsed_seconds.count();
            break;
        }
    }
    end = chrono::steady_clock::now();
    elapsed_seconds = end - start;
    time_for_solve = elapsed_seconds.count();
    obj_val = data.number_customers * alpha - original_obj;
    //obj_val = original_obj;
    
    cout << "Objective value: " << setprecision(5) << obj_val << endl;

    cout << "Solving, it took " << time_for_solve << " seconds" << endl;
    fstream report_results;
    report_results.open(out_res_csv, ios::app);
    
    report_results.close();
    env.end();

    return check_solve;
}