#include "CuttingPlaneSolver.h"
#include<chrono>

CuttingPlaneSolver::CuttingPlaneSolver() {

}

CuttingPlaneSolver::CuttingPlaneSolver(Data data, double time_limit, string outfile) {
    this->data = data;
    this->time_limit = time_limit;
    this->out_res_csv = outfile;
}

vector<vector<double>> CuttingPlaneSolver::create_optimal_sub_intervals(Data data, int number_itervals) {
    double alpha = -1;
    for (int j = 0; j < data.number_products; ++j)
        if (data.revenue[j] > alpha)
            alpha = data.revenue[j];

    vector<vector<double>> c(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        c[i].resize(number_itervals + 1, 0);
    for (int i = 0; i < data.number_customers; ++i) {
        c[i][0] = log(alpha * data.no_purchase[i]);
        double upper = log(alpha * data.no_purchase[i]);
        for (int j = 0; j < data.number_products; ++j)
            upper += (alpha - data.revenue[j]) * data.utilities[i][j];
        c[i][c[i].size() - 1] = log(upper);
        for (int k = 1; k < c[i].size() - 1; ++k)
            c[i][k] = c[i][0] + k * (c[i][c[i].size() - 1] - c[i][0]) / number_itervals;
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
    //cout << "y = " << endl;
    vector<double> y(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i) {
        double tmp_y = alpha * data.no_purchase[i];
        for (int j = 0; j < data.number_products; ++j)
            tmp_y += (alpha - data.revenue[j]) * x[j] * data.utilities[i][j];
        y[i] = log(tmp_y);
        //cout << y[i] << " ";
    }
    //cout << endl;
    return y;
}

vector<double> CuttingPlaneSolver::calculate_z(Data data, vector<int> x, double alpha) {
    //cout << "z = " << endl;
    vector<double> z(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i) {
        double tmp_z = data.no_purchase[i];
        for (int j = 0; j < data.number_products; ++j)
            tmp_z += x[j] * data.utilities[i][j];
        z[i] = -log(tmp_z);
        //cout << z[i] << " ";
    }
    //cout << endl;
    return z;
}

double  CuttingPlaneSolver::calculate_original_obj(Data data, vector<int> x, double alpha) {
    double obj = 0;
    for (int i = 0; i < data.number_customers; ++i) {
        double ts = alpha * data.no_purchase[i], ms = data.no_purchase[i];
        for (int j = 0; j < data.number_products; ++j) {
            ts += (alpha - data.revenue[j]) * x[j] * data.utilities[i][j];
            ms += x[j] * data.utilities[i][j];
        }
        obj += ts / ms;
    }
    return obj;
}

double  CuttingPlaneSolver::calculate_master_obj(Data data, vector<int> x) {
    double obj = 0;
    for (int i = 0; i < data.number_customers; ++i) {
        double ts = 0, ms = data.no_purchase[i];
        for (int j = 0; j < data.number_products; ++j) {
            ts += data.revenue[j] * x[j] * data.utilities[i][j];
            ms += x[j] * data.utilities[i][j];
        }
        obj += ts / ms;
    }
    return obj;
}

void CuttingPlaneSolver::solve(Data data) {
    double alpha = -1;
    for (int j = 0; j < data.number_products; ++j)
        if (data.revenue[j] > alpha)
            alpha = data.revenue[j];

    //Calculate initial_x, initial_y, initial_z
    vector<int> initial_x(data.number_products, 1);
    vector<double> initial_y = calculate_y(data, initial_x, alpha);
    vector<double> initial_z = calculate_z(data, initial_x, alpha);

    //create bounds c^i_k for e^{y_i}
    vector<vector<double>> c = create_optimal_sub_intervals(data, 20);
    vector<int> number_sub_intervals(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        number_sub_intervals[i] = c[i].size() - 1;

    //create sigma^i_k
    //vector<vector<double>> sigma (data.number_customers);
    //for (int i = 0; i < data.number_customers; ++i)
    //    for (int k = 0; k < number_sub_intervals[i]; ++k)
    //        sigma[i].push_back((exp(c[i][k + 1]) - exp(c[i][k])) / (c[i][k + 1] - c[i][k]));

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
    //exp(c[i][0]) = alpha * no_purchase[i] => remove form both sides
    for (int i = 0; i < data.number_customers; ++i) {
        IloExpr sum_r(env), sum_x(env);
        for (int k = 0; k < number_sub_intervals[i]; ++k) {
            //sum_r += sigma[i][k] * (c[i][k + 1] - c[i][k]) * r[i][k];
            sum_r += (exp(c[i][k + 1]) - exp(c[i][k])) * r[i][k];
        }
        for (int j = 0; j < data.number_products; ++j)
            sum_x += (alpha - data.revenue[j]) * x[j] * data.utilities[i][j];

        IloConstraint constraint;
        constraint = IloConstraint(sum_r >= sum_x);
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
    cplex.setParam(IloCplex::Threads, 8);
    cplex.exportModel("cpao.lp");
    string log_file;
    ofstream logfile(log_file);
    cplex.setOut(logfile);

    auto start = chrono::steady_clock::now(); //get start time
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;

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
        cplex.setParam(IloCplex::Threads, 8);
        //cplex.exportModel("cpao_.lp");

        cout << "\nInitial product list: " << endl;
        for (int j = 0; j < data.number_products; ++j)
            if (initial_x[j] == 1)
                cout << j << " ";
        cout << endl;

        if (cplex.solve()) {
            cout << "\nIteration " << num_iterative << ": Solved!" << endl;
            //update obj, variables
            obj_val_cplex = cplex.getObjValue();
            cout << "\nResult product list: " << endl;
            for (int j = 0; j < data.number_products; ++j)
                if (cplex.getValue(x[j]) > 0.5) {
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
            cout << "origin obj = " << std::setprecision(5) << fixed << calculate_original_obj(data, initial_x, alpha) << endl;
            cout << "master obj = " << std::setprecision(5) << fixed << calculate_master_obj(data, initial_x) << endl;

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
            end = chrono::steady_clock::now();
            elapsed_seconds = end - start;
            time_for_solve = elapsed_seconds.count();
            break;
        }
    }
    end = chrono::steady_clock::now();
    elapsed_seconds = end - start;
    time_for_solve = elapsed_seconds.count();
    master_obj_val = calculate_master_obj(data, initial_x);
    
    cout << "\nObjective value: " << setprecision(5) << master_obj_val << endl;
    cout << "Solving, it took " << time_for_solve << " seconds" << endl;

    ofstream report_results(out_res_csv, ofstream::out);
    report_results.precision(10);
    report_results << master_obj_val << " " << time_for_solve << endl;
    for (int j = 0; j < data.number_products; ++j)
        if(initial_x[j] == 1)
            report_results << j << " ";
    report_results.close();
    env.end();
}