#include "CuttingPlaneSolver.h"

CuttingPlaneSolver::CuttingPlaneSolver() {

}

CuttingPlaneSolver::CuttingPlaneSolver(Data data, const char* output, string log, string out_res_csv, double time_limit) {
    this->data = data;
    this->time_limit = time_limit;
    this->output = output;
    this->log = log;
    this->out_res_csv = out_res_csv;
    this->time_limit = time_limit;
}

vector<vector<double>> CuttingPlaneSolver::create_optimal_sub_intervals(Data data, vector<int> initial_x) {
    vector<vector<double>> c(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        c[i].resize(10, 0);
    for (int i = 0; i < data.number_customers; ++i)
        for (int k = 0; k < c[i].size(); ++k)
            c[i][k] = k;
    return c;
}

bool CuttingPlaneSolver::solve(Data data, vector<int> initial_x) {
    double alpha = -1;
    for (int j = 0; j < data.number_products; ++j)
        if (data.revenue[j] > alpha)
            alpha = data.revenue[j];

    //create bounds c^i_k for e^{y_i}
    vector<vector<double>> c = create_optimal_sub_intervals(data, initial_x);
    vector<int> number_sub_intervals(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        number_sub_intervals[i] = c[i].size() - 1;
    //create sigma^i_k
    vector<vector<double>> sigma (data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        for (int k = 0; k < number_sub_intervals[i]; ++k)
            sigma[i].push_back((exp(c[i][k + 1]) - exp(c[i][k])) / (c[i][k + 1] - c[i][k]));

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
                sprintf_s(var_name, "s(%d,%d)", i, k);
                constraint_s.setName(var_name);
                model.add(constraint_s);
                //cout << "ct_s_" << i << "_" << k << endl;

                IloConstraint constraint_sr;
                constraint_sr = IloConstraint(s[i][k] >= r[i][k + 1]);
                sprintf_s(var_name, "sr(%d,%d)", i, k);
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
    model.add(IloMaximize(env, obj));

    IloCplex cplex(model);
    cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 1e-8);
    double runtime = 3600;
    cplex.setParam(IloCplex::TiLim, 3600);
    cplex.exportModel("cpao.lp");
    //IloNum tol = cplex.getParam(IloCplex::EpInt);
    //ofstream logfile(log);
    //cplex.setOut(logfile);
}