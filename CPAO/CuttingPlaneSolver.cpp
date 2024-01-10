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

bool CuttingPlaneSolver::solve(Data data, vector<int> initial_x) {
    //create bounds c^{i}_k for e^{y_i}
    vector<vector<double>> c = create_optimal_sub_intervals(data, initial_x);

    ILOSTLBEGIN
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
        r[i] = IloIntVarArray(env, c[i].size());
        for (int k = 0; k < c[i].size(); ++k) {
            sprintf_s(var_name, "r(%d,%d)", i, k);
            r[i][k] = IloIntVar(env, 0, 1, var_name);
        }
    }

    //Decision variables: s_{ik}
    IloArray<IloNumVarArray> s(env, data.number_customers);
    for (int i = 0; i < data.number_customers; ++i) {
        s[i] = IloNumVarArray(env, c[i].size());
        for (int k = 0; k < c[i].size(); ++k) {
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
}