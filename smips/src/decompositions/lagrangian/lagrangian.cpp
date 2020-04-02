#include "decompositions/lagrangian.h"

Lagrangian::Lagrangian(GRBEnv &env, Problem const &problem) :
    Relaxation(env, problem)
{
    // adding first-stage variables (z)
    char zTypes[d_problem.d_n1];
    std::fill_n(zTypes, problem.nFirstStageIntVars(), GRB_INTEGER);
    std::fill_n(zTypes + problem.nFirstStageIntVars(),
                d_problem.d_n1 - problem.nFirstStageIntVars(),
                GRB_CONTINUOUS);

    d_z_vars = d_model.addVars(problem.d_l1.memptr(),
                               problem.d_u1.memptr(),
                               nullptr,
                               zTypes,
                               nullptr,
                               d_problem.d_n1);  // cost coeffs set by update()

    // TODO: include first-stage constraints

    // adding second-stage variables (y)
    // variable types
    char yTypes[d_problem.d_n2];
    std::fill_n(yTypes, problem.nSecondStageIntVars(), GRB_INTEGER);
    std::fill_n(yTypes + problem.nSecondStageIntVars(),
                d_problem.d_n2 - problem.nSecondStageIntVars(),
                GRB_CONTINUOUS);

    GRBVar *y_vars = d_model.addVars(problem.d_l2.memptr(),
                                     problem.d_u2.memptr(),
                                     problem.d_q.memptr(),
                                     yTypes,
                                     nullptr,
                                     d_problem.d_n2);

    size_t ss_leq = problem.d_ss_leq;
    size_t ss_geq = problem.d_ss_geq;

    // constraint senses
    char senses[d_problem.d_Tmat.n_cols];
    std::fill(senses, senses + ss_leq, GRB_LESS_EQUAL);
    std::fill(senses + ss_leq, senses + ss_leq + ss_geq, GRB_GREATER_EQUAL);
    std::fill(senses + ss_leq + ss_geq, senses + d_problem.d_Tmat.n_cols, GRB_EQUAL);

    // constraint rhs
    double rhs[d_problem.d_Tmat.n_cols];
    std::fill(rhs, rhs + d_problem.d_Tmat.n_cols, 0.0);

    // constraint lhs
    GRBLinExpr TxWy[d_problem.d_Tmat.n_cols];

    for (size_t conIdx = 0; conIdx != d_problem.d_Tmat.n_cols; ++conIdx)
    {
        TxWy[conIdx].addTerms(problem.d_Tmat.colptr(conIdx),
                              d_z_vars,
                              d_problem.d_Tmat.n_rows);

        TxWy[conIdx].addTerms(problem.d_Wmat.colptr(conIdx),
                              y_vars,
                              d_problem.d_Wmat.n_rows);
    }

    // add constraints
    d_constrs = d_model.addConstrs(TxWy, senses, rhs, nullptr, d_problem.d_Tmat.n_cols);
    d_model.update();

    delete[] y_vars;
}
