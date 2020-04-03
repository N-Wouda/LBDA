#include "cuts/strongbenders.h"


StrongBenders::StrongBenders(GRBEnv &env, Problem const &problem) :
    Cut(env, problem)
{
    auto const &Amat = problem.Amat();

    // adding first-stage variables (z)
    char zTypes[Amat.n_rows];
    std::fill_n(zTypes, problem.nFirstStageIntVars(), GRB_INTEGER);
    std::fill_n(zTypes + problem.nFirstStageIntVars(),
                Amat.n_rows - problem.nFirstStageIntVars(),
                GRB_CONTINUOUS);

    d_z_vars = d_model.addVars(problem.d_l1.memptr(),
                               problem.d_u1.memptr(),
                               nullptr,
                               zTypes,
                               nullptr,
                               Amat.n_rows);

    // TODO: include first-stage constraints

    // adding second-stage variables (y)
    auto const &Wmat = problem.Wmat();
    auto const &Tmat = d_problem.Tmat();

    // variable types
    char yTypes[Wmat.n_rows];
    std::fill_n(yTypes, problem.nSecondStageIntVars(), GRB_INTEGER);
    std::fill_n(yTypes + problem.nSecondStageIntVars(),
                Wmat.n_rows - problem.nSecondStageIntVars(),
                GRB_CONTINUOUS);

    GRBVar *y_vars = d_model.addVars(problem.d_l2.memptr(),
                                     problem.d_u2.memptr(),
                                     problem.d_q.memptr(),
                                     yTypes,
                                     nullptr,
                                     Wmat.n_rows);

    size_t ss_leq = problem.d_ss_leq;
    size_t ss_geq = problem.d_ss_geq;

    // constraint senses
    char senses[Tmat.n_cols];
    std::fill(senses, senses + ss_leq, GRB_LESS_EQUAL);
    std::fill(senses + ss_leq, senses + ss_leq + ss_geq, GRB_GREATER_EQUAL);
    std::fill(senses + ss_leq + ss_geq, senses + Tmat.n_cols, GRB_EQUAL);

    // constraint rhs
    double rhs[Tmat.n_cols];
    std::fill(rhs, rhs + Tmat.n_cols, 0.0);

    // constraint lhs
    GRBLinExpr TxWy[Tmat.n_cols];

    for (size_t conIdx = 0; conIdx != Tmat.n_cols; ++conIdx)
    {
        TxWy[conIdx].addTerms(Tmat.colptr(conIdx), d_z_vars, Tmat.n_rows);
        TxWy[conIdx].addTerms(Wmat.colptr(conIdx), y_vars, Wmat.n_rows);
    }

    // add constraints
    d_constrs = d_model.addConstrs(TxWy, senses, rhs, nullptr, Tmat.n_cols);
    d_model.update();

    delete[] y_vars;
}
