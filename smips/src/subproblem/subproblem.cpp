#include "subproblem.h"

SubProblem::SubProblem(GRBEnv &env, Problem const &problem) :
    d_model(env),
    d_problem(problem)
{
    // variable types
    char vTypes[d_problem.d_n2];
    std::fill(vTypes, vTypes + d_problem.d_n2, GRB_CONTINUOUS);

    d_vars = d_model.addVars(problem.d_l2.memptr(),
                             problem.d_u2.memptr(),
                             problem.d_q.memptr(),
                             vTypes,
                             nullptr,
                             d_problem.d_n2);

    size_t const ss_leq = problem.d_ss_leq;
    size_t const ss_geq = problem.d_ss_geq;

    // constraint senses
    char senses[d_problem.d_Wmat.n_cols];
    std::fill(senses, senses + ss_leq, GRB_LESS_EQUAL);
    std::fill(senses + ss_leq, senses + ss_leq + ss_geq, GRB_GREATER_EQUAL);
    std::fill(senses + ss_leq + ss_geq,
              senses + d_problem.d_Wmat.n_cols,
              GRB_EQUAL);

    // constraint rhs
    double rhs[d_problem.d_Wmat.n_cols];
    std::fill(rhs, rhs + d_problem.d_Wmat.n_cols, 0.0);

    // constraint lhs
    auto const &Wmat = problem.d_Wmat;
    GRBLinExpr Wy[d_problem.d_Wmat.n_cols];

    for (size_t conIdx = 0; conIdx != d_problem.d_Wmat.n_cols; ++conIdx)
        Wy[conIdx].addTerms(Wmat.colptr(conIdx),
                            d_vars,
                            d_problem.d_Wmat.n_rows);

    // add constraints
    d_constrs = d_model.addConstrs(Wy,
                                   senses,
                                   rhs,
                                   nullptr,
                                   d_problem.d_Wmat.n_cols);
    d_model.update();
}
