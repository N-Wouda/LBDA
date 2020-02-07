#include "decompositions/gomory.h"

Gomory::Gomory(GRBEnv &env, Problem const &problem) :
    Relaxation(env),
    d_m2(problem.d_m2),
    d_n2(problem.d_n2),
    d_ss_leq(problem.d_ss_leq),
    d_ss_geq(problem.d_ss_geq),
    d_l2(problem.d_l2.memptr()),
    d_u2(problem.d_u2.memptr())
{
    arma::Col<char> vTypes(d_n2);
    vTypes.head(problem.d_p2).fill(GRB_INTEGER);
    vTypes.tail(d_n2 - problem.d_p2).fill(GRB_CONTINUOUS);

    d_vars = d_model.addVars(d_l2,
                             d_u2,
                             problem.d_q.memptr(),
                             vTypes.memptr(),
                             nullptr,
                             d_n2);

    arma::Col<char> senses(d_m2);
    senses.fill(GRB_GREATER_EQUAL);
    senses.head(d_ss_leq).fill(GRB_LESS_EQUAL);
    senses.tail(d_m2 - d_ss_leq - d_ss_geq).fill(GRB_EQUAL);

    GRBLinExpr lhs[d_m2];

    for (size_t conIdx = 0; conIdx != d_m2; ++conIdx)
        lhs[conIdx].addTerms(problem.d_Wmat.colptr(conIdx), d_vars, d_n2);

    arma::vec rhs = arma::zeros(d_m2);

    d_constrs = d_model.addConstrs(lhs,
                                   senses.memptr(),
                                   rhs.memptr(),
                                   nullptr,
                                   d_m2);

    d_model.update();
}
