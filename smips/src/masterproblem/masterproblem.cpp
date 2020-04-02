#include "masterproblem.h"

MasterProblem::MasterProblem(GRBEnv &env, GRBenv *c_env, Problem &problem) :
    d_n1(problem.d_n1),
    d_nSlacks(problem.d_fs_leq + problem.d_fs_geq)
{
    GRBnewmodel(c_env,  // The C API gives access to advanced simplex routines.
                &d_cmodel,
                nullptr,
                0,
                nullptr,
                nullptr,
                nullptr,
                nullptr,
                nullptr);

    GRBaddvar(d_cmodel,  // theta
              0,
              nullptr,
              nullptr,
              1.0,
              problem.d_L,
              arma::datum::inf,
              GRB_CONTINUOUS,
              nullptr);

    char vtypes[d_n1];
    std::fill_n(vtypes, problem.d_p1, GRB_INTEGER);
    std::fill(vtypes + problem.d_p1, vtypes + d_n1, GRB_CONTINUOUS);

    GRBaddvars(d_cmodel,  // first-stage (x) variables.
               d_n1,
               0,
               nullptr,
               nullptr,
               nullptr,
               problem.d_c.memptr(),
               problem.d_l1.memptr(),
               problem.d_u1.memptr(),
               vtypes,
               nullptr);

    arma::Col<int> cind = arma::ones<arma::Col<int>>(d_n1);

    for (size_t con = 0; con != problem.d_Amat.n_cols; ++con)  // constraints
        GRBaddconstr(d_cmodel,
                     d_n1,
                     cind.memptr(),
                     problem.d_Amat.colptr(con),
                     GRB_EQUAL,
                     problem.d_b(con),
                     nullptr);

    // Add slack variables.
    arma::Col<int> vbeg = arma::linspace<arma::Col<int>>(0,
                                                         d_nSlacks,
                                                         d_nSlacks);

    arma::vec vval(d_nSlacks);
    vval.head(problem.d_fs_leq).fill(1);
    vval.tail(d_nSlacks - problem.d_fs_leq).fill(-1);

    GRBaddvars(d_cmodel,
               d_nSlacks,
               d_nSlacks,
               vbeg.memptr(),
               vbeg.memptr(),
               vval.memptr(),
               nullptr,
               nullptr,
               nullptr,
               nullptr,
               nullptr);

    GRBupdatemodel(d_cmodel);
}
