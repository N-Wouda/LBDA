#include "deterministicequivalent.h"

DeterministicEquivalent::DeterministicEquivalent(GRBEnv &env,
                                                 Problem const &problem) :
    d_problem(problem),
    d_model(GRBModel(env)),
    d_status(status::UNSOLVED)
{
    initFirstStage(d_problem.Amat().n_rows,
                   problem.nFirstStageIntVars(),
                   problem.d_fs_leq,
                   problem.d_fs_geq,
                   problem.d_l1.memptr(),
                   problem.d_u1.memptr(),
                   problem.d_c.memptr(),
                   problem.d_b.memptr(),
                   problem.d_Amat);

    initSecondStage(d_problem.Amat().n_rows,
                    problem.d_n2,
                    problem.nSecondStageIntVars(),
                    problem.d_Wmat.n_cols,
                    problem.nScenarios(),
                    problem.d_ss_leq,
                    problem.d_ss_geq,
                    problem.d_l2.memptr(),
                    problem.d_u2.memptr(),
                    problem.d_probs.memptr(),
                    problem.d_q.memptr(),
                    problem.d_Tmat,
                    problem.d_Wmat,
                    problem.d_omega);
}
