#include "deterministicequivalent.h"

DeterministicEquivalent::DeterministicEquivalent(GRBEnv &env, Problem &problem) :
    d_n1(problem.d_n1),
    d_model(GRBModel(env)),
    d_status(status::UNSOLVED),
    d_isMip(problem.isMixedIntegerProblem())
{
    initFirstStage(d_n1,
                   problem.nFirstStageIntVars(),
                   problem.d_fs_leq,
                   problem.d_fs_geq,
                   problem.d_l1.memptr(),
                   problem.d_u1.memptr(),
                   problem.d_c.memptr(),
                   problem.d_b.memptr(),
                   problem.d_Amat);

    initSecondStage(d_n1,
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
