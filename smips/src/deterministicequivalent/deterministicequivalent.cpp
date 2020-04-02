#include "deterministicequivalent.h"

DeterministicEquivalent::DeterministicEquivalent(GRBEnv &env, Problem &problem) :
    d_n1(problem.d_n1),
    d_model(GRBModel(env)),
    d_status(status::UNSOLVED),
    d_isMip(problem.d_p1 != 0 && problem.d_p2 != 0)
{
    initFirstStage(d_n1,
                   problem.d_p1,
                   problem.d_fs_leq,
                   problem.d_fs_geq,
                   problem.d_l1.memptr(),
                   problem.d_u1.memptr(),
                   problem.d_c.memptr(),
                   problem.d_b.memptr(),
                   problem.d_Amat);

    initSecondStage(d_n1,
                    problem.d_n2,
                    problem.d_p2,
                    problem.d_m2,
                    problem.d_S,
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
