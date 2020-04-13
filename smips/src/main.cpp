#include "main.h"

#include <cassert>

int main()
{
    Data rand(31415);

    GRBEnv env;
    env.set(GRB_IntParam_OutputFlag, 0);
    env.set(GRB_IntParam_Threads, 1);

    GRBenv *c_env;
    GRBloadenv(&c_env, nullptr);
    GRBsetintparam(c_env, "OutputFlag", 0);
    GRBsetintparam(c_env, "Threads", 1);

    Problem problem(rand, env);
    problem.ssv95(11, true, true, true);

    DeterministicEquivalent detEqv(env, problem);
    auto ptr = detEqv.solve(900.0);
    auto res = *ptr;

    std::cout << res;
    std::cout << "\ncx + Q(x) = " << problem.evaluate(res) << '\n';

    MasterProblem master{env, c_env, problem};

    Benders lshaped(master);
    LpDual lpCut{env, problem};
    ptr = lshaped.solve(lpCut);
    res = *ptr;

    std::cout << res;
    std::cout << "\ncx + Q(x) = " << problem.evaluate(res) << '\n';
    assert(problem.evaluate(res) + 50.814 <= 0.001);

    arma::vec alpha = arma::zeros(problem.Tmat().n_cols);

    Benders lbda = lshaped;
    LooseBenders lbdaCut{env, problem, alpha, 1.0};
    ptr = lbda.solve(lbdaCut);
    res = *ptr;

    std::cout << res;
    std::cout << "\ncx + Q(x) = " << problem.evaluate(res) << '\n';
    assert(problem.evaluate(res) + 59.8893 <= 0.001);

    Benders sb = lshaped;
    StrongBenders sbCut{env, problem};
    ptr = sb.solve(sbCut);
    res = *ptr;

    std::cout << res;
    std::cout << "\ncx + Q(x) = " << problem.evaluate(res) << '\n';
    assert(problem.evaluate(res) + 59.8893 <= 0.001);

    GRBfreeenv(c_env);
}
