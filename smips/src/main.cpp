#include "main.h"


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
    size_t n1 = problem.d_n1;

    double *x;
    DeqForm DEF(env, problem);
    DEF.solve(900.0);
    x = DEF.d_xVals;
    for (size_t var = 0; var != n1; ++var)
        std::cout << x[var] << ' ';
    std::cout << "\ncx + Q(x) = " << problem.evaluate(x) << '\n';

    MasterProblem master{env, c_env, problem};

    Benders lshaped(master);
    LpDual lpCut{env, problem};
    auto ptr = lshaped.solve(lpCut);
    auto res = *ptr;

    std::cout << res;
    std::cout << "\ncx + Q(x) = " << problem.evaluate(res.memptr()) << '\n';

    arma::vec alpha = arma::zeros(problem.d_m2);

    Benders lbda = lshaped;
    LooseBenders lbdaCut{env, problem, alpha, 1.0};
    ptr = lbda.solve(lbdaCut);
    res = *ptr;

    std::cout << res;
    std::cout << "\ncx + Q(x) = " << problem.evaluate(res.memptr()) << '\n';

    // TODO fix valgrind here (it's probably a small issue).
    Benders sb = lshaped;
    StrongBenders sbCut{env, problem};
    ptr = sb.solve(sbCut);
    res = *ptr;

    std::cout << res;
    std::cout << "\ncx + Q(x) = " << problem.evaluate(res.memptr()) << '\n';

    GRBfreeenv(c_env);
}
