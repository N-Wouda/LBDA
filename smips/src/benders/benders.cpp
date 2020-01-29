#include "benders.h"

Benders::Benders(GRBEnv &env, GRBenv *c_env, Problem &problem) :
    d_nCuts(0),
    d_runTime(0),
    d_env(env),
    d_problem(problem),
    d_master(env, c_env, problem)
{
}
