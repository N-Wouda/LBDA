#include "benders.h"

Benders::Benders(GRBEnv &env, GRBenv *c_env, Problem &problem) :
    d_nCuts(0),
    d_runTime(0),
    d_env(env),
    d_problem(problem),
    d_master(env, c_env, problem),
    d_lr(env, problem),
    d_gomory(env, problem),
    d_visited(problem.d_S),
    d_objectives(problem.d_S)
{
}
