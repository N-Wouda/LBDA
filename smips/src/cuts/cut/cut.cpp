#include "cuts/cut.h"


Cut::Cut(GRBEnv &env, Problem const &problem) :
    d_env(env),
    d_model(env),
    d_problem(problem)
{
}
