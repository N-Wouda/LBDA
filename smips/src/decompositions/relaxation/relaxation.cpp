#include "decompositions/relaxation.h"


Relaxation::Relaxation(GRBEnv &env, Problem const &problem) :
    d_model(env),
    d_problem(problem)
{
}
