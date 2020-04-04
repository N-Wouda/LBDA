#include "problem.h"

Problem::Problem(Data &generator, GRBEnv &env) :
    d_gen(generator),
    d_sub(env),
    d_isSubProblemInitialised(false),
    d_L(0)
{
}