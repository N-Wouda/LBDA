#include "cuts/loosebenders.h"


LooseBenders::LooseBenders(GRBEnv &env,
                           Problem const &problem,
                           arma::vec const &alpha,
                           double timeLimit) :
    Cut(env, problem),
    d_alpha(alpha),
    d_gomory(env, problem),
    d_visited(problem.nScenarios()),
    d_objectives(problem.nScenarios())
{
    d_gomory.setTimeLimit(timeLimit);
}
