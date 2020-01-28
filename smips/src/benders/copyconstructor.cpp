#include "benders.h"

Benders::Benders(Benders const &other) :
    d_nCuts(other.d_nCuts),
    d_runTime(other.d_runTime),
    d_env(other.d_env),
    d_problem(other.d_problem),
    d_master(other.d_master),
    d_lr(other.d_lr),
    d_gomory(other.d_gomory),
    d_visited(other.d_visited),
    d_objectives(other.d_objectives)
{
}
