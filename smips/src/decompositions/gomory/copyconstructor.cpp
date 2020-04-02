#include "decompositions/gomory.h"

Gomory::Gomory(Gomory const &other) :
    Relaxation(other),
    d_constrs(d_model.getConstrs()),
    d_vars(d_model.getVars())
{
}
