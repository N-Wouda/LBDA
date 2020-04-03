#include "cuts/loosebenders.h"

LooseBenders::~LooseBenders()
{
    delete[] d_vars;
    delete[] d_constrs;
}
