#ifndef BEND_H
#define BEND_H
#include <cplex.h>
#include "exact.h"
#include "tsp.h"

void ben_add_sec(CPXENVptr env, CPXLPptr lp, multitour_sol* mlt_sol, const instance* inst);
	void ben_reduce_comp(CPXENVptr env, CPXLPptr lp, instance * inst, multitour_sol * mlt_sol);
void ben_solve(instance * inst);

#endif
