#ifndef BEND_H
#define BEND_H
#include <cplex.h>
#include "exact.h"
#include "tsp.h"
#include "plot.h"

void copy_mlt_sol(multitour_sol* dest_sol, const multitour_sol* source_sol);
void ben_add_sec(CPXENVptr env, CPXLPptr lp, multitour_sol* mlt_sol, const instance* inst);
void ben_reduce_comp(char patching, CPXENVptr env, CPXLPptr lp, instance * inst, multitour_sol * mlt_sol);
void update_best_delta(float* best_delta, int* best_i, int* best_j, int* comp_to_kill, float delta_ij, int i, int j, int k2);
void ben_solve(char patching, instance * inst);
void ben_patching(const multitour_sol * curr_sol, multitour_sol * patched_sol, const instance * inst);


#endif
