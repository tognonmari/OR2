#ifndef EXACT_H
#define EXACT_H

#include "utils.h"
#include <cplex.h>  
#include "tsp.h"
#include <stdio.h>
#include <math.h>

#define EPS 1e-7
double dist(int i, int j, instance* inst);
void build_model(instance* inst, CPXENVptr env, CPXLPptr lp);
int xpos(int i, int j, instance* inst);
void set_init_param(CPXENVptr env, const instance* inst);
void build_sol(const double* xstar, instance* inst, int* succ, int* comp, int* ncomp);
int TSPopt(instance* inst);

#endif