#ifndef TSPMODEL_H
#define TSPMODEL_H

#include <cplex.h>  
#include "tsp.h"

void print_error(const char *err);
double dist(int i, int j, instance* inst);
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
int xpos(int i, int j, instance *inst);
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);
int TSPopt(instance *inst);

#endif