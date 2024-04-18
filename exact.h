#ifndef EXACT_H
#define EXACT_H

#include "utils.h"
#include <cplex.h>  
#include "tsp.h"
#include <stdio.h>
#include <math.h>
#include "vns.h"

typedef struct {
	int* succ; // It contains the sequence of indices of the nodes involved in $xstar, observe that succ always starts with index 0. (succ
	int ncomp; // Number of connected components.
	int* comp; //comp[i] indicate the index of the connected component of node with index i;
	double z;
} multitour_sol;
#define EPS 1e-7
void handleCPXResult(int flag, int result, char* format);
void cpx_convert_succ_in_path(const multitour_sol* mlt, int* path, int n);
void cpx_convert_path_in_succ(const int* path, multitour_sol* mlt, int n);
int cpx_update_best(char flag, instance* inst, CPXENVptr env, CPXLPptr lp, const multitour_sol* sol);
double dist(int i, int j, instance* inst);
void build_model(instance* inst, CPXENVptr env, CPXLPptr lp);
int xpos(int i, int j, instance* inst);
void set_init_param(CPXENVptr env, const instance* inst, char* log_name, size_t size);
void init_multitour_sol(multitour_sol* sol, int n);
void free_multitour_sol(multitour_sol* sol);
void build_sol(const double* xstar, instance* inst, int* succ, int* comp, int* ncomp);
void print_selected_arcs(char flag, const double* xstar, const instance* inst);
//void plot_multitour(char figure_flag, char debug_flag, const char figure_name[], const multitour_sol* sol, const point* points);
int TSPopt(instance* inst);

#endif