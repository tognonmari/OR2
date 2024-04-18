#ifndef PLOT_H
#define PLOT_H

#include "tsp.h"
#include "exact.h"
#include "utils.h"

void plot_edge(FILE* gnuplot_pipe, const point* p1, const point* p2);
void plot_tour(FILE* gnuplot_pipe, const int* succ, int start, const point* nodes);
void plot_multitour(char flag, const multitour_sol* mlt, int n, const point* nodes, const char figure_name[]);
void plot_path(char flag, const int* path, int n, double cost, const point* nodes, const char figure_name[]);


#endif
