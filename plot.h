#ifndef PLOT_H
#define PLOT_H

#include "tsp.h"
#include "cpx.h"
#include "utils.h"
#include "malloc.h"
#include "string.h"

void plot_edge(FILE* gnuplot_pipe, const point* p1, const point* p2);
void plot_tour(FILE* gnuplot_pipe, const int* succ, int start, const point* nodes);
void plot_multitour(char flag, const multitour_sol* mlt, int n, const point* nodes, const char figure_name[]);
void plot_path(char flag, const int* path, int n, double cost, const point* nodes, const char figure_name[]);
void fprintf_center(FILE* file, const char text[], int field_width);
void make_first_row(char flag, FILE* file, int num_cols, const char* text);
void make_last_row(char flag, FILE* file, int num_cols);
void make_table_row(char flag, FILE* file, int num_cols, char* text[]);



#endif
