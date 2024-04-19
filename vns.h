#ifndef VNS_H
#define VNS_H
#include "tsp.h"
#include "utils.h"
#include <stdlib.h>
#include "plot.h"

typedef struct{
    int max_kicks;
    int min_kicks;
    FILE* gnuplot_pipe;
} vns_params;

void vns(instance* inst);

vns_params parse_and_init_vns_params();

void kick(instance* inst, int* sol_to_kick); 

void vns_close_plot_pipe(char flag, FILE* gnuplotPipe);

void vns_init_plot_iter_cost(char flag, FILE* gnuplotPipe, instance* inst);
#endif

