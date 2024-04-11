#include "tsp.h"
#include "utils.h"

typedef struct{
    int max_kicks;
    int min_kicks;
    FILE* gnuplot_pipe;
} vns_params;

char opt2_move(instance* inst, int* incumbent_sol, double* incumbent_cost);

void vns(instance* inst);

vns_params parse_and_init_vns_params();

void kick(instance* inst, int* sol_to_kick); 

void vns_close_plot_pipe(char flag, FILE* gnuplotPipe);

void vns_init_plot_iter_cost(char flag, FILE* gnuplotPipe, instance* inst);

