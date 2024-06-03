#ifndef TABU_H
#define TABU_H

#include "tsp.h"
#include "greedy.h"
#include "limits.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

typedef struct{
    int vertex_to_swap_1;
    int vertex_to_swap_2;
    char improvement;
    double delta;
}move;

typedef struct{
    int* best_sol;
    double zbest;
    double tbest;
    int* curr_sol;
    double zcurr;
    int iter;
    int tenure;
    char tenure_is_variable;
    int* tabu_list;
    int printing_period;
    char figure_cost_flag;
    FILE* pipe;
    FILE* data_iter_and_cost;
    move best_admissible_move;
}tabu;

int tabu_get_tenure(instance* inst,int num_iterations, int nnodes);

int* init_tabu_list(size_t n);

int* compute_tabu_init_sol(instance* inst, tabu* tabu);

void tabu_update_best(tabu* tabu, int n, int verbose);

void tabu_init(char tenure_is_variable, tabu* tabu, instance* inst);

char isTabu(tabu* tabu, int vertex);

void update_move(move* mov, int i1, int i2, double delta);

char tabu_find_best_admissible_move(tabu* tabu, instance* inst);

void tabu_update_current(tabu* tabu, int verbose);

void tabu_update_list(tabu* tabu, int verbose);

void tabu_update_tenure(tabu* tabu,instance* inst);

void tabu_update(tabu* tabu, instance* inst);

void tabu_close(char flag, FILE* gnuplotPipe, FILE* data_file);

void tabu_free(tabu* tabu);

void tabu_search(char tenure_is_variable, instance* inst);

#endif