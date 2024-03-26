#ifndef TABU_H
#define TABU_H

#include "tsp.h"
#include "limits.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#define MAX_ITER 500
#define Y_RANGE_MIN 125000
#define Y_RANGE_MAX 140000
#define TENURE 30

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
    int* tabu_list;
    int iter_stop;
    FILE* pipe;
    move best_admissible_move;
}tabu;

int* init_tabu_list(size_t n);

int* compute_tabu_init_sol(instance* inst, tabu* tabu);

void tabu_update_best(tabu* tabu, int n);

void tabu_init(tabu* tabu, instance* inst);

char isTabu(tabu* tabu, int vertex);

void update_move(move* mov, int i1, int i2, double delta);

char tabu_find_best_admissible_move(tabu* tabu, instance* inst);

void tabu_update_current(tabu* tabu);

void tabu_update_list(tabu* tabu);

void tabu_update(tabu* tabu, int n);

void tabu_free(tabu* tabu);

void tabu_search(instance* inst);

#endif