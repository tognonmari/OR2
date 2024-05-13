#ifndef GRE_H
#define GRE_H
#include "tsp.h"
#include "utils.h"
#include "plot.h"

typedef struct {
    int* best_sol;
    int* curr_sol;
    double zbest;
    double tbest;
    double zcurr;
    int best_start;
    char check_feasibility;
    char plots_on_screen;
    char table_flag;
    char apply_opt2;
    float* gre_dist_matrix;
    FILE* pipe;
}greedy;

void gre_update_sol(instance* inst, greedy* gre, int new_best_start);
int* gre_search_min(const int* p, const int* end, const float* dist_matrix, double* min);
int* gre_compute_path(int index_first, const instance * inst, greedy * gre);
void gre_init_table(char table_flag);
void gre_fill_table(char table_flag, int start, double cost, char is_new_best);
void gre_init(instance * inst, greedy * gre);
void gre_close(greedy * gre);
void gre_solve(instance * inst, char apply_opt2);

#endif