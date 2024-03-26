#include "tsp.h"
#include "utils.h"

typedef struct{
    int max_kicks;
    int min_kicks;
} vns_params;

char opt2_move(instance* inst, int* incumbent_sol, double* incumbent_cost);

void vns(instance* inst);

vns_params parse_vns_params();

void kick(instance* inst, int* sol_to_kick);