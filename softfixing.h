#ifndef SOFTFIX
#define SOFTFIX
#include "cpx.h"
#include "tsp.h"
#include "tabu.h"
typedef struct{
    char check_feasibility;
    char plots_on_screen;
    char table_flag;
    char apply_opt2; //
    double tl_mipcall;
    int status_mipcall;
    int k;
    int k_scaling;
    int nr_call;
    char improvement;
    FILE* table_out;
    FILE* pipe;
}softfixing;

void sf_init_table(softfixing* sf);
void sf_init(softfixing* sf, instance* inst);
void sf_fill_table(softfixing* sf, int nr_call, double cost, char is_new_best);
void sf_close(softfixing* sf, instance* inst);
void sf_fast_mip_solve(softfixing* sf, instance* inst, CPXENVptr env, CPXLPptr lp);
void sf_remove_constraint(instance* inst, CPXENVptr env, CPXLPptr lp);
void sf_add_constraint(int* xheu_path, softfixing* sf, instance* inst, CPXENVptr env, CPXLPptr lp);
void sf_solve(instance* inst);

#endif