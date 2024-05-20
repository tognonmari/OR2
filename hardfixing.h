#ifndef HARDFIX
#define HARDFIX
#include "cpx.h"
#include "tsp.h"

typedef enum {
	RND
}hf_policy;

typedef struct {
    int policy;
    char check_feasibility;
    char plots_on_screen;
    char table_flag;
    double tl_mipcall;
    int status_mipcall;
    double p_fix;
    int nr_call;
    FILE* table_out;
    FILE* pipe;
}hardfixing;

void hf_init_table(hardfixing* hf);
void hf_fill_table(hardfixing* hf, int nr_call, double cost, char is_new_best);
void hf_init(hardfixing* hf, instance* inst);
void hf_close(hardfixing* hf, instance* inst);
void hf_fixing_rnd(int* xheu_path, instance* inst, double p_fix, CPXENVptr env, CPXLPptr lp);
void hf_fixing(int* xheu_path, hardfixing* hf, instance* inst, CPXENVptr env, CPXLPptr lp);
void hf_unfix(instance* inst, CPXENVptr env, CPXLPptr lp);
void hf_fast_mip_solve(hardfixing* hf, instance* inst, CPXENVptr env, CPXLPptr lp);
void hf_solve(instance* inst);
#endif