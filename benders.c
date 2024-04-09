#include "benders.h"



void ben_add_sec(CPXENVptr env, CPXLPptr lp, multitour_sol* mlt_sol, const instance* inst) {
	int ncols = CPXgetnumcols(env, lp);
	int n = inst->nnodes;
	int izero = 0;
	int* index = (int*)malloc(ncols * sizeof(int));
	double* value = (double*)malloc(ncols * sizeof(double));
	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));
	for (int k = 1; k <= mlt_sol->ncomp; k++) {
		int nnz = 0;
		sprintf(cname[0], "SEC for %d comp", k);
		char sense = 'L';
		double rhs = 0.0;
		for (int i = 0; i < n; i++) {
			if (mlt_sol->comp[i] != k) {
				continue;
			}
			for (int j = i + 1; j < n; j++) {
				if (mlt_sol->comp[j] != k) {
					continue;
				}
				index[nnz] = xpos(i, j, inst);
				value[nnz] = 1.0;
				nnz++;
				rhs += 1.0;
			}
		}
		if(CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0])) { exit(main_error_text(-9, "CPXgetx() error")); }
	}

}
void ben_reduce_comp(CPXENVptr env, CPXLPptr lp, instance* inst, multitour_sol* mlt_sol) {
	tsp_debug(inst->verbose >= 100, 1, "Initial solution has %d connected components", mlt_sol->ncomp);
	if (mlt_sol->ncomp == 1) { return; }
	while (mlt_sol->ncomp > 1) {
		CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit - get_timer());
		ben_add_sec(env, lp, mlt_sol, inst);
		if (CPXmipopt(env, lp)) { exit(main_error_text(-9, "CPXmipopt() error")); }
		int ncols = CPXgetnumcols(env, lp);
		double* xstar = (double*)calloc(ncols, sizeof(double));
		{
			{
				if (CPXgetx(env, lp, xstar, 0, ncols - 1)) { exit(main_error_text(-9, "CPXgetx() error")); }
				tsp_debug(inst->verbose >= 100, 1, "CPXgetx SUCCESSFUL\n");
			}
		}
		build_sol((const double*)xstar, (const instance*)inst, mlt_sol->succ, mlt_sol->comp, &(mlt_sol->ncomp));
		tsp_debug(inst->verbose >= 100, 1, "Current solution has %d connected components", mlt_sol->ncomp);
		free(xstar);

	}
}
void ben_solve(instance* inst) {
	// open CPLEX model
	int error;
	char log_name[64];
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP_Problem");
	multitour_sol curr_sol;
	build_model(inst, env, lp);
	tsp_debug(inst->verbose >= 100, 1, "build_model SUCCESSFUL");
	// Cplex's parameter setting
	set_init_param(env, (const instance*)inst, log_name, sizeof(log_name));
	tsp_debug(inst->verbose >= 100, 1, "set_init_param SUCCESSFUL");
	if (CPXmipopt(env, lp)) { exit(main_error_text(-9, "CPXmipopt() error")); }
	tsp_debug(inst->verbose >= 100, 1, "CPXmipopt SUCCESSFUL");
	// use the optimal solution found by CPLEX

	int ncols = CPXgetnumcols(env, lp);
	double* xstar = (double*)calloc(ncols, sizeof(double));
	{
		{
			if (CPXgetx(env, lp, xstar, 0, ncols - 1)) { exit(main_error_text(-9, "CPXgetx() error")); }
			tsp_debug(inst->verbose >= 100, 1, "CPXgetx SUCCESSFUL\n");
		}
	}
	print_selected_arcs(inst->verbose >= 100, (const double*)xstar, (const instance*)inst);
	tsp_debug(inst->verbose >= 100, 1, "print_seleceted_arcs SUCCESSFUL");
	init_multitour_sol(&curr_sol, inst->nnodes);
	build_sol((const double*)xstar, (const instance*)inst, curr_sol.succ, curr_sol.comp, &curr_sol.ncomp);
	tsp_debug(inst->verbose >= 100, 1, "build_sol SUCCESSFUL: there are %d connected components", curr_sol.ncomp);
	printf("\n\n NCOMP = %d \n\n ", curr_sol.ncomp);
	char figure_name[64];
	generate_name(figure_name, sizeof(figure_name), "figures/ben_%d_%d_multitour.png", inst->nnodes, inst->randomseed);
	plot_multitour(inst->verbose >= 1, inst->verbose >= 100, figure_name, (const multitour_sol*)&curr_sol, inst->nodes);
	ben_reduce_comp(env,lp, inst, &curr_sol);
	free_multitour_sol(&curr_sol);
	// free and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
	tsp_debug(inst->verbose >= 100, 1, "ben_solve SUCCESSFUL");
	return 0; // or an appropriate nonzero error code
}