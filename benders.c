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
		double rhs = -1.0;
		for (int i = 0; i < n; i++) {
			if (mlt_sol->comp[i] != k) {
				continue;
			}
			rhs += 1.0;
			for (int j = i + 1; j < n; j++) {
				if (mlt_sol->comp[j] != k) {
					continue;
				}
				index[nnz] = xpos(i, j, inst);
				value[nnz] = 1.0;
				nnz++;
			}
		}
		if(CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0])) { exit(main_error_text(-9, "CPXaddrows() error")); }
	}
	free(cname[0]);
	free(cname);
}
/*
* Add SEC constraints to the model for each connected component, solve the model
* and update the solution. The previous procedure is repeated until time limit is exceeded
* or the nr of connected component is reduced to one.
* If the function ends for exceeding time limit -> OV message with the current LB of objval
*/
void ben_reduce_comp(CPXENVptr env, CPXLPptr lp, instance* inst, multitour_sol* mlt_sol) {
	int n_call_ben_add_sec = 0;
	int tot_add_sec = 0;
	if (mlt_sol->ncomp == 1) { return; }
	while (mlt_sol->ncomp > 1) {
		if (is_time_limit_exceeded(inst->timelimit)) {
			tsp_debug(inst->verbose >= 1, 1, "Ending ben_reduce_comp for exceeding time limit");
			break;
		}
		char text[256];
		int init_nr_cons = CPXgetnumrows(env, lp);
		tsp_debug(inst->verbose >= 100, 1, "------- Try #%d to reduce with Benders -------", ++n_call_ben_add_sec);
		tsp_debug(inst->verbose >= 100, 1, "Initial solution has %d connected components", mlt_sol->ncomp);
		CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit - get_timer());
		ben_add_sec(env, lp, mlt_sol, inst);
		int curr_add_sec = CPXgetnumrows(env, lp) - init_nr_cons;
		tot_add_sec += curr_add_sec;
		tsp_debug(inst->verbose >= 100, 1, "#%d SEC cons has been added to the problem", curr_add_sec);
		if (CPXmipopt(env, lp)) { exit(main_error_text(-9, "CPXmipopt() error")); }
		sprintf(text, "CPXResult for CPXmipopt (with #%d SEC)", tot_add_sec);
		handleCPXResult(inst->verbose > 1, CPXgetstat(env, lp), text);
		int ncols = CPXgetnumcols(env, lp);
		double* xstar = (double*)calloc(ncols, sizeof(double));
		{
			{
				if (CPXgetx(env, lp, xstar, 0, ncols - 1)) { exit(main_error_text(-9, "CPXgetx() error: Impossible to get x because \nno feasible solution to the problem has been found")); }
			}
		}
		build_sol((const double*)xstar, (const instance*)inst, mlt_sol->succ, mlt_sol->comp, &(mlt_sol->ncomp));
		if (CPXgetobjval(env, lp, &(mlt_sol->z))) { exit(main_error_text(-9, "CPXgetobjval() error")); }
		tsp_debug(inst->verbose >= 100, 1, "Current solution has %d connected components", mlt_sol->ncomp);
		tsp_debug(inst->verbose >= 100, 1, "---------------- End Try #%d ----------------\n", n_call_ben_add_sec);
		free(xstar);

	}
	handleCPXResult(inst->verbose > 1,CPXgetstat(env, lp), "CPXResult for ben_reduce_comp:");
}
void ben_solve(instance* inst) {
	// open CPLEX model
	int error;
	char log_name[64];
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP_Problem");
	multitour_sol curr_sol;
	// Cplex's parameter setting
	set_init_param(env, (const instance*)inst, log_name, sizeof(log_name));
	tsp_debug(inst->verbose >= 200, 1, "set_init_param SUCCESSFUL");
	build_model(inst, env, lp);
	tsp_debug(inst->verbose >= 100, 1, "build_model SUCCESSFUL");
	if (CPXmipopt(env, lp)) { exit(main_error_text(-9, "CPXmipopt() error")); }
	handleCPXResult(inst->verbose > 1, CPXgetstat(env, lp), "CPXResult for initial CPXmipopt (no SEC):");
	// use the optimal solution found by CPLEX

	int ncols = CPXgetnumcols(env, lp);
	double* xstar = (double*)calloc(ncols, sizeof(double));
	{
		{
			if (CPXgetx(env, lp, xstar, 0, ncols - 1)) { exit(main_error_text(-9, "CPXgetx() error")); }
		}
	}
	print_selected_arcs(inst->verbose >= 200, (const double*)xstar, (const instance*)inst);
	tsp_debug(inst->verbose >= 200, 1, "print_seleceted_arcs SUCCESSFUL");
	init_multitour_sol(&curr_sol, inst->nnodes);
	build_sol((const double*)xstar, (const instance*)inst, curr_sol.succ, curr_sol.comp, &curr_sol.ncomp);
	if (CPXgetobjval(env, lp, &(curr_sol.z))) { exit(main_error_text(-9, "CPXgetobjval() error")); }
	tsp_debug(inst->verbose >= 100, 1, "build_sol SUCCESSFUL: there are %d connected components", curr_sol.ncomp);
	char figure_name[64];
	generate_name(figure_name, sizeof(figure_name), "figures/ben_%d_%d_multitour.png", inst->nnodes, inst->randomseed);
	plot_multitour(inst->verbose >= 1, inst->verbose >= 200, figure_name, (const multitour_sol*)&curr_sol, inst->nodes);
	ben_reduce_comp(env,lp, inst, &curr_sol);
	handleCPXResult(inst->verbose > 1, CPXgetstat(env, lp), "Final CPXResult:");
	if (cpx_update_best(inst->verbose >= 1, inst, env, lp, &curr_sol)) {
		generate_name(figure_name, sizeof(figure_name), "figures/ben_%d_%d_path.png", inst->nnodes, inst->randomseed);
		plot_path(inst->verbose >= 1, figure_name, inst->best_sol, inst->nodes, inst->nnodes);
		//plot_multitour(inst->verbose >= 1, inst->verbose >= 200, figure_name, (const multitour_sol*)&curr_sol, inst->nodes);
	}
	free_instance(inst);
	free_multitour_sol(&curr_sol);
	// free and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
	tsp_debug(inst->verbose >= 100, 1, "ben_solve SUCCESSFUL");
	return 0; // or an appropriate nonzero error code
}