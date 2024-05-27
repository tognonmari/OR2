#include "softfixing.h"
#include "plot.h"

void sf_init_table(softfixing* sf) {
	if (!sf->table_flag) { return; }
	char* table_fields[] = { "#Nr. Call", "Best Cost", "New best", "k" };
	int num_cols_table = sizeof(table_fields) / sizeof(table_fields[0]);
	make_first_row(sf->table_flag, sf->table_out, num_cols_table, "SOFTFIXING");
	make_table_row(sf->table_flag, sf->table_out, num_cols_table, table_fields);
}
void sf_init(softfixing* sf, instance* inst) {
	sf->check_feasibility = 1;
	sf->plots_on_screen = 1;
	sf->pipe = _popen("gnuplot -persist", "w");
	sf->table_flag = (inst->verbose >= 1);
	sf->tl_mipcall = inst->timelimit/5; //da studiare

	// decommenta se devi fare perf plot
	/*
	sf->k = inst->sf_k_start;
	sf->k_scaling = inst->sf_k_scaling;*/
	// usa queste se non devi fare perf plot
	sf->apply_opt2 = 0;
	sf->k = 20;
	sf->k_scaling = 5;
	sf->nr_call = 0;
	char table_file_name[64];
	generate_name(table_file_name, sizeof(table_file_name), "figures/sf_table_%d_%d.txt", inst->nnodes, inst->randomseed);
	sf->table_out = fopen(table_file_name, "w");
	sf_init_table(sf);
	tsp_debug(inst->verbose >= 300, 1, "sf_init SUCCESSFUL");
}

// Fill table with "#Nr Call", "Best Cost", "New best"
void sf_fill_table(softfixing* sf, int nr_call, double cost, char is_new_best) {
	char table_0[32];
	char table_1[32];
	char table_2[32];
	char table_3[32];
	sprintf(table_0, "%d", nr_call);
	sprintf(table_1, "%.3f", cost);
	if (is_new_best) { sprintf(table_2, "Best"); }
	else { sprintf(table_2, " "); }
	if (sf->nr_call == 0) { sprintf(table_3, "HEUR"); }
	else { sprintf(table_3, "%d", sf->k); }
	char* table_fields[] = { table_0, table_1, table_2, table_3 };
	int num_cols_table = sizeof(table_fields) / sizeof(table_fields[0]);
	make_table_row(sf->table_flag, sf->table_out, num_cols_table, table_fields);
}

void sf_close(softfixing* sf, instance* inst) {
	fclose(sf->pipe);
	make_last_row(sf->table_flag, sf->table_out, 4);
	fclose(sf->table_out);
	tsp_debug(inst->verbose >= 300, 1, "sf_close SUCCESSFUL");
}

void sf_fast_mip_solve(softfixing* sf, instance* inst, CPXENVptr env, CPXLPptr lp) {
	sf->improvement = 0;
	int ncols = inst->ncols;
	inst->posting = 0;
	if (ncols == 0) { printf("ncols not initialized!"); exit(main_error(-7)); }
	multitour_sol curr_sol;
	//set the time to give to the mip_solver
	if ((get_timer() + sf->tl_mipcall) > inst->timelimit) { sf->tl_mipcall = inst->timelimit - get_timer(); }
	CPXsetdblparam(env, CPX_PARAM_TILIM, sf->tl_mipcall);
	//install callback
	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION;
	if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst)) { exit(main_error_text(-9, "CPXcallbacksetfunc() error")); }
	int status = CPXmipopt(env, lp);
	if (status) {
		tsp_debug(1, 1, "The status is %d.\n", status);
		exit(main_error_text(-9, "CPXmipopt() error "));
	}
	tsp_debug(inst->verbose >= 300, 1, "CPXmipopt SUCCESSFUL");
	double* xstar = (double*)calloc(ncols, sizeof(double));
	{
		{
			if (CPXgetx(env, lp, xstar, 0, ncols - 1)) { exit(main_error_text(-9, "CPXgetx() error")); }
			tsp_debug(inst->verbose >= 300, 1, "CPXgetx SUCCESSFUL\n");
		}
	}

	init_multitour_sol(&curr_sol, inst->nnodes);
	build_sol((const double*)xstar, (const instance*)inst, curr_sol.succ, curr_sol.comp, &curr_sol.ncomp);
	if (CPXgetobjval(env, lp, &(curr_sol.z))) { exit(main_error_text(-9, "CPXgetobjval() error")); }
	sf->status_mipcall = CPXgetstat(env, lp);
	double old_zbest = inst->zbest;
	handleCPXResult(inst->verbose > 1, sf->status_mipcall, "Final CPXResult:");
	if (cpx_update_best(inst->verbose >= 1, inst, env, lp, &curr_sol)) {
		plot_path_on_screen(sf->pipe, sf->plots_on_screen, inst->best_sol, inst->nnodes, inst->zbest, inst->nodes);
		if (sf->apply_opt2) { opt2(inst, inst->best_sol, &(inst->zbest), 0); }
	}
	sf->improvement = inst->zbest < old_zbest;
	sf_fill_table(sf, sf->nr_call, inst->zbest, sf->improvement);
	(sf->nr_call)++;
	free(xstar);
	free_multitour_sol(&curr_sol);
}

void sf_remove_constraint(instance* inst,CPXENVptr env,CPXLPptr lp) {
	int num_rows = CPXgetnumrows(env, lp);
	if(CPXdelrows(env, lp, num_rows - 1, num_rows -1) ) { exit(main_error_text(-9, "Wrong CPXdelrows[softfixing]")); }
	tsp_debug(inst->verbose >= 300, 1, "sf_remove_constraint SUCCESSFUL");
}
void sf_add_constraint(int* xheu_path, softfixing* sf, instance* inst, CPXENVptr env, CPXLPptr lp) {
	int ncols = inst->ncols;
	int n = inst->nnodes;
	int k = sf->k;
	int* index = (int*)malloc(ncols * sizeof(int));
	double* value = (double*)malloc(ncols * sizeof(double));
	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));
	// add the k constraints
	double rhs = (double)(n - k);
	char sense = 'G';                     // 'G' for greater equal constraint 
	int izero = 0;
	sprintf(cname[0], "softfix(%d)", sf->nr_call);
	int nnz = 0;
	for (int i = 0; i < n - 1; i++)
	{
		index[nnz] = xpos(xheu_path[i], xheu_path[i+1], inst);
		value[nnz] = 1.0;
		nnz++;
	}
	index[nnz] = xpos(xheu_path[n-1], xheu_path[0], inst);
	value[nnz] = 1.0;
	nnz++;
	tsp_debug(inst->verbose >= 300, 1, "Num Rows constraint before addrows = %d",CPXgetnumrows(env,lp));
	if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0])) { exit(main_error_text(-9, "Wrong CPXaddrows[softfixing]"));}
	tsp_debug(inst->verbose >= 300, 1, "Num Rows constraint after addrows = %d", CPXgetnumrows(env, lp));
	free(value);
	free(index);
	free(cname[0]);
	free(cname);
	tsp_debug(inst->verbose >= 300, 1, "sf_add_constraint SUCCESSFUL");
}

void sf_solve(instance* inst) {
	softfixing sf;
	// open CPLEX model
	int error;
	char log_name[64];
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP_Problem");
	build_model(inst, env, lp);
	// Cplex's parameter setting
	set_init_param(env, (const instance*)inst, log_name, sizeof(log_name));
	inst->ncols = CPXgetnumcols(env, lp);
	sf_init(&sf, inst); //TODO
	//1: initialize with an euristich method xheu
	gre_partial_solve(inst, 1, inst->nnodes / 100);
	/* better but require more time
	inst->timelimit *= 0.1;
	tabu_search(1, inst);
	inst->timelimit *= 9;
	*/
	int* xheu_path = inst->best_sol;
	sf_fill_table(&sf, sf.nr_call, inst->zbest, 1);
	(sf.nr_call)++;
	//2:fixing until time limit is not exceeded
	while (!is_time_limit_exceeded(inst->timelimit)) {
		//2.1: add mipstart
		cpx_add_mipstart(xheu_path, inst, env, lp);
		//2.2: add softfixing constraint
		sf_add_constraint(xheu_path, &sf, inst, env, lp);
		//2.3: pass the restricted MIP model to the MIP solver
		sf_fast_mip_solve(&sf, inst, env, lp);
		//2.4: remove softfixing constraint
		sf_remove_constraint(inst, env, lp);
		//2.5: check the actual best sol is feasible
		check_sol_is_feasible(sf.check_feasibility, inst, inst->best_sol, inst->zbest);
		//2.6: update k
		sf.k += sf.k_scaling;
		/*
		if (!sf.improvement) {
			sf.k += sf.k_scaling;
		}
		*/

	}
	sf_close(&sf, inst);
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
	tsp_debug(inst->verbose >= 300, 1, "sf_solve SUCCESSFUL");
}