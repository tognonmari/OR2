#include "hardfixing.h"
#include "plot.h"

void hf_init_table(hardfixing* hf) {
	if (!hf->table_flag) { return; }
	char* table_fields[] = { "#Nr. Call", "Best Cost", "New best", "p_fix" };
	int num_cols_table = sizeof(table_fields) / sizeof(table_fields[0]);
	make_first_row(hf->table_flag, hf->table_out, num_cols_table, "HARDFIXING");
	make_table_row(hf->table_flag, hf->table_out, num_cols_table, table_fields);
}

// Fill table with "#Nr Call", "Best Cost", "New best"
void hf_fill_table(hardfixing* hf, int nr_call, double cost, char is_new_best) {
	char table_0[32];
	char table_1[32];
	char table_2[32];
	char table_3[32];
	sprintf(table_0, "%d", nr_call);
	sprintf(table_1, "%.3f", cost);
	if (is_new_best) { sprintf(table_2, "Best"); }
	else { sprintf(table_2, " "); }
	if (hf->nr_call == 0) { sprintf(table_3, "HEUR"); }
	else { sprintf(table_3, "%.2f", hf->p_fix); }
	char* table_fields[] = { table_0, table_1, table_2, table_3 };
	int num_cols_table = sizeof(table_fields) / sizeof(table_fields[0]);
	make_table_row(hf->table_flag, hf->table_out, num_cols_table, table_fields);
}

void hf_init(hardfixing* hf, instance* inst) {
	hf->check_feasibility = 1;
	hf->plots_on_screen = 0; //it slower a lot the program
	hf->pipe = _popen("gnuplot -persist", "w");
	hf->table_flag = (inst->verbose >= 1);
	hf->tl_mipcall = inst->timelimit / 5.0;

	hf->apply_opt2 = 1;
	hf->p_fix = 0.7;
	hf->p_fix_scaling = 0.03;
	// decommenta se devi fare perf plot
	//hf->apply_opt2 = inst->hf_opt2;
	hf->p_fix = inst->hf_pfix_start;
	hf->p_fix_scaling = inst->hf_pfix_scaling;

	hf->policy = RND;
	hf->nr_call = 0;
	char table_file_name[64];
	generate_name(table_file_name, sizeof(table_file_name), "figures/hf_table_%d_%d.txt", inst->nnodes, inst->randomseed);
	hf->table_out = fopen(table_file_name, "w");
	hf_init_table(hf);
	tsp_debug(inst->verbose >= 300, 1, "hf_init SUCCESSFUL");
}

void hf_close(hardfixing* hf, instance* inst) {
	fclose(hf->pipe);
	make_last_row(hf->table_flag, hf->table_out, 4);
	fclose(hf->table_out);
	tsp_debug(inst->verbose >= 300, 1, "hf_close SUCCESSFUL");
}

void hf_fixing_rnd(int* xheu_path, instance* inst, double p_fix, CPXENVptr env, CPXLPptr lp) {
	char bd = 'L';
	double lb = 1.0;
	for (int i = 0; i < inst->nnodes - 1; i++) {
		double random_01 = rand_01();
		tsp_debug(inst->verbose >= 300, 0, "rand = %.2f, p_fix = %.2f", random_01, p_fix);
		if ((random_01 < p_fix)) {
			int index = xpos(xheu_path[i], xheu_path[i + 1], inst);
			if (CPXchgbds(env, lp, 1, &index, &bd, &lb)) { exit(main_error_text(-9, "CPXchgbds() error.\n")); }
			tsp_debug(inst->verbose >= 300, 0, "fixed edge (%d,%d)", xheu_path[i], xheu_path[i + 1]);
		}
		else {
			tsp_debug(inst->verbose >= 300, 0, "edge (%d,%d) has been left unfixed", xheu_path[i], xheu_path[i + 1]);
		}
	}
	tsp_debug(inst->verbose >= 300, 1, "hf_fixing_rnd SUCCESSFUL");
}

void hf_fixing(int* xheu_path, hardfixing* hf, instance* inst, CPXENVptr env, CPXLPptr lp) {
	switch (hf->policy) {
	case RND:
		hf_fixing_rnd(xheu_path, inst, hf->p_fix, env, lp);
		return;
	default:
		exit(main_error(-11));
		return;
	}
	return;
	tsp_debug(inst->verbose >= 300, 1, "hf_fixing SUCCESSFUL");
}

void hf_unfix(instance* inst, CPXENVptr env, CPXLPptr lp) {
	/*
	int cnt = inst->ncols;
	char* bds = (char*)malloc(cnt * sizeof(char));
	double* lb = (double*)malloc(cnt * sizeof(double));
	int* index = (int*)malloc(cnt * sizeof(int));
	char* ptr_bds = bds;
	double* ptr_lb = lb;
	int* ptr_index = index;
	for (int i = 0; i < cnt; i++) {
		(*ptr_bds) = 'L';
		(*ptr_lb) = 0.0;
		(*ptr_index) = 1;
		ptr_lb++;
		ptr_bds++;
		ptr_index++;
	}
	if (CPXchgbds(env, lp, cnt, index, bds, lb)) { exit(main_error_text(-9, "CPXchgbds() error.\n")); }
	free(index);
	free(lb);
	free(bds);
	*/
	char bd = 'L';
	double lb = 0.0;
	for (int i = 0; i < inst->ncols; i++) {
		int pos = i;
		if (CPXchgbds(env, lp, 1, &pos, &bd, &lb)) { exit(main_error_text(-9, "CPXchgbds() error.\n")); }
	}
	tsp_debug(inst->verbose >= 10, 1, "hf_unfix SUCCESSFUL");

}

//It does CPXmipopt() with the callback installed, this is the function to solve mip in the fastest known way.
void hf_fast_mip_solve(hardfixing* hf, instance* inst, CPXENVptr env, CPXLPptr lp) {
	char improvement = 0;
	int ncols = inst->ncols;
	inst->posting = 0;
	if (ncols == 0) { printf("ncols not initialized!"); exit(main_error(-7)); }
	multitour_sol curr_sol;
	//set the time to give to the mip_solver
	if ((get_timer() + hf->tl_mipcall) > inst->timelimit) { hf->tl_mipcall = inst->timelimit - get_timer(); }
	CPXsetdblparam(env, CPX_PARAM_TILIM, hf->tl_mipcall);
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
	hf->status_mipcall = CPXgetstat(env, lp);
	double old_zbest = inst->zbest;
	handleCPXResult(inst->verbose > 1, hf->status_mipcall, "Final CPXResult:");
	if (cpx_update_best(inst->verbose >= 1, inst, env, lp, &curr_sol)) {
		plot_path_on_screen(hf->pipe, hf->plots_on_screen, inst->best_sol, inst->nnodes, inst->zbest, inst->nodes);
		if (hf->apply_opt2) { opt2(inst, inst->best_sol, &(inst->zbest), 0); }
	}
	improvement = inst->zbest < old_zbest;
	hf_fill_table(hf, hf->nr_call, inst->zbest, improvement);
	(hf->nr_call)++;
	free(xstar);
	free_multitour_sol(&curr_sol);
}

void hf_solve(instance* inst) {
	hardfixing hf;
	// open CPLEX model
	int error;
	char log_name[64];
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP_Problem");
	build_model(inst, env, lp);
	// Cplex's parameter setting
	set_init_param(env, (const instance*)inst, log_name, sizeof(log_name));
	inst->ncols = CPXgetnumcols(env, lp);
	hf_init(&hf, inst); //TODO
	//1: initialize with an euristich method xheu
	gre_partial_solve(inst, 1, inst->nnodes / 100);
	int* xheu_path = inst->best_sol;
	hf_fill_table(&hf, hf.nr_call, inst->zbest, 1);
	(hf.nr_call)++;
	//2:fixing until time limit is not exceeded
	while (!is_time_limit_exceeded(inst->timelimit)) {
		//2.1: add mipstart
		cpx_add_mipstart(xheu_path, inst, env, lp);
		//2.2: fix variables to 1 (hardfixing)
		hf_fixing(xheu_path, &hf, inst, env, lp);
		//2.3: pass the restricted MIP model to the MIP solver
		hf_fast_mip_solve(&hf, inst, env, lp);
		//2.4: unfix variables
		hf_unfix(inst, env, lp);
		//2.5: check the actual best sol is feasible
		check_sol_is_feasible(hf.check_feasibility, inst, inst->best_sol, inst->zbest);
		//2.6: update p_fix
		int status = hf.status_mipcall;
		if ((status == CPXMIP_OPTIMAL_TOL) || (status == CPXMIP_OPTIMAL)) {
			hf.p_fix -= hf.p_fix_scaling;
			if (hf.p_fix < 0) { break; }
		}
		else { hf.tl_mipcall *= 2; }
	}
	hf_close(&hf, inst);
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
	tsp_debug(inst->verbose >= 300, 1, "hf_solve SUCCESSFUL");
}