#include "cpx.h"
#include "benders.h"
#include "plot.h"
#include "mincut.h"
//Function that handles CPX function results
void handleCPXResult(int flag, int result, char* format) {
	tsp_debug(flag,1, "%s", format);
	tsp_debug_inline(flag, " CPXstatus = %d", result);
	switch (result) {
	case CPX_STAT_OPTIMAL:
		tsp_debug_inline(flag, " The optimal solution has been found");
		break;
	case CPX_STAT_INFEASIBLE:
		tsp_debug_inline(flag, " Infeasible problem");
		break;
	case CPX_STAT_UNBOUNDED:
		tsp_debug_inline(flag, " Unbounded problem");
		break;
	case CPX_STAT_ABORT_TIME_LIM:
		tsp_debug(flag,1,"The resolution of the problem is ended due by exceeding internal CPX_PARAM_TILIM.\n");
		break;
	case CPXMIP_OPTIMAL:
		tsp_debug_inline(flag, " An optimal solution has been found, without tolerance");
		break;
	case CPXMIP_TIME_LIM_FEAS:
		tsp_debug_inline(flag, " CPXmipopt stopped because timelimit has been exceeded.");
		tsp_debug(flag,0, "An integer feasible solution is been found for the problem.");
		break;
	case CPXMIP_OPTIMAL_TOL:
		tsp_debug_inline(flag, " An optimal solution has been found, with tolerance.");
		break;
	case CPXMIP_TIME_LIM_INFEAS:
		tsp_debug_inline(flag, " CPXmipopt stopped because timelimit has been exceeded.");
		tsp_debug(flag, 0, "No feasible solution has been found for the problem.");
		break;
	default:
		printf("\n\nWARNING it has been returned result: %d, which is not handled. Please add the case %d in exact.c -> handleCPXResult \n\n",result, result);
	}
}


static int CPXPUBLIC my_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle){
	instance* inst = (instance*)userhandle;
	// select what to do according to the context
	if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE) {
		return my_callback_candidate(context, inst);
	}

	if (contextid == CPX_CALLBACKCONTEXT_RELAXATION) {
		return my_callback_relaxation(context,contextid, inst);
	}

	exit(main_error_text(-9, "Called the callback with the wrong contextid.\n"));
	
}


static int CPXPUBLIC my_callback_relaxation(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle) {
	
	int mynode = -1; 
	
	CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode);
	if(mynode!=0){
	//if (mynode % 10) {
		return 0;
	}
	//Step 0: initializations and allocations
	instance* inst = (instance*)userhandle;
	double* xstar = calloc(inst->ncols, sizeof(double));
	int* elist = calloc(2 * inst->ncols, sizeof(int));
	int ncomp=0;
	int* compscount = NULL;
	int* comps = NULL;
	double relaxation_cost;
	
	//Step 1: get the solution to the relaxation from the context

	if (CPXcallbackgetrelaxationpoint(context, xstar, 0, inst->ncols-1, &relaxation_cost)) {
	
		exit(main_error_text(-9, "Failed at retrieving the relaxation.\n"));

	}
	//Step 2: turn it into a concorde-like solution
	int num_edges = 0;
	//int t = 0;
	int* pointer_to_elist = elist;
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = i + 1; j < inst->nnodes; j++) {
			
			(*pointer_to_elist) = i;
			pointer_to_elist++;
			(*pointer_to_elist) = j;
			pointer_to_elist++;
			
			/*
			elist[t++] = i;
			elist[t++] = j;
			*/
			num_edges++;
		}

	}
	//Step 3: get the number of connected components

	if (CCcut_connect_components(inst->nnodes, num_edges, elist, xstar, &ncomp, &compscount, &comps)) {

		exit(main_error_text(-11, "Error in CCcut_connect_components.\n"));

	}

	tsp_debug(inst->verbose >= 1, 1, "There are %d components at node %d.\n",ncomp, mynode);
	//Step 4: add the cuts according to the number of components

	if (ncomp == 1) {

		violated_cuts_params params = { .inst = inst, .context = context };
		//concorde method to find the cuts, given that the graph is connected: the method will itself use a callback -> a new structure of parameters must be used
		if (CCcut_violated_cuts(inst->nnodes, inst->ncols, elist, xstar, 1.9,cpx_add_cut_single_comp, &params )) {

			exit(main_error_text(-11, "Error in CCcut_violated_cuts().\n"));

		}
		

	}
	else {
		//manually add the cuts from compscount e comps

		int* my_comp = (int*) calloc(inst->nnodes, sizeof(int));
		int start_component = 0; 

		//for each component of the solution, update nodes'component in my_comp

		for (int i = 0; i < ncomp; i++) {
			
			for (int j = start_component; j < start_component + compscount[i]; j++) {
				
				my_comp[comps[j]] = i + 1;
			}
			
			start_component += compscount[i];
		}
		//for each connected components we add the violated SEC

		for (int i = 0; i < ncomp; i++) {


			cpx_add_violated_SECs_fractional(context, inst, my_comp, i+1);

		}

		free(my_comp);
	}
	
	free(xstar);
	free(elist);
	return 0;
}

static int cpx_add_cut_single_comp(double cutval, int num_nodes_in_cut, int* nodes_in_cut, void* userhandle) {

	// Step 0: cast the userhandle
	violated_cuts_params* params = (violated_cuts_params*)userhandle;
	instance* inst = params->inst;
	CPXCALLBACKCONTEXTptr context = params->context;

	int edges_in_cut =(num_nodes_in_cut * (num_nodes_in_cut - 1)/2);
	int* index = (int*)malloc(edges_in_cut * sizeof(int));
	double* value = (double*)malloc(edges_in_cut * sizeof(double));
	int matbeg = 0;
	int purgeable = CPX_USECUT_FILTER;
	int local = 0;
	int nnz = 0;
	char sense = 'L';
	double rhs = num_nodes_in_cut-1;

	//scorro nodes_in_cut e metto in index[nnz] and 

	for (int i = 0; i < num_nodes_in_cut; i++) {
		
		for (int j = 0; j < num_nodes_in_cut; j++) {
			
			if (nodes_in_cut[i] >= nodes_in_cut[j])
				continue;
			
			index[nnz] = xpos(nodes_in_cut[i], nodes_in_cut[j], inst);
			value[nnz] = 1.0;
			nnz++;
		}
	}

	if (CPXcallbackaddusercuts(context, 1, edges_in_cut, &rhs, &sense, &matbeg, index, value, &purgeable, &local)) {
		exit(main_error_text(-11, "Could not add fractional cut for 1 connected component. Error within the callback of concorde.\n"));
	}

	free(index);
	free(value);

	return 0;
}

static int cpx_add_violated_SECs_fractional(CPXCALLBACKCONTEXTptr context, instance* inst, int* my_comp, int connected_component_id ) {

	int* index = (int*)malloc(inst->ncols * sizeof(int));
	double* value = (double*)malloc(inst->ncols * sizeof(double));
	int matbeg = 0;
	int purgeable = CPX_USECUT_FILTER; //cosa cambia tra filter e purge? non mi ricordo
	int local = 0;
	int nnz = 0;
	char sense = 'L';
	double rhs;
	int num_nodes_in_cut = 0;
	for (int i = 0; i < inst->nnodes; i++) {
		if (my_comp[i] != connected_component_id) {
			continue;
		}
		num_nodes_in_cut++;
		for (int j = i + 1; j < inst->nnodes; j++) {
			if (my_comp[j] != connected_component_id) {
				continue;
			}
			index[nnz] = xpos(i, j, inst);
			value[nnz] = 1.0;
			nnz++;
		}
	}
	rhs = (double)(num_nodes_in_cut - 1);
	if (CPXcallbackaddusercuts(context, 1, nnz, &rhs, &sense, &matbeg, index, value, &purgeable, &local)) {
		exit(main_error_text(-11, "Could not add fractional cut for more than 1 connected component.\n"));
	}

	free(index);
	free(value);

}

static int CPXPUBLIC my_callback_candidate(CPXCALLBACKCONTEXTptr context, instance* inst) {
	//Add cuts from an integer multitour solution
	multitour_sol mlt;
	double* xstar = (double*)malloc(inst->ncols * sizeof(double));
	double objval = CPX_INFBOUND;
	// 1: ottieni la soluzione candidata
	if (CPXcallbackgetcandidatepoint(context, xstar, 0, inst->ncols - 1, &objval)) { exit(main_error_text(-9, "CPXcallbackgetcandidatepoint error")); }
	// 2: ottieni informazioni sui nodi
	int mythread = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread);
	int mynode = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode);
	double incumbent = CPX_INFBOUND; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent);
	tsp_debug(inst->verbose >= 10,1,"callback at node %5d thread %2d incumbent %.12g, candidate value %10.2lf\n", mynode, mythread, incumbent, objval);
	//3: inizializzazione soluzione multitour
	init_multitour_sol(&mlt, inst->nnodes);
	//4: costruzione della soluzione multitour form input solution $xstar is built a solution for the tsp problem
	build_sol((const double*)xstar, (const instance*)inst, mlt.succ, mlt.comp, &mlt.ncomp);
	//4b: acquisizione del costo della soluzione multitour
	
	//5: se ci sono più componente connesse creo vincolo di taglio da aggiungere al modello e rifiuto il candidato
	if (mlt.ncomp > 1) {
		int izero = 0;
		int* index = (int*)malloc(inst->ncols * sizeof(int));
		double* value = (double*)malloc(inst->ncols * sizeof(double));
		for (int k = 1; k <= mlt.ncomp; k++) {
			int nnz = 0;
			char sense = 'L';
			double rhs = -1.0;
			for (int i = 0; i < inst->nnodes; i++) {
				if (mlt.comp[i] != k) {
					continue;
				}
				rhs += 1.0;
				for (int j = i + 1; j < inst->nnodes; j++) {
					if (mlt.comp[j] != k) {
						continue;
					}
					index[nnz] = xpos(i, j, inst);
					value[nnz] = 1.0;
					nnz++;
				}
			}
			if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, value)) { exit(main_error_text(-9, "CPXcallbackrejectcandidate() error")); }// reject the solution and adds one cut
		}
		free(index);
		free(value);
		//6: esegui gluing soluzione e postala
		if (inst->posting) {
			multitour_sol patched_sol;
			init_multitour_sol(&patched_sol, inst->nnodes);
			ben_patching(&mlt, &patched_sol, inst);
			int* path = (int*)calloc(inst->nnodes, sizeof(int));
			cpx_convert_succ_in_path(&patched_sol, path, inst->nnodes);
			double zheu = compute_path_length(path, inst->nnodes, inst->nodes);
			opt2(inst, path, &zheu, 0);
			double* cplex_format_xheu = (double*)calloc(inst->ncols, sizeof(double));
			int* ind = (int*)malloc(inst->ncols * sizeof(int));
			cpx_convert_path_to_cplex(path, cplex_format_xheu, inst);
			for (int j = 0; j < inst->ncols; j++) ind[j] = j;
			CPXcallbackpostheursoln(context, inst->ncols, ind, cplex_format_xheu, zheu, CPXCALLBACKSOLUTION_CHECKFEAS);
			tsp_debug((inst->verbose >= 10), 1, "Posting Heuristic (Gluing) with Cost = %.3f", zheu);
			free_multitour_sol(&patched_sol);
			free(path);
			free(ind);
			free(cplex_format_xheu);
		}
	}
	free(xstar);
	return 0;
}

void cpx_branch_and_cut(char relaxation, char mipstart, char posting, instance* inst) {
	inst->posting = posting;
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
	// install a "lazyconstraint" callback to cut infeasible integer sol.s (found e.g. by heuristics) 
	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE ; 
	if (relaxation) {
		contextid |=  CPX_CALLBACKCONTEXT_RELAXATION;
	}
	int ncols = CPXgetnumcols(env, lp);

	inst->ncols = ncols;
	if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst)) { exit(main_error_text(-9, "CPXcallbacksetfunc() error") ); }
	if (mipstart) {
		//Generate a greedy solution and initialize the incumbent with it
		gre_partial_solve(inst, 1, inst->nnodes/30);
		double* cplex_format_xheu = (double*)calloc(inst->ncols, sizeof(double));
		int* ind = (int*)malloc(inst->ncols * sizeof(int));
		cpx_convert_path_to_cplex(inst->best_sol, cplex_format_xheu, inst);
		for (int j = 0; j < inst->ncols; j++) ind[j] = j;
		int effortlevel = CPX_MIPSTART_CHECKFEAS;
		int beg = 0;
		if (CPXaddmipstarts(env, lp, 1, inst->ncols, &beg, ind, cplex_format_xheu, &effortlevel, NULL)) {
			exit(main_error_text(-9, "Could not add MipStart.\n"));
		}

		tsp_debug(1, 1, "Updated warm start successfully.\n");
		free(ind);
		free(cplex_format_xheu);

	}
	CPXsetintparam(env, CPX_PARAM_THREADS, 1); 	// just for debugging
	int status = CPXmipopt(env, lp);
	if (status) { 
		tsp_debug(1, 1, "The status is %d.\n",status);
		exit(main_error_text(-9, "CPXmipopt() error ")); }
	tsp_debug(inst->verbose >= 100, 1, "CPXmipopt SUCCESSFUL");
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
	if (CPXgetobjval(env, lp, &(curr_sol.z) ) ) { exit(main_error_text(-9, "CPXgetobjval() error")); }
	tsp_debug(inst->verbose >= 100, 1, "build_sol SUCCESSFUL: there are %d connected components", curr_sol.ncomp);
	handleCPXResult(inst->verbose > 1, CPXgetstat(env, lp), "Final CPXResult:");

	if (cpx_update_best(inst->verbose >= 1, inst, env, lp, &curr_sol)) {
		char figure_name[64];
		generate_name(figure_name, sizeof(figure_name), "figures/branch_%d_%d_path.png", inst->nnodes, inst->randomseed);
		plot_path(inst->verbose >= 1, inst->best_sol, inst->nnodes, inst->zbest, inst->nodes, figure_name);
	}


	free(xstar);
	free_multitour_sol(&curr_sol);


	// free and close cplex model 
	  
	CPXfreeprob(env, &lp);

	CPXcloseCPLEX(&env);
	tsp_debug(inst->verbose >= 100, 1, "cpx_branch_and_cut SUCCESSFUL");
	return 0; // or an appropriate nonzero error code
}
void cpx_convert_succ_in_path(const multitour_sol* mlt, int* path, int n) {
	if (mlt->ncomp > 1) { 
		exit(main_error_text(-10, "Cannot convert multitour in a full connected path\nYou must provide a single tour!")); }
	int curr = 0;
	path[0] = curr;
	for (int i = 1; i < n; i++) {
		curr = (mlt->succ)[curr];
		path[i] = curr;
	}
}
void cpx_convert_path_in_succ(const int* path, multitour_sol* mlt, int n) {
	if (mlt->ncomp > 1) { exit(main_error_text(-10, "Cannot convert path into multitour.\nYou must provide a single tour!")); }
	for (int t = 0; t < n-1; t++) {
		// t moves along path array 
		mlt->succ[path[t]]=path[t+1];
	}
	mlt->succ[path[n - 1]] = path[0];
}
void cpx_convert_path_to_cplex(const int* xheu_path, double* cplex_like_format, instance* inst) {

	int n = inst->nnodes;
	int* succ = (int*)calloc(n, sizeof(int));
	//Step 0: Convert path in succ
	for (int t = 0; t < n - 1; t++) {
		succ[xheu_path[t]] = xheu_path[t + 1];
	}
	succ[xheu_path[n - 1]] = xheu_path[0];

	//Step 1: convert succ in cplex 
	cpx_convert_succ_to_cplex(succ, cplex_like_format, inst);
	free(succ);

}

void cpx_convert_succ_to_cplex(const int* succ, double* cplex_like_format, instance* inst) {

	int n = inst->nnodes;
	
	//Step 1: convert succ in cplex 
	for (int i = 0; i < n; i++) {
		cplex_like_format[xpos(i, succ[i], inst)] = 1.0;
	}
}
int cpx_update_best(char flag, instance* inst, CPXENVptr env, CPXLPptr lp, const multitour_sol* sol) {
	double z;
	double lb;
	int result = CPXgetstat(env, lp);
	int* path;
	char figure_name[128];
	int n = inst->nnodes;
	switch (result) {
	case CPXMIP_OPTIMAL: case CPXMIP_OPTIMAL_TOL:
		tsp_debug(flag, 1, " Updating solution for the problem with an optimal one");
		CPXgetbestobjval(env, lp, &z);
		CPXgetobjval(env, lp, &lb);
		path = (int*)calloc(n, sizeof(int));
		cpx_convert_succ_in_path(sol, path, n);
		update_best(inst, z, get_timer(), path);
		tsp_debug(flag, 1, "LB = %.2lf ZBEST = %.2lf REALCOST = %.2lf", z, lb, compute_path_length(path, n, inst->nodes));
		break;
	case CPXMIP_TIME_LIM_FEAS:
		if (sol->ncomp == 1) {
			tsp_debug(flag, 1, "Updating solution for the problem with a feasible one");
			z = sol->z;
			path = (int*)calloc(n, sizeof(int));
			cpx_convert_succ_in_path(sol, path, n);
			if (inst->solver == GLU) {
				generate_name(figure_name, sizeof(figure_name), "figures/ben_%d_%d_%d_pre2opt.png", inst->nnodes, inst->randomseed);
				plot_path(inst->verbose >= 1, (const int*)path, inst->nnodes, z, inst->nodes, figure_name);
				//opt 2 optimization for patched solution 
				opt2(inst, path, &z, inst->verbose >=100);
				//From path to succ: from path, update succ of patched sol
				generate_name(figure_name, sizeof(figure_name), "figures/ben_%d_%d_%d_postpatch.png", inst->nnodes, inst->randomseed);
				plot_path(inst->verbose >= 1, (const int*)path, inst->nnodes, z, inst->nodes, figure_name);
			}
			update_best(inst, sol->z, get_timer(), path);
			break;
		}
		else {
			CPXgetbestobjval(env, lp, &z);
			tsp_debug(flag, 1, "Try to updating solution, but tsp require single tour solution");
			tsp_debug(flag, 0, "While the current solution has a multitour solution");
			tsp_debug(flag, 0, "Probably you need to increase the time limit for this instance");
			tsp_debug(flag, 0, "The current cost of the solution is: %.2f", z);
			return 0;
		}
	case CPXMIP_TIME_LIM_INFEAS:
		tsp_debug(flag, 1, "Try to updating solution, but no feasible solution has been found.");
		tsp_debug(flag, 0, "Probably you need to increase the time limit for this instance");
		return 0;
	default:
		printf("\n\nWARNING it has been returned result: %d, which is not handled. Please add the case %d in handleCPXResult and cpx_update_best \n\n",result, result);
		return 0;
	}
	return 1;
}

double dist(int i, int j, instance* inst) {
	float* dist_matr = inst->dist_matrix;
	if (!inst->integer_costs) return (double)(get_dist_matrix((const float*)(dist_matr), i, j));
	int dis = (double)(get_dist_matrix((const float*)(dist_matr), i, j)) + 0.499999999;
	return dis + 0.0;
}
/*
 * Calculates the index of the variable x(i,j) in the model.
 *
 * IP i First node index (0-based).
 * IP j Second node index (0-based).
 * IP inst Pointer to the instance data.
 * OR Index of the variable x(i,j) in the model.
 */
int xpos(int i, int j,const instance* inst)                                      
{
	if (i == j) exit(main_error_text(-9," i == j in xpos , error because x_ii is not a variable of the model"));
	if (i > j) return xpos(j, i, inst);
	int pos = i * inst->nnodes + j - ((i + 1) * (i + 2)) / 2;
	return pos;
}

/*
 * Builds the mathematical model for solving basic TSP problem (no subtour constraints).
 *
 * IP inst Pointer to the instance data (assumed to be defined elsewhere).
 * IOP env CPLEX environment pointer.
 * OP lp CPLEX linear programming pointer.
*/
void build_model(const instance* inst, CPXENVptr env, CPXLPptr lp)
{
	const int n = inst->nnodes;
	int izero = 0;
	char binary = 'B';

	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));
	// add binary var.s x(i,j) for i < j 
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);  // x(1,2), x(1,3) ....
			double obj = dist(i, j, inst); // cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) { printf(main_error_text(-9, "Wrong CPXnewcols on x var.s")); }
			//observe that the number of added cols - 1 is equal to the actual index in the model of x_ij
			if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, inst)) { printf(main_error_text(-9, "Wrong position for x var.s")); }
		}
	}

	// add degree constr.s 

	int* index = (int*)malloc(n * sizeof(int));
	double* value = (double*)malloc(n * sizeof(double));

	// add the degree constraints
	for (int h = 0; h < n; h++)  // degree constraints
	{
		double rhs = 2.0;
		char sense = 'E';                     // 'E' for equality constraint 
		sprintf(cname[0], "degree(%d)", h + 1);
		int nnz = 0;
		for (int i = 0; i < n; i++)
		{
			if (i == h) continue;
			index[nnz] = xpos(i, h, inst);
			value[nnz] = 1.0;
			nnz++;
		}

		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0])) { exit(main_error_text(-9, "Wrong CPXaddrows[degree]")); }
	}
	free(value);
	free(index);

	if (inst->verbose >= -100) CPXwriteprob(env, lp, "model.lp", NULL);
	tsp_debug(inst->verbose >= 100, 1, "CPXwriteprob SUCCESSFUL");
	free(cname[0]);
	free(cname);
}

void set_init_param(CPXENVptr env, const instance* inst, char* log_name, size_t size_log_name){
	// increased precision for big-M models
	CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);		// very important if big-M is present
	CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9);
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit - get_timer());
	CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4);
	if (inst->verbose >= 101) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	//CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);

	//CPXsetdblparam(env, CPX_PARAM_TILIM, CPX_INFBOUND + 0.0);
	MKDIR("logs");
	generate_name(log_name, size_log_name, "logs/log_%d_%d_%d.txt", inst->solver,  inst->nnodes, inst->randomseed);
	CPXsetlogfilename(env, (const char*) log_name, "w");
	//CPXsetintparam(env, CPX_PARAM_NODELIM, 0); 		// abort Cplex after the root node
	//CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);	// abort Cplex after the first incumbent update
	//CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-4);  	// abort Cplex when gap below 10%
}
void init_multitour_sol(multitour_sol* sol, int n) {
	sol->succ = (int*)calloc(n, sizeof(int));
	sol->ncomp = 0;
	sol->comp = (int*)calloc(n, sizeof(int));
	sol->z = -1;
}
void free_multitour_sol(multitour_sol* sol) {
	free(sol->comp);
	free(sol->succ);
}

void update_mlt_cost(multitour_sol* mlt, instance* inst) {
	mlt->z = compute_mlt_cost(mlt, inst);
}
double compute_mlt_cost(const multitour_sol* mlt, instance* inst) {
	int* succ = mlt->succ;
	int* comp = mlt->comp;
	int ncomp = mlt->ncomp;
	int n = inst->nnodes;
	float** dist_m = inst->dist_matrix;
	double cost = 0;
	for (int k = 1; k <= mlt->ncomp; k++) {
		for (int i = 0; i < n; i++) {
			//select a node in component k and set it to start
			if (comp[i] != k) { continue; }
			int start = i;
			int curr = start;
			int next = succ[start];
			//compute component cost
			while (next != start) {
				cost += get_dist_matrix(dist_m, curr, next);
				curr = next;
				next = succ[next];
			}
			cost += get_dist_matrix(dist_m, curr, next);
			break;
		}
	}
	return cost;
}
/*From input solution $xstar is built a solution for tsp problem, the solution is structured as follows:
* OP succ It contains the sequence of indices of the nodes involved in $xstar, observe that initially subtour constraints are not in the model,
*		  hence we can have several tours.
* OP ncomp Number of connected components.
* OP comp Array specifyng the connected components of the nodes, thus comp[i] indicate the index of the
*		  connected component of node with index i;
* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
* 
* WARNING! THIS FUNCTION DOES NOT UPDATE THE COST OF THE MULTITOUR SOLUTION 
* 
* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
* The function does the following:
* 1)Node Degree Calculation:
*		The array degree is initialized with zeros and represents the degree of each node.
*		It iterates over all pairs of nodes (i, j), and if the edge [i, j] is selected based on the values in xstar,
*		the degree of both nodes is incremented.
* 2)Node Degree Verification:
*		Subsequently, it checks that the degree of each node is exactly 2.
*		If not, an error is generated.
* 3)Connected Component Search:
*		The function examines nodes and creates connected components.
*		For each unvisited node, a new component is created.
*		A path is followed through the edges selected in xstar until the component is complete.
*/
void build_sol(const double* xstar, const instance* inst, int* succ, int* comp, int* ncomp) {
	const int n = inst->nnodes;
	int* degree = (int*)calloc(n, sizeof(int));
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			int k = xpos(i, j, inst);
			if (	( fabs(xstar[k]) > EPS ) && ( fabs(xstar[k] - 1.0) > EPS)		)	
				{ exit(main_error_text(-9, "Wrong xstar in build_sol()")); }
			if (xstar[k] > 0.5)
			{
				++degree[i];
				++degree[j];
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		if (degree[i] != 2){ exit(main_error_text(-9, "Wrong degree in build_sol()"));
	}
	}
	free(degree);
	tsp_debug(inst->verbose >= 200, 1, "checking degree constraints SUCCESSFUL");
	*ncomp = 0;
	for (int i = 0; i < n; i++)
	{
		succ[i] = -1;
		comp[i] = -1;
	}

	for (int start = 0; start < n; start++)
	{
		if (comp[start] >= 0) continue;  // node "start" was already visited, just skip it

		// a new component is found
		(*ncomp)++;
		int i = start;
		int done = 0;
		while (!done)  // go and visit the current component
		{
			comp[i] = *ncomp;
			done = 1;
			for (int j = 0; j < n; j++)
			{
				if (i != j && xstar[xpos(i, j, inst)] > 0.5 && comp[j] == -1) // the edge [i,j] is selected in xstar and j was not visited before 
				{
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}
		succ[i] = start;  // last arc to close the cycle
		tsp_debug(inst->verbose >= 200, 1, "build_sol for component:%d SUCCESSFUL",comp[start]);
		// go to the next component...
	}
}
/*Print all the edges selected for the solution (edges corresponding to a variable equals to 1 in the solution)
* IP flag Print flag 
* IP xstar Solution variables value
* IP inst Instance pointer
* OV if $flag print the selected edges in the solution, otherwise print nothing
*/
void print_selected_arcs(char flag, const double* xstar,const instance* inst) {

	int n = inst->nnodes;
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			if (xstar[xpos(i, j, inst)] > 0.5) tsp_debug_inline(flag,"  ... x(%3d,%3d) = 1\n", i + 1, j + 1);
		}
	}
}
/*
void plot_multitour(char figure_flag,char debug_flag, const char figure_name[], const multitour_sol* sol, const point* points){
	FILE* gnuplotPipe = _popen("gnuplot -persist", "w");

	if (gnuplotPipe == NULL) {
		fclose(gnuplotPipe);
		exit(main_error_text(-1, "Failed to open the pipeline to gnuplot"));
	}

	fprintf(gnuplotPipe, "set output '%s'\n", figure_name); //set output name
	fprintf(gnuplotPipe, "set terminal png\n"); //set extension
	fprintf(gnuplotPipe, "set title 'Solution without SEC'\n");
	fprintf(gnuplotPipe, "set key off\n"); //if key on the name of the file is printed on the plot
	fprintf(gnuplotPipe, "set pointsize 0.5\n");
	fprintf(gnuplotPipe, "set grid \n");
	fprintf(gnuplotPipe, "plot '-' with linespoints pt 1 lc rgb %s\n", generate_color(rand() % 101) );
	int comp_counter = 1;
	int i = 0; //index_node counter
	int ncomp = sol->ncomp;
	int* succ = sol->succ;
	int* comp = sol->comp;
	tsp_debug(debug_flag, 1, "START PLOTTING");
	while (comp_counter <= ncomp) {
		int start = 0;
		while (comp[start] != comp_counter) {
			tsp_debug(debug_flag, 1, "START = %d, comp[start] = %d ",start, comp[start]); //
			start++;
		}
		i = start;
		fprintf(gnuplotPipe, "%lf %lf\n", points[i].x, points[i].y);
		i = succ[i];
		while(i!=start){
			//tsp_debug(debug_flag, 1, "PLOT COMP %d, i = %d, start = %d", comp_counter, i, start); //
			fprintf(gnuplotPipe, "%lf %lf\n", points[i].x, points[i].y);
			i = succ[i];
		}
		fprintf(gnuplotPipe, "%lf %lf\n", points[i].x, points[i].y);
		tsp_debug(debug_flag, 1, "Plot of component %d SUCCESSFUL",comp_counter);
		comp_counter++;
		// Usa replot per sovrapporre i plot successivi
		if (comp_counter <= ncomp) {
			fprintf(gnuplotPipe, "replot\n");
		}
	}
	fprintf(gnuplotPipe, "e\n"); //This line serves to terminate the input of data for the plot command in gnuplot.
	fclose(gnuplotPipe);
	tsp_debug(debug_flag, 1, "Plot of multitour_sol SUCCESSFUL");
}
*/ // TO ERASE
int TSPopt(instance* inst)
{

	// open CPLEX model
	int error;
	char log_name[64];
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP_Problem");
	multitour_sol curr_sol;
	build_model(inst, env, lp);
	tsp_debug(inst->verbose>=100,1,"build_model SUCCESSFUL");
	// Cplex's parameter setting
	set_init_param(env, (const instance*) inst, log_name, sizeof(log_name));
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
	print_selected_arcs(inst->verbose >= 100, (const double*) xstar, (const instance*) inst);
	tsp_debug(inst->verbose >= 100, 1, "print_seleceted_arcs SUCCESSFUL");
	init_multitour_sol(&curr_sol,inst->nnodes);
	build_sol((const double*) xstar, (const instance*) inst, curr_sol.succ, curr_sol.comp, &curr_sol.ncomp);
	tsp_debug(inst->verbose >= 100, 1, "build_sol SUCCESSFUL: there are %d connected components", curr_sol.ncomp);
	char figure_name[64];
	generate_name(figure_name, sizeof(figure_name), "figures/exact_%d_%d.png", inst->nnodes, inst->randomseed);
	plot_multitour(inst->verbose >= 1, &curr_sol, inst->nnodes, inst->nodes, figure_name);
	free(xstar);
	free_multitour_sol(&curr_sol);
	// free and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
	tsp_debug(inst->verbose >= 100, 1, "TSPopt SUCCESSFUL");
	return 0; // or an appropriate nonzero error code

}