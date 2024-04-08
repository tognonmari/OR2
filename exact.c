#include "exact.h"
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
	tsp_debug((inst->verbose == 101),0, "LIST OF SUCCESSFUL ITERATIONS");
	for (int i = 0; i < n; i++)
	{
		tsp_debug((inst->verbose == 101), 0, "	i=%d", i);
		for (int j = i + 1; j < n; j++)
		{
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);  // x(1,2), x(1,3) ....
			double obj = dist(i, j, inst); // cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			tsp_debug((inst->verbose == 101), 0, "		j=%d");
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

	if (VERBOSE >= -100) CPXwriteprob(env, lp, "model.lp", NULL);
	tsp_debug(inst->verbose >= 100, 1, "CPXwriteprob SUCCESSFUL");
	free(cname[0]);
	free(cname);
	tsp_debug(inst->verbose >= 100, 1, "free cname SUCCESSFUL");
}

void set_init_param(CPXENVptr env, const instance* inst){
	// increased precision for big-M models
	CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);		// very important if big-M is present
	CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9);

	CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4);
	if (inst->verbose >= 200) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);

	CPXsetdblparam(env, CPX_PARAM_TILIM, CPX_INFBOUND + 0.0);

	//CPXsetintparam(env, CPX_PARAM_NODELIM, 0); 		// abort Cplex after the root node
	//CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);	// abort Cplex after the first incumbent update
	//CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-4);  	// abort Cplex when gap below 10%
}
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

		// go to the next component...
	}
}

int TSPopt(instance* inst)
{

	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP_Problem");

	build_model(inst, env, lp);
	tsp_debug(inst->verbose>=100,1,"build_model SUCCESSFUL");
	// Cplex's parameter setting
	set_init_param(env, (const instance*) inst);
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
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = i + 1; j < inst->nnodes; j++)
		{
			if (xstar[xpos(i, j, inst)] > 0.5) printf("  ... x(%3d,%3d) = 1\n", i + 1, j + 1);
		}
	}
	tsp_debug(inst->verbose >= 100, 1, "print solution SUCCESSFUL");
	free(xstar);

	// free and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
	tsp_debug(inst->verbose >= 100, 1, "TSPopt SUCCESSFUL");
	return 0; // or an appropriate nonzero error code

}