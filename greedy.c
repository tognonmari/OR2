#include "greedy.h"
#include "plot.h"

void gre_update_sol(instance* inst, greedy* gre, int new_best_start) {
	gre->best_start = new_best_start;
	gre->zbest = gre->zcurr;
	gre->best_sol = gre->curr_sol;
	update_best(inst, gre->zbest, get_timer(), gre->best_sol);
}
/*
* Searches for the node with the minimum distance to $p from the points allocated in [p+1 ; end] (address space).
*
* This function iterates through the nodes starting from the specified current node
* and calculates the distance between each node and the reference node (p).
* It updates the minimum distance and the pointer to the closest node accordingly.
* The computed minimum distance is stored in the variable pointed to by current_cost.
*
* IP p Pointer to the reference node for distance calculation.
* IP end Pointer to the last node in the search range.
* OP current_cost Pointer to a variable to store the computed minimum distance.
* OR Pointer to the node with the minimum distance.
* @note A node is represented by a int value which is its label, thus its position in the instance point array $(inst->nodes).
*/
int* gre_search_min(const int* p, const int* end, const float* dist_matrix, double* min) {
	double min_distance = DBL_MAX;
	int* closest_node = NULL;
	int* current = (int*)p + 1;
	while (current <= end) {
		double distance = get_dist_matrix(dist_matrix, *p, *current);
		if (distance < min_distance) {
			min_distance = distance;
			closest_node = current;
		}
		current++;
	}
	(*min) = min_distance;
	return closest_node;
}
/*
* Computes a greedy path starting from a specified index and calculates its cost.
* This function constructs a greedy path starting from the node at the given index.
* It iteratively selects the closest node to the last visited node, swaps it
* into the path, and calculates the cost of the path.
*
* IP index_first Index of the starting node for the greedy path.
* IP inst Pointer to the instance containing node information.
* OP gre Pointer to the structure containing parameters for greedy
*/
int* gre_compute_path(int index_first, const instance* inst, greedy* gre) {
	int n = inst->nnodes;
	int* path = (int*)calloc(n, sizeof(int));
	assert(path != NULL);
	init_path(path, n);
	swap(&(path[0]), &(path[index_first]));
	int* next = &(path[1]); //pointer to the next node to visit in the path (at the start path[0] is already correct)
	int* end = &(path[n - 1]); //pointer to the ending node of the path
	int* min; //pointer to the closest node of the last node visited;
	double current_cost = 0;
	double aggregate_cost = 0;
	while (next < end) {
		min = gre_search_min((next - 1), end, (const float*)(gre->gre_dist_matrix), &current_cost);
		swap(next, min);
		next++;
		aggregate_cost += current_cost;
	}
	current_cost = get_dist_matrix((const float*)(gre->gre_dist_matrix), *(end - 1), *end);
	aggregate_cost += (current_cost + get_dist_matrix((const float*)(gre->gre_dist_matrix), path[0], *end));
	gre->zcurr = aggregate_cost;
	plot_path_on_screen(gre->pipe, gre->plots_on_screen, path, inst->nnodes, gre->zcurr, inst->nodes);
	if (gre->apply_opt2) { 
		opt2(inst, path, &(gre->zcurr), 0);
		plot_path_on_screen(gre->pipe, gre->plots_on_screen, path, inst->nnodes, gre->zcurr, inst->nodes);
	}
	gre_fill_table(gre->table_flag, index_first, gre->zcurr, (gre->zcurr < gre->zbest) );
	return path;
}
// Fill table with "Start", "Cost", "New best"
void gre_fill_table(char table_flag, int start, double cost, char is_new_best) {
	char table_0[32];
	char table_1[32];
	char table_2[32];
	sprintf(table_0, "%d", start);
	sprintf(table_1, "%.3f", cost);
	if (is_new_best) { sprintf(table_2, "Best"); } else { sprintf(table_2, " "); }
	char* table_fields[] = { table_0, table_1, table_2 };
	int num_cols_table = sizeof(table_fields) / sizeof(table_fields[0]);
	make_table_row(table_flag, stdout, num_cols_table, table_fields);
}
void gre_init_table(char table_flag) {
	if (!table_flag) { return; }
	char* table_fields[] = { "Start", "Cost", "New best" };
	int num_cols_table = sizeof(table_fields) / sizeof(table_fields[0]);
	make_first_row(table_flag, stdout, num_cols_table, "GREEDY");
	make_table_row(table_flag, stdout, num_cols_table, table_fields);
}
void gre_init(instance* inst, greedy* gre) {
	gre->check_feasibility = 1;
	gre->plots_on_screen = 0; //it slower a lot the program
	gre->table_flag = (inst->verbose) >= 1;
	gre->pipe = _popen("gnuplot -persist", "w");
	gre->gre_dist_matrix = inst->dist_matrix;
	gre->best_start = 0;
	gre->zbest = DBL_MAX;
	gre_init_table(gre->table_flag);
	gre->curr_sol = gre_compute_path(0, inst, gre);
	check_sol_is_feasible(gre->check_feasibility, inst, gre->curr_sol, gre->zcurr);
	gre_update_sol(inst, gre, 0);
	update_best(inst, gre->zbest, get_timer(), gre->best_sol);
}
void gre_close(greedy* gre) {
	fclose(gre->pipe);
	make_last_row(gre->table_flag, stdout, 3);
}
void gre_solve(instance* inst, char apply_opt2) {
	greedy gre;
	gre.apply_opt2 = apply_opt2;
	gre_init(inst, &gre);
	for (int i = 1; i < inst->nnodes; i++) {
		if (is_time_limit_exceeded(inst->timelimit)) {
			gre_close(&gre);
			tsp_debug(inst->verbose, 1, "Could not finish greedy NN due to time constraints: visited up until node %d, best start with %d", i, gre.best_start);
			return;
		}
		gre.curr_sol = gre_compute_path(i,inst,&gre);
		check_sol_is_feasible(gre.check_feasibility, inst, gre.curr_sol, gre.zcurr);
		if (gre.zcurr < gre.zbest) {
			gre_update_sol(inst, &gre, i);
		}
		else {
			free(gre.curr_sol); //sol is not the best path, then i can free it
		}

	}
	gre_close(&gre);	
}