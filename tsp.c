#include "tsp.h"
#include "utils.h"
#include "cpx.h"
#include "tabu.h"
#include "vns.h"
#include "greedy.h"
#include "hardfixing.h"
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include "plot.h"

/*
* Conditional debugging function that prints
* a formatted message to the console if the given $flag is true,
* the time of the execution is printed together with the message
* if the $time_flag is true.
* The printed text always start with \n.
*/
void tsp_debug(char flag, int time_flag, char* format, ...)
{
	if (flag) {
		va_list args;
		va_start(args, format);
		if (time_flag) { printf("\n%12.6f|", get_timer()); }
		else { printf("\n%12s|", ""); }
		vprintf(format, args);
		va_end(args);
		return;
	}
}
/*
* Works similarly to @func tsp_debug, but it is used to print the formatted
* message on the same line and it doesn't print the execution time.
*/
void tsp_debug_inline(char flag, char* format, ...)
{
	if (flag) {
		va_list args;
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
		return;
	}
}

/*
* Function that prints the best_solution of $inst if the given $flag is true.
*/
void print_best_sol(char flag, instance* inst) {
	if (!(inst->is_best_sol_avail)) {
		tsp_debug_inline(flag, "\n Warning. You try to print best solution, but it is not available\n");
		return;
	}
	check_sol_is_feasible(inst->check_feasibility, inst, inst->best_sol, inst->zbest);
	tsp_debug_inline(flag, "\n -------------- Best solution: --------------\n");
	tsp_debug_inline(flag, "Best Distance (zbest): %.2lf\n", inst->zbest);
	tsp_debug_inline(flag, "Best Solution found at time: %12.6f", inst->tbest);
	tsp_debug_inline(flag && (inst->verbose >= 500), "The following sequence of nodes is the best solution to the problem:");
	print_path(flag && (inst->verbose >= 500), "", inst->best_sol, inst->nodes, inst->nnodes);
	tsp_debug_inline(flag, "\n ------------ End best solution: ------------\n");
}

/**
 * Prints the elements of a triangular matrix up to the main diagonal.
 * IP flag Print flag
 * IP matrix A pointer to a 2D array representing the triangular matrix saved in memory with a contiguos 1D array.
 * IP nrows The number of rows in the triangular matrix.
 * OV The triangular matrix represented by $matrix if $flag and if ($nrows<=10)
 */
void print_triangular_matrix_as_array(char flag, const char text_to_print[], const float* matrix, int nrows) {
	if (nrows > 10) { return; }
	if (text_to_print[0] != '\0') {
		tsp_debug(flag, 0, "%s", text_to_print);
	}
	for (int row = 0; row < nrows; row++) {
		tsp_debug_inline(flag, "\n");
		for (int col = 0; col <= row; col++) {
			tsp_debug_inline(flag, "%8.1f ", get_dist_matrix(matrix, row, col));
		}
	}
}
/* Print the coordinates of a point
* IP flag Print flag
* IP text_to_print
* IP p Point to print
* OV $text_to print followed by $point if $flag is true.
*/
void print_point(char flag, const char text_to_print[], const point* p) {
	int x = (int)(p->x);
	int y = (int)(p->y);
	if (text_to_print[0] != '\0') {
		tsp_debug(flag, 0, "%s", text_to_print);
	}
	tsp_debug_inline(flag, "(%d, %d)", x, y);
}
/* Print the coordinates of the nodes of the graph
* IP flag Print flag
* IP text_to_print.
* IP inst Graph to print
* IP n Size of the graph
* OV nodes of the graph associated to &inst if $flag.
*/
void print_nodes(char flag, const char text_to_print[], const instance* inst, int n) {
	tsp_debug(flag, 1, "%s", text_to_print);
	point* current_point = inst->nodes;
	for (int i = 0; i < n; i++) {
		tsp_debug(flag, 0, "node[%d]: ", i);
		print_point(flag, "", current_point);
		current_point++;
	}
}
/* Print the coordinates of the point of a path
* IP flag Print flag
* IP text_to_print.
* IP path Path to print (sequence of indices of array $nodes)
* IP nodes Array of nodes on which the path is defined
* IP n Size of the path
* OV Path $path if $flag
*/
void print_path(char flag, const char text_to_print[], const int* path, const point* nodes, int n) {
	if (!flag) { return; }
	if (text_to_print[0] != '\0') {
		tsp_debug(flag, 1, "%s", text_to_print);
	}
	for (int i = 0; i < n; i++) {
		tsp_debug(flag, 0, "node[%d]: ", i);
		print_point(flag, "", &nodes[path[i]]);
	}
}

/*
* The points in an array are transformed s.t. their coordinates assume random values.
* x field assume random values belonging to [0,$max_x]
* y field assume random values belonging to [0,$max_y]
* IP n Length of the array
* OP nodes Array of points to be modified with random values
* IP max_x, max_y Maximum value for x and y.
*/
void generate_nodes(int n, point* nodes, int max_x, int max_y) {
	point* pointer_nodes = nodes;
	for (int i = 0; i < n; i++) {
		pointer_nodes->x = rand() % (max_x + 1);
		pointer_nodes->y = rand() % (max_y + 1);
		pointer_nodes++;
	}
}
void generate_name(char buffer[], size_t bufferSize, const char* format, ...) {
	va_list args;
	va_start(args, format);
	int length = vsnprintf(buffer, bufferSize, format, args);
	if (length >= bufferSize) {
		exit(main_error_text(-4, "%d %d", bufferSize, length));
	}
	va_end(args);
}

/*
* A complete graph is implicited defined by its nodes -> this function generates $n nodes
* in a random way, thus each generated pair of coordinates (x_i,y_i) is s.t.:
* x_i belongs to [0,$max_x]
* x_i belongs to [0,$max_y]
* IP n Number of nodes to generate
* IP seed Seed of the pseudorandom sequence
* IOP inst Instance to be initialized with random nodes using $(inst->randomseed) as rarndom seed
*/
void generate_instance(instance* inst) {

	if (strcmp(inst->input_file, "0") == 0) {
		//usual generation

		depolarize_pseudornd_seq();

		//only generate random instance if file name is not provided, so only if no tsplib inst is specified

		inst->nodes = (point*)malloc(inst->nnodes * sizeof(point));
		generate_nodes(inst->nnodes, inst->nodes, MAX_X, MAX_Y);
		inst->best_ub = DBL_MAX;
		inst->best_lb = DBL_MIN;
		inst->is_best_sol_avail = 0;
		compute_dist_matrix(inst);
	}
	else {
		generate_instance_from_tsplib(inst);
		compute_dist_matrix(inst);
	}
	//print_triangular_matrix((inst->verbose>0),"The cost matrix is: \n", (const float**)inst->dist_matrix, inst->nnodes);
}
/*
 * Calculates the Euclidean distance between two points in a two-dimensional space.
 * This function uses the Euclidean distance formula:
 * distance = sqrt((x2 - x1)^2 + (y2 - y1)^2)
 * IP p1 Point 1.
 * IP p2 Point 2.
 * OR Euclidean distance between $p1 and $p2
 */
double get_distance(const point* p1, const point* p2) {
	double x1, x2, y1, y2;
	x1 = p1->x; x2 = p2->x; y1 = p1->y; y2 = p2->y;
	double dist = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
	return dist;
}
/*
* Compute the dist matrix, thus the matrix containing in each entry (i,j) the
* distance between node i and node j.
* IOP inst The $inst field $cost_matrix is computed with the distance of the nodes.
*/
void compute_dist_matrix(instance* inst) {
	int nrow = inst->nnodes;
	float* matrix = (float*)alloc_triangular_matrix_as_array(nrow, sizeof(float));
	float* curr = matrix;
	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j <= i; j++) {
			(*curr) = get_distance(&(inst->nodes[i]), &(inst->nodes[j]));
			curr++;
		}
	}
	inst->dist_matrix = matrix;
}

/*
 * Ottiene il valore di un elemento della matrice simmetrica $matrix
 * IP matrix La matrice simmetrica.
 * IP a L'indice di riga (0-based).
 * IP b L'indice di colonna (0-based).
 * OR Il valore dell'entry associata ad a e b.
 */
float get_dist_matrix(const float* matrix, int row, int col) {
	if (col <= row) { return matrix[(row * (row + 1) / 2) + col]; }
	return matrix[(col * (col + 1) / 2) + row];
}

/*
* The file data_file is modified in the following way:
* each line of the file contains the coordinates of each point
* separeted by a space.
* IP inst Instance from which the data are taken
* OF data_file File to modify
*/
void make_datafile(instance* inst, FILE* data_file) {
	point* current_point = inst->nodes;
	for (int i = 0; i < inst->nnodes; i++) {
		int x = (int)(current_point->x);
		int y = (int)(current_point->y);
		fprintf(data_file, "%d %d\n", x, y);
		current_point++;
	}
}

void parse_command_line(int argc, char** argv, instance* inst) {

	//Default parameters

	// inst->model_type = 0;
	// inst->old_benders = 0;
	strcpy(inst->input_file, "0");
	inst->randomseed = 0;
	// inst->num_threads = 0;
	inst->timelimit = 10; //seconds
	// inst->cutoff = -1; 
	inst->integer_costs = 0;
	inst->tenure_is_variable = 1;
	inst->verbose = VERBOSE;
	inst->available_memory = 12000;
	inst->zbest = DBL_MAX;
	// inst->max_nodes = -1; 	
	inst->nnodes = 0;
	inst->ncols = 0;
	inst->check_feasibility = 1;
	inst->solver = NOT_DEF;
	inst->sf_k_start = 20;
	inst->sf_k_scaling = 10;

	int help = 0; if (argc < 1) help = 1;
	for (int i = 1; i < argc; i++)
	{
		//hf_param
		if (strcmp(argv[i], "-hf2opt") == 0) { inst->hf_opt2 = atoi(argv[++i]); continue; }
		if (strcmp(argv[i], "-hfpstart") == 0) { inst->hf_pfix_start = atof(argv[++i]); continue; }
		if (strcmp(argv[i], "-hfpscal") == 0) { inst->hf_pfix_scaling = atof(argv[++i]); continue; }
		//end hf_param
		//sf_param
		if (strcmp(argv[i], "-sfkstart") == 0) { inst->sf_k_start = atoi(argv[++i]); continue; }
		if (strcmp(argv[i], "-sfkscal") == 0) { inst->sf_k_scaling = atoi(argv[++i]); continue; }
		//end sf_param
		//gre_param
		if (strcmp(argv[i], "-greperc") == 0) { inst->gre_perc_start = atof(argv[++i]); continue; }
		//end gre_param
		//tabu param
		if (strcmp(argv[i], "-varten") == 0) { inst->tenure_is_variable = atoi(argv[++i]); continue; }
		if (strcmp(argv[i], "-tabufreq") == 0) { inst->tabu_freq = atof(argv[++i]); continue; }
		if (strcmp(argv[i], "-tabuavg") == 0) { inst->tabu_avg = atof(argv[++i]); continue; }
		if (strcmp(argv[i], "-tabuamp") == 0) { inst->tabu_amp = atof(argv[++i]); continue; }
		//end tabu param
		//vns param
		if (strcmp(argv[i], "-vns_min_kicks") == 0) { inst->vns_min_kicks = atof(argv[++i]); continue; }
		if (strcmp(argv[i], "-vns_max_kicks") == 0) { inst->vns_max_kicks = atof(argv[++i]); continue; }
		//end vns param
		if ( strcmp(argv[i], "-nnodes") == 0) { inst->nnodes = atoi(argv[++i]); continue; } 			// input file
		if (strcmp(argv[i], "-file") == 0) { strcpy(inst->input_file, argv[++i]); printf("I have just wriitten into instance this string %s", inst->input_file); continue; } 			// input file
		if (strcmp(argv[i], "-input") == 0) { strcpy(inst->input_file, argv[++i]); continue; } 			// input file
		if (strcmp(argv[i], "-f") == 0) { strcpy(inst->input_file, argv[++i]); continue; } 				// input file
		if (strcmp(argv[i], "-tl") == 0) { inst->timelimit = atof(argv[++i]); continue; }		// total time limit
		if (strcmp(argv[i], "-solver") == 0) { inst->solver = parse_solver(argv[++i]); continue; }		//solver
		if (strcmp(argv[i], "-s") == 0) { inst->solver = parse_solver(argv[++i]); continue; }			//solver
		// if ( strcmp(argv[i],"-model_type") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 	// model type
		// if ( strcmp(argv[i],"-old_benders") == 0 ) { inst->old_benders = atoi(argv[++i]); continue; } 	// old benders
		// if ( strcmp(argv[i],"-model") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 			// model type
		if (strcmp(argv[i], "-seed") == 0) { inst->randomseed = abs(atoi(argv[++i])); continue; } 		// random seed
		// if ( strcmp(argv[i],"-threads") == 0 ) { inst->num_threads = atoi(argv[++i]); continue; } 		// n. threads
		if (strcmp(argv[i], "-memory") == 0) { inst->available_memory = atoi(argv[++i]); continue; }	// available memory (in MB)
		if (strcmp(argv[i], "-verbose") == 0) { inst->verbose = atoi(argv[++i]); continue; }				// verbosity
		if (strcmp(argv[i], "-v") == 0) { inst->verbose = atoi(argv[++i]); continue; }					// verbosity
		if (strcmp(argv[i], "-csv_column_name") == 0) { strcpy(inst->csv_column_name, argv[++i]); continue; }
		//if ( strcmp(argv[i],"-node_file") == 0 ) { strcpy(inst->node_file,argv[++i]); continue; }		// cplex's node file
		// if ( strcmp(argv[i],"-max_nodes") == 0 ) { inst->max_nodes = atoi(argv[++i]); continue; } 		// max n. of nodes
		// if ( strcmp(argv[i],"-cutoff") == 0 ) { inst->cutoff = atof(argv[++i]); continue; }				// master cutoff
		// if ( strcmp(argv[i],"-int") == 0 ) { inst->integer_costs = 1; continue; } 						// inteher costs
		if (strcmp(argv[i], "-test_bed_size") == 0) { break; }//ignore the parameter }
		if (strcmp(argv[i], "-help") == 0) { help = 1; continue; } 									// help
		if (strcmp(argv[i], "--help") == 0) { help = 1; continue; } 									// help
		help = 1;
	}

	if (help) exit(1);
}

void print_instance_parameters(instance* inst) {
	printf("-------- Selected input parameters: -------\n");
	printf("File: %s\n", inst->input_file);
	if (inst->timelimit > 0) { printf("Time Limit: %lf\n", inst->timelimit); }
	else { printf("Time Limit: Unsetted\n"); }
	// printf("Model Type: %d\n", inst.model_type); 
	// printf("Old Benders: %d\n", inst.old_benders); 
	printf("Number of nodes: %d\n", inst->nnodes);
	printf("Seed: %d\n", inst->randomseed);
	printf("Level of verbosity: %d\n", inst->verbose);
	// printf("Threads: %d\n", inst.num_threads);  
	// printf("Max Nodes: %d\n", inst.max_nodes); 
	printf("Memory: %d\n", inst->available_memory);
	// printf("Integer Costs: %d\n", inst.integer_costs); 
	// printf("Node File: %s\n", inst.node_file);
	// printf("Cutoff: %lf\n", inst.cutoff); 
	printf("-------------------------------------------\n");
}

void free_matrix(void** matrix, int rows) {
	for (int i = 0; i < rows; i++) {
		free(matrix[i]);
	}
	free(matrix);

}
void free_instance(instance* inst) {
	//TODO: free memory accroding to how instance is allocated
	// free(inst->demand);
	free(inst->nodes);
	if (inst->is_best_sol_avail) {
		free(inst->best_sol);
		inst->is_best_sol_avail = 0;
	}
	free(inst->dist_matrix);
}

double compute_path_length(int* path, int nodes_number, point* nodes) {
	int a, b;

	double path_length = 0;
	for (int i = 0; i < nodes_number - 1; i++) {
		a = path[i];
		b = path[i + 1];
		path_length += (float)get_distance(&nodes[a], &nodes[b]);
	}
	//last edge
	a = path[0];

	return (path_length + (float)get_distance(&nodes[a], &nodes[b]));
}

/**
 * Reverses a sequence within limits min and max, not included, in the array path.
 * IP: array of integers path, limits to the segment to reverse ( min and max )
 * OP: none
*/
void reverse_sequence(int* path, int min, int max) {

	int* s = &path[min + 1];
	int* t = &path[max - 1];
	while (s < t) {

		swap(s, t);
		s++;
		t--;

	}

}

/**
	Given a solution (i.e. a path), it swaps the positions of two nodes if an improvement is found.
	IP: instance inst passed by reference
	OP: formally none, but inst's best_sol is modified
*/
void swap_2_opt(int* path, int i, int j) {

	//swapping the extremes
	int temp = path[i];
	path[i] = path[j];
	path[j] = temp;
	int min = i;
	int max = j;
	if (i > j) {
		min = j;
		max = i;
	}
	//reverse the sequence

	reverse_sequence(path, min, max);

}
/*
 * Swaps the values of two integers.
 * IOP a Pointer to the first integer.
 * IOP b Pointer to the second integer.
 */
void swap(int* a, int* b) {
	int temp = *a;
	*a = *b;
	*b = temp;
}
void opt2_init_table(char table_flag) {
	if (!table_flag) { return; }
	char* table_fields[] = { "Nr swap", "node1", "node2", "Delta", "New cost" };
	int num_cols_table = sizeof(table_fields) / sizeof(table_fields[0]);
	make_first_row(table_flag, stdout, num_cols_table, "2OPT_SWAPS_TABLE");
	make_table_row(table_flag, stdout, num_cols_table, table_fields);
}
void opt2_fill_table(char table_flag, int* nr_swap, int best_i, int best_j, double best_delta, double incumbent_cost) {
	if (!table_flag) { return; }
	char table_0[32];
	char table_1[32];
	char table_2[32];
	char table_3[32];
	char table_4[32];
	sprintf(table_0, "%d", ++(*nr_swap));
	sprintf(table_1, "%d", best_i + 1);
	sprintf(table_2, "%d", best_j);
	sprintf(table_3, "%.1f", best_delta);
	sprintf(table_4, "%.1f", incumbent_cost);
	char* table_fields[] = { table_0, table_1, table_2, table_3, table_4 };
	int num_cols_table = sizeof(table_fields) / sizeof(table_fields[0]);
	make_table_row(table_flag, stdout, num_cols_table, table_fields);
}
void opt2(instance* inst, int* incumbent_sol, double* incumbent_cost, char table_flag) {
	opt2_init_table(table_flag);
	int nr_swap = 0;
	while (opt2_move(table_flag, inst, incumbent_sol, incumbent_cost, &nr_swap))
	{
	}
	make_last_row(table_flag, stdout, 5);
}

/**
 * Returns 0 if the path cannot be optimized any further with 2 opt moves
*/
char opt2_move(char table_flag, instance* inst, int* incumbent_sol, double* incumbent_cost, int* nr_swap) {
	int num_cols_table = 5;
	int n = inst->nnodes;
	point* nodes_list = inst->nodes;
	double best_delta = 0;
	char improvement = 0;
	int best_i = -1;
	int best_j = -1;
	float* dist_matrix = inst->dist_matrix;
	for (int i = 0; i <= n - 3; i++) //cambiato da n-2!
	{
		for (int j = i + 2; j <= n - 1; j++)
		{
			double current_dist = (double)get_dist_matrix((const float*)(dist_matrix), incumbent_sol[i], incumbent_sol[(i + 1)]) + (double)get_dist_matrix((const float*)(dist_matrix), incumbent_sol[j], incumbent_sol[(j + 1) % n]);
			double changed_dist = (double)get_dist_matrix((const float*)(dist_matrix), incumbent_sol[i], incumbent_sol[(j)]) + (double)get_dist_matrix((const float*)(dist_matrix), incumbent_sol[i + 1], incumbent_sol[(j + 1) % n]);

			//double current_dist = get_distance(&nodes_list[incumbent_sol[i % n]], &nodes_list[incumbent_sol[(i + 1) % n]]) + get_distance(&nodes_list[incumbent_sol[j % n]], &nodes_list[incumbent_sol[(j + 1) % n]]);
			//double changed_dist = get_distance(&nodes_list[incumbent_sol[i % n]], &nodes_list[incumbent_sol[j % n]]) + get_distance(&nodes_list[incumbent_sol[(i + 1) % n]], &nodes_list[incumbent_sol[(j + 1) % n]]);

			double delta = changed_dist - current_dist;

			if (delta < best_delta)
			{
				improvement = 1;
				//printf("found a new best delta\n");
				best_i = i % n;
				best_j = j % n;
				best_delta = delta;
				//printf("Best delta : %lf\n",best_delta);
			}
		}
	}
	if (improvement) {
		tsp_debug(0, 0, "I am swapping nodes %d and %d", best_i + 1, best_j);
		swap_2_opt(incumbent_sol, (best_i + 1) % n, (best_j) % n);
		(*incumbent_cost) += best_delta;
		opt2_fill_table(table_flag, nr_swap, best_i + 1, best_j, best_delta, *incumbent_cost);
		/*
		char table_0[32];
		char table_1[32];
		char table_2[32];
		char table_3[32];
		char table_4[32];
		sprintf(table_0, "%d", ++(*nr_swap));
		sprintf(table_1, "%d", best_i + 1);
		sprintf(table_2, "%d", best_j);
		sprintf(table_3, "%.1f", best_delta);
		sprintf(table_4, "%.1f", *incumbent_cost);
		char* table_fields[] = { table_0, table_1, table_2, table_3, table_4 };
		make_table_row(table_flag, stdout, num_cols_table, table_fields);
		*/
		return 1;
	}
	else {
		//printf("no further improvement with 2opt.\n");
		return 0;
	}
}

/*
 * Generic swap function for swapping data of arbitrary types.
 * This function takes two pointers to data ($a and $b) and the size of each data element.
 * IP a,b Pointers to the first and the second data elements.
 * IP size Size of each data element in bytes.
 * ACTUALLY NOT USED!
 */
void swap_space(void* a, void* b, size_t size) {
	void* temp = malloc(size);
	assert(temp != NULL);
	memcpy(temp, a, size);
	memcpy(a, b, size);
	memcpy(b, temp, size);
	free(temp);
}
/*
 * Copies data from one dynamically allocated array to another.
 * This function checks if the sizes of the input arrays in the heap are equal.
 * If the sizes match, it performs a memory copy from the source array $a2 to
 * the destination array $a1.
 * IP a1 Pointer to the destination array.
 * IP a2 Pointer to the source array.
 */
void copy_array(void* a1, const void* a2) {
	size_t size_a1 = _msize(a1);
	/*
	* the(void*)a2 cast is used to reassure the compiler that we are intentionally treating a2 as a non-constant pointer
	* for the purposes of this specific function call.
	*/
	size_t size_a2 = _msize((void*)a2);
	if (size_a1 == size_a2) {
		memcpy(a1, a2, size_a1);
	}
	else {
		fprintf(stderr, "Error: The program tried to copy arrays that don't have the same space available in the heap.\n");
		fprintf(stderr, "Size a1 = %d , size a2 = %d\n", (int)size_a1, (int)size_a2);
	}
}
/*
 * Copies data from one dynamically allocated array to another.
 * This function checks if the sizes of the input arrays in the heap are equal.
 * If the sizes match, it performs a memory copy from the source array a2 to
 * the destination array a1.
 *   OP a1 Pointer to the destination array.
 *   IP a2 Pointer to the source array.
 *   IP elem_size Size of each element in bytes.
 *   IP num_elems Number of elements in the arrays.
 */
void copy_din_array(void* a1, const void* a2, size_t elem_size, size_t num_elems) {
	if (a1 == NULL || a2 == NULL) {
		fprintf(stderr, "Error: null pointers.\n");
		exit(main_error(-9));
	}
	size_t total_size = elem_size * num_elems;
	memcpy(a1, a2, total_size);
}

/**
 * Updates the best solution information in the instance.
 * OP inst Pointer to the instance structure to modify.
 * IP z New best objective value.
 * IP t New best execution time.
 * IP sol Pointer to the best solution (array of integers).
 */
void update_best(instance* inst, double z, double t, int* sol) {
	char flag = inst->verbose >= 1;
	tsp_debug_inline(flag, "\n -------------- Update Best : --------------\n");
	tsp_debug_inline(flag, "Old Best : %.20g\n", inst->zbest);
	tsp_debug_inline(flag, "New Best : %.20g\n", z);
	tsp_debug_inline(flag, "\n --------------- End Update : --------------\n");
	inst->zbest = z;
	inst->tbest = t - inst->tstart;
	inst->is_best_sol_avail = 1;
	inst->best_sol = sol;
	check_sol_is_feasible(inst->check_feasibility, inst, inst->best_sol, inst->zbest);
}
void update_lb(instance* inst, double lb) {
	if (lb > inst->best_lb) {
		inst->best_lb = lb;
	}
}
void update_ub(instance* inst, double ub) {
	if (ub < inst->best_ub) {
		inst->best_ub = ub;
	}
}
/*
 * Initializes an array $path with values from 0 to n-1.
 * OP path Pointer to the array to be initialized.
 * IP n Number of elements in the array.
 */
void init_path(int* path, size_t n) {
	for (int i = 0; i < n; i++) {
		(*path) = i;
		path++;
	}
}
void init_data_file(char flag, FILE* data_file, instance* inst) {
	char figure_name[64];
	if (flag) {
		generate_name(figure_name, sizeof(figure_name), "data/tabu_%d_%d_cost.dat", inst->nnodes, inst->randomseed);
		data_file = fopen(figure_name, "w");
		if (data_file == NULL) {
			fclose(data_file);
			exit(main_error_text(-2, "Failed to open the file %s", figure_name));
		}
	}
}
/*
* Checks if the input solution $sol is feasible. If it is not feasible returns error.
*/
void check_sol_is_feasible(char check_feasibility, instance* inst, int* sol, double zsol) {
	if (!check_feasibility) { return; }
	if (!is_feasible_solution(inst, sol, zsol)) {
		exit(main_error_text(-99, "Stop! Check the code! Current solution is not feasible."));
	}
}
/**
 * Checks whether sol_path and sol_cost consitute a feasible solution for inst
 * IP: instance i
 * IP: sol_path , candidate path
 * IP: sol_cost, cost of the candidate solution
 * OP: 1 if sol_path and sol_Cost are coherent with one another and sol_path consitutes a valid path
*/
char is_feasible_solution(instance* inst, int* sol_path, double sol_cost) {
	char flag = (inst->verbose >= 200);
	int n = inst->nnodes;
	int* counter = (int*)calloc(n, sizeof(int));
	//print_path(1, "path: \n", sol_path, inst->nodes,n);
	//check for repeated vertices or out of range content of path
	for (int i = 0; i <= n - 1; i++) {

		counter[sol_path[i]]++;

	}
	for (int i = 0; i <= n - 1; i++) {

		if (counter[sol_path[i]] != 1) {
			tsp_debug(flag, 0, "Node %d is visited %d", i, counter[sol_path[i]]);
			free(counter);
			return 0;

		}
	}

	free(counter);

	//check cost coincides with the path length 

	double path_length = compute_path_length(sol_path, n, inst->nodes);
	tsp_debug(flag, 0, "Real Path length : %.6f", path_length);
	tsp_debug(flag, 0, "Sol passed length : %.6f", sol_cost);
	tsp_debug(flag, 0, "is_feasible = %d", is_equal_double(path_length, sol_cost, EPSILON));
	return is_equal_double(path_length, sol_cost, EPSILON);
}

/**
 * Utility method to return a solver given a command line input
 * IP: string from command line
 * OP: appropriate solver_id, NOT_DEF if incorrectly specified
*/

solver_id parse_solver(char* solver_input) {
	if (strcmp(solver_input, "greedy") == 0 || strcmp(solver_input, "nn") == 0) {
		return NN;
	}
	else if (strcmp(solver_input, "random") == 0 || strcmp(solver_input, "rand") == 0) {
		return RANDOM_SOL;
	}
	else if (strcmp(solver_input, "2opt") == 0 || strcmp(solver_input, "opt2") == 0) {
		return OPT_2;
	}
	else if (strcmp(solver_input, "tabu") == 0) {
		return TABU;
	}
	else if (strcmp(solver_input, "vns") == 0) {
		return VNS;
	}
	else if (strcmp(solver_input, "exact_cplex_only") == 0) {
		return EX;
	}
	else if (strcmp(solver_input, "benders") == 0) {
		return BEN;
	}
	else if (strcmp(solver_input, "gluing") == 0) {
		return GLU;
	}
	else if (strcmp(solver_input, "bcp") == 0) {
		return BCP;
	}
	else if (strcmp(solver_input, "bc") == 0) {

		return BC;
	}
	else if (strcmp(solver_input, "bcf") == 0) {
		return BCF;
	}
	else if (strcmp(solver_input, "bcfm") == 0) {
		return BCFM;
	}
	else if (strcmp(solver_input, "bcm") == 0) {
		return BCM;
	}
	else if (strcmp(solver_input, "bcmp") == 0) {
		return BCMP;
	}
	else if (strcmp(solver_input, "bcp") == 0) {
		return BCP;
	}
	else if (strcmp(solver_input, "bcfm") == 0) {
		return BCMP;
	}
	else if (strcmp(solver_input, "bcfp") == 0) {
		return BCFP;
	}
	else if (strcmp(solver_input, "bcp") == 0) {
		return BCP;
	}
	else if (strcmp(solver_input, "bcfmp") == 0) {
		return BCFMP;
	}
	else if (strcmp(solver_input, "hf") == 0) {
		return HF;
	}
	else if (strcmp(solver_input, "sf") == 0) {
		return SF;
	}
	exit(main_error_text(-8, "%s", "solver"));
}


void tsp_solve(instance* inst) {
	inst->tstart = get_timer();
	inst->timelimit += get_timer();
	MKDIR("figures");
	MKDIR("data");
	char figure_name[64];
	int check_truncation = 0;
	switch (inst->solver) {
	case NN:
		gre_partial_solve(inst, 0, (int)(inst->nnodes * inst->gre_perc_start));
		print_best_sol((inst->verbose >= 1), inst);
		//generate_name(figure_name, sizeof(figure_name), "figures/greedy_%d_%d.png", inst->nnodes, inst->randomseed);
		//plot_path(inst->verbose >= 1000, inst->best_sol, inst->nnodes, inst->zbest, inst->nodes, figure_name);
		//init_data_file((inst->verbose>-1),(inst->best_sol_data), inst);
		//plot_generator(inst);
		break;
	case OPT_2:
		gre_partial_solve(inst, 1, (int)(inst->nnodes * inst->gre_perc_start));
		generate_name(figure_name, sizeof(figure_name), "figures/greedy+opt2_%d_%d.png", inst->nnodes, inst->randomseed);
		plot_path(inst->verbose >= 1000, inst->best_sol, inst->nnodes, inst->zbest, inst->nodes, figure_name);
		print_best_sol((inst->verbose >= 1), inst);
		break;
	case TABU:
		tabu_search(inst->tenure_is_variable, inst);
		print_best_sol((inst->verbose >= 1), inst);
		generate_name(figure_name, sizeof(figure_name), "figures/tabu_%d_%d.png", inst->nnodes, inst->randomseed);
		plot_path(inst->verbose >= 1, inst->best_sol, inst->nnodes, inst->zbest, inst->nodes, figure_name);
		//init_data_file((inst->verbose>-1),(inst->best_sol_data), inst);
		break;
	case VNS:
		
		
    	gre_partial_solve(inst,1,1); 
		
		print_best_sol((inst->verbose>=5), inst);
		generate_name(figure_name, sizeof(figure_name), "figures/greedy_%d_%d.png", inst->nnodes, inst->randomseed);
		vns(inst);
		print_best_sol((inst->verbose>=5), inst);
		generate_name(figure_name, sizeof(figure_name), "figures/vns_%d_%d.png", inst->nnodes, inst->randomseed);
		//plot_path(inst->verbose >= 1000, inst->best_sol, inst->nnodes, inst->zbest, inst->nodes, figure_name);
		//init_data_file((inst->verbose>-1),(inst->best_sol_data), inst);
		break;
	case EX:
		TSPopt(inst);
		break;
	case BEN:
		ben_solve(0, inst);
		print_best_sol((inst->verbose >= 0), inst);
		break;
	case GLU:
		ben_solve(1, inst);
		print_best_sol((inst->verbose >= 0), inst);
		break;
	case BC:
		cpx_branch_and_cut(0, 0, 0, inst);
		print_best_sol((inst->verbose >= 0), inst);
		break;
	case BCM:
		cpx_branch_and_cut(0, 1, 0, inst);
		print_best_sol((inst->verbose >= 0), inst);
		break;
	case BCP:
		cpx_branch_and_cut(0, 0, 1, inst);
		print_best_sol((inst->verbose >= 0), inst);
		break;
	case BCMP:
		cpx_branch_and_cut(0, 1, 1, inst);
		print_best_sol((inst->verbose >= 0), inst);
		break;
	case BCF:
		cpx_branch_and_cut(1, 0, 0, inst);
		print_best_sol((inst->verbose >= 0), inst);
		break;
	case BCFM:
		cpx_branch_and_cut(1, 1, 0, inst);
		print_best_sol((inst->verbose >= 0), inst);
		break;
	case BCFP:
		cpx_branch_and_cut(1, 0, 1, inst);
		print_best_sol((inst->verbose >= 0), inst);
		break;
	case BCFMP:
		cpx_branch_and_cut(1, 1, 1, inst);
		print_best_sol((inst->verbose >= 0), inst);
		break;
	case HF:
		hf_solve(inst);
		print_best_sol((inst->verbose >= 0), inst);
		generate_name(figure_name, sizeof(figure_name), "figures/hf_%d_%d.png", inst->nnodes, inst->randomseed);
		plot_path(inst->verbose >= 1, inst->best_sol, inst->nnodes, inst->zbest, inst->nodes, figure_name);
		break;
	case SF:
		sf_solve(inst);
		print_best_sol((inst->verbose >= 0), inst);
		generate_name(figure_name, sizeof(figure_name), "figures/sf_%d_%d.png", inst->nnodes, inst->randomseed);
		plot_path(inst->verbose >= 1, inst->best_sol, inst->nnodes, inst->zbest, inst->nodes, figure_name);
		break;
	default:
		exit(main_error(-7));
	}
}

void update_solver(instance* inst) {

	int selection;
	char buf[8];
	printf("----Choose a solver, from these options:----\n");
	printf("0: Nearest Neighbor\n");
	printf("1: Nearest Neighbor and OPT2\n");
	printf("2: Nearest Neighbor and TABU search\n");
	printf("3: Nearest Neighbor and VNS\n");
	printf("4: Exact Method with CPLEX, no subtour constraints\n");
	printf("5: Exact Method with CPLEX, with SECs, Benders' method\n");
	printf("6: Exact Method with CPLEX, with SECs, Benders' method with patching\n");
	printf("7: Exact Method with CPLEX, Branch and Cut Method, no frational cut\n");
	printf("8: Exact Method with CPLEX, Branch and Cut Method, no fractional cut, MIPSTART enabled\n");
	printf("9: Exact Method with CPLEX, Branch and Cut Method, no fractional cut with posting\n");
	printf("10: Exact Method with CPLEX, Branch and Cut Method, no fractional cut with posting, MIPSTART enabled\n");
	printf("11: Exact Method with CPLEX, Branch and Cut Method, with fractional cut\n");
	printf("12: Exact Method with CPLEX, Branch and Cut Method, with fractional cut, MIPSTART enabled\n");
	printf("13: Exact Method with CPLEX, Branch and Cut Method, with fractional cut and posting\n");
	printf("14: Exact Method with CPLEX, Branch and Cut Method, with fractional cut and posting, MIPSTART enabled\n");
	printf("15: Matheuristic with CPLEX, Hardfixing Method\n");
	printf("16: Matheuristic with CPLEX, Softfixing Method\n");
	printf("---------------------------------------------\n");
	fgets(buf, 8, stdin);
	selection = atoi(buf);

	switch (selection) {
	case 0:
	{
		inst->solver = NN;

		printf("successful update.\n");
		break;
	}
	case 1:
	{
		inst->solver = OPT_2;
		printf("successful update.\n");
		break;
	}
	case 2:
	{
		inst->solver = TABU;
		printf("successful update.\n");
		break;
	}
	case 3:
	{
		inst->solver = VNS;
		printf("successful update. \n");
		break;
	}
	case 4:
	{
		inst->solver = EX;
		printf("successful update. \n");
		break;
	}
	case 5:
	{
		inst->solver = BEN;
		printf("successful update. \n");
		break;
	}
	case 6:
	{
		inst->solver = GLU;
		printf("successful update. \n");
		break;
	}
	case 7:
	{
		inst->solver = BC;
		printf("successful update. \n");
		break;
	}
	case 8:
	{
		inst->solver = BCM;
		printf("successful update. \n");
		break;
	}
	case 9:
	{
		inst->solver = BCP;
		printf("successful update. \n");
		break;
	}
	case 10:
	{
		inst->solver = BCMP;
		printf("successful update. \n");
		break;
	}
	case 11:
	{
		inst->solver = BCF;
		printf("successful update. \n");
		break;
	}
	case 12:
	{
		inst->solver = BCFM;
		printf("successful update. \n");
		break;
	}
	case 13:
	{
		inst->solver = BCFP;
		printf("successful update. \n");
		break;
	}
	case 14:
	{
		inst->solver = BCFMP;
		printf("successful update. \n");
		break;
	}
	case 15:
	{
		inst->solver = HF;
		printf("successful update. \n");
		break;
	}
	case 16:
	{
		inst->solver = SF;
		printf("successful update. \n");
		break;
	}
	default:
	{
		inst->solver = NN;
		printf("successful update.\n");
		break;
	}
	}


}
void generate_test_bed(int size_test_bed, int argc, char** argv, instance* test_bed) {

	/// inizializza la prima istanza e parsa il test bed size

	//con i parametri della prima istanza istanzia le altre 

	for (int i = 0; i < size_test_bed; i++) {
		//initialize instance
		parse_command_line(argc, argv, &test_bed[i]);
		printf("successful parsing \n");
	}
	srand(test_bed[0].randomseed);

	for (int i = 0; i < size_test_bed; i++) {
		generate_instance(&test_bed[i]);
	}
}
void read_test_bed_size(int* test_bed_size, int argc, char** argv) {

	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-test_bed_size") == 0) {
			*(test_bed_size) = abs(atoi(argv[++i]));
			return;
		}
	}

}
void generate_csv_file(int size_test_bed, instance* test_bed) {

	int n = test_bed[0].nnodes;
	int seed = test_bed[0].randomseed;
	MKDIR("runs");
	char csv_name[128];
	generate_name(csv_name, sizeof(csv_name), "runs/%s_%d_%d_%d_cost.csv", test_bed[0].csv_column_name, n, seed, size_test_bed);

	//Print a single-column .csv file
	FILE* file_handler = fopen(csv_name, "w");
	int error;
	//Print the first row
	error = fprintf(file_handler, "1, %s\n", test_bed[0].csv_column_name);
	fflush(file_handler);

	//Print the remaining rows
	for (int i = 0; i < size_test_bed; i++) {
		fprintf(file_handler, "%d_%d_[%d], %lf\n", n, seed, i, test_bed[i].zbest);
		fflush(file_handler);
	}
	//

	error = fclose(file_handler);
	if (error) {
		printf("error in closing csv\n");
	}
	char csv_time_name[128];
	generate_name(csv_time_name, sizeof(csv_time_name), "runs/%s_%d_%d_%d_time.csv", test_bed[0].csv_column_name, n, seed, size_test_bed);
	FILE* time_file_handler = fopen(csv_time_name, "w");
	//Print a single-column .csv file

	//Print the first row
	fprintf(time_file_handler, "1, %s\n", test_bed[0].csv_column_name);
	fflush(time_file_handler);
	//Print the remaining rows
	for (int i = 0; i < size_test_bed; i++) {
		fprintf(time_file_handler, "%d_%d_[%d], %lf\n", n, seed, i, test_bed[i].tbest);
		fflush(time_file_handler);
	}
	//
	fclose(time_file_handler);

}

void generate_instance_from_tsplib(instance* inst) {

	if (strcmp(inst->input_file, "0") == 0) {

		printf("invalid input file\n");
		exit(main_error(-1));
	}

	FILE* input_file = fopen(inst->input_file, "r");

	if (input_file == NULL) {
		exit(main_error(-1));
	}
	char line[100];

	while (fgets(line, sizeof(line), input_file) != NULL) {
		//Read nnodes
		char* id = strtok(line, ": ");
		if (strcmp(id, "DIMENSION") == 0) {
			inst->nnodes = atoi(strtok(NULL, ": "));
			//tsp_debug(1, 1, "nodes numerr: %d\n", inst->nnodes);
			break;
		}
	}

	inst->nodes = (point*)malloc(inst->nnodes * sizeof(point));

	//Read coordinates
	while (fgets(line, sizeof(line), input_file) != NULL) {
		char* id;
		id = strtok(line, "\n");
		if (strcmp(id, "NODE_COORD_SECTION") == 0) {
			//printf("Found node coord section");
			//parse the nodes and add them as points
			break;
		}
	}
	int i = 0;
	while (fgets(line, sizeof(line), input_file) != NULL) {

		if (strncmp(line, "EOF", 3) == 0) {
			//printf("breaking");
			break;
		}

		char* index = strtok(line, " ");// skip it
		char* x = strtok(NULL, " ");
		char* y = strtok(NULL, " ");
		inst->nodes[i].x = atof(x);
		inst->nodes[i].y = atof(y);
		//printf("I have just added node (%.2lf, %.2lf)", inst->nodes[i].x, inst->nodes[i].y);
		i++;
	}



}