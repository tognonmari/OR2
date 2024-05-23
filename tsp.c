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
void tsp_debug(char flag,int time_flag, char* format, ...)
{
	if (flag) {
		va_list args;
		va_start(args, format);
		if (time_flag) { printf("\n%12.6f|", get_timer()); }
        else { printf("\n%12s|",""); }
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
	tsp_debug_inline(flag, "\n -------------- Best solution: --------------\n");
	tsp_debug_inline(flag, "Best Distance (zbest): %.2lf\n",inst->zbest);
	tsp_debug_inline(flag, "Best Solution found at time: %12.6f",inst->tbest);
	tsp_debug_inline(flag && ( inst->verbose >= 500 ), "The following sequence of nodes is the best solution to the problem:");
	print_path(flag && ( inst->verbose >= 500 ),"", inst->best_sol, inst->nodes, inst->nnodes);
	tsp_debug_inline(flag, "\n ------------ End best solution: ------------\n");
}

//ACTUALLY NOT USED BECAUSE DIST_MATRIX IS SAVED AS CONTIGUOS ARRAY
/**
 * Prints the elements of a triangular matrix up to the main diagonal.
 * IP flag Print flag
 * IP matrix A pointer to a 2D array representing the triangular matrix.
 * IP nrows The number of rows in the triangular matrix.
 * OV The triangular matrix represented by $matrix if $flag and if ($nrows<=10)
 */
void print_triangular_matrix(char flag, const char text_to_print[], const float** matrix, int nrows) {
	if(nrows>10){ return; }
	if(text_to_print[0]!= '\0'){
		tsp_debug(flag,0,"%s", text_to_print);
	}
	for (int row = 0; row < nrows; row++) {
		tsp_debug_inline(flag,"\n");
		for (int col = 0; col <= row; col++) {
			tsp_debug_inline(flag,"%8.1f ", get_cost_matrix(matrix,row,col));
		}
	}
}
/**
 * Prints the elements of a triangular matrix up to the main diagonal.
 * IP flag Print flag
 * IP matrix A pointer to a 2D array representing the triangular matrix saved in memory with a contiguos 1D array.
 * IP nrows The number of rows in the triangular matrix.
 * OV The triangular matrix represented by $matrix if $flag and if ($nrows<=10)
 */
void print_triangular_matrix_as_array(char flag, const char text_to_print[], const float* matrix, int nrows) {
	if(nrows>10){ return; }
	if(text_to_print[0]!= '\0'){
		tsp_debug(flag,0,"%s", text_to_print);
	}
	for (int row = 0; row < nrows; row++) {
		tsp_debug_inline(flag,"\n");
		for (int col = 0; col <= row; col++) {
			tsp_debug_inline(flag,"%8.1f ", get_dist_matrix(matrix,row,col));
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
	if(text_to_print[0]!= '\0'){
		tsp_debug(flag,0,"%s", text_to_print);
	}
	tsp_debug_inline(flag,"(%d, %d)", x, y);
}
/* Print the coordinates of the nodes of the graph
* IP flag Print flag
* IP text_to_print.
* IP inst Graph to print
* IP n Size of the graph
* OV nodes of the graph associated to &inst if $flag.
*/
void print_nodes(char flag, const char text_to_print[], const instance *inst, int n) {
	tsp_debug(flag ,1,"%s", text_to_print);
	point* current_point = inst->nodes;
	for (int i = 0; i < n; i++) {
		tsp_debug(flag ,0,"node[%d]: ", i);
		print_point(flag ,"",current_point);
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
    if(text_to_print[0]!= '\0'){
        tsp_debug(flag, 1, "%s", text_to_print);
    }
	for (int i = 0; i < n; i++) {
		tsp_debug(flag, 0, "node[%d]: ", i);
		print_point(flag,"", &nodes[path[i]] );
	}
}
/*
* The values of an array are transformed in random values belonging to [0,$max_value]
* IP n Length of the array
* OP array Array to be modified with random values
* IP max_value Maximum value of the generated values
*/
//ACTUALLY NOT USED TO ERASE
void generate_array(int n, double *array, int max_value) {
	for (int i = 0; i < n; i++) {
		*array = rand() % (max_value+1);
		array++;
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
void generate_name(char buffer[], size_t bufferSize, const char *format, ...) {
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
void generate_instance(instance *inst){
	
	
	depolarize_pseudornd_seq();

	//only generate random instance if file name is not provided, so only if no tsplib inst is specified
	
	inst->nodes = (point*)malloc(inst->nnodes * sizeof(point));
	generate_nodes(inst->nnodes, inst->nodes, MAX_X, MAX_Y);
	inst->best_ub = DBL_MAX;
	inst->best_lb = DBL_MIN;
	inst->is_best_sol_avail = 0;
	compute_dist_matrix(inst);
	
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
void compute_dist_matrix(instance* inst){
	int nrow = inst->nnodes;
	float* matrix = (float*)alloc_triangular_matrix_as_array(nrow, sizeof(float));
	float* curr = matrix;
	for (int i = 0; i < nrow; i++ ) {
		for (int j = 0; j <= i; j++) {
			(*curr) = get_distance(&(inst->nodes[i]),&(inst->nodes[j]));
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
float get_dist_matrix(const float* matrix, int row, int col){
	if(col<=row) {return matrix[(row * (row+1) / 2) + col];}
	return matrix[(col * (col+1) / 2) + row];
}

////ACTUALLY NOT USED BECAUSE DIST_MATRIX IS SAVED AS CONTIGUOS ARRAY AND THIS FUNCTION WORKS FOR 2D ARRAYS
/*
 * Ottiene il valore di un elemento della matrice simmetrica $matrix
 * IP matrix La matrice simmetrica.
 * IP a L'indice di riga (0-based).
 * IP b L'indice di colonna (0-based).
 * OR Il valore dell'entry associata ad a e b.
 */
double get_cost_matrix(const float** matrix, int a, int b){
	if(b<=a) { return matrix[a][b];}
	return matrix[b][a];
}
//TO ERASE

/*
* The file data_file is modified in the following way:
* each line of the file contains the coordinates of each point
* separeted by a space.
* IP inst Instance from which the data are taken
* OF data_file File to modify
*/
void make_datafile(instance *inst, FILE* data_file) {
	point* current_point = inst->nodes;
	for (int i = 0; i < inst->nnodes; i++) {
		int x =(int) (current_point->x);
		int y =(int) (current_point->y);
		fprintf(data_file, "%d %d\n", x, y);
		current_point++;
	}
}

void parse_command_line(int argc, char** argv, instance *inst){

    //Default parameters

    // inst->model_type = 0;
	// inst->old_benders = 0;
	strcpy(inst->input_file, "0");
	inst->randomseed = 0; 
	// inst->num_threads = 0;
	inst->timelimit = 10; //seconds
	// inst->cutoff = -1; 
	inst->integer_costs = 0;
    inst->verbose = VERBOSE;
	inst->available_memory = 12000;   
	inst->zbest = DBL_MAX;
	// inst->max_nodes = -1; 	
	inst->nnodes =0;
	inst->ncols = 0;
	inst->solver = NOT_DEF;
    int help = 0; if ( argc < 1 ) help = 1;	
	for ( int i = 1; i < argc; i++ ) 
	{	
		//hf_param
		if (strcmp(argv[i], "-hf2opt") == 0) { inst->hf_opt2 = atoi(argv[++i]); continue; }
		if (strcmp(argv[i], "-hfpstart") == 0) { inst->hf_pfix_start = atof(argv[++i]); continue; }
		if (strcmp(argv[i], "-hfpscal") == 0) { inst->hf_pfix_scaling = atof(argv[++i]); continue; }
		//end hf_param

		if ( strcmp(argv[i], "-nnodes") == 0) { inst->nnodes = atoi(argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-input") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-f") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 				// input file
		if ( strcmp(argv[i],"-tl") == 0 ) { inst->timelimit = atof(argv[++i]); continue; }		// total time limit
		if ( strcmp(argv[i], "-solver") == 0){ inst->solver = parse_solver(argv[++i]); continue; }		//solver
		if ( strcmp(argv[i], "-s") == 0){ inst->solver = parse_solver(argv[++i]); continue; }			//solver
		// if ( strcmp(argv[i],"-model_type") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 	// model type
		// if ( strcmp(argv[i],"-old_benders") == 0 ) { inst->old_benders = atoi(argv[++i]); continue; } 	// old benders
		// if ( strcmp(argv[i],"-model") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 			// model type
		if ( strcmp(argv[i],"-seed") == 0 ) { inst->randomseed = abs(atoi(argv[++i])); continue; } 		// random seed
		// if ( strcmp(argv[i],"-threads") == 0 ) { inst->num_threads = atoi(argv[++i]); continue; } 		// n. threads
		if ( strcmp(argv[i],"-memory") == 0 ) { inst->available_memory = atoi(argv[++i]); continue; }	// available memory (in MB)
		if ( strcmp(argv[i], "-verbose")==0) { inst->verbose = atoi(argv[++i]); continue;}				// verbosity
		if ( strcmp(argv[i], "-v") ==0 ) {inst->verbose = atoi(argv[++i]);continue;}					// verbosity
		if (strcmp(argv[i], "-csv_column_name") == 0) { strcpy(inst->csv_column_name, argv[++i]); printf("succesful columnname\n"); continue; }
		//if ( strcmp(argv[i],"-node_file") == 0 ) { strcpy(inst->node_file,argv[++i]); continue; }		// cplex's node file
		// if ( strcmp(argv[i],"-max_nodes") == 0 ) { inst->max_nodes = atoi(argv[++i]); continue; } 		// max n. of nodes
		// if ( strcmp(argv[i],"-cutoff") == 0 ) { inst->cutoff = atof(argv[++i]); continue; }				// master cutoff
		// if ( strcmp(argv[i],"-int") == 0 ) { inst->integer_costs = 1; continue; } 						// inteher costs
		if (strcmp(argv[i], "-test_bed_size") == 0) { break; }//ignore the parameter }
		if ( strcmp(argv[i],"-help") == 0 ) { help = 1; continue; } 									// help
		if ( strcmp(argv[i],"--help") == 0 ) { help = 1; continue; } 									// help
		help = 1;
    }      
	
	if ( help ) exit(1);
}

void print_instance_parameters(instance* inst){
    printf("-------- Selected input parameters: -------\n");
	printf("File: %s\n", inst->input_file);
	if(inst->timelimit > 0){ printf("Time Limit: %lf\n", inst->timelimit); }	
	else{ printf("Time Limit: Unsetted\n"); }
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
void free_instance(instance *inst){
    //TODO: free memory accroding to how instance is allocated
    // free(inst->demand);
    free(inst->nodes);
	if (inst->is_best_sol_avail) {
		free(inst->best_sol);
		inst->is_best_sol_avail = 0;
	}
	free(inst->dist_matrix);
}

double compute_path_length(int* path, int nodes_number, point* nodes){
	int a, b;
	
	double path_length = 0;
	for (int i = 0; i < nodes_number - 1; i++){
		a = path[i];
		b = path[i+1];
		path_length += (float) get_distance(&nodes[a],&nodes[b]);
	}
	//last edge
	a = path[0];

	return (path_length + (float) get_distance(&nodes[a],&nodes[b]));
}

/**
 * Reverses a sequence within limits min and max, not included, in the array path. 
 * IP: array of integers path, limits to the segment to reverse ( min and max )
 * OP: none
*/
void reverse_sequence(int* path, int min, int max){

	int* s = &path[min+1];
	int* t = &path[max-1];
	while (s<t){

		swap(s,t);
		s++;
		t--;

	}

}

/**
	Given a solution (i.e. a path), it swaps the positions of two nodes if an improvement is found.
	IP: instance inst passed by reference
	OP: formally none, but inst's best_sol is modified
*/
void swap_2_opt(int* path, int i, int j){

	//swapping the extremes
	int temp = path[i];
	path[i] = path[j];
	path[j] = temp;
	int min = i;
	int max = j;
	if(i>j){
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
			double current_dist= (double)get_dist_matrix((const float*) (dist_matrix),incumbent_sol[i], incumbent_sol[(i+1)]) + (double)get_dist_matrix((const float*) (dist_matrix),incumbent_sol[j], incumbent_sol[(j+1)%n]);
			double changed_dist = (double)get_dist_matrix((const float*) (dist_matrix),incumbent_sol[i], incumbent_sol[(j)])+ (double)get_dist_matrix((const float*) (dist_matrix),incumbent_sol[i+1], incumbent_sol[(j+1)%n]);

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
		tsp_debug(0, 0, "I am swapping nodes %d and %d", best_i+1,best_j);
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

//NINO: I think it is unuseful and it can be erased, it is sufficient to use classic opt2 passing as parameter inst->best_sol
/**
 * 	Given an instance inst, it performs opt-2 refinement.
 * IP: instance inst passed by reference
 * OP: formally none, but inst's best_sol and zbest are updated if an improvement is found
*/
void opt2_optimize_best_sol(instance* inst) {

	int n = inst->nnodes;
	int* path = inst->best_sol;
	point* nodes_list = inst->nodes;
	double path_length = inst->zbest;
	tsp_debug((inst->verbose > 49), 0,"Initial real path length: %lf", compute_path_length(path,n,nodes_list));
	tsp_debug((inst->verbose > 49), 0,"initial zbest : %lf", path_length);
	double best_delta = -1;
	
    while (best_delta < 0 && !(is_time_limit_exceeded(inst->timelimit)))
    {	
		char improvement = 0;		
        best_delta = 0;
        int best_i = -1;
        int best_j = -1;
        for (int i = 0; i <= n -3; i++) //cambiato da n-2!
        {
            for (int j = i + 2; j <= n-1; j++)
            {
                double current_dist= get_distance(&nodes_list[inst->best_sol[i%n]], &nodes_list[inst->best_sol[(i+1)%n]]) + get_distance(&nodes_list[inst->best_sol[j%n]], &nodes_list[inst->best_sol[(j+1)%n]]);
                double changed_dist = get_distance(&nodes_list[inst->best_sol[i%n]], &nodes_list[inst->best_sol[j%n]]) + get_distance(&nodes_list[inst->best_sol[(i+1)%n]], &nodes_list[inst->best_sol[(j+1)%n]]);
                double delta = changed_dist - current_dist;

                if (delta < best_delta)
                {
					improvement = 1;
                    best_i = i%n;
                    best_j = j%n;
                    best_delta = delta;
                }
            }
        }
        if (improvement)
        {
			tsp_debug((inst->verbose > 49), 0, "I am swapping nodes %d and %d", best_i+1,best_j);
            swap_2_opt(inst->best_sol, (best_i + 1)%n, (best_j)%n);
            path_length += best_delta;
			inst->tbest = get_timer();
        }
    }
	if(is_time_limit_exceeded(inst->timelimit)){

		tsp_debug((inst->verbose>=5),0, "Could not finish 2 opt before time limit: current time is %lf",get_timer());
	}
	if(is_feasible_solution(inst,path,path_length)){

		printf("Solution is feasible!\n");

	}
	inst->zbest = path_length;

}
//TO ERASE
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
//IT IS USED ONLY IN VNS, IM NOT SURE IT IS RESILIENT FROM ERRORS, BUT I THINK IT IS OK
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
//I THINK IT IS OK
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
void copy_din_array(void *a1, const void *a2, size_t elem_size, size_t num_elems) {
    if (a1 == NULL || a2 == NULL) {
        fprintf(stderr, "Error: null pointers.\n");
        exit(main_error(-9));
    }
    size_t total_size = elem_size * num_elems;
    memcpy(a1, a2, total_size);
}
//TO ERASE, WE HAVE IMPLEMENTED THE NEW VERSION IN greedy.c
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
int* search_min(const int* p, const int* end, const float* dist_matrix, double* min){
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
//TO ERASE, WE HAVE IMPLEMENTED THE NEW VERSION IN greedy.c
/*
* Computes a greedy path starting from a specified index and calculates its cost.
* This function constructs a greedy path starting from the node at the given index.
* It iteratively selects the closest node to the last visited node, swaps it
* into the path, and calculates the cost of the path.
* 
* IP index_first Index of the starting node for the greedy path.
* IP inst Pointer to the instance containing node information.
* OP path_cost Pointer to a variable to store the computed path cost.
* OR A pointer to the constructed greedy path.
*/
int* compute_greedy_path(int index_first, instance* inst, double* path_cost) {
	int n = inst->nnodes;
	int* path =(int*) calloc(n , sizeof(int)); //manca un assert, non sarebbe meglio metterlo dentro init_path?
	init_path(path, n);
	swap(&(path[0]),&(path[index_first]));
	int* next = &(path[1]); //pointer to the next node to visit in the path (at the start path[0] is already correct)
	int* end = &(path[n - 1]); //pointer to the ending node of the path
	int* min; //pointer to the closest node of the last node visited;
	double current_cost = 0;
	double aggregate_cost = 0;
	while (next < end) {
		min = search_min((next - 1), end, (const float*)(inst->dist_matrix), &current_cost);
		swap(next, min);
		next++;
		aggregate_cost += current_cost;
	}
	current_cost = get_dist_matrix((const float*)(inst->dist_matrix), *(end-1), *end);
	aggregate_cost += (current_cost + get_dist_matrix((const float*)(inst->dist_matrix), path[0], *end));
	(*path_cost) = aggregate_cost;
	return path;
}
//TO ERASE, WE HAVE IMPLEMENTED THE NEW VERSION IN greedy.c
/*
* This function must update the data of $inst related to the best solution using
* a greedy approach that provides a good solution for the tsp problem.
* The idea is to build a path in which the initial node 0 (position 0 in the vector representing the path) 
* is arbitrarily chosen among the nodes of the graph, 
* while node i is chosen by selecting the node at minimum distance from node i-1.
* INVARIANT FOR THE PATH BUILDED IN GREEDY TSP.
* Let PATH = [v1, v2, ..., vK] the sequence of visited node we have that
* for each (1 <= i <= k-1) : d(v_i, v_i+1) <= d(v_i, u) for each u not yet visited.
* IOP inst Pointer to the instance containing node information. The following field are updated:
*	tbest Execution time to find the best solution.
*	zbest Cost of the best solution.
*	best_sol Sequence of indices representing the best solution.
*/
void greedy_tsp(instance *inst){
	int n = inst->nnodes;
	int verbose = inst->verbose;
	double current_cost;
	double min_cost;
	int* sol;
	int* min_sol = compute_greedy_path(0, inst, &min_cost);
	update_best(inst, min_cost, get_timer(), min_sol);
	int i;
	for (i = 1; i < n; i++) {
		sol = compute_greedy_path(i,inst,&current_cost);
		if (current_cost < min_cost) {
			free(min_sol); //I am going to update min_sol, then I can free the actual min_sol
			min_cost = current_cost; //update cost
			min_sol = sol; //update min_sol
			update_best(inst, min_cost, get_timer(), min_sol);
		}
		else{
			free(sol); //sol is not the best path, then i can free it
		}
		tsp_debug((inst->verbose>99), 0, "iter %d ", i);
		tsp_debug((inst->verbose>99), 0, "cost = %.2f", current_cost);
		tsp_debug((inst->verbose>99), 0, "zbest = %.2f", min_cost);
	}


	//If we exited the cycle due to time constraints, print an error message saying the reached point
	if(i<n){
		tsp_debug(inst->verbose,0,"Could not finish greedy NN due to time constraints: visited up until node %d", i);
	}
}

/**
 * Updates the best solution information in the instance.
 * OP inst Pointer to the instance structure to modify.
 * IP z New best objective value.
 * IP t New best execution time.
 * IP sol Pointer to the best solution (array of integers).
 */
void update_best(instance* inst, double z, double t, int* sol){
	char flag = inst->verbose >= 1;
	tsp_debug_inline(flag, "\n -------------- Update Best : --------------\n");
	tsp_debug_inline(flag, "Old Best : %.20g\n", inst->zbest);
	tsp_debug_inline(flag, "New Best : %.20g\n", z);
	tsp_debug_inline(flag, "\n --------------- End Update : --------------\n");
	inst->zbest = z;
	inst->tbest = t - inst->tstart;
	inst->is_best_sol_avail = 1;
	inst->best_sol = sol;
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
void init_path(int* path, size_t n){
	for(int i=0; i<n; i++){
		(*path) = i;
		path++;
	}
}
void init_data_file(char flag, FILE* data_file, instance* inst){
    char figure_name[64];
    if(flag){
         generate_name(figure_name, sizeof(figure_name), "data/tabu_%d_%d_cost.dat", inst->nnodes, inst->randomseed);
         data_file = fopen(figure_name, "w");
         if(data_file == NULL){
            fclose(data_file);
            exit(main_error_text(-2,"Failed to open the file %s", figure_name));
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
char is_feasible_solution(instance* inst, int* sol_path, double sol_cost){

	int n = inst->nnodes;
	int* counter = (int*) calloc(n, sizeof(int));
	//print_path(1, "path: \n", sol_path, inst->nodes,n);
	//check for repeated vertices or out of range content of path
	for (int i =0; i<=n-1; i++){

		counter[sol_path[i]]++;

	}
	for(int i= 0; i<=n-1; i++){

		if(counter[sol_path[i]]!=1){
			tsp_debug((inst->verbose >= 200), 0, "Node %d is visited %d", i, counter[sol_path[i]]);
			free(counter);
			return 0;

		}
	}

	free(counter);

	//check cost coincides with the path length 

	double path_length = compute_path_length(sol_path, n, inst->nodes);
	tsp_debug((inst->verbose >= 200), 0, "Real Path length : %.6f", path_length);
	tsp_debug((inst->verbose >= 200), 0, "Sol passed length : %.6f",sol_cost);
	tsp_debug((inst->verbose >= 200), 0, "is_feasible = %d", is_equal_double(path_length, sol_cost, EPSILON));
	return is_equal_double(path_length, sol_cost, EPSILON);
}

/**
 * Utility method to return a solver given a command line input
 * IP: string from command line
 * OP: appropriate solver_id, NOT_DEF if incorrectly specified
*/

solver_id parse_solver(char* solver_input){
	if ( strcmp(solver_input, "greedy") == 0 || strcmp(solver_input, "nn") == 0){
		return NN;
	}
	else if ( strcmp(solver_input, "random") == 0 || strcmp(solver_input, "rand") == 0 ){
		return RANDOM_SOL;
	}
	else if ( strcmp(solver_input, "2opt") == 0|| strcmp(solver_input, "opt2") == 0 ){
		return OPT_2;
	}
	else if ( strcmp(solver_input, "tabu") == 0){
		return TABU;
	}
	else if( strcmp(solver_input, "vns")==0){
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
	else if (strcmp(solver_input, "bcf") == 0) {

		return BCF;
	}
	else if (strcmp(solver_input, "bcm") == 0) {
		return BCM;
	}
	else if (strcmp(solver_input, "bcfm") == 0) {
		return BCFM;
	}
	else if (strcmp(solver_input, "bc") == 0) {
		return BC;
	}
	exit(main_error_text(-8, "%s","solver"));
}


void tsp_solve(instance* inst){
	inst->tstart = get_timer();
	inst->timelimit += get_timer();
	MKDIR("figures");
	MKDIR("data");
	char figure_name[64];
	int check_truncation = 0;
	switch (inst->solver) {
	case NN:
		gre_solve(inst, 0);
		print_best_sol((inst->verbose>=1), inst);
		generate_name(figure_name, sizeof(figure_name), "figures/greedy_%d_%d.png", inst->nnodes, inst->randomseed);
		plot_path(inst->verbose >= 1000, inst->best_sol, inst->nnodes, inst->zbest, inst->nodes, figure_name);
		//init_data_file((inst->verbose>-1),(inst->best_sol_data), inst);
		//plot_generator(inst);
		break;
	case OPT_2:
		gre_solve(inst, 1);
		generate_name(figure_name, sizeof(figure_name), "figures/greedy+opt2_%d_%d.png", inst->nnodes, inst->randomseed);
		plot_path(inst->verbose >= 1000, inst->best_sol, inst->nnodes, inst->zbest, inst->nodes, figure_name);
		print_best_sol((inst->verbose>=1), inst);
		break;
	case TABU:
		tabu_search(1,inst);
		print_best_sol((inst->verbose>=1), inst);
		generate_name(figure_name, sizeof(figure_name), "figures/tabu_%d_%d.png", inst->nnodes, inst->randomseed);
		plot_path(inst->verbose >= 1, inst->best_sol, inst->nnodes, inst->zbest, inst->nodes, figure_name);
		//init_data_file((inst->verbose>-1),(inst->best_sol_data), inst);
		break;
	case VNS:
		printf("Initializing a greedy solution\n");
		//WARNING greedy_tsp will be erased
    	greedy_tsp(inst); //to change and use gre_solve instead!!
		//WARNING
		print_best_sol((inst->verbose>=5), inst);
		generate_name(figure_name, sizeof(figure_name), "figures/greedy_%d_%d.png", inst->nnodes, inst->randomseed);
		plot_path(inst->verbose >= 1, inst->best_sol, inst->nnodes, inst->zbest, inst->nodes, figure_name);
		init_data_file((inst->verbose>-1),(inst->best_sol_data), inst);
		vns(inst);
		print_best_sol((inst->verbose>=5), inst);
		generate_name(figure_name, sizeof(figure_name), "figures/vns_%d_%d.png", inst->nnodes, inst->randomseed);
		plot_path(inst->verbose >= 1, inst->best_sol, inst->nnodes, inst->zbest, inst->nodes, figure_name);
		//init_data_file((inst->verbose>-1),(inst->best_sol_data), inst);
		break;
	case EX:
		TSPopt(inst);
		break;
	case BEN:
		ben_solve(0, inst);
		break;
	case GLU:
		ben_solve(1, inst);
		break;
	case BC:
		cpx_branch_and_cut(0, 0, 0, inst);
		break;
	case BCM:
		cpx_branch_and_cut(0, 1, 0, inst);
		break;
	case BCP:
		cpx_branch_and_cut(0, 0, 1, inst);
		break;
	case BCMP:
		cpx_branch_and_cut(0, 1, 1, inst);
		break;
	case BCF:
		cpx_branch_and_cut(1,0,0,inst);
		break;
	case BCFM:
		cpx_branch_and_cut(1,1,0, inst);
		break;
	case BCFP:
		cpx_branch_and_cut(1, 0, 1, inst);
		break;
	case BCFMP:
		cpx_branch_and_cut(1, 1, 1, inst);
		break;
	case HF:
		hf_solve(inst);
		print_best_sol((inst->verbose >= 0), inst);
		generate_name(figure_name, sizeof(figure_name), "figures/hf_%d_%d.png", inst->nnodes, inst->randomseed);
		plot_path(inst->verbose >= 1, inst->best_sol, inst->nnodes, inst->zbest, inst->nodes, figure_name);
		break;
	default:
		exit(main_error(-7));
	}
}

void update_solver(instance* inst){

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
	printf("---------------------------------------------\n");
	fgets(buf, 8, stdin);
    selection = atoi(buf);
	
	switch(selection){
		case 0:
			{inst->solver = NN;
			
			printf("successful update.\n");
			break;}
		case 1:
			{inst->solver = OPT_2;
			printf("successful update.\n");
			break;}
		case 2:
			{inst->solver = TABU;
			printf("successful update.\n");
			break;}
		case 3:
			{inst->solver = VNS;
			printf("successful update. \n");
			break;}
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
		default:
			{inst->solver = NN;
			printf("successful update.\n");
			break;}
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
	printf("heree\n");
	srand(test_bed[0].randomseed);

	for (int i = 0; i < size_test_bed; i++) {
		generate_instance(&test_bed[i]);
	}
}
void read_test_bed_size(int* test_bed_size, int argc, char** argv) {
	
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-test_bed_size") == 0) {
			*(test_bed_size) =abs(atoi(argv[++i]));
			return;
		}
	}

}
void generate_csv_file(int size_test_bed, instance* test_bed) {

	int n = test_bed[0].nnodes;
	int seed = test_bed[0].randomseed;

	MKDIR("runs");
	char csv_name[128];
	generate_name(csv_name, sizeof(csv_name), "runs/%s_%d_%d_cost.csv", test_bed[0].csv_column_name, n, seed);
	
	//Print a single-column .csv file
	FILE* file_handler = fopen(csv_name, "w");
	//Print the first row
	fprintf(file_handler, "1, %s\n", test_bed[0].csv_column_name);

	//Print the remaining rows
	for (int i = 0; i < size_test_bed; i++) {
		fprintf(file_handler, "%d_%d_[%d], %lf\n", n, seed, i, test_bed[i].zbest);
	}
	//
	fclose(file_handler);

	generate_name(csv_name, sizeof(csv_name), "runs/%s_%d_%d_time.csv", test_bed[0].csv_column_name, n, seed);

	//Print a single-column .csv file
	file_handler = fopen(csv_name, "w");
	//Print the first row
	fprintf(file_handler, "1, %s\n", test_bed[0].csv_column_name);

	//Print the remaining rows
	for (int i = 0; i < size_test_bed; i++) {
		fprintf(file_handler, "%d_%d_[%d], %lf\n", n, seed, i, test_bed[i].tbest);
	}
	//
	fclose(file_handler);

}
