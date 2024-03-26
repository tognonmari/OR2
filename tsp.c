#include "tsp.h"
#include "utils.h"
#include "tabu.h"
#include "vns.h"
#include <math.h>
#include <malloc.h>
#include <string.h>

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
	tsp_debug_inline(flag, "\n -------------- Best solution: --------------\n");
	tsp_debug_inline(flag, "Best Distance (zbest): %.2lf\n",inst->zbest);
	tsp_debug_inline(flag, "Best Solution found at time: %12.6f\n",inst->tbest);
	tsp_debug_inline(flag && ( inst->verbose > 49 ), "The following sequence of nodes is the best solution to the problem:");
	print_path(flag && ( inst->verbose > 49 ),"", inst->best_sol, inst->nodes, inst->nnodes);
	tsp_debug_inline(flag, "\n ------------ End best solution: ------------\n");
}
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
//ACTUALLY NOT USED
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
void generate_figure_name(char buffer[], size_t bufferSize, const char *format, ...) {
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
	
	srand(inst->randomseed);
	depolarize_pseudornd_seq();

	//only generate random instance if file name is not provided, so only if no tsplib inst is specified
	
	inst->nodes = (point*)malloc(inst->nnodes * sizeof(point));
	generate_nodes(inst->nnodes, inst->nodes, MAX_X, MAX_Y);
	
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
/*
* Plot the input datafile using GNUPLOT
* IP graph_data Path of the data file to be plotted
* IP graph Path of the figure
* @note ALWAYS CHECK THAT THE DATA PASSED BY PARAMETER IS CLOSED,
* 	passing a file open for the writing can be problematic and make the plot a blank plot.
*/
int plot_graph(const char graph_data[], const char graph[]) {
	FILE* gnuplotPipe = popen("gnuplot -persist", "w");

	if (gnuplotPipe == NULL) {
		fclose(gnuplotPipe);
		exit(main_error_text(-1,"Failed to open the pipeline to gnuplot"));
	}

	//Send GNUPLOT commands through the pipe
	fprintf(gnuplotPipe, "set output '%s'\n", graph); //set output name
	fprintf(gnuplotPipe, "set terminal png\n"); //set extension
	fprintf(gnuplotPipe, "set title 'Graph'\n"); 
	fprintf(gnuplotPipe, "set key off\n"); //if key on the name of the file is printed on the plot
	fprintf(gnuplotPipe, "set xrange [0:%d]\n",MAX_X);
	fprintf(gnuplotPipe, "set yrange [0:%d]\n",MAX_Y);
	fprintf(gnuplotPipe, "set pointsize 0.5\n");
	fprintf(gnuplotPipe, "plot '%s'\n", graph_data);

	fclose(gnuplotPipe);
}
int plot_path(char flag, const char figure_name[], const int* indices, const point* points, int num_points) {
    FILE* gnuplotPipe = popen("gnuplot -persist", "w");

    if (gnuplotPipe == NULL) {
        fclose(gnuplotPipe);
		exit(main_error_text(-1,"Failed to open the pipeline to gnuplot"));
    }

    fprintf(gnuplotPipe, "set output '%s'\n", figure_name); //set output name
	fprintf(gnuplotPipe, "set terminal png\n"); //set extension
	fprintf(gnuplotPipe, "set title 'Greedy Path'\n"); 
	fprintf(gnuplotPipe, "set key off\n"); //if key on the name of the file is printed on the plot
	fprintf(gnuplotPipe, "set xrange [0:%d]\n",MAX_X);
	fprintf(gnuplotPipe, "set yrange [0:%d]\n",MAX_Y);
	fprintf(gnuplotPipe, "set pointsize 0.5\n");
	fprintf(gnuplotPipe, "set grid \n");
	fprintf(gnuplotPipe, "plot '-' with linespoints pt 1 lc rgb '#800080'\n");

    for (int i = 0; i < num_points; ++i) {
        fprintf(gnuplotPipe, "%lf %lf\n", points[indices[i]].x, points[indices[i]].y);
    }
	fprintf(gnuplotPipe, "%lf %lf\n", points[indices[0]].x, points[indices[0]].y);
    fprintf(gnuplotPipe, "e\n"); //This line serves to terminate the input of data for the plot command in gnuplot.

    fclose(gnuplotPipe);
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
	// inst->integer_costs = 0;
    inst->verbose = VERBOSE;
	inst->available_memory = 12000;   			
	// inst->max_nodes = -1; 	
	inst->nnodes =0;					       
	inst->solver = NOT_DEF;
    int help = 0; if ( argc < 1 ) help = 1;	
	for ( int i = 1; i < argc; i++ ) 
	{ 
		if ( strcmp(argv[i], "-nnodes") == 0) { inst->nnodes = atoi(argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-input") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-f") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 				// input file
		if ( strcmp(argv[i],"-time_limit") == 0 ) { inst->timelimit = atof(argv[++i]); continue; }		// total time limit
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
		//if ( strcmp(argv[i],"-node_file") == 0 ) { strcpy(inst->node_file,argv[++i]); continue; }		// cplex's node file
		// if ( strcmp(argv[i],"-max_nodes") == 0 ) { inst->max_nodes = atoi(argv[++i]); continue; } 		// max n. of nodes
		// if ( strcmp(argv[i],"-cutoff") == 0 ) { inst->cutoff = atof(argv[++i]); continue; }				// master cutoff
		// if ( strcmp(argv[i],"-int") == 0 ) { inst->integer_costs = 1; continue; } 						// inteher costs
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
	free(inst->best_sol);
	free(inst->dist_matrix);
}

double compute_path_length(int* path, int nodes_number, point* nodes){
	int a, b;
	
	double path_length = 0;
	for (int i = 0; i < nodes_number - 1; i++){
		a = path[i];
		b = path[i+1];
		path_length += get_distance(&nodes[a],&nodes[b]);
	}
	//last edge
	a = path[0];

	return (path_length + get_distance(&nodes[a],&nodes[b]));
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
	
    while (best_delta < 0)
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
        }
    }

	if(is_feasible_solution(inst,path,path_length)){

		printf("Solution is feasible!\n");

	}
	inst->zbest = path_length;

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
void copy_din_array(void *a1, const void *a2, size_t elem_size, size_t num_elems) {
    if (a1 == NULL || a2 == NULL) {
        fprintf(stderr, "Error: null pointers.\n");
        exit(main_error(-9));
    }
    size_t total_size = elem_size * num_elems;
    memcpy(a1, a2, total_size);
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
	for (int i = 1; i < n; i++) {
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
}

/**
 * Updates the best solution information in the instance.
 *
 * OP inst Pointer to the instance structure to modify.
 * IP z New best objective value.
 * IP t New best execution time.
 * IP sol Pointer to the best solution (array of integers).
 */
void update_best(instance* inst, double z, double t, int* sol){
	inst->zbest = z;
	inst->tbest = t;
	inst->best_sol = sol;
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

			free(counter);
			return 0;

		}
	}

	free(counter);

	//check cost coincides with the path length 

	double path_length = compute_path_length(sol_path, n, inst->nodes);
	//tsp_debug((inst->verbose > 49), 0, "Real Path length : %.2f", path_length);
	//tsp_debug((inst->verbose > 49), 0, "Sol passed length : %.2f",sol_cost);
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
	exit(main_error_text(-8, "%s","solver"));
}


void tsp_solve(instance* inst){
	get_timer();
	MKDIR("figures");
	char figure_name[64];
	int check_truncation;
	switch (inst->solver) {
	case NN:
		greedy_tsp(inst);
		print_best_sol((inst->verbose>-1), inst);
		generate_figure_name(figure_name, sizeof(figure_name), "figures/greedy_%d_%d.png", inst->nnodes, inst->randomseed);
		plot_path((inst->verbose>-1),figure_name, inst->best_sol, inst->nodes, inst->nnodes);
		break;
	case OPT_2:
		greedy_tsp(inst);
		print_best_sol((inst->verbose>=5), inst);
		generate_figure_name(figure_name, sizeof(figure_name), "figures/greedy_%d_%d.png", inst->nnodes, inst->randomseed);
		plot_path((inst->verbose>-1),figure_name, inst->best_sol, inst->nodes, inst->nnodes);
		//optimize 
		opt2_optimize_best_sol(inst);
		generate_figure_name(figure_name, sizeof(figure_name), "figures/greedy_%d_%d_opt2.png", inst->nnodes, inst->randomseed);
		plot_path((inst->verbose>-1),figure_name, inst->best_sol, inst->nodes, inst->nnodes);
		print_best_sol((inst->verbose>=5), inst);
		break;
	case TABU:
		tabu_search(inst);
		print_best_sol((inst->verbose>=1), inst);
		generate_figure_name(figure_name, sizeof(figure_name), "figures/tabu_%d_%d.png", inst->nnodes, inst->randomseed);
		plot_path((inst->verbose>-1),figure_name, inst->best_sol, inst->nodes, inst->nnodes);
		break;
	case VNS:
		printf("Initializing a greedy solution\n");
    	greedy_tsp(inst);
		print_best_sol((inst->verbose>=5), inst);
		generate_figure_name(figure_name, sizeof(figure_name), "figures/greedy_%d_%d.png", inst->nnodes, inst->randomseed);
		plot_path((inst->verbose>-1),figure_name, inst->best_sol, inst->nodes, inst->nnodes);
		vns(inst);
		print_best_sol((inst->verbose>=5), inst);
		generate_figure_name(figure_name, sizeof(figure_name), "figures/vns_%d_%d.png", inst->nnodes, inst->randomseed);
		plot_path((inst->verbose>-1),figure_name, inst->best_sol, inst->nodes, inst->nnodes);
	default:
		exit(-7);
	}
}

void update_solver(instance* inst){

	int selection;
	printf("----Choose a solver, from these options:----\n");
	printf("0: Nearest Neighbor\n");
	printf("1: Nearest Neighbor and OPT2\n");
	printf("2: Nearest Neighbor and TABU search\n");
	printf("3: Nearest Neighbor and VNS\n");
	printf("---------------------------------------------\n");
	scanf("%d",&selection);
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
			printf("successful update.\n");
			break;}
		default:
			{inst->solver = NN;
			printf("successful update.\n");
			break;}
	}


}