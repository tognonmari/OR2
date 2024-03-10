#include "tsp1.h"
#include "utils.h"
#include <math.h>
/* 
* Conditional debugging function that prints
* a formatted message to the console if the given $flag is true,
* the time of the execution is printed together with the message
* if the $time_flag is true.
*/
void tsp_debug(int flag,int time_flag, char* format, ...)
{
	if (flag) {
		va_list args;
		va_start(args, format);
		if (time_flag) { printf("\n%12.6f|", get_timer()); }
		vprintf(format, args);
		va_end(args);
		return;
	}
}
/*
*/
void print_triangular_matrix(const double** matrix, int nrows) {
	for (int row = 0; row < nrows; row++) {
		for (int col = 0; col <= row; col++) {
			tsp_debug(1,0,"%.3f ", matrix[row][col]);
		}
		tsp_debug(1,0,"\n");
	}
}
/* Print the coordinates of a point
* IP text_to_print
* IP p Point to print
*/
void print_point(const char text_to_print[], const point* p) {
	int x = (int)(p->x);
	int y = (int)(p->y);
	tsp_debug(1,1,"%s", text_to_print);
	tsp_debug(1,1,"(%d, %d)\n", x, y);
}
/* Print the coordinates of the nodes of the graph
* IP text_to_print.
* IP inst Graph to print
* IP n Size of the graph
*/
void print_nodes(const char text_to_print[], const instance *inst, int n) {
	int i;
	tsp_debug(1,1,"%s", text_to_print);
	point* current_point = inst->nodes;
	for (i = 0; i < n; i++) {
		tsp_debug(1,1,"node[%d]: ", i);
		print_point("",current_point);
		current_point++;
	}
}
/*
* The values of an array are transformed in random values belonging to [0,$max_value]
* IP n Length of the array
* OP array Array to be modified with random values
* IP max_value Maximum value of the generated values
*/
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

/*
* A complete graph is implicited defined by its nodes -> this function generates $n nodes
* in a random way, thus each generated pair of coordinates (x_i,y_i) is s.t.:
* x_i belongs to [0,$max_x]
* x_i belongs to [0,$max_y]
* IP n Number of nodes to generate
* IP seed Seed of the pseudorandom sequence
* IOP inst Instance to be initialized with random nodes
*/
void generate_instance(instance *inst){
	inst->nodes = (point*)malloc(inst->nnodes * sizeof(point));
	srand(inst->randomseed);
	depolarize_pseudornd_seq();
	generate_nodes(inst->nnodes, inst->nodes, MAX_X, MAX_Y);
}
/*
 * Calculates the Euclidean distance between two points in a two-dimensional space.
 * This function uses the Euclidean distance formula:
 * distance = sqrt((x2 - x1)^2 + (y2 - y1)^2)
 * IP p1 Point 1.
 * IP p2 Point 2.
 * OR Euclidean distance between $p1 and $p2
 */
double get_distance(point* p1, point* p2) {
	double x1, x2, y1, y2;
	x1 = p1->x; x2 = p2->x; y1 = p1->y; y2 = p2->y;
	double dist = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
	return dist;
}
/*
* Compute the cost matrix, thus the matrix containing in each entry (i,j) the
* distance between node i and node j.
* IOP inst The $inst field $cost_matrix is computed with the distance of the nodes.
*/
void compute_cost_matrix(instance* inst) {
	int nrow = inst->nnodes;
	double** matrix = (double**)alloc_triangular_matrix(nrow, sizeof(double));
	for (int i = 0; i < nrow; i++ ) {
		for (int j = 0; j <= i; j++) {
			matrix[i][j] = get_distance(&(inst->nodes[i]),&(inst->nodes[j]));
		}
	}
	inst->cost_matrix = matrix;
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
		printf("ciao");
		current_point++;
	}
}
/*
* Plot the input datafile using GNUPLOT
* IP graph_data Path of the data file to be plotted
* IP graph Path of the figure 
*/
int plot_graph(const char graph_data[], const char graph[]) {
	FILE* gnuplotPipe = popen("gnuplot -persistent", "w");

	if (gnuplotPipe == NULL) {
		fprintf(stderr, "Error starting GNUPLOT.\n");
		exit(-1);
	}

	// Send GNUPLOT commands through the pipe
	fprintf(gnuplotPipe, "set terminal png\n");
	fprintf(gnuplotPipe, "set output '%s'\n", graph);
	fprintf(gnuplotPipe, "set key off\n");
	fprintf(gnuplotPipe, "set title 'Graph'\n");
	fprintf(gnuplotPipe, "set pointsize 0.5\n");
	fprintf(gnuplotPipe, "plot '%s' with points pointtype 7\n", graph_data);

	fclose(gnuplotPipe);
}

void parse_command_line(int argc, char** argv, instance *inst){

    //Default parameters

    inst->model_type = 0;
	inst->old_benders = 0;
	strcpy(inst->input_file, "NULL");
	inst->randomseed = 0; 
	inst->num_threads = 0;
	inst->timelimit = -1; 
	inst->cutoff = -1; 
	inst->integer_costs = 0;
    inst->verbose = 10;
	inst->available_memory = 12000;   			
	inst->max_nodes = -1; 	
	inst->nnodes =0;					       

    int help = 0; if ( argc < 1 ) help = 1;	
	for ( int i = 1; i < argc; i++ ) 
	{ 
		if ( strcmp(argv[i], "-nnodes") == 0) { inst->nnodes = atoi(argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-input") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-f") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 				// input file
		if ( strcmp(argv[i],"-time_limit") == 0 ) { inst->timelimit = atof(argv[++i]); continue; }		// total time limit
		if ( strcmp(argv[i],"-model_type") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 	// model type
		if ( strcmp(argv[i],"-old_benders") == 0 ) { inst->old_benders = atoi(argv[++i]); continue; } 	// old benders
		if ( strcmp(argv[i],"-model") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 			// model type
		if ( strcmp(argv[i],"-seed") == 0 ) { inst->randomseed = abs(atoi(argv[++i])); continue; } 		// random seed
		if ( strcmp(argv[i],"-threads") == 0 ) { inst->num_threads = atoi(argv[++i]); continue; } 		// n. threads
		if ( strcmp(argv[i],"-memory") == 0 ) { inst->available_memory = atoi(argv[++i]); continue; }	// available memory (in MB)
		//if ( strcmp(argv[i],"-node_file") == 0 ) { strcpy(inst->node_file,argv[++i]); continue; }		// cplex's node file
		if ( strcmp(argv[i],"-max_nodes") == 0 ) { inst->max_nodes = atoi(argv[++i]); continue; } 		// max n. of nodes
		if ( strcmp(argv[i],"-cutoff") == 0 ) { inst->cutoff = atof(argv[++i]); continue; }				// master cutoff
		if ( strcmp(argv[i],"-int") == 0 ) { inst->integer_costs = 1; continue; } 						// inteher costs
		if ( strcmp(argv[i],"-help") == 0 ) { help = 1; continue; } 									// help
		if ( strcmp(argv[i],"--help") == 0 ) { help = 1; continue; } 									// help
		help = 1;
    }      
	
	if ( help ) exit(1);
}

void print_instance_parameters(instance inst){
    printf("-------- Selected input parameters: -------\n");
	printf("File: %s\n", inst.input_file); 
	printf("Time Limit: %lf\n", inst.timelimit); 
	printf("Model Type: %d\n", inst.model_type); 
	printf("Old Benders: %d\n", inst.old_benders); 
	printf("Seed: %d\n", inst.randomseed); 
	printf("Threads: %d\n", inst.num_threads);  
	printf("Max Nodes: %d\n", inst.max_nodes); 
	printf("Memory: %d\n", inst.available_memory); 
	printf("Integer Costs: %d\n", inst.integer_costs); 
	printf("Node File: %s\n", inst.node_file);
	printf("Cutoff: %lf\n", inst.cutoff); 
	printf("-------------------------------------------\n");
}

void free_instance(instance *inst){
    //TODO: free memory accroding to how instance is allocated
    free(inst->demand);
    free(inst->nodes);
	free(inst->best_sol);

}
/*
	Returns euclidean distance between two points.
	IP: two points passed by reference
	OP: the Euclidean Distance between the two points
*/
double euclidean_dist(point* p1, point* p2){
	
	double dx = p1->x -p2->x;
	double dy = p1->y-p2->y;
	return sqrt(pow(dx,2) + pow(dy,2));
}


double compute_path_length(point* path, int nodes_number){
	point p1, p2;
	
	double path_length = 0;
	for (int i = 0; i < nodes_number; i++){
		p1 = path[i];
		p2 = path[i+1];
		path_length += euclidean_dist(&p1,&p2);
	}
	//last edge
	p1 = path[0];

	return (path_length+ euclidean_dist(&p1,&p2));
}


/*
	Given a solution (i.e. a path), it swaps the positions of two nodes if an improvement is found.
	IP: instance inst passed by reference
	OP: formally none, but inst's best_sol is modified
*/
void swap_2_opt(int* path, int i, int j){

	int temp = path[i+1];
	path[i+1] = path[j];
	path[j] = temp;

}


/*
	Given an instance inst, it performs opt-2 refinement.
	IP: instance inst passed by reference
	OP: formally none, but inst's best_sol and zbest are updated if an improvement is found
*/
void opt2_optimize_best_sol(instance *inst){

	//Preliminary steps: parameters and pointers initialization

	int nodes_number = inst->nnodes;
	double path_length = inst->zbest;
	int* path = (inst->best_sol);
	point* nodes_list = inst->nodes;

	//Opt-2 algorithm

	int improvement = 1;

	while(improvement ==1){

		improvement =0;

		for (int i =0; i< nodes_number-1; i++){

			for (int j= i+1; j<nodes_number; j++){

				double delta = -euclidean_dist(&nodes_list[path[i]],&nodes_list[path[i+1]])-euclidean_dist(&nodes_list[path[j]],&nodes_list[path[j+1]])+euclidean_dist(&nodes_list[path[i]], &nodes_list[path[j]])+euclidean_dist(&nodes_list[path[i+1]],&nodes_list[path[j+1]]);
				if(delta <0){

					swap_2_opt(path, i, j); //swap operations if an improvement is found
					path_length += delta;
					improvement = 1;

				}
			}
		}
	}

	//updating new best cost

	inst->zbest = path_length;

}



