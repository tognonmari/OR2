#include "tsp1.h"
#include "utils.h"
/* 
* Conditional debugging function that prints
* a formatted message to the console if the given flag is true.
*/
void tsp_debug(int flag, char* format, ...)
{
	if (flag) {
		va_list args;
		va_start(args, format);
		printf("\n%12.6f|", get_timer());
		vprintf(format, args);
		va_end(args);
		return;
	}
}

/* Print the coordinates of a point
* IP text_to_print
* IP p Point to print
*/
void print_point(const char text_to_print[], const point* p) {
	int x = (int)(p->x);
	int y = (int)(p->y);
	printf("%s", text_to_print);
	printf("(%d, %d)\n", x, y);
}
/* Print the coordinates of the nodes of the graph
* IP text_to_print.
* IP inst Graph to print
* IP n Size of the graph
*/
void print_nodes(const char text_to_print[], const instance *inst, int n) {
	int i;
	printf("%s", text_to_print);
	point* current_point = inst->nodes;
	for (i = 0; i < n; i++) {
		printf("node[%d]: ", i);
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
* OP inst Instance to be initialized with random nodes
*/
void generate_instance(instance *inst){

	inst->nodes = (point*)malloc(inst->nnodes * sizeof(point));
	srand(inst->randomseed);
	depolarize_pseudornd_seq();
	generate_nodes(inst->nnodes, inst->nodes, MAX_X, MAX_Y);
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
//Utility methods for input parsing

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

}
/*
OR Esito(
	0: elaborazione riuscita;
	-1: apertura fallita di una pipe;
	-2: apertura fallita di un file;
	-3: numero di parametri inseriti da linea di comando non corretto).
	-4: buffer truncation;
*/ 
