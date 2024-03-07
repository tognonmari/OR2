#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#if defined(_WIN32) || defined(_WIN64)
#include <direct.h>
#define MKDIR(dir) _mkdir(dir)
#else
#include <sys/stat.h>
#define MKDIR(dir) mkdir(dir, 0755)
#endif

/*
This function generates and discards 100 pseudorandom numbers using the rand() function.
* The purpose is to mitigate polarization effects in the pseudorandom sequence.
* @note The generated random numbers are not used or returned by this function.
*/

void depolarize_pseudornd_seq() {
	for (int i = 1; i <= 100; i++) {
		rand();
	}
}

/* 
* Function to initialize a timer and get the current time from the first call of the function.
*/
double get_timer() {
	static clock_t start_time = 0;  // Memorizza il tempo di avvio in modo persistente
	if (start_time == 0) {
		start_time = clock();  // Inizializza il tempo di avvio alla prima chiamata
	}

	clock_t current_time = clock();
	double elapsed_time = (double)(current_time - start_time) / CLOCKS_PER_SEC;
	return elapsed_time;
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