//Definition of Data Structures: for now we only keep in the struct information about the single problem intance and the input parameters, 
//since we are not solving anything yet.

typedef struct{
    double x;
    double y;
} point;

typedef struct {

    //input data
	int nnodes; 	
	double *demand;   
    point *nodes;
	int depot;
	double capacity; 
	int nveh;

	// parameters 
	int model_type; 
	int old_benders;
	int randomseed;
	int num_threads;
	double timelimit;						// overall time limit, in sec.s
	char input_file[1000];		  			// input file
	char node_file[1000];		  			// cplex node file
	int available_memory;
	int max_nodes; 							// max n. of branching nodes in the final run (-1 unlimited)
	double cutoff; 							// cutoff (upper bound) for master
	int integer_costs;
	int verbose;


} instance;