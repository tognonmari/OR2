#ifndef TSP_H_

#define TSP_H_
//parameters for the random values generation
#define SEED 14
#define MAX_X 10000  //maximum value for ascisse of generated points
#define MAX_Y 10000  //maximum value for ordinate of generated points

typedef struct {
	double x;
	double y;
} point;

typedef struct {

	//input data
	int nnodes;
	double* demand;
	point* nodes;
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

	//global data
	double	tstart;
	double zbest;							// best sol. available  
	double tbest;							// time for the best sol. available  
	double* best_sol;						// best sol. available    
	double	best_lb;						// best lower bound available  
	double* load_min;						// minimum load when leaving a node
	double* load_max;						// maximum load when leaving a node

	// model;     
	int xstart;
	int qstart;
	int bigqstart;
	int sstart;
	int bigsstart;
	int ystart;
	int fstart;
	int zstart;
} instance;

#endif   /* TSP_H_ */ 