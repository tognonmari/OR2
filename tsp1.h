#ifndef TSP_H
#define TSP_H
#include<stdarg.h>
#include <stdio.h>
#include <stdlib.h>
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
	double** cost_matrix;					//cost_matrix[i][j] represent the distance between node i and node j.
	double	tstart;
	double zbest;							// best sol. available  
	double tbest;							// time for the best sol. available  
	int* best_sol;						// best sol. available    
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

int plot_graph(const char graph_data[], const char graph[]);

void make_datafile(instance *inst, FILE* data_file);

void generate_instance(instance *inst);

void generate_nodes(int n, point* nodes, int max_x, int max_y);

void generate_array(int n, double *array, int max_value);

double get_distance(point* p1, point* p2);

void print_nodes(const char text_to_print[], const instance *inst, int n);

void print_point(const char text_to_print[], const point* p);

void tsp_debug(int flag,int flag_time, char* format, ...);

void free_instance(instance *inst);

void print_instance_parameters(instance inst);

void parse_command_line(int argc, char** argv, instance *inst);

double euclidean_dist(point* p1, point* p2);

void opt2_optimize_best_sol(instance *inst);

double compute_path_length(point* path, int nodes_number);

void print_triangular_matrix(const double** matrix, int nrows);

void compute_cost_matrix(instance* inst);
#endif