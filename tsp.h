#ifndef TSP_H
#define TSP_H
#include<stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
//parameters for the random values generation
#define MAX_X 1000  //maximum value for ascisse of generated points
#define MAX_Y 1000//maximum value for ordinate of generated points

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
	int* best_sol;							// best sol. available    
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

int plot_path(char flag, const char figure_name[], const int* indices, const point* points, int num_points);

void make_datafile(instance *inst, FILE* data_file);

void generate_instance(instance *inst);

void generate_nodes(int n, point* nodes, int max_x, int max_y);

void generate_array(int n, double *array, int max_value);

double get_distance(const point* p1,const point* p2);

void print_nodes(char flag, const char text_to_print[], const instance *inst, int n);

void print_path(char flag, const char text_to_print[], const int* path, const point* nodes, int n);

void print_best_sol(char flag, instance* inst);

void print_point(char flag, const char text_to_print[], const point* p);

void print_triangular_matrix(char flag,const char text[], const double** matrix, int nrows);

void tsp_debug(char flag,int flag_time, char* format, ...);

void tsp_debug_inline(char flag, char* format, ...);

void free_matrix(void** matrix, int rows);

void free_instance(instance *inst);

void print_instance_parameters(instance inst);

void parse_command_line(int argc, char** argv, instance *inst);

void opt2_optimize_best_sol(instance *inst);

double compute_path_length(point* path, int nodes_number);

void compute_cost_matrix(instance* inst);

double get_cost_matrix(const double** matrix, int a, int b);

void swap_space(void* a, void* b, size_t size);

void swap(int* a, int* b);

void copy_array(void* a1, const void* a2);

int* search_min(const int* p, const int* end, const double** cost_matrix, double* current_cost);

int* compute_greedy_path(int index_first, instance* inst, double* path_cost);

void greedy_tsp(char flag, instance* inst);

void update_best(instance* inst, double z, double t, int* sol);

void init_path(int* path, size_t n);
#endif