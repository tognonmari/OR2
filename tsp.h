#ifndef TSP_H
#define TSP_H
#include<stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
//parameters for the random values generation
#define MAX_X 10000  //maximum value for ascisse of generated points
#define MAX_Y 10000//maximum value for ordinate of generated points
//default verbose
<<<<<<< HEAD
#define VERBOSE 50
=======
#define VERBOSE 0

typedef enum {

	NOT_DEF,
	RANDOM_SOL,
	NN,
	OPT_2,
	TABU

} solver_id;

>>>>>>> a388ecb2eba8c50bb0897d35427cb28547bed972
typedef struct {
	double x;
	double y;
} point;

typedef struct {

	//input data
	int nnodes;
	point* nodes;

	// parameters 
	solver_id solver;
	int randomseed;
	double timelimit;						// overall time limit, in sec.s
	char input_file[1000];		  			// input file
	int available_memory;
	int verbose;

	//global data
	float* dist_matrix;
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

void print_triangular_matrix_as_array(char flag, const char text_to_print[], const float* matrix, int nrows);

void tsp_debug(char flag,int flag_time, char* format, ...);

void tsp_debug_inline(char flag, char* format, ...);

void free_matrix(void** matrix, int rows);

void free_instance(instance *inst);

void print_instance_parameters(instance* inst);

void parse_command_line(int argc, char** argv, instance *inst);

void opt2_optimize_best_sol(instance *inst);

void reverse_sequence(int* path, int min, int max);

void compute_dist_matrix(instance* inst);

double get_cost_matrix(const double** matrix, int a, int b);

float get_dist_matrix(const float* matrix, int row, int col);

void swap_space(void* a, void* b, size_t size);

void swap(int* a, int* b);

void copy_array(void* a1, const void* a2);

int* search_min(const int* p, const int* end, const float* cost_matrix, double* current_cost);

int* compute_greedy_path(int index_first, instance* inst, double* path_cost);

void greedy_tsp(instance* inst);

void update_best(instance* inst, double z, double t, int* sol);

void init_path(int* path, size_t n);

double compute_path_length(int* path, int nodes_number, point* nodes);

char is_feasible_solution(instance* inst, int* sol_path, double sol_cost);

solver_id parse_solver(char* solver_input);

void tsp_solve(instance* inst);

void update_solver(instance* inst);

#endif 