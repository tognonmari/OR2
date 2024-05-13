#ifndef TSP_H
#define TSP_H
#include<stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "assert.h"
//parameters for the random values generation
#define MAX_X 10000  //maximum value for ascisse of generated points
#define MAX_Y 10000//maximum value for ordinate of generated points
//default verbose
#define VERBOSE 1

typedef enum {

	NOT_DEF,
	BAD_DEF,
	RANDOM_SOL,
	NN,
	OPT_2,
	TABU,
	VNS,
	EX,
	BEN,
	GLU,
	BC,
	BCM,
	BCP,
	BCMP,
	BCF,
	BCFM,
	BCFP,
	BCFMP

} solver_id;

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
	char* csv_column_name[128];
	char* legend;
	int randomseed;
	double timelimit;						// overall time limit, in sec.s
	char input_file[1000];		  			// input file
	int available_memory;
	int verbose;
	int integer_costs;
	FILE* best_sol_data;
	//global data
	float* dist_matrix;
	double	tstart;
	double zbest;							// best sol. available  
	double tbest;							// time for the best sol. available  
	int* best_sol;							// best sol. available
	char is_best_sol_avail;					// flag that tells if best sol. is available in the heap
	double	best_lb;						// best lower bound available  
	double best_ub;
	double* load_min;						// minimum load when leaving a node
	double* load_max;						// maximum load when leaving a node
	char posting;
	// model;
	int ncols;
	/*
	int xstart;
	int qstart;
	int bigqstart;
	int sstart;
	int bigsstart;
	int ystart;
	int fstart;
	int zstart;
	*/
} instance;


void make_datafile(instance *inst, FILE* data_file);

void generate_instance(instance *inst);

void generate_name(char buffer[], size_t bufferSize, const char *format, ...);

void generate_nodes(int n, point* nodes, int max_x, int max_y);

void generate_array(int n, double *array, int max_value);

double get_distance(const point* p1,const point* p2);

void print_nodes(char flag, const char text_to_print[], const instance *inst, int n);

void print_path(char flag, const char text_to_print[], const int* path, const point* nodes, int n);

void print_best_sol(char flag, instance* inst);

void print_point(char flag, const char text_to_print[], const point* p);

void print_triangular_matrix(char flag,const char text[], const float** matrix, int nrows);

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

double get_cost_matrix(const float** matrix, int a, int b);

float get_dist_matrix(const float* matrix, int row, int col);

void swap_space(void* a, void* b, size_t size);

void swap(int* a, int* b);

void swap_2_opt(int* path, int i, int j);

char opt2_move(char table_flag, instance* inst, int* incumbent_sol, double* incumbent_cost, int* nr_swap);

void opt2_init_table(char table_flag);

void opt2(instance* inst, int* incumbent_sol, double* incumbent_cost, char table_flag);

void copy_array(void* a1, const void* a2);

void copy_din_array(void *a1, const void *a2, size_t elem_size, size_t num_elems);

int* search_min(const int* p, const int* end, const float* cost_matrix, double* current_cost);

int* compute_greedy_path(int index_first, instance* inst, double* path_cost);

void greedy_tsp(instance* inst);

void update_best(instance* inst, double z, double t, int* sol);

void init_path(int* path, size_t n);

void init_data_file(char flag, FILE* data_file, instance* inst);

double compute_path_length(int* path, int nodes_number, point* nodes);

void check_sol_is_feasible(char check_feasibility, instance* inst, int* sol, double zsol);

char is_feasible_solution(instance* inst, int* sol_path, double sol_cost);

solver_id parse_solver(char* solver_input);

void tsp_solve(instance* inst);

void update_solver(instance* inst);

void update_best(instance* inst, double z, double t, int* sol);

void update_lb(instance* inst, double lb);

void update_ub(instance* inst, double ub);

void generate_test_bed(int size_test_bed, int argc, char** argv, instance* test_bed);

void generate_csv_file(int size_test_bed, instance* test_bed);

void read_test_bed_size(int* test_bed_size, int argc, char** argv);

#endif