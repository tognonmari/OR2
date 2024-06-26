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
	BCFMP,
	HF,
	SF

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
	char check_feasibility;
	// model;
	int ncols;

	//hf param
	char hf_opt2;
	double hf_pfix_start;
	double hf_pfix_scaling;

	//sf param
	int sf_k_start;
	int sf_k_scaling;

	//greedy param
	double gre_perc_start;

	//tabu param
	char tenure_is_variable;
	double tabu_amp;
	double tabu_avg;
	double tabu_freq;

	//vns paaram
	int vns_min_kicks;
	int vns_max_kicks;

} instance;


void make_datafile(instance* inst, FILE* data_file);

void generate_instance(instance* inst);

void generate_name(char buffer[], size_t bufferSize, const char* format, ...);

void generate_nodes(int n, point* nodes, int max_x, int max_y);

double get_distance(const point* p1, const point* p2);

void print_nodes(char flag, const char text_to_print[], const instance* inst, int n);

void print_path(char flag, const char text_to_print[], const int* path, const point* nodes, int n);

void print_best_sol(char flag, instance* inst);

void print_point(char flag, const char text_to_print[], const point* p);

void print_triangular_matrix_as_array(char flag, const char text_to_print[], const float* matrix, int nrows);

void tsp_debug(char flag, int flag_time, char* format, ...);

void tsp_debug_inline(char flag, char* format, ...);

void free_matrix(void** matrix, int rows);

void free_instance(instance* inst);

void print_instance_parameters(instance* inst);

void parse_command_line(int argc, char** argv, instance* inst);

void reverse_sequence(int* path, int min, int max);

void compute_dist_matrix(instance* inst);

float get_dist_matrix(const float* matrix, int row, int col);

void swap_space(void* a, void* b, size_t size);

void swap(int* a, int* b);

void swap_2_opt(int* path, int i, int j);

char opt2_move(char table_flag, instance* inst, int* incumbent_sol, double* incumbent_cost, int* nr_swap);

void opt2_init_table(char table_flag);

void opt2(instance* inst, int* incumbent_sol, double* incumbent_cost, char table_flag);

void copy_array(void* a1, const void* a2);

void copy_din_array(void* a1, const void* a2, size_t elem_size, size_t num_elems);

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

void generate_instance_from_tsplib(instance* inst);

#endif