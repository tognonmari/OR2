#include <stdio.h>
#include <stdlib.h>
#include "tsp1.h"
#include "utils.h"
#include<stdarg.h>
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

/*
OR Esito(
	0: elaborazione riuscita;
	-1: apertura fallita di una pipe;
	-2: apertura fallita di un file;
	-3: numero di parametri inseriti da linea di comando non corretto).
	-4: buffer truncation;
*/ 
//TODO Creare una funzione per ogni tipo di errore che lanci un messaggio su stderror adeguato e restituisca l'int associato all'errore.
int main(int argc, char** argv) {
	/*
	*/
	instance inst;
	FILE (*data_file);
	char data_file_name[256];
	char figure_name[256];
	int check_truncation;
	
	parse_command_line(argc, argv, &inst);
	if (inst.verbose >= 10) {

		print_instance_parameters(inst);

	}
	get_timer();
	MKDIR("figures");
	MKDIR("data");
	generate_instance(&inst); //TODO n e SEED devono essere letti da linea di tastiera
	// The name of the file have the form "graph_data_$nnodes_$randomseed.txt" where $randomseed is the seed used for the random generation
	// and $nnodes is the size of the graph. The file is created in the directory "data".
	check_truncation = snprintf(data_file_name, sizeof(data_file_name), "data/graph_data_%d_%d.txt", inst.nnodes, inst.randomseed);
	if (check_truncation < 0 || check_truncation >= sizeof(data_file_name)) {
		fprintf(stderr, "A buffer has been truncated.\n");
		return -4;
	}
	data_file = fopen(data_file_name, "w");
	if (data_file == NULL) {
		fclose(data_file);
		return -2;
	}

	// The name of the figure have the form "graph_$nnodes_$randomseed.png" where $randomseed is the seed used for the random generation
	// and $nnodes is the size of the graph. The file is created in the directory "figures".
	check_truncation = snprintf(figure_name, sizeof(figure_name), "figures/graph_%d_%d.png", inst.nnodes, inst.randomseed);
	if (check_truncation < 0 || check_truncation >= sizeof(data_file_name)) {
		fprintf(stderr, "A buffer has been truncated.\n");
		return -4;
	}
	make_datafile(&inst, data_file);
	plot_graph(data_file_name,figure_name);
	//JUST FOR TESTING
	print_nodes("The nodes of the graph are\n", &inst, inst.nnodes);
	fclose(data_file);
	free_instance(&inst);
	return 0;
}