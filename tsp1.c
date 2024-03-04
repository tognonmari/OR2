#include <stdio.h>
#include <stdlib.h>
#include "tsp1.h"
//SOLO PER TESTING
void printArrayDouble(const char etiche[], const double a[], int n) {
	int i;
	printf("%s", etiche);
	for (i = 0; i < n; i++)
		printf("a[%d]: %f\n", i, a[i]);
} /* printArrDouble */

/*
* The values of an array are transformed in random values belonging to [0,$max_value]
* IP n Length of the array
* P array Array to be modified with random values
* IP max_value Maximum value of the generated values
*/
void generate_array(int n, double *array,int max_value) {
	for (int i = 0; i < n; i++) {
		*array = rand() % (max_value+1);
		array++;
	}
}
/*
* A complete graph is implicited defined by its nodes -> this function generates $n nodes
* in a random way, thus each generated pair of coordinates (x_i,y_i) is s.t.:
* x_i belongs to [0,$max_x]
* x_i belongs to [0,$max_y]
* IP n Number of nodes to generate
* OP inst Instance to be initialized with random nodes
*/
void generate_instance(int n, instance *inst){
	inst->nnodes = n; 
	inst->xcoord = (double*)calloc(inst->nnodes, sizeof(double));
	inst->ycoord = (double*)calloc(inst->nnodes, sizeof(double));
	generate_array(n, inst->xcoord, MAX_X);
	generate_array(n, inst->ycoord, MAX_Y);
}

/*
* A datafile "graph_data_$(inst->nnodes)_$SEED.txt" containing the coordinates of the nodes
* is created in the "data" directory. Each line of the file contains the coordinates of each point
* separeted by a space.
* IP inst Instance from which the data are taken
*/
void make_datafile(instance *inst) {
	// Open the file in writing mode ("w")
	char file_name[256];
	// The name of the file have the form "graph_data_$nnodes_SEED.txt" where SEED is the seed used for the random generation
	// and $nnodes is the size of the graph. The file is created in the directory "data".
	snprintf(file_name, sizeof(file_name), "data/graph_data_%d_%d.txt", inst->nnodes, SEED);
	FILE* file = fopen(file_name, "w");
	if (file == NULL) {
		fprintf(stderr, "The file could not be opened for writing.\n");
		return;
	}

	double* x_pointer = inst->xcoord;
	double* y_pointer = inst->ycoord;

	for (int i = 0; i < inst->nnodes; i++) {
		fprintf(file, "%d %d\n", (int)(*x_pointer), (int)(*y_pointer));
		x_pointer++;
		y_pointer++;
	}
	fclose(file);
}
/*
* Plot the input datafile using GNUPLOT
* IP graph_data Path of the data file to be plotted
* IP graph Path of the figure 
*/
void plot_graph(const char graph_data[], const char graph[]) {
	FILE* gnuplotPipe = popen("gnuplot -persistent", "w");

	if (gnuplotPipe == NULL) {
		fprintf(stderr, "Error starting GNUPLOT.\n");
		return;
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

int main(void) {
	/*
	*/
	instance inst;
	srand(SEED);
	generate_instance(12, &inst);
	make_datafile(&inst);
	plot_graph("data/graph_data_12_13.txt","figures/graph.png");
	printArrayDouble("Stampo xcoord\n", inst.xcoord, inst.nnodes);
	printArrayDouble("Stampo ycoord\n", inst.ycoord, inst.nnodes);
	return 0;
}