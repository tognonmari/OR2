#include <stdio.h>
#include "tsp1.h"
#include "utils.h"


/*
OR Esito(
	0: elaborazione riuscita;
	-1: apertura fallita di una pipe;
	-2: apertura fallita di un file;
	-3: numero di parametri inseriti da linea di comando non corretto).
	-4: buffer truncation;
	-5: finito spazio nell'heap (fallisce allocazione);
*/
//TODO Creare una funzione per ogni tipo di errore che lanci un messaggio su stderror adeguato e restituisca l'int associato all'errore.

int main(int argc, char** argv) {
	
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
	/**/
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
	//plot_graph(data_file_name,figure_name);
	//JUST FOR TESTING
	print_nodes("The nodes of the graph are\n", &inst, inst.nnodes);
	compute_cost_matrix(&inst);
	print_triangular_matrix((const double**)inst.cost_matrix, inst.nnodes);
	greedy_tsp(&inst);
	print_best_sol(1, &inst);
	fclose(data_file);
	free_instance(&inst);
	return 0;
}

