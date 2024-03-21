#include <stdio.h>
#include "tsp.h"
#include "utils.h"

/*
@ho scoperto che in verità %f è per i float e per i double dovremmo usare %lf.
*/
/*
OR Esito(
	 0: elaborazione riuscita;
	-1: apertura fallita di una pipe;
	-2: apertura fallita di un file;
	-3: numero di parametri inseriti da linea di comando non corretto).
	-4: buffer truncation;
	-5: finito spazio nell'heap (fallisce allocazione);
	-6: Time limit exceeded
	-7: 
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

		print_instance_parameters(&inst);

	}
	get_timer();
	MKDIR("figures");
	MKDIR("data");
	generate_instance(&inst); 
	// The name of the file have the form "graph_data_$nnodes_$randomseed.txt" where $randomseed is the seed used for the random generation
	// and $nnodes is the size of the graph. The file is created in the directory "data".
	/**/
	check_truncation = snprintf(data_file_name, sizeof(data_file_name), "data/graph_data_%d_%d.txt", inst.nnodes, inst.randomseed);
	if (check_truncation < 0 || check_truncation >= sizeof(data_file_name)) {
		return main_error_text(-4, "Buffer: %d char\nText: %d char", sizeof(data_file_name), check_truncation);
	}
	data_file = fopen(data_file_name, "w");
	if (data_file == NULL) {
		fclose(data_file);
		return main_error_text(-2,"File: %s",data_file_name);
	}

	// The name of the figure have the form "graph_$nnodes_$randomseed.png" where $randomseed is the seed used for the random generation
	// and $nnodes is the size of the graph. The file is created in the directory "figures".
	check_truncation = snprintf(figure_name, sizeof(figure_name), "figures/greedy_%d_%d.png", inst.nnodes, inst.randomseed);
	if (check_truncation < 0 || check_truncation >= sizeof(figure_name)) {
		return main_error_text(-4, "Buffer: %d char, Text: %d char", sizeof(figure_name), check_truncation);
	}
	make_datafile(&inst, data_file);
	fclose(data_file);
	//JUST FOR TESTING
	//print_nodes(1, "The nodes of the graph are\n", &inst, inst.nnodes);
	compute_cost_matrix(&inst);
	print_triangular_matrix((inst.verbose>99),"The cost matrix is: \n", (const double**)inst.cost_matrix, inst.nnodes);
	greedy_tsp(&inst);
	print_best_sol((inst.verbose>2),&inst);
	plot_path((inst.verbose>-1), figure_name, inst.best_sol, inst.nodes, inst.nnodes);
	opt2_optimize_best_sol(&inst);
	check_truncation = snprintf(figure_name, sizeof(figure_name), "figures/greedy_%d_%d_2opt.png", inst.nnodes, inst.randomseed);
	if (check_truncation < 0 || check_truncation >= sizeof(figure_name)) {
		return main_error_text(-4, "Buffer: %d char, Text: %d char", sizeof(figure_name), check_truncation);
	}
	plot_path((inst.verbose>-1), figure_name, inst.best_sol, inst.nodes, inst.nnodes);
	print_best_sol((inst.verbose>2), &inst);
	free_instance(&inst);
	return 0;
}
