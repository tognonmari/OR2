
#include "tabu.h"
#include "plot.h"

int tabu_get_tenure(instance* inst, int num_iterations, int nnodes) {
    //Parameters
    double amplitude = 0.1 * nnodes;
    double frequency = 0.1;
    double phase_shift = 0.0; //do not change
    double average = 2.0;

    //for perf prof
    //amplitude = inst->tabu_amp * nnodes;
    frequency = inst->tabu_freq;
    //average = inst->tabu_avg;


    int tenure =(int) amplitude *(average + sin(frequency * num_iterations + phase_shift));
    return tenure;
}
void tabu_init_data_iter_and_cost(char flag, FILE* data_file, int y_range_min, int y_range_max, instance* inst){
    char figure_name[64];
    if(flag){
         generate_name(figure_name, sizeof(figure_name), "data/tabu_%d_%d_cost.txt", inst->nnodes, inst->randomseed);
         data_file = fopen(figure_name, "w");
         if(data_file == NULL){
            fclose(data_file);
            exit(main_error_text(-2,"Failed to open the file %s", figure_name));
         }
    }
}
void tabu_init_plot_iter_and_cost(char flag, FILE* gnuplotPipe, int y_range_min, int y_range_max, instance* inst){
    char figure_name[64];
    if (gnuplotPipe == NULL) {
        fclose(gnuplotPipe);
		exit(main_error_text(-1,"Failed to open the pipeline to gnuplot"));
    }
    if(flag){
        generate_name(figure_name, sizeof(figure_name), "figures/tabu_%d_%d_cost.png", inst->nnodes, inst->randomseed);
        fprintf(gnuplotPipe, "set output '%s'\n", figure_name); //set output name
        fprintf(gnuplotPipe, "set terminal png\n"); //set extension
        fprintf(gnuplotPipe, "set title 'Tabu Search cost'\n"); 
        fprintf(gnuplotPipe, "set key off\n"); //if key on the name of the file is printed on the plot
        fprintf(gnuplotPipe, "set xlabel 'iter'\n");
        fprintf(gnuplotPipe, "set ylabel 'cost'\n");
        //fprintf(gnuplotPipe, "set xrange [%d:%d]\n", 0, 1000);
        //fprintf(gnuplotPipe, "set yrange [%d:%d]\n",y_range_min, y_range_max);
        fprintf(gnuplotPipe, "set pointsize 0.5\n");
        fprintf(gnuplotPipe, "set grid \n");
        //fprintf(gnuplotPipe, "plot '-' with lines pt 1 lc rgb '#800080'\n");
        fprintf(gnuplotPipe, "plot '-' with points lc rgb '#800080'\n");
    }
}
/*
Initialize each value of tabu_list to the minimum int
*/
int* init_tabu_list(size_t n){
    int* tabu_list = (int*) malloc(n*sizeof(int));
    assert(tabu_list!=NULL);
    int* pointer = tabu_list;
    for(int i=0; i<n; i++){
		(*pointer) = INT_MIN;
		pointer++;
	} 
    return tabu_list;
}

/*
Computes initial solution by constructing a greedy solution starting from node 0,
updates the field zcurr.

* IP inst Pointer to the instance containing node information
* IOP tabu Pointer to the tabu structure 
*       OP tabu->zcurr The cost of best sol is updated with current sol
* OR A pointer to the constructed greedy path chose as initial solution
*/

int* compute_tabu_init_sol(instance* inst, tabu* tabu){
    /*
    int* init_sol = compute_greedy_path(0, inst, &(tabu->zcurr));
    return init_sol;
    */
    greedy gre;
    gre_init(inst, &gre);
    tabu->zcurr = gre.zcurr;
    tabu->zbest = gre.zcurr;
    return gre.curr_sol;
}
/*
Updates best solution with current solution.
*/
void tabu_update_best(tabu* tabu, int n, int verbose){
    char flag = verbose >= 10;
    copy_din_array(tabu->best_sol, tabu->curr_sol, sizeof(int), n);
    tsp_debug_inline(flag, "\n -------- Update Tabu Best at iter %d : ---------\n", tabu->iter);
    tsp_debug_inline(flag, "Old Best : %.20g\n", tabu->zbest);
    tsp_debug_inline(flag, "New Best : %.20g\n", tabu->zcurr);
    tsp_debug_inline(flag, " --------------- End Tabu Update : --------------\n");
    tabu->zbest = tabu->zcurr;
    tabu->tbest = get_timer();
}

/*
Initialize the tabu structure:
Step 1) initialize tabu parameters (iter, tenure)
Step 2) Initialize best_sol and curr_sol with a greedy path
*/
void tabu_init(char tenure_is_variable, tabu* tabu, instance* inst){
    tabu->iter = 0; 
    tabu->tenure_is_variable = tenure_is_variable;
    tabu->tenure = tabu_get_tenure(inst, 0, inst->nnodes);
    tabu->figure_cost_flag = (inst->verbose) >= 1;
    tabu->tabu_list = init_tabu_list(inst->nnodes);
    tabu->printing_period = 1;
    tabu->curr_sol = compute_tabu_init_sol(inst, tabu);
    tabu->best_sol = (int*) calloc(inst->nnodes, sizeof(int));
    assert(tabu->best_sol!=NULL);
    tabu->pipe = _popen("gnuplot -persist", "w");
    tabu_init_plot_iter_and_cost( tabu->figure_cost_flag, tabu->pipe,  (int)( (tabu->zcurr) * 0.8), (int)( (tabu->zcurr*0.85) ), inst);
    //init_data_file(tabu -> figure_cost_flag, tabu->data_iter_and_cost, inst);
    //tabu_init_plot_iter_and_cost( tabu->figure_cost_flag, tabu->pipe,  215000, 220000, inst);

    tabu_update_best(tabu, inst->nnodes, inst->verbose);
    tsp_debug(inst->verbose >= 100, 1, "tabu_init SUCCESSFULL");
}
char isTabu(tabu* tabu, int vertex){
    return tabu->iter < tabu->tabu_list[tabu->curr_sol[vertex]] + tabu->tenure;
}
void update_move(move* mov, int i1, int i2, double delta){
    mov->vertex_to_swap_1 = i1;
    mov->vertex_to_swap_2 = i2;
    mov->delta = delta;
    mov->improvement = ( delta < 0 );
}
/*
Update $(tabu->best_admissible_move) with the current best move which is not a tabu move
(hence, neither i,i+1,j and j+1 are tabu), the function returns true if an
admissable move (which is not tabu) exists
IP inst
IOP tabu
    OP tabu->best_admissible move updated with the value of the current move if $(tabu->exist_admissibile_move)
OR exist_admissible_move True if exists a no tabu move false otherwise.
*/
char tabu_find_best_admissible_move(tabu* tabu, instance* inst){
    char exist_admissible_move = 0;
    float* dist_matrix = inst->dist_matrix;
    int n = inst->nnodes;
    int* curr = tabu->curr_sol;
    point* nodes_list = inst->nodes;
    double best_delta = +DBL_MAX;
    int best_i = -1;
    int best_j = -1;
    for (int i = 0; i <= n-3 ; i++){
        char jump_to_next_iter;
        jump_to_next_iter = ( isTabu(tabu, i) || isTabu(tabu, i+1) );
        if(jump_to_next_iter) { continue;}
        for (int j = i+2; j <= n-1; j++){
            jump_to_next_iter = ( isTabu(tabu, j) || isTabu(tabu, (j+1)%n) ) || ( (i==0) && (j+1 == n) );
            // se i = 0 e j+1 = n quindi j+1 %n = 0 allora sto considerando due lati adiacenti, mentre 2 opt deve considerare sempre lati non adaicenti 
            if(jump_to_next_iter) { continue;}
            //length (i,i+1) and (j,j+1)
            double current_dist =(double)get_dist_matrix((const float*)(dist_matrix), curr[i], curr[i+1])  + (double)get_dist_matrix((const float*)(dist_matrix), curr[j], curr[(j+1)%n] );
            //length (i,j) and (i+1,j+1)
            double changed_dist =(double)get_dist_matrix((const float*)(dist_matrix), curr[i], curr[j]) + (double)get_dist_matrix((const float*)(dist_matrix), curr[(i+1)], curr[(j+1)%n]);
            double delta = changed_dist - current_dist;
            if( delta < best_delta){
                best_i = i;
                best_j = j;
                best_delta = delta;
                exist_admissible_move = 1;
            }
        }
    }
    if(exist_admissible_move){
        update_move(&(tabu->best_admissible_move), best_i + 1, best_j, best_delta);
    }
    tsp_debug(inst->verbose >= 100, 1, "tabu_find_best_move at iter %d SUCCESSFULL", tabu->iter);
    return exist_admissible_move;
}

void tabu_update_current(tabu* tabu, int verbose){
    int v1 = (tabu->best_admissible_move).vertex_to_swap_1;
    int v2 = (tabu->best_admissible_move).vertex_to_swap_2;
    swap_2_opt(tabu->curr_sol, v1, v2);
    tabu->zcurr += (tabu->best_admissible_move).delta;
    if((tabu->figure_cost_flag) && (tabu->iter % tabu->printing_period == 0)){
        fprintf(tabu->pipe, "%d %lf\n", tabu->iter, tabu->zcurr);
    }
    tsp_debug(verbose >= 100, 0, "  tabu_update_current SUCCESSFULL");
}

void tabu_update_list(tabu* tabu, int verbose){
    int j = (tabu->best_admissible_move).vertex_to_swap_2;
    (tabu->tabu_list)[tabu->curr_sol[j]] = tabu->iter;
    tsp_debug(verbose >= 100, 0, "  tabu_update_list SUCCESSFULL");
}

void tabu_update_tenure(tabu* tabu, instance* inst) {
    int n = inst->nnodes;
    int verbose = inst->verbose;
    tabu->tenure = tabu_get_tenure(inst, tabu->iter, n);
    tsp_debug(verbose >= 100, 0, "  tabu_update_tenure SUCCESSFULL");
}
void tabu_update(tabu* tabu, instance* inst){
    tabu_update_list(tabu, inst->verbose); 
    tabu_update_current(tabu, inst->verbose);
    if(tabu->zcurr < tabu-> zbest){
        tabu_update_best(tabu, inst->nnodes, inst->verbose);
    }
    if (tabu->tenure_is_variable) {
        tabu_update_tenure(tabu, inst);
    }
    tsp_debug(inst->verbose >= 100, 1, "tabu_update SUCCESSFULL");
}

void tabu_debug(char flag, tabu* tabu, instance* inst){
    move mov = tabu->best_admissible_move;
    tsp_debug_inline(flag, "\nITER = %d\nSWAP BETWEEN %d %d\nDELTA =%.1f\n", tabu->iter, mov.vertex_to_swap_1, mov.vertex_to_swap_2, mov.delta);
    tsp_debug_inline(flag, "TABU = ");
    for(int i = 0; i<inst->nnodes; i++){
        if(isTabu(tabu, i)){
            tsp_debug_inline(flag, "%d ",i);
        }
    }
    
}

void tabu_close(char flag, FILE* gnuplotPipe, FILE* data_file){
    if(flag){
        fprintf(gnuplotPipe, "e\n"); //This line serves to terminate the input of data for the plot command in gnuplot.
    }
    fclose(gnuplotPipe);
}

void tabu_free(tabu* tabu){
    //free(tabu->best_sol); NON DEVI FARE QUESTO PERCHe inst best sol viene aggiornata a tabu best sol alla fine, quindi devi solo liberare inst best sol e lo fai gia in tsp_free
    free(tabu->curr_sol);
    free(tabu->tabu_list);
}

void tabu_search(char tabu_is_variable, instance* inst){
    tabu tab;
    tabu_init(tabu_is_variable, &tab, inst);
    while(!(is_time_limit_exceeded(inst->timelimit)) ){
        tsp_debug(inst->verbose >= 100, 1, "---------------- Iter #%d ---------------- ", tab.iter);
        if(tabu_find_best_admissible_move(&tab, inst)){
            tabu_update(&tab, inst);
        }
        else {
            tsp_debug((inst->verbose>0),0,"There aren't admissible moves at iter %d", tab.iter);
        }
        tabu_debug((inst->verbose)>=200, &tab, inst);
        (tab.iter)++;
    }
    tsp_debug(inst->verbose >= 100, 1, "time limit has been exceeded");
    update_best(inst, tab.zbest, tab.tbest, tab.best_sol);
    opt2(inst, inst->best_sol, &(inst->zbest), 0); // Attenzione la figura iter_and_cost non contiene il costo finale migliorato da opt2!
    tabu_close( tab.figure_cost_flag , tab.pipe, tab.data_iter_and_cost);
    tabu_free(&tab);
}
