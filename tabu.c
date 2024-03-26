
#include "tabu.h"
void printIntArray(const int *arr, int size) {
    printf("[");
    for (int i = 0; i < size; i++) {
        printf("%d", arr[i]);
        if (i < size - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}

void tabu_init_plot_iter_and_cost(char flag, FILE* gnuplotPipe, int y_range_min, int y_range_max, instance* inst){
    char figure_name[64];
    if (gnuplotPipe == NULL) {
        fclose(gnuplotPipe);
		exit(main_error_text(-1,"Failed to open the pipeline to gnuplot"));
    }
    if(flag){
        generate_figure_name(figure_name, sizeof(figure_name), "figures/tabu_%d_%d_cost.png", inst->nnodes, inst->randomseed);
        fprintf(gnuplotPipe, "set output '%s'\n", figure_name); //set output name
        fprintf(gnuplotPipe, "set terminal png\n"); //set extension
        fprintf(gnuplotPipe, "set title 'Tabu Search cost'\n"); 
        fprintf(gnuplotPipe, "set key off\n"); //if key on the name of the file is printed on the plot
        fprintf(gnuplotPipe, "set xlabel 'iter'\n");
        fprintf(gnuplotPipe, "set ylabel 'cost'\n");
        fprintf(gnuplotPipe, "set yrange [%d:%d]\n",y_range_min, y_range_max);
        fprintf(gnuplotPipe, "set pointsize 0.5\n");
        fprintf(gnuplotPipe, "set grid \n");
        fprintf(gnuplotPipe, "plot '-' with linespoints pt 1 lc rgb '#800080'\n");
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
    int* init_sol = compute_greedy_path(0, inst, &(tabu->zcurr) );
    return init_sol;
}
/*
Updates best solution with current solution.
*/
void tabu_update_best(tabu* tabu, int n){
    copy_din_array(tabu->best_sol, tabu->curr_sol, sizeof(int), n);
    tabu->zbest = tabu->zcurr;
    tabu->tbest = get_timer();
}

/*
Initialize the tabu structure:
Step 1) initialize tabu parameters (iter_stop, iter, tenure)
Step 2) Initialize best_sol and curr_sol with a greedy path
*/
void tabu_init(tabu* tabu, instance* inst){
    tabu->iter_stop = MAX_ITER; 
    tabu->iter = 0; 
    tabu->tenure = TENURE;
    tabu->tabu_list = init_tabu_list(inst->nnodes);
    tabu->curr_sol = compute_tabu_init_sol(inst, tabu);
    tabu->best_sol = (int*) calloc(inst->nnodes, sizeof(int));
    assert(tabu->best_sol!=NULL);
    tabu->pipe = popen("gnuplot -persist", "w");
    tabu_init_plot_iter_and_cost( (inst->verbose)>=1, tabu->pipe,  (int)( (tabu->zcurr) * 0.75), (int)( (tabu->zcurr) * 0.8), inst);
    tabu_update_best(tabu, inst->nnodes);
}
char isTabu(tabu* tabu, int vertex){
    return tabu->iter < tabu->tabu_list[vertex] + tabu->tenure;
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
            jump_to_next_iter = ( isTabu(tabu, j) || isTabu(tabu, (j+1)%n) );
            if(jump_to_next_iter) { continue;}
            //length (i,i+1) and (j,j+1)
            double current_dist = (double)(get_dist_matrix((const float*)(dist_matrix), curr[i], curr[i+1] ) + get_dist_matrix((const float*)(dist_matrix), curr[j], curr[(j+1)%n] ));
            //length (i,j) and (i+1,j+1)
            double changed_dist = (double)(get_dist_matrix((const float*)(dist_matrix), curr[i], curr[j] )+ get_dist_matrix((const float*)(dist_matrix), curr[(i+1)], curr[(j+1)%n] ));
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
    return exist_admissible_move;
}

void tabu_update_current(tabu* tabu){
    int v1 = (tabu->best_admissible_move).vertex_to_swap_1;
    int v2 = (tabu->best_admissible_move).vertex_to_swap_2;
    swap_2_opt(tabu->curr_sol, v1, v2);
    tabu->zcurr += (tabu->best_admissible_move).delta;
    fprintf(tabu->pipe, "%d %lf\n", tabu->iter, tabu->zcurr);
}

void tabu_update_list(tabu* tabu){
    int j = (tabu->best_admissible_move).vertex_to_swap_2;
    (tabu->tabu_list)[j] = tabu->iter;
}

void tabu_update(tabu* tabu, int n){
    tabu_update_list(tabu);
    tabu_update_current(tabu);
    if(tabu->zcurr < tabu-> zbest){
        tabu_update_best(tabu, n);
    }
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

void tabu_close_plot(FILE* gnuplotPipe){
    fprintf(gnuplotPipe, "e\n"); //This line serves to terminate the input of data for the plot command in gnuplot.
    fclose(gnuplotPipe);
}
void tabu_free(tabu* tabu){
    //free(tabu->best_sol); NON DEVI FARE QUESTO PERCHe inst best sol viene aggiornata a tabu best sol alla fine, quindi devi solo liberare inst best sol e lo fai gia in tsp_free
    free(tabu->curr_sol);
    free(tabu->tabu_list);
}

void tabu_search(instance* inst){
    tabu tab;
    tabu_init(&tab, inst);
    while(tab.iter <= tab.iter_stop && !(is_time_limit_exceeded(inst->timelimit)) ){
        if(tabu_find_best_admissible_move(&tab, inst)){
            tabu_update(&tab, inst->nnodes);
        }
        else {
            tsp_debug((inst->verbose>0),0,"There aren't admissible moves at iter %d", tab.iter);
        }
        tabu_debug(1, &tab, inst);
        (tab.iter)++;
    }
    tabu_close_plot(tab.pipe);
    update_best(inst, tab.zbest, tab.tbest, tab.best_sol);
    tabu_free(&tab);
}