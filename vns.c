#include "tsp.h"
#include "utils.h"
#include "vns.h"
#include <stdlib.h>


/**
 * Gathers parameter information for vns from cmd line
*/

vns_params parse_and_init_vns_params(){

    //Initialization

    vns_params pars;
    char buf[2];
    char buf2[2];
    int min;
    int max;

    //Scanning for user input: 

    printf("----Choose a minimum number of kicks:----\n");
    fgets(buf, 2, stdin);
    min = atoi(buf);
    printf("you chose as minimum %d\n", min);
    printf("----Choose a maximium number of kicks:----\n");
    fgets(buf2, 2, stdin);
    max = atoi(buf2);
    printf("you chose as maximum %d\n", max);
    assert(min<=max);
    printf("------------------------------------------\n");
    pars.min_kicks = min;
    pars.max_kicks = max;
    pars.gnuplot_pipe = _popen("gnuplot -persist", "w");


    return pars;
}

/**
 * Solve with vns
*/
void vns(instance* inst){

    
    char plot_desired = inst->verbose >= 10;
    
    //Step 1: Parse vns parameters + copy the current best_sol to optimize + initializing gnuplot pipe with correct axis
    vns_params params = parse_and_init_vns_params();

    // copying the current best sol to optimize with opt2 and then to kick
    int* incumbent_sol = (int*) calloc(inst->nnodes, sizeof(int));

    copy_array(incumbent_sol, inst->best_sol);

    double incumbent_cost = inst->zbest;

    //initializing gnuplot

    vns_init_plot_iter_cost(plot_desired, params.gnuplot_pipe, inst);

    //Step 2: Solve with vns upon the incumbent: need to rewrite opt 2
    printf("-----Starting vns heuristics:-------\n");
    int max_iter = 1500;
    int t=1;
    srand(time(NULL));
    //opt2_optimize_best_sol(inst);
    //generate_name(figure_name, sizeof(figure_name), "figures/after2opt_%d_%d.png", inst->nnodes, inst->randomseed);
	//plot_path((inst->verbose>-1),figure_name,incumbent_sol, inst->nodes, inst->nnodes);
    while(!(is_time_limit_exceeded(inst->timelimit))){

        //while its possible do a single 2 opt move, modifying the correct values
        int swaps =0;
        char improvement = 1;
        while(improvement){
            improvement = opt2_move(inst, incumbent_sol, &incumbent_cost);
            swaps++;
           
            t++;
        }

        if (plot_desired) {
            fprintf(params.gnuplot_pipe, "%d %lf\n", t, incumbent_cost);
        }
        
	    //generate_name(figure_name, sizeof(figure_name), "figures/after2opt_%d_%d.png", inst->nnodes, inst->randomseed);
	    //plot_path((inst->verbose>-1),figure_name,incumbent_sol, inst->nodes, inst->nnodes);
        //printf("Finished intensification\n");
        //update the solution if we found a better one
        if (incumbent_cost <inst->zbest){

            printf("Better solution found: updating the instance\n");
            printf("zbest: %lf\n",incumbent_cost);
            copy_array(inst->best_sol, incumbent_sol);
            inst->zbest = incumbent_cost;
            inst->tbest = get_timer();  

        }
        
        //params.min_kicks = 3;
        //params.max_kicks = 3;
        //printf("I am kicking with min : %d, max: %d", params.min_kicks,params.max_kicks);
        int kicks = 4;
        //printf("I want to kick %d times\n", kicks);
        for (int jj = 0; jj<kicks; jj++){
            //printf("kicking!\n");
            kick(inst, incumbent_sol);

            if (plot_desired) {

                fprintf(params.gnuplot_pipe, "%d %lf\n", t, compute_path_length(incumbent_sol, inst->nnodes, inst->nodes));
            }

            t++;
        }
        
       
        incumbent_cost = compute_path_length(incumbent_sol, inst->nnodes, inst->nodes);

        //t += (swaps+kicks);
       
    }

    free(incumbent_sol);
    vns_close_plot_pipe(plot_desired, params.gnuplot_pipe);

}
/**
 * Returns 0 if the path cannot be optimized any further with 2 opt moves
*/
char opt2_move(instance* inst, int* incumbent_sol, double* incumbent_cost){
    
    int n = inst->nnodes;
    point* nodes_list = inst->nodes;
    double best_delta = 0;
    char improvement = 0;
    int best_i = -1;
    int best_j = -1;
    float* dist_matrix = inst->dist_matrix;
    for (int i = 0; i <= n -3; i++) //cambiato da n-2!
    {
        for (int j = i + 2; j <= n-1; j++)
        {
            //double current_dist= (double)(get_dist_matrix((const float*) (dist_matrix),incumbent_sol[i], incumbent_sol[(i+1)]) +get_dist_matrix((const float*) (dist_matrix),incumbent_sol[j], incumbent_sol[(j+1)%n]));
            //double changed_dist = (double)(get_dist_matrix((const float*) (dist_matrix),incumbent_sol[i], incumbent_sol[(j)])+ get_dist_matrix((const float*) (dist_matrix),incumbent_sol[i+1], incumbent_sol[(j+1)]%n));
            
            double current_dist= get_distance(&nodes_list[incumbent_sol[i%n]], &nodes_list[incumbent_sol[(i+1)%n]]) + get_distance(&nodes_list[incumbent_sol[j%n]], &nodes_list[incumbent_sol[(j+1)%n]]);
            double changed_dist = get_distance(&nodes_list[incumbent_sol[i%n]], &nodes_list[incumbent_sol[j%n]]) + get_distance(&nodes_list[incumbent_sol[(i+1)%n]], &nodes_list[incumbent_sol[(j+1)%n]]);
            
            double delta = changed_dist - current_dist;

            if (delta < best_delta)
            {
                improvement =1;
                //printf("found a new best delta\n");
                best_i = i%n;
                best_j = j%n;
                best_delta = delta;
                //printf("Best delta : %lf\n",best_delta);
            }
        }
    }
    if (improvement){
		//tsp_debug((inst->verbose > 49), 0, "I am swapping nodes %d and %d", best_i+1,best_j);
        swap_2_opt(incumbent_sol, (best_i + 1)%n, (best_j)%n);
        (*incumbent_cost) += best_delta;
        return 1;
    }
    else{
        //printf("no further imprevement with 2opt.\n");
        return 0;
    }
}


void kick(instance* inst, int* sol_to_kick){

    int n = inst->nnodes;
    //srand(time(NULL));
    int* random_indexes = (int*) calloc(3, sizeof(int));
    //Step 1: Extract 3 random edges and check their validity
    random_indexes[0] = rand() % n;
    random_indexes[1] = rand() % n;
    random_indexes[2] = rand() % n;
    //check they are distinct
    char success = 0;
    while(!success){
        success =1;
        for (int i=0; i<3; i++){
            for (int j = i+1; j<3; j++ ){
                if(random_indexes[i]==random_indexes[j]){
                    success = 0;
                }
            }
        }
        if (!success){
            random_indexes[0] = rand() % n;
            random_indexes[1] = rand() % n;
            random_indexes[2] = rand() % n;
        }

    }
    //Sort the array
    qsort(random_indexes, 3, sizeof(int), cmp_int_increasing);
    //sort_int_array(random_indexes, 3);
    //printf("Splitting at i = %d,, j= %d, k= %d\n", random_indexes[0], random_indexes[1], random_indexes[2]);
    //Step 2: Reconnect the edges - TODO: implement more than 1 way for reconnections

    int* kicked_sol = (int*) calloc(n, sizeof(int));
    int r =0;
    //Copy up to index random_indexes[0]
    for (int t =0; t<=random_indexes[0]; t++){
        kicked_sol[r] = sol_to_kick[t];
        r++;
    }
    //copy segment from j+1 to k, both included
    for (int t = random_indexes[1]+1; t<=random_indexes[2]; t++){
        kicked_sol[r] = sol_to_kick[t];
        r++;
    }
    
    //copy segment from i+1 to j both included
    for (int t = random_indexes[0]+1; t<=random_indexes[1]; t++){
        kicked_sol[r] = sol_to_kick[t];
        r++;
    }
    
    //copy remaining part of the array from k+1 to 0-1
    for (int t = (random_indexes[2]+1); (t%inst->nnodes) !=0; t++){
        kicked_sol[r] = sol_to_kick[t];
        r++;
    }
    /*
    if (!is_feasible_solution(inst, kicked_sol,compute_path_length(kicked_sol,n,inst->nodes))){
        printf("Wrong kicking: revise implementation.\n");
    }
    */
    //If everything is correct, then write kicked sol into sol to kick
    copy_array(sol_to_kick, kicked_sol);
    free(kicked_sol);
}


void vns_init_plot_iter_cost(char flag, FILE* gnuplotPipe, instance* inst) {

    char figure_name[64];
    if (gnuplotPipe == NULL) {
        fclose(gnuplotPipe);
        exit(main_error_text(-1, "Failed to open the pipeline to gnuplot"));
    }
    if (flag) {
        generate_name(figure_name, sizeof(figure_name), "figures/vns_%d_%d_cost.png", inst->nnodes, inst->randomseed);
        fprintf(gnuplotPipe, "set output '%s'\n", figure_name); //set output name
        fprintf(gnuplotPipe, "set terminal png \n"); //set extension
        fprintf(gnuplotPipe, "set title 'VNS cost'\n");
        fprintf(gnuplotPipe, "set key off\n"); //if key on the name of the file is printed on the plot
        fprintf(gnuplotPipe, "set xlabel 'iter'\n");
        fprintf(gnuplotPipe, "set ylabel 'cost'\n");
        //fprintf(gnuplotPipe, "set yrange [%d:%d]\n",y_range_min, y_range_max);
        fprintf(gnuplotPipe, "set pointsize 0.5\n");
        fprintf(gnuplotPipe, "set grid \n");
        fprintf(gnuplotPipe, "plot '-' with linespoints pt 0.3 lc rgb '#800080'\n");
    }
}

void vns_close_plot_pipe(char flag, FILE* gnuplotPipe) {
    if (flag) {
        fprintf(gnuplotPipe, "e\n"); //This line serves to terminate the input of data for the plot command in gnuplot.
    }
    fclose(gnuplotPipe);
}