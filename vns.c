
#include "vns.h"

/**
 * Gathers parameter information for vns from cmd line
*/

void init_vns_params(vns_params* pars) {
    pars->gnuplot_pipe = _popen("gnuplot -persist", "w");
    pars->min_kicks = 1;
    pars->max_kicks = 1;
}
/**
 * Solve with vns
*/
void vns(instance* inst) {


    char plot_desired = inst->verbose >= 1;


    //Step 1: Parse vns parameters + copy the current best_sol to optimize + initializing gnuplot pipe with correct axis
    vns_params params;

    init_vns_params(&params);

    //to be commented if not passing cmd line arguments
    //params.min_kicks = inst->vns_min_kicks;
    //params.max_kicks = inst->vns_max_kicks;

    // copying the current best sol to optimize with opt2 and then to kick
    int* incumbent_sol = (int*) calloc(inst->nnodes, sizeof(int));
    int* best_vns_sol = (int*)calloc(inst->nnodes, sizeof(int));
    copy_array(incumbent_sol, inst->best_sol);
    copy_array(best_vns_sol, inst->best_sol);
    double incumbent_cost = inst->zbest;
    double best_vns_cost = inst->zbest;
    double best_time = inst->tbest;
    //initializing gnuplot

    vns_init_plot_iter_cost(plot_desired, params.gnuplot_pipe, inst);

    //Step 2: Solve with vns upon the incumbent: need to rewrite opt 2
    tsp_debug(inst->verbose >=10, 1,"-----Starting vns heuristics:-------\n");

    int t=1;
    srand(time(NULL));
    //generate_name(figure_name, sizeof(figure_name), "figures/after2opt_%d_%d.png", inst->nnodes, inst->randomseed);
	//plot_path((inst->verbose>-1),figure_name,incumbent_sol, inst->nodes, inst->nnodes);
    while(!(is_time_limit_exceeded(inst->timelimit))){

        //while its possible do a single 2 opt move, modifying the correct values
        int swaps =0;
        char improvement = 1;
        while(improvement){
            improvement = opt2_move(0, inst, incumbent_sol, &incumbent_cost, &swaps);
           
            
        }
        t++;

        if (plot_desired) {
            fprintf(params.gnuplot_pipe, "%d %lf\n", t, incumbent_cost);
        }
        
	    //generate_name(figure_name, sizeof(figure_name), "figures/after2opt_%d_%d.png", inst->nnodes, inst->randomseed);
	    //plot_path((inst->verbose>-1),figure_name,incumbent_sol, inst->nodes, inst->nnodes);
        //printf("Finished intensification\n");
        //update the solution if we found a better one
        if (incumbent_cost <best_vns_cost){

            tsp_debug(inst->verbose >=100, 1,"Better solution found: updating the instance\n");
            
            //update_best(inst, incumbent_cost, get_timer(), incumbent_sol);
            //update the best vns sol
            copy_array(best_vns_sol, incumbent_sol);
            best_vns_cost = incumbent_cost;
            best_time = get_timer();

        }
        
        //params.min_kicks = 3;
        //params.max_kicks = 3;
        //printf("I am kicking with min : %d, max: %d", params.min_kicks,params.max_kicks);
        int kicks = (int) (rand_01()*(params.max_kicks- params.min_kicks))+params.min_kicks;
        //printf("I want to kick %d times\n", kicks);
        for (int jj = 0; jj<kicks; jj++){
            //printf("kicking!\n");
            kick(inst, incumbent_sol);

            if (plot_desired) {

                fprintf(params.gnuplot_pipe, "%d %lf\n", t, compute_path_length(incumbent_sol, inst->nnodes, inst->nodes));
            }

            
        }
        t++;
       
        incumbent_cost = compute_path_length(incumbent_sol, inst->nnodes, inst->nodes);

        
       
    }

    vns_close_plot_pipe(plot_desired, params.gnuplot_pipe);

    if (best_vns_cost < inst->zbest) {

        update_best(inst, best_vns_cost, best_time, best_vns_sol);
    }

    free(incumbent_sol);

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
        fprintf(gnuplotPipe, "plot '-' with points lc rgb '#800080'\n");
    }
}

void vns_close_plot_pipe(char flag, FILE* gnuplotPipe) {
    if (flag) {
        fprintf(gnuplotPipe, "e\n"); //This line serves to terminate the input of data for the plot command in gnuplot.
    }
    fclose(gnuplotPipe);
}