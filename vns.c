#include "tsp.h"
#include "utils.h"
#include "vns.h"

/**
 * Gathers parameter information for vns from cmd line
*/

vns_params parse_vns_params(){

    int min;
    int max;
    printf("----Choose a minimum number of kicks:----\n");
    scanf("%d",&min);
    printf("----Choose a maximium number of kicks:----\n");
    scanf("%d",&max);
    assert(min<=max);
    printf("------------------------------------------\n");
}

/**
 * Solve with vns
*/
void vns(instance* inst){


    //Step 1: Parse vns parameters + copy the current best_sol to optimize
    vns_params params = parse_vns_params();

    // copying the current best sol to optimize with opt2 and then to kick
    int* incumbent_sol = (int*) calloc(inst->nnodes, sizeof(int));

    copy_array(incumbent_sol, inst->best_sol);

    double incumbent_cost = inst->zbest;
    //Step 2: Solve with vns upon the incumbent: need to rewrite opt 2
    printf("-----Starting vns heuristics:-------\n");
    int max_iter = 1500;
    int t=0;
    //opt2_optimize_best_sol(inst);
    while(t< max_iter){

        //while its possible do a single 2 opt move, modifying the correct values
        int swaps =0;
        char improvement = 1;
        while(improvement){
            improvement = opt2_move(inst, incumbent_sol, &incumbent_cost);
            swaps++;
        }
        
        char figure_name[64];
	    generate_figure_name(figure_name, sizeof(figure_name), "figures/after2opt_%d_%d.png", inst->nnodes, inst->randomseed);
	    plot_path((inst->verbose>-1),figure_name,incumbent_sol, inst->nodes, inst->nnodes);
        //printf("Finished intensification\n");
        //update the solution if we found a better one
        if (is_feasible_solution(inst, incumbent_sol, incumbent_cost) && incumbent_cost <inst->zbest){

            printf("Better solution found: updating the instance\n");
            printf("zbest: %lf\n",incumbent_cost);
            copy_array(inst->best_sol, incumbent_sol);
            inst->zbest = incumbent_cost;
            inst->tbest = get_timer();  

        }
        
        srand(time(NULL));
        params.min_kicks = 3;
        params.max_kicks = 3;
        //printf("I am kicking with min : %d, max: %d", params.min_kicks,params.max_kicks);
        int kicks = (rand() % (params.max_kicks - params.min_kicks+1)) + params.min_kicks;
        //printf("I want to kick %d times\n", kicks);
        for (int jj = 0; jj<kicks; jj++){
            kick(inst, incumbent_sol);
            //printf("kicking!\n");
        }
        
       
        incumbent_cost = compute_path_length(incumbent_sol, inst->nnodes, inst->nodes);

        t += (swaps+kicks);
       
    }

    free(incumbent_sol);
    //opt2_optimize_best_sol(inst);
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
    for (int i = 0; i <= n -3; i++) //cambiato da n-2!
    {
        for (int j = i + 2; j <= n-1; j++)
        {
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
    srand(time(NULL));
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
    
    if (!is_feasible_solution(inst, kicked_sol,compute_path_length(kicked_sol,n,inst->nodes))){
        printf("Wrong kicking: revise implementation.\n");
    }

    //If everything is correct, then write kicked sol into sol to kick
    copy_array(sol_to_kick, kicked_sol);
    free(kicked_sol);
}

