#include <stdio.h>
#include "tsp.h"
#include "utils.h"
//#include "TSPmodel.h"


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
    -7: Wrong solver definition
    -8: parametri inseriti da linea di comando errati.
    -9: error in CPLEX enviroment
   -10: Incompatible type error
*/
//TODO Creare una funzione per ogni tipo di errore che lanci un messaggio su stderror adeguato e restituisca l'int associato all'errore.

int main(int argc, char** argv) {

    //Step 0: hardcoded instance test_bed_size, default should be 1
    int test_bed_size = 1;
    // Step 0.1: parse the input to get test bed size
    read_test_bed_size(&test_bed_size, argc, argv);

    //Step 1: Allocate the array of instances
    instance* inst_array = (instance*) calloc(test_bed_size, sizeof(instance));

    //Step 2: Initialize all instances of the test bed
    generate_test_bed(test_bed_size, argc, argv, inst_array);

    //Step 3.5 : print instance information up to now, according to verbosity

    print_instance_parameters(&inst_array[0]);

    //TSPopt(&inst);
    //Step 4 : choose a solver, if not done already

    if (inst_array[0].solver == NOT_DEF) {

        //Ask the user for a solver to use and update it upon all instances
        update_solver(&inst_array[0]);
        solver_id s = inst_array[0].solver;
        for (int i = 1; i < test_bed_size; i++) {
            inst_array[i].solver = s;
        }

    }

    //STEP 5: solve + plots : TODO: Name the plot also with the index of the instance
    for (int i = 0; i < test_bed_size; i++) {
        tsp_solve(&inst_array[i]);
    }

    //Step 6 : generate file for perf prof
    
    generate_csv_file(test_bed_size, inst_array);

    //Step 7 : free memory
    for (int j = 0; j < test_bed_size; j++) {
        free_instance(&inst_array[j]);
    }

    //Step 8 : return
 
    return main_error(0);



    

}
