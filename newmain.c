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
*/
//TODO Creare una funzione per ogni tipo di errore che lanci un messaggio su stderror adeguato e restituisca l'int associato all'errore.

int main(int argc, char** argv) {

    // Step 1: declare the instance
    instance inst;
    //Step 2: parse the cmd line 
    parse_command_line(argc, argv, &inst);

    //Step 3 : initialize the instance

    generate_instance(&inst); // generates the random instance or scans a file of tsplib

    //Step 3.5 : print instance information up to now, according to verbosity

    print_instance_parameters(&inst);

    //TSPopt(&inst);
    //Step 4 : choose a solver, if not done already

    if (inst.solver == NOT_DEF){

        //Ask the user for a solver to use
        update_solver(&inst);

    }

    //STEP 5: solve + plots

    tsp_solve(&inst);
    //Step 7 : free memory
    
    free_instance(&inst);

    //Step 8 : return
 
    return main_error(0);



    

}
