#include <stdio.h>
#include "TSP.h"
#include "utils.h"

int main(int argc,char **argv){

    instance inst;

    parse_command_line(argc,argv,&inst);

    if(inst.verbose>= 10){

        print_instance_parameters(inst);

    }
    //TODO: read input instance to build the effective instance 
    free_instance(&inst);
    return 0;
    
}