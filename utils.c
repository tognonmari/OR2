#include "utils.h"
#include <stdio.h>
/*
This function generates and discards 100 pseudorandom numbers using the rand() function.
* The purpose is to mitigate polarization effects in the pseudorandom sequence.
* @note The generated random numbers are not used or returned by this function.
*/

void depolarize_pseudornd_seq() {
	for (int i = 1; i <= 100; i++) {
		rand();
	}
}

/* 
* Function to initialize a timer and get the current time from the first call of the function.
*/
double get_timer() {
	static clock_t start_time = 0;  // Memorizza il tempo di avvio in modo persistente
	if (start_time == 0) {
		start_time = clock();  // Inizializza il tempo di avvio alla prima chiamata
	}

	clock_t current_time = clock();
	double elapsed_time = (double)(current_time - start_time) / CLOCKS_PER_SEC;
	return elapsed_time;
}
/* Returns 1 if two double differ by less than a parameter EPSILON, 0 otherwise
*/
int is_equal_double(double d1, double d2){
	if (fabs(d1-d2)>=EPSILON){
		return 0;
	}
	return 1;
}
