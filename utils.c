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
/*
* Method for allocating matrices of size $nrow x $ncol where each element have size $size_type.
* IP nrow, ncol Dimensions of the matrix
* IP size_type Size of the type of the matrix elements
* @howtouse for example to make a 10x10 matrix of double:
*    double** matrix = (double**) alloc_matrix(10,10, sizeof(double))
*/
void** alloc_matrix(int nrow, int ncol, size_t size_type) {
	void** matrix = malloc(nrow * sizeof(void*));

	if (matrix == NULL) {
		exit(-5);
	}
	for (int i = 0; i < nrow; i++) {
		matrix[i] = malloc(ncol * size_type);
		if (matrix[i] == NULL) {
			exit(-5);
		}
	}
	return matrix;
}

/*
* Method for allocating triangular matrices with $nrow rows where each element have size $size_type.
* IP nrow Number of rows
* IP size_type Size of the type of the matrix elements
* @howtouse for example to make a triangular matrix of double of 3 rows:
*    double** matrix = (double**) alloc_matrix(10,10, sizeof(double))
* @note with triangular matrix is meant the following structure (each x represent an entry):
* (example for nrows = 4)
*   0 1 2 3
* 0 x 
* 1 x x
* 2 x x x
* 3 x x x x
*/
void** alloc_triangular_matrix(int nrow, size_t size_type) {
	void** matrix = malloc(nrow * sizeof(void*));

	if (matrix == NULL) {
		exit(-5);
	}
	for (int i = 0; i < nrow; i++) {
		matrix[i] = malloc( (i+1) * size_type);
		if (matrix[i] == NULL) {
			exit(-5);
		}
	}
	return matrix;
}