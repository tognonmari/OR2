#ifndef UTILS_H
#define UTILS_H
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>


#if defined(_WIN32) || defined(_WIN64)
#include <direct.h>
#define MKDIR(dir) _mkdir(dir)
#else
#include <sys/stat.h>
#define MKDIR(dir) mkdir(dir, 0755)
#endif

void depolarize_pseudornd_seq();

double get_timer();

void** alloc_matrix(int nrow, int ncol, size_t size_type);

void** alloc_triangular_matrix(int nrow, size_t size_type);

#endif