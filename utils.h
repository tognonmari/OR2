#ifndef UTILS_H
#define UTILS_H
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include<stdarg.h>


#if defined(_WIN32) || defined(_WIN64)
#include <direct.h>
#define MKDIR(dir) _mkdir(dir)
#elif defined(__linux__) || defined(__APPLE__)
#include <sys/stat.h>
#define MKDIR(dir) mkdir(dir, 0755)
#else
#error ""Unsupported operating system. This program is currently available only for Windows, Linux, and macOS.""
#endif

#define EPSILON 1e-5

int main_error_text(int error, char* format, ...);

int main_error(int error);

void depolarize_pseudornd_seq();

double get_timer();

char is_time_limit_exceeded(double time_limit);

char is_equal_double(double d1, double d2, double epsilon);

void** alloc_matrix(int nrow, int ncol, size_t size_type);

void** alloc_triangular_matrix(int nrow, size_t size_type);

void* alloc_triangular_matrix_as_array(int nrow, size_t size_type);

int cmp_int_increasing(const void* a, const void* b);

char* concatenate_strings(const char* str1, const char* str2);


#endif