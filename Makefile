main : main.c tsp.o utils.o
	gcc main.c utils.o tsp.o tsp.h -o main

tsp.o : tsp.c
	gcc -c tsp.c -o tsp.o

utils.o : utils.c 
	gcc -c utils.c -o utils.o