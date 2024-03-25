main : newmain.c tsp.o utils.o tabu.o
	gcc newmain.c utils.o tsp.o tabu.o -o newmain

tsp.o : tsp.c
	gcc -c tsp.c -o tsp.o

utils.o : utils.c 
	gcc -c utils.c -o utils.o

tabu.o : tabu.c
	gcc -c tabu.c -o tabu.o
