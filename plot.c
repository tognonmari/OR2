#include "plot.h"

void plot_edge(FILE* gnuplot_pipe, const point* p1, const point* p2) {
    fprintf(gnuplot_pipe, "%lf %lf\n", p1->x, p1->y);
    fprintf(gnuplot_pipe, "%lf %lf\n", p2->x, p2->y);
}

void plot_tour(FILE* gnuplot_pipe, const int* succ, int start, const point* nodes) {
    fprintf(gnuplot_pipe, "\n");
    int i1 = start;
    int i2 = succ[i1];
    while (i2 != start) {
        plot_edge(gnuplot_pipe, &nodes[i1], &nodes[i2]);
        i1 = i2;
        i2 = succ[i1];
    }
    plot_edge(gnuplot_pipe, &nodes[i1], &nodes[i2]); //last edge that close the tour
    fflush(gnuplot_pipe);
}

void plot_multitour(const multitour_sol* mlt, int n, const point* nodes, const char figure_name[]) {
    FILE* gnuplot_pipe = _popen("gnuplot -persist", "w");
    fprintf(gnuplot_pipe, "set output '%s'\n", figure_name); //set output name
    fprintf(gnuplot_pipe, "set terminal png\n"); //set extension
    fprintf(gnuplot_pipe, "set title 'Multitour'\n");
    fprintf(gnuplot_pipe, "set key off\n"); //if key on the name of the file is printed on the plot
    fprintf(gnuplot_pipe, "set grid \n");
    fprintf(gnuplot_pipe, "plot '-' with lines lc rgb '#800080'\n");
    int* start = (int*)calloc(mlt->ncomp + 1, sizeof(int));
    for (int i = 0; i < n; i++) {
        start[mlt->comp[i]] = i;
    }
    plot_tour(gnuplot_pipe, mlt->succ, start[1], nodes);
    for (int k = 2; k <= mlt->ncomp; k++) {
        plot_tour(gnuplot_pipe, mlt->succ, start[k], nodes);
    }
    free(start);
    fprintf(gnuplot_pipe, "e\n"); //This line serves to terminate the input of data for the plot command in gnuplot.
    fflush(gnuplot_pipe);
    fclose(gnuplot_pipe);
}

/*
* void plot_figures(solver_id solv_id, const instance* inst, const multitour_sol* mlt) {
    FILE* gnuplotPipe = _popen("gnuplot -persist", "w");
    char figure_name[64];
    switch (solv_id) {
    case BEN:
        generate_name(figure_name, sizeof(figure_name), "figures/ben_%d_%d_multitour.png", inst->nnodes, inst->randomseed);


    }
    fprintf(gnuplotPipe, "set output '%s'\n", figure_name); //set output name
    fprintf(gnuplotPipe, "set terminal png\n"); //set extension
}

*/
