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

void plot_multitour(char flag, const multitour_sol* mlt, int n, const point* nodes, const char figure_name[]) {
    if (!flag) { return; }
    FILE* gnuplot_pipe = _popen("gnuplot -persist", "w");
    fprintf(gnuplot_pipe, "set output '%s'\n", figure_name); //set output name
    fprintf(gnuplot_pipe, "set terminal png\n"); //set extension
    fprintf(gnuplot_pipe, "set title 'Multitour (Cost  = %.1lf'\n", mlt->z);
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

void plot_path(char flag, const int* path, int n, double cost, const point* nodes, const char figure_name[]) {
    if (!flag) { return; }
    FILE* gnuplot_pipe = _popen("gnuplot -persist", "w");
    fprintf(gnuplot_pipe, "set output '%s'\n", figure_name); //set output name
    fprintf(gnuplot_pipe, "set terminal png\n"); //set extension
    fprintf(gnuplot_pipe, "set title 'Path (Cost  = %.1lf'\n", cost);
    fprintf(gnuplot_pipe, "set key off\n"); //if key on the name of the file is printed on the plot
    fprintf(gnuplot_pipe, "set grid \n");
    fprintf(gnuplot_pipe, "plot '-' with lines lc rgb '#800080'\n");
    for (int i = 0; i < n - 1 ; i++) {
        plot_edge(gnuplot_pipe, &nodes[path[i]], &nodes[path[i + 1]]);
    }
    plot_edge(gnuplot_pipe, &nodes[path[n-1]], &nodes[path[0]]);
    fprintf(gnuplot_pipe, "e\n"); //This line serves to terminate the input of data for the plot command in gnuplot.
    fflush(gnuplot_pipe);
    fclose(gnuplot_pipe);
}
// Funzione per centrare una stringa all'interno di un campo di larghezza specificata
void fprintf_center(FILE* file, const char text[], int field_width) {
    int text_length = strlen(text);
    int padding = (field_width - text_length) / 2;

    // Stampa spazi a sinistra per centrare il testo
    for (int i = 0; i < padding; ++i) {
        fprintf(file, " ");
    }

    // Stampa il testo
    fprintf(file, "%s", text);

    // Stampa spazi a destra (se necessario)
    for (int i = 0; i < padding; ++i) {
        fprintf(file, " ");
    }

    // Se la lunghezza del testo è dispari, aggiungi uno spazio finale
    if (text_length % 2 != 0) {
        fprintf(file, " ");
    }
}


void make_first_row(char flag, FILE* file, int num_cols, const char* text) {
    fprintf(file, "\n+");
    for (int i = 0; i < num_cols - 1; i++) {
        fprintf(file, "%s", "---------------------");
    }
    fprintf(file, "--------------------+\n");
    fprintf(file, "|");
    fprintf_center(file, text, 21 * num_cols - 2);
    if (num_cols % 2 != 0) {
        fprintf(file, " ");
    }
    fprintf(file, " ");
    fprintf(file, "|\n");
}
void make_last_row(char flag, FILE* file, int num_cols) {
    if (!flag) { return; }
    for (int i = 0; i < num_cols; i++) {
        fprintf(file, "+%s", "--------------------");
    }
    fprintf(file, "+\n");
}
void make_table_row(char flag, FILE* file, int num_cols, char* text[]) {
    if (!flag) { return; }
    fprintf(file, "|%s", "--------------------");
    for (int i = 1; i < num_cols; i++) {
        fprintf(file, "+%s", "--------------------");
    }
    fprintf(file, "|\n");
    for (int i = 0; i < num_cols; i++) {
        fprintf(file, "|");
        fprintf_center(file, text[i], 20);
    }
    fprintf(file, "|\n");
}

