#ifndef GRIDPROCESSING_H
#define GRIDPROCESSING_H

int check_protein_neighbours(int ***A, int i, int j, int k, int m, int n, int o);
int define_surface_points(int ***A, int i, int j, int k, int m, int n, int o);
void SAS(int ***A, int m, int n, int o, double h, double probe,
         double X1, double Y1, double Z1);

#endif