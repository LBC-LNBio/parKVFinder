#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "utils.h"

/*
 * Function: check_protein_neighbours
 * ----------------------------------
 *
 * Checks if a cavity point on the grid is next to a protein point (0 or -2)
 *
 * A: 3D grid
 * i: x coordinate of cavity point
 * j: y coordinate of cavity point
 * k: z coordinate of cavity point
 * m: x grid units
 * n: y grid units
 * o: z grid units
 *
 * returns: true (int 1) or false (int 0)
 */
int check_protein_neighbours(int ***A, int i, int j, int k, int m, int n,
                             int o) {
  int a, b, c;

  /* Loop around one position in each direction from starting point */
  for (a = i - 1; a <= i + 1; a++)
    for (b = j - 1; b <= j + 1; b++)
      for (c = k - 1; c <= k + 1; c++) {

        /* Check if point is inside 3D grid */
        if (a < 0 || b < 0 || c < 0 || a > m - 1 || b > n - 1 || c > o - 1)
          ;
        else
          /* If point next to a protein point, return True */
          if (A[a][b][c] == 0 || A[a][b][c] == -2)
            return 1;
      }

  return 0;
}

/*
 * Function: define_surface_points
 * -------------------------------
 *
 * Identify surface points based on neighboring points
 *
 * cavities: cavities 3D grid
 * i: x coordinate of cavity point
 * j: y coordinate of cavity point
 * k: z coordinate of cavity point
 * m: x grid units
 * n: y grid units
 * o: z grid units
 *
 * returns: cavity identifier (>1) or medium point (-1)
 */
int define_surface_points(int ***A, int i, int j, int k, int m, int n, int o) {

  /* Check if a protein point(0) is next to a cavity point(>=1) */
  if (i - 1 >= 0)
    if (A[i - 1][j][k] == 0)
      return A[i][j][k];
  if (i + 1 < m)
    if (A[i + 1][j][k] == 0)
      return A[i][j][k];
  if (j - 1 >= 0)
    if (A[i][j - 1][k] == 0)
      return A[i][j][k];
  if (j + 1 < n)
    if (A[i][j + 1][k] == 0)
      return A[i][j][k];
  if (k - 1 >= 0)
    if (A[i][j][k - 1] == 0)
      return A[i][j][k];
  if (k + 1 < o)
    if (A[i][j][k + 1] == 0)
      return A[i][j][k];

  return -1;
}

/* Mark the grid leaving a probe size around the protein */
void SAS(int ***A, int m, int n, int o, double h, double probe, double X1,
         double Y1, double Z1) {

  /* Declare variables */
  int i, j, k;
  double distance, x, y, z, xaux, yaux, zaux, H, x1, y1, z1;
  atom *p;

  /* Loop around PDB linked list */
  for (p = v; p != NULL; p = p->next) {

    /* Standardize each position */
    x1 = (p->x - X1) / h;
    y1 = (p->y - Y1) / h;
    z1 = (p->z - Z1) / h;
    xaux = x1 * cosb + z1 * sinb;
    yaux = y1;
    zaux = -x1 * sinb + z1 * cosb;
    x1 = xaux;
    y1 = yaux * cosa - zaux * sina;
    z1 = yaux * sina + zaux * cosa;

    /* Create a variable for space occupied by probe and radius of atom */
    H = (probe + p->radius) / h;

    /* Loop around space occupied by probe and radius of atom from atom position
     */
    int imax, jmax, kmax;
    imax = ceil(x1 + H);
    jmax = ceil(y1 + H);
    kmax = ceil(z1 + H);
#pragma omp parallel default(none),                                            \
    shared(H, x1, y1, z1, m, n, o, A, imax, jmax, kmax),                       \
    private(i, j, k, distance)
#pragma omp for collapse(3)
    for (i = floor(x1 - H); i <= imax; i++)
      for (j = floor(y1 - H); j <= jmax; j++)
        for (k = floor(z1 - H); k <= kmax; k++) {
          /* Get absolute distance between protein and grid point inside box */
          distance = sqrt(pow(i - x1, 2) + pow(j - y1, 2) + pow(k - z1, 2));
          /* Mark the grid with 0, leaving a probe size around the protein */
          if (distance < H)
            if (i >= 0 && i < m && j >= 0 && j < n && k >= 0 && k < o)
              A[i][j][k] = 0;
        }
  }
}
