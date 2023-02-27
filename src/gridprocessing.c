#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "utils.h"

/* Grid initialization */

/*
 * Function: igrid
 * ---------------
 *
 * Fill integer grid with 1
 *
 * grid: empty 3D grid
 * m: x grid units
 * n: y grid units
 * o: z grid units
 *
 */
int ***igrid(int m, int n, int o) {
  int i, j, k, ***A;

  A = (int ***)calloc(m, sizeof(int **));
  for (i = 0; i < m; i++) {
    A[i] = (int **)calloc(n, sizeof(int *));
    for (j = 0; j < n; j++) {
      A[i][j] = (int *)calloc(o, sizeof(int));
      for (k = 0; k < o; k++) {
        A[i][j][k] = 1;
      }
    }
  }
  return A;
}

/*
 * Function: dgrid
 * ---------------
 *
 * Fill double grid with 0.0
 *
 * grid: empty 3D grid
 * m: x grid units
 * n: y grid units
 * o: z grid units
 *
 */
double ***dgrid(int m, int n, int o) {
  int i, j, k;
  double ***M;

  M = (double ***)calloc(m, sizeof(double **));
  for (i = 0; i < m; i++) {
    M[i] = (double **)calloc(n, sizeof(double *));
    for (j = 0; j < n; j++) {
      M[i][j] = (double *)calloc(o, sizeof(double));
      for (k = 0; k < o; k++) {
        M[i][j][k] = 0.0;
      }
    }
  }
  return M;
}

/* Molecular representation */

/*
 * Function: check_protein_neighbours
 * ----------------------------------
 *
 * Checks if a cavity point on the grid is next to a protein point (0 or -2).
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
 *
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
 * Function: SAS
 * -------------
 *
 * Insert atoms with a probe addition inside a 3D grid, producing a Solvent
 * Accessible Surface (SAS).
 *
 * A: 3D grid
 * m: x grid units
 * n: y grid units
 * o: z grid units
 * h: 3D grid spacing (A)
 * probe: Probe size (A)
 * X1: x coordinate of P1
 * Y1: y coordinate of P1
 * Z1: z coordinate of P1
 *
 */
void SAS(int ***A, int m, int n, int o, double h, double probe, double X1,
         double Y1, double Z1) {

  /* Declare variables */
  int i, j, k, imax, jmax, kmax;
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

/*
 * Function: SES
 * -------------
 *
 * Adjust surface representation to Solvent Excluded Surface (SES).
 *
 * A: 3D grid
 * m: x grid units
 * n: y grid units
 * o: z grid units
 * h: 3D grid spacing (A)
 * probe: Probe size (A)
 *
 */
void SES(int ***A, int m, int n, int o, double h, double probe) {
  int i, j, k, i2, j2, k2, aux;
  double distance;

  /* Calculate sas limit in 3D grid units */
  aux = ceil(probe / h);

  /* Set number of processes in OpenMP */
  int ncores = omp_get_num_procs() - 1;
  omp_set_num_threads(ncores);
  omp_set_nested(1);

/* Create a parallel region */
#pragma omp parallel default(none), shared(A, m, n, o, aux, probe, h),         \
    private(i, j, k, i2, j2, k2, distance)
  {
/* Create a parallel loop, collapsing 3 loops inside 1, which will send values
 * to next loop before it ends */
#pragma omp for schedule(dynamic) collapse(3)
    /* Loop around the search box */
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        for (k = 0; k < o; k++) {

          /* If a given cavity point on the grid is next to a protein point, do
           * ... */
          if (A[i][j][k] == 1)
            if (check_protein_neighbours(A, i, j, k, m, n, o)) {
              /* Loop around space occupied by radius of atom from atom position
               */
              for (i2 = i - aux; i2 <= i + aux; i2++)
                for (j2 = j - aux; j2 <= j + aux; j2++)
                  for (k2 = k - aux; k2 <= k + aux; k2++) {

                    /* If point inside analysis box, do ... */
                    if (i2 > 0 && j2 > 0 && k2 > 0 && i2 < m && j2 < n &&
                        k2 < o) {

                      /* Get absolute distance between point inside radius from
                       * atom and atom point */
                      distance = sqrt(pow(i - i2, 2) + pow(j - j2, 2) +
                                      pow(k - k2, 2));
                      /* If distance inside radius and point is a cavity, do ...
                       */
                      if (distance < (probe / h))
                        if (A[i2][j2][k2] == 0)
                          /* Mark space occupied by a big probe size from
                           * protein surface */
                          A[i2][j2][k2] = -2;
                    }
                  }
            }
        }

/* Create a parallel loop, collapsing 3 loops inside 1, which will receive
 * values from the above loop before it ends */
#pragma omp for collapse(3)
    /* Loop around analysis box */
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        for (k = 0; k < o; k++) {
          /* Mark space occupied by a big probe size from protein surface */
          if (A[i][j][k] == -2)
            A[i][j][k] = 1;
        }
  }
}

/* Cavity detection (Probe In - Probe Out) */

/*
 * Function: subtract
 * ------------------
 *
 * Compare Probe In and Probe Out 3D grids to define biomolecular cavities,
 * selecting points where Probe In reached and Probe Out did not.
 *
 * A: Probe In 3D grid
 * S: Probe Out 3D grid
 * m: x grid units
 * n: y grid units
 * o: z grid units
 * h: 3D grid spacing (A)
 * removal_distance: Length to be removed from the cavity-bulk frontier (A)
 *
 */
void subtract(int ***A, int ***S, int m, int n, int o, double h,
              double removal_distance) {
  /* Declare variables */
  int i, j, k, i2, j2, k2, rd;

  rd = ceil(removal_distance / h);

  /* Set number of processes in OpenMP */
  int ncores = omp_get_num_procs() - 1;
  omp_set_num_threads(ncores);
  omp_set_nested(1);

/* Create a parallel region */
#pragma omp parallel default(none), shared(A, S, rd, m, n, o, i, j, k),        \
    private(j2, i2, k2)
  {
/* Create a parallel loop, collapsing 3 loops inside 1 and schedule dynamic
 * allocation of threads */
#pragma omp for schedule(dynamic) collapse(3)
    /* Loop around the search box */
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        for (k = 0; k < o; k++) {

          /* If point is a cavity, do ... */
          if (S[i][j][k]) {

            /*Loops around space occupied by probe from atom position*/
            // #pragma omp taskloop
            for (i2 = i - rd; i2 <= i + rd; i2++)
              for (j2 = j - rd; j2 <= j + rd; j2++)
                for (k2 = k - rd; k2 <= k + rd; k2++)
                  /*If inside box, do... */
                  if (i2 >= 0 && i2 < m && j2 >= 0 && j2 < n && k2 >= 0 &&
                      k2 < o)
                    if (A[i2][j2][k2] == 1)
                      /* Mark points where big probe passed in cavities in A */
                      A[i2][j2][k2] = -1;
          }
        }
  }
}

/*
 * Function: filter_noise
 * ----------------------
 *
 * Removes cavities points (1) surrounded by biomolecule (0) or medium (-1)
 * points.
 *
 * A: cavities 3D grid
 * m: x grid units
 * n: y grid units
 * o: z grid units
 *
 */
void filter_noise(int ***A, int m, int n, int o) {
  int i, j, k, contacts;

  /* Set number of processes in OpenMP */
  int ncores = omp_get_num_procs() - 1;
  omp_set_num_threads(ncores);
  omp_set_nested(1);

#pragma omp parallel default(none), shared(A, m, n, o),                        \
    private(contacts, i, j, k)
#pragma omp for collapse(3) schedule(static)
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < o; k++) {

        if (A[i][j][k] == 1) {

          /* Initialize counter */
          contacts = 0;

          /* Check if a protein point (0) or a medium point (-1) is next to a
           * cavity point (>=1) */
          if (i - 1 >= 0)
            if (A[i - 1][j][k] == 0 || A[i - 1][j][k] == -1)
              contacts++;
          if (i + 1 < m)
            if (A[i + 1][j][k] == 0 || A[i + 1][j][k] == -1)
              contacts++;
          if (j - 1 >= 0)
            if (A[i][j - 1][k] == 0 || A[i][j - 1][k] == -1)
              contacts++;
          if (j + 1 < n)
            if (A[i][j + 1][k] == 0 || A[i][j + 1][k] == -1)
              contacts++;
          if (k - 1 >= 0)
            if (A[i][j][k - 1] == 0 || A[i][j][k - 1] == -1)
              contacts++;
          if (k + 1 < o)
            if (A[i][j][k + 1] == 0 || A[i][j][k + 1] == -1)
              contacts++;

          /* Cavity point is a medium point */
          if (contacts == 6)
            A[i][j][k] = -1;
        }
      }
}

/* Ligand adjustment */

/*
 * Function: adjust2ligand
 * -----------------------
 *
 * Adjust cavities to a radius around atoms of a target ligand.
 *
 * A: cavities 3D grid
 * m: x grid units
 * n: y grid units
 * o: z grid units
 * h: Grid spacing (A)
 * limit: Radius value to limit a space around a ligand (A)
 * X1: x coordinate of P1
 * Y1: y coordinate of P1
 * Z1: z coordinate of P1
 *
 */
void adjust2ligand(int ***A, int m, int n, int o, double h, double limit,
                   double X1, double Y1, double Z1) {
  /* Declare variables */
  int i, j, k, inside, aux;
  double distance, x, y, z, xaux, yaux, zaux;
  atom *p;

  /* Loop around analysis box */
#pragma omp parallel default(none),                                            \
    shared(A, m, n, o, h, sina, sinb, cosa, cosb, limit, X1, Y1, Z1, p, v),    \
    private(inside, i, j, k, x, y, z, xaux, yaux, zaux, distance)
  {
#pragma omp for collapse(3) schedule(static) // nowait
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        for (k = 0; k < o; k++) {
          inside = 0;

          /* Loop around ligand information linked list */
          for (p = v; p != NULL && !inside; p = p->next) {

            /* Get ligand atom coordinates */
            x = i * h;
            y = j * h;
            z = k * h;
            xaux = x * cosb + y * sina * sinb - z * cosa * sinb;
            yaux = y * cosa + z * sina;
            zaux = x * sinb - y * sina * cosb + z * cosa * cosb;
            xaux += X1;
            yaux += Y1;
            zaux += Z1;

            /* Get distance between protein and grid point inside box */
            distance = sqrt(pow(xaux - p->x, 2) + pow(yaux - p->y, 2) +
                            pow(zaux - p->z, 2));

            /* Mark Point (i,j,k) is inside limited region */
            if (distance < limit)
              inside = 1;
          }

          /* Cavity point is not inside ligand search space */
          if (inside == 0 && A[i][j][k])
            A[i][j][k] = -1;
        }
  }
}

/* Box adjustment */

/*
 * Function: filter2box
 * --------------------
 *
 * Adjust cavities to a search box.
 *
 * A: cavities 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * P1: xyz coordinates of 3D grid origin
 * ndims: number of coordinates (3: xyz)
 * P2: xyz coordinates of x-axis vertice
 * nndims: number of coordinates (3: xyz)
 * sincos: sin and cos of 3D grid angles
 * nvalues: number of sin and cos (sina, cosa, sinb, cosb)
 * step: 3D grid spacing (A)
 * probe_out: Radius value to limit a space around a ligand (A)
 * nthreads: number of threads for OpenMP
 *
 */
/* Filter search space based on box adjusment mode.
Analyze if points are outside the user defined search space and points outside
it are excluded */
void filter2box(int ***A, int m, int n, int o, double h, double bX1, double bY1,
                double bZ1, double bX2, double bY2, double bZ2, double norm1) {
  /* Declare variables */
  int i, j, k;
  double aux, normB;

  /* Set number of processes in OpenMP */
  int ncores = omp_get_num_procs() - 1;
  omp_set_num_threads(ncores);
  omp_set_nested(1);

  normB = sqrt(pow(bX2 - bX1, 2) + pow(bY2 - bY1, 2) + pow(bZ2 - bZ1, 2));
  aux = floor(norm1 - normB) / (2 * h);

/* Create a parallel region */
#pragma omp parallel default(shared), private(i, j, k)
  {

    for (i = 0; i <= aux; i++)
/* Create a parallel loop, collapsing 2 loops inside 1, which will send values
 * to next loop before it ends */
#pragma omp for collapse(2) nowait
      for (j = 0; j < n; j++)
        for (k = 0; k < o; k++)
          A[i][j][k] = -1;

    for (i = m - 1; i >= m - aux - 1; i--)
/* Create a parallel loop, collapsing 2 loops inside 1, which will send values
 * to next loop before it ends */
#pragma omp for collapse(2) nowait
      for (j = 0; j < n; j++)
        for (k = 0; k < o; k++)
          A[i][j][k] = -1;

    for (j = 0; j <= aux; j++)
/* Create a parallel loop, collapsing 2 loops inside 1, which will send values
 * to next loop before it ends */
#pragma omp for collapse(2) nowait
      for (i = 0; i < m; i++)
        for (k = 0; k < o; k++)
          A[i][j][k] = -1;

    for (j = n - 1; j >= n - aux - 1; j--)
/* Create a parallel loop, collapsing 2 loops inside 1, which will send values
 * to next loop before it ends */
#pragma omp for collapse(2) nowait
      for (i = 0; i < m; i++)
        for (k = 0; k < o; k++)
          A[i][j][k] = -1;

    for (k = 0; k <= aux; k++)
/* Create a parallel loop, collapsing 2 loops inside 1, which will send values
 * to next loop before it ends */
#pragma omp for collapse(2) nowait
      for (j = 0; j < n; j++)
        for (i = 0; i < m; i++)
          A[i][j][k] = -1;

    for (k = o - 1; k >= o - aux - 1; k--)
/* Create a parallel loop, collapsing 2 loops inside 1, which will send values
 * to next loop before it ends */
#pragma omp for collapse(2) nowait
      for (j = 0; j < n; j++)
        for (i = 0; i < m; i++)
          A[i][j][k] = -1;
  }
}

/* Cavity clustering and volume estimation */

/*
 * Function: check_unclustered_neighbours
 * --------------------------------------
 *
 * Checks if a cavity point on the grid is next to a unclustered cavity point
 * (1)
 *
 * A: 3D grid
 * dx: x grid units
 * dy: y grid units
 * dz: z grid units
 * i: x coordinate of cavity point
 * j: y coordinate of cavity point
 * k: z coordinate of cavity point
 *
 * returns: true (int 1) or false (int 0)
 */
int check_unclustered_neighbours(int ***A, int m, int n, int o, int i, int j,
                                 int k) {
  int x, y, z;

  // Loop around neighboring points
  for (x = i - 1; x <= i + 1; x++)
    for (y = j - 1; y <= j + 1; y++)
      for (z = k - 1; z <= k + 1; z++) {
        // Check if point is inside 3D grid
        if (x < 0 || y < 0 || z < 0 || x > m - 1 || y > n - 1 || z > o - 1)
          ;
        else if (A[x][y][z] > 1)
          return A[x][y][z];
      }

  return 0;
}

/*
 * Function: remove_cavity
 * -----------------------
 *
 * Untag cavity that does not reach volume cutoff.
 *
 * A: cavities 3D grid
 * m: x grid units
 * n: y grid units
 * o: z grid units
 * tag: cavity integer identifier
 *
 */
void remove_cavity(int ***A, int m, int n, int o, int tag) {
  int i, j, k;

  /* Set number of processes in OpenMP */
  int ncores = omp_get_num_procs() - 1;
  omp_set_num_threads(ncores);
  omp_set_nested(1);

/*Create a parallel region*/
#pragma omp parallel default(shared)
#pragma omp for collapse(3)
  /* Loop around the search box */
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < o; k++)
        /*Remove tag*/
        if (A[i][j][k] == tag)
          A[i][j][k] = -1;
}

/*
 * Function: DFS
 * -------------
 *
 * Recursive Depth-First Search (DFS) algorithm.
 *
 * A: cavities 3D grid
 * m: x grid units
 * n: y grid units
 * o: z grid units
 * i: x coordinate of cavity point
 * j: y coordinate of cavity point
 * k: z coordinate of cavity point
 * tag: cavity integer identifier
 *
 */
void DFS(int ***A, int m, int n, int o, int i, int j, int k, int tag) {
  int x, y, z;

  /* Ignore points in border */
  if (i == 0 || i == m - 1 || j == 0 || j == n - 1 || k == 0 || k == o - 1)
    return;

  /* If point is a cavity point, do ... */
  if (A[i][j][k] == 1 && !big) {
    A[i][j][k] = tag;
    volume++;

    /* Split big cavities */
    if (volume == 10000)
      big = 1;

    /* Loop around one position in each direction from starting point */
    if (!big) {
      for (x = i - 1; x <= i + 1; x++)
        for (y = j - 1; y <= j + 1; y++)
          for (z = k - 1; z <= k + 1; z++)
            /* Recursive call */
            DFS(A, m, n, o, x, y, z, tag);
    }
  }
}

/*
 * Function: cluster
 * -----------------
 *
 * Cluster consecutive cavity points together, by applying Depth-First Search in
 * accessible cavity points (1; nodes). During clustering, it calculates volume
 * based on cavity points with the same numeric tag.
 *
 * NOTE: Due to memory restrictions, the recursion is divided for big cavities.
 *
 * A: cavities 3D grid
 * m: x grid units
 * n: y grid units
 * o: z grid units
 * h: 3D grid spacing (A)
 * volume_cutoff: Cavities volume filter (A3)
 *
 */
int clustering(int ***A, int m, int n, int o, double h, double volume_cutoff) {
  /* Declare variables */
  int i, j, k, i2, j2, k2, tag, volume_aux;
  node *p;

  /* Initialize variables*/
  big = 0;
  volume_aux = 0;
  tag = 1;

  /* NULL pointer to last item in volume linked list */
  V = NULL;

  /*Loops around analysis box*/
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < o; k++)

        /* If point is a cavity point, do ... */
        if (A[i][j][k] == 1) {
          tag++;
          volume = 0;

          /* Call DFS algorithm */
          DFS(A, m, n, o, i, j, k, tag);
          volume_aux = volume;

          /* Loop for big cavities */
          while (big) {

            volume_aux = 0;
            /* Loop around search box */
            for (i2 = 0; i2 < m; i2++)
              for (j2 = 0; j2 < n; j2++)
                for (k2 = 0; k2 < o; k2++) {
                  big = 0;
                  volume_aux += volume;
                  volume = 0;
                  /* For a given identified cavity point, check if there is a
                  unidentified cavity point around it */
                  if (A[i2][j2][k2] == 1 && check_unclustered_neighbours(
                                                A, m, n, o, i2, j2, k2) == tag)
                    /* Call DFS algorithm */
                    DFS(A, m, n, o, i2, j2, k2, tag);
                }
          }
          /* Cavity volume */
          volume = volume_aux;

          /* If volume is less than threshold, remove cavity */
          if ((double)volume * pow(h, 3) < volume_cutoff) {
            remove_cavity(A, m, n, o, tag);
            tag--;
          } else {
            /* Append item to volume linked list */
            p = malloc(sizeof(node));
            p->volume = (double)volume * pow(h, 3);
            p->pos = tag - 2;

            if (V != NULL)
              p->next = V;
            else
              p->next = NULL;
            V = p;
          }
        }

  /* Return number of cavities */
  return tag - 1;
}

/* Cavity surface and area estimation */

/*
 * Function: define_surface_points
 * -------------------------------
 *
 * Identify surface points based on neighboring points.
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
 *
 */
int define_surface_points(int ***A, int m, int n, int o, int i, int j, int k) {

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

/*
 * Function: filter_surface
 * ------------------------
 *
 * Inspect cavities 3D grid and mark detected surface points on a surface 3D
 * grid.
 *
 * A: cavities 3D grid
 * S: surface points 3D grid
 * m: x grid units
 * n: y grid units
 * o: z grid units
 *
 */
void filter_surface(int ***A, int ***S, int m, int n, int o) {
  int i, j, k;

  /* Set number of processes in OpenMP */
  int ncores = omp_get_num_procs() - 1;
  omp_set_num_threads(ncores);
  omp_set_nested(1);

#pragma omp parallel default(none), shared(A, S, m, n, o), private(i, j, k)
  {
#pragma omp for collapse(3) schedule(static)
    /* Loop around the search box */
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        for (k = 0; k < o; k++)
          if (A[i][j][k] > 1) {
            /* Define surface cavity points (tag) */
            S[i][j][k] = define_surface_points(A, m, n, o, i, j, k);
          } else {
            /* Define protein point (0) */
            if (A[i][j][k] == 0)
              S[i][j][k] = 0;
            /* Define medium point (-1) */
            else
              S[i][j][k] = -1;
          }
  }
}

/*
 * Function: check_voxel_class
 * ---------------------------
 *
 * Identify voxel class of surface voxel and return class weight.
 *
 * S: surface points 3D grid
 * i: x coordinate of cavity point
 * j: y coordinate of cavity point
 * k: z coordinate of cavity point
 *
 * returns: voxel class weight (double)
 */
double check_voxel_class(int ***S, int i, int j, int k) {
  int contacts = 0;
  double weight = 1.0;

  /* If face accessible to protein point, increment contacts */
  if (S[i - 1][j][k] == 0)
    contacts++;
  if (S[i + 1][j][k] == 0)
    contacts++;
  if (S[i][j - 1][k] == 0)
    contacts++;
  if (S[i][j + 1][k] == 0)
    contacts++;
  if (S[i][j][k - 1] == 0)
    contacts++;
  if (S[i][j][k + 1] == 0)
    contacts++;

  /* Attribute weight based on voxel class */
  switch (contacts) {
    /*One face accessible to protein*/
  case 1:
    weight = 0.894;
    break;

  /*Two faces accessible to protein*/
  case 2:
    weight = 1.3409;
    break;

  /*Three non-consecutive faces accessible to protein*/
  case 3:
    if ((S[i + 1][j][k] == 0 && S[i - 1][j][k] == 0) ||
        (S[i][j + 1][k] == 0 && S[i][j - 1][k] == 0) ||
        (S[i][j][k + 1] == 0 && S[i][j][k - 1] == 0)) {
      weight = 2;
    } else {
      weight = 1.5879;
    }
    break;

    /* Four consecutive faces accessible to protein */
  case 4:
    weight = 2.6667;
    break;

    /* Five consecutive faces accessible to protein */
  case 5:
    weight = 3.3333;
    break;
  }

  return weight;
}

/*
 * Function: area
 * --------------
 *
 * Calculate area of cavities, using Mullikin and Verbeek method.
 *
 * S: surface points 3D grid
 * m: x grid units
 * n: y grid units
 * o: z grid units
 * h: 3D grid spacing (A)
 * ncav: number of cavities
 *
 */
void area(int ***S, int m, int n, int o, double h, int ncav) {
  /* Declare variables */
  int i, j, k;
  double *area;

  /* Set number of processes in OpenMP */
  int ncores = omp_get_num_procs() - 1;
  omp_set_num_threads(ncores);
  omp_set_nested(1);

  /* Initialize area object */
  area = (double *)calloc(ncav, sizeof(double));
  for (i = 0; i < ncav; i++)
    area[i] = 0.0;

/* Create a parallel loop and schedule dynamic allocation of threads */
#pragma omp parallel for shared (S, i, j, k, m, n, o, h) schedule(dynamic) reduction (+: area[:ncav])
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < o; k++) {
        if (S[i][j][k] > 1)
          area[S[i][j][k] - 2] += check_voxel_class(S, i, j, k) * pow(h, 2);
      }

  /* Save area in KVFinder results struct */
  for (i = 0; i < ncav; i++) {
    KVFinder_results[i].area = area[i];
  }

  /* Free area object from memory */
  free(area);
}

/* Constituional characterization */

/*
 * Function: create
 * ----------------
 *
 * Create a res node
 *
 * resnumber: residue number
 * resname: residue name
 * chain: chain identifier
 *
 * returns: residues_info node with residue information (residue number, residue
 * name, chain identier)
 */
residues_info *_create_residue(int resnumber, char resname, char chain) {
  residues_info *new = (residues_info *)malloc(sizeof(residues_info));

  new->resnumber = resnumber;
  new->resname = resname;
  new->chain = chain;
  new->next = NULL;

  return new;
}

/*
 * Function: _insert_residue
 * -------------------------
 *
 * Insert residues_info node in linked list
 *
 * head: pointer to linked list head
 * new: residues_info node
 *
 */
void _insert_residue(residues_info **head, residues_info *new) {
  residues_info *current;

  if (*head == NULL || (*head)->resnumber > new->resnumber) {
    new->next = *head;
    *head = new;
  } else {
    current = *head;
    while (current->next != NULL && current->next->resnumber < new->resnumber) {
      current = current->next;
    }
    new->next = current->next;
    current->next = new;
  }
}

/*
 * Function: _remove_duplicate_residue
 * -----------------------------------
 *
 * Remove duplicated residues_info nodes from linked list
 *
 * head: pointer to linked list head
 *
 */
void _remove_duplicate_residue(residues_info *head) {
  /* Pointer to traverse the linked list */
  residues_info *current = head;

  /* Pointer to store the next pointer of a node to be deleted*/
  residues_info *next_next;

  /* do nothing if the list is empty */
  if (current == NULL)
    return;

  /* Traverse the list till last node */
  while (current->next != NULL) {
    /* Compare current node with next node */
    if (current->resnumber == current->next->resnumber &&
        current->chain == current->next->chain) {
      /* The sequence of steps is important*/
      next_next = current->next->next;
      free(current->next);
      current->next = next_next;
    } else /* This is tricky: only advance if no deletion */
    {
      current = current->next;
    }
  }
}

/*
 * Function: interface
 * -------------------
 *
 * Retrieve interface residues surrounding cavities.
 *
 * A: cavities 3D grid
 * m: x grid units
 * n: y grid units
 * o: z grid units
 * h: 3D grid spacing (A)
 * probe_in: Probe In size (A)
 * ncav: number of cavities
 * X1: x coordinate of P1
 * Y1: y coordinate of P1
 * Z1: z coordinate of P1
 *
 */
void interface(int ***A, int m, int n, int o, double h, double probe, int ncav,
               double X1, double Y1, double Z1) {
  int i, j, k, imax, jmax, kmax, tag, old_num = -1, old_tag = -1;
  double x, y, z, xaux, yaux, zaux, distance, H;
  atom *p;
  residues_info *new;

  /* Loop around PDB linked list */
  for (p = v; p != NULL; p = p->next) {

    /* Standardize each position */
    x = (p->x - X1) / h;
    y = (p->y - Y1) / h;
    z = (p->z - Z1) / h;
    xaux = x * cosb + z * sinb;
    yaux = y;
    zaux = (-x) * sinb + z * cosb;
    x = xaux;
    y = yaux * cosa - zaux * sina;
    z = yaux * sina + zaux * cosa;

    /* Create a variable for space occupied by probe and radius of atom */
    H = (probe + p->radius) / h;

    /* Loop around space occupied by probe and radius of atom from atom
     * position */
    imax = ceil(x + H);
    jmax = ceil(y + H);
    kmax = ceil(z + H);
#pragma omp parallel default(none),                                            \
    shared(p, KVFinder_results, H, x, y, z, m, n, o, A, imax, jmax, kmax),     \
    private(i, j, k, distance, tag, new, old_tag, old_num)
#pragma omp for collapse(3)
    for (i = floor(x - H); i <= imax; i++)
      for (j = floor(y - H); j <= jmax; j++)
        for (k = floor(z - H); k <= kmax; k++)
          /* If inside box, do ... */
          if (i < m && i > 0 && j < n && j > 0 && k < o && k > 0) {
            if (abs(A[i][j][k]) > 1) {
              tag = A[i][j][k] - 2;
              distance = sqrt(pow(i - x, 2) + pow(j - y, 2) + pow(k - z, 2));
              if (distance <= H) {
                if (old_num != p->resnumber || old_tag != tag) {
                  new = _create_residue(p->resnumber, p->resname, p->chain);
                  _insert_residue(&KVFinder_results[tag].res_info, new);
                }
                old_num = p->resnumber;
                old_tag = tag;
              }
            }
          }
  }

  /* Remove duplicates */
  for (i = 0; i < ncav; i++) {
    _remove_duplicate_residue(KVFinder_results[i].res_info);
  }
}

/* Cavity boundary and depth estimation */

int define_boundary_points(int ***A, int m, int n, int o, int i, int j, int k) {
  if (i - 1 >= 0)
    if (A[i - 1][j][k] == -1)
      return -(A[i][j][k]);
  if (i + 1 < m)
    if (A[i + 1][j][k] == -1)
      return -(A[i][j][k]);
  if (j - 1 >= 0)
    if (A[i][j - 1][k] == -1)
      return -(A[i][j][k]);
  if (j + 1 < n)
    if (A[i][j + 1][k] == -1)
      return -(A[i][j][k]);
  if (k - 1 >= 0)
    if (A[i][j][k - 1] == -1)
      return -(A[i][j][k]);
  if (k + 1 < o)
    if (A[i][j][k + 1] == -1)
      return -(A[i][j][k]);

  return A[i][j][k];
}

void filter_boundary(int ***A, int m, int n, int o, int ncav) {
  int i, j, k, tag;

  // Set number of threads in OpenMP
  int ncores = omp_get_num_procs();
  omp_set_num_threads(ncores);
  omp_set_nested(1);

  for (i = 0; i < ncav; i++) {
    cavity[i].Xmin = m;
    cavity[i].Xmax = 0;
    cavity[i].Ymin = n;
    cavity[i].Ymax = 0;
    cavity[i].Zmin = o;
    cavity[i].Zmax = 0;

    /*Start min and max frontier coordinates*/
    boundary[i].Xmin = m;
    boundary[i].Xmax = 0;
    boundary[i].Ymin = n;
    boundary[i].Ymax = 0;
    boundary[i].Zmin = o;
    boundary[i].Zmax = 0;
  }

#pragma omp parallel default(none), shared(A, m, n, o, cavity, boundary),      \
    private(i, j, k, tag)
  {
#pragma omp for collapse(3) schedule(static)
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        for (k = 0; k < o; k++)
#pragma omp critical
          if (A[i][j][k] > 1) {
            // Get cavity identifier
            tag = A[i][j][k] - 2;

            // Get min and max coordinates of each cavity
            cavity[tag].Xmin = min(cavity[tag].Xmin, i);
            cavity[tag].Ymin = min(cavity[tag].Ymin, j);
            cavity[tag].Zmin = min(cavity[tag].Zmin, k);
            cavity[tag].Xmax = max(cavity[tag].Xmax, i);
            cavity[tag].Ymax = max(cavity[tag].Ymax, j);
            cavity[tag].Zmax = max(cavity[tag].Zmax, k);

            // Define cavity-bulk boundary points
            A[i][j][k] = define_boundary_points(A, m, n, o, i, j, k);

            // Get min and max coordinates of each cavity-bulk boundary
            if (A[i][j][k] < -1) {
              boundary[tag].Xmin = min(boundary[tag].Xmin, i);
              boundary[tag].Ymin = min(boundary[tag].Ymin, j);
              boundary[tag].Zmin = min(boundary[tag].Zmin, k);
              boundary[tag].Xmax = max(boundary[tag].Xmax, i);
              boundary[tag].Ymax = max(boundary[tag].Ymax, j);
              boundary[tag].Zmax = max(boundary[tag].Zmax, k);
            }
          }
  }
}

void remove_boundary(int ***A, int m, int n, int o, int ncav) {
  int i, j, k, tag;

  // Set number of threads in OpenMP
  int ncores = omp_get_num_procs();
  omp_set_num_threads(ncores);
  omp_set_nested(1);

#pragma omp parallel default(none), shared(A, boundary, ncav, m, n, o),        \
    private(tag, i, j, k)
#pragma omp for schedule(dynamic)
  for (tag = 0; tag < ncav; tag++)
    for (i = boundary[tag].Xmin; i <= boundary[tag].Xmax; i++)
      for (j = boundary[tag].Ymin; j <= boundary[tag].Ymax; j++)
        for (k = boundary[tag].Zmin; k <= boundary[tag].Zmax; k++)
          if (A[i][j][k] < -1)
            // Untag cavity-bulk boundary points
            A[i][j][k] = abs(A[i][j][k]);
}

void depth(int ***A, double ***M, int m, int n, int o, double h, int ncav) {
  int i, j, k, i2, j2, k2, count, tag;
  double distance, tmp;

  // Set number of threads in OpenMP
  int ncores = omp_get_num_procs();
  omp_set_num_threads(ncores);
  omp_set_nested(1);

#pragma omp parallel default(none),                                            \
    shared(A, M, m, n, o, h, ncav, cavity, boundary, KVFinder_results),        \
    private(tmp, tag, i, j, k, i2, j2, k2, distance, count)
  {
#pragma omp for schedule(dynamic)
    for (tag = 0; tag < ncav; tag++) {
      KVFinder_results[tag].max_depth = 0.0;
      KVFinder_results[tag].avg_depth = 0.0;
      count = 0;

      for (i = cavity[tag].Xmin; i <= cavity[tag].Xmax; i++)
        for (j = cavity[tag].Ymin; j <= cavity[tag].Ymax; j++)
          for (k = cavity[tag].Zmin; k <= cavity[tag].Zmax; k++)
            if (abs(A[i][j][k]) == (tag + 2)) {
              tmp = sqrt(pow(m, 2) + pow(n, 2) + pow(o, 2)) * h;
              count++;

              if (boundary[tag].Xmin == m && boundary[tag].Ymin == n &&
                  boundary[tag].Zmin == o && boundary[tag].Xmax == 0.0 &&
                  boundary[tag].Ymax == 0.0 && boundary[tag].Zmax == 0.0) {
                // Cavity without boundary (void)
                tmp = 0.0;
              } else {
                for (i2 = boundary[tag].Xmin; i2 <= boundary[tag].Xmax; i2++)
                  for (j2 = boundary[tag].Ymin; j2 <= boundary[tag].Ymax; j2++)
                    for (k2 = boundary[tag].Zmin; k2 <= boundary[tag].Zmax;
                         k2++)
                      if (A[i2][j2][k2] == -(tag + 2)) {
                        distance = sqrt(pow(i2 - i, 2) + pow(j2 - j, 2) +
                                        pow(k2 - k, 2)) *
                                   h;
                        if (distance < tmp)
                          tmp = distance;
                      }
              }

              // Save depth for cavity point
              M[i][j][k] = tmp;

              // Save maximum depth for cavity tag
              if (tmp > KVFinder_results[tag].max_depth)
                KVFinder_results[tag].max_depth = tmp;

              // Add cavity point depth to average depth for cavity tag
              KVFinder_results[tag].avg_depth += tmp;
            }
      // Divide sum of depths by number of cavity points for cavity tag
      KVFinder_results[tag].avg_depth /= count;
    }
  }

  /* Remove boundary */
  remove_boundary(A, m, n, o, ncav);
}

/* Export cavity PDB file */

/*
 * Function: _filter_cavity
 * ------------------------
 *
 * Filter cavity points to cavity exportion.
 *
 * A: cavities 3D grid
 * m: x grid units (cavities)
 * n: y grid units (cavities)
 * o: z grid units (cavities)
 * i: x coordinate of cavity point
 * j: y coordinate of cavity point
 * k: z coordinate of cavity point
 *
 */
int _filter_cavity(int ***A, int m, int n, int o, int i, int j, int k) {
  /* Declare variables */
  int a, b, c;

  /* Loop around one position in each direction from starting point */
  for (a = i - 1; a <= i + 1; a++)
    for (b = j - 1; b <= j + 1; b++)
      for (c = k - 1; c <= k + 1; c++) {

        /*If point inside analysis box, do...*/
        if (a < 0 || b < 0 || c < 0 || a > m - 1 || b > n - 1 || c > o - 1)
          ;
        else
          /*If point next to a protein point, return tag number*/
          if (A[a][b][c] == 0)
            return A[i][j][k];
      }

  return 0;
}

/*
 * Function: export
 * ----------------
 *
 * Export cavities to PDB file.
 *
 * output_pdb: cavity PDB filename
 * A: cavities 3D grid
 * S: surface points 3D grid
 * M: b-factor 3D grid (depths)
 * m: x grid units (cavities)
 * n: y grid units (cavities)
 * o: z grid units (cavities)
 * h: Grid spacing (A)
 * ncav: number of cavities
 * X1: x coordinate of P1
 * Y1: y coordinate of P1
 * Z1: z coordinate of P1
 *
 */
void export(char *output_pdb, int ***A, int ***S, double ***M, double ***HP,
            int kvp_mode, int m, int n, int o, double h, int ncav, double X1,
            double Y1, double Z1) {
  /* Declare variables */
  int i, j, k, count, tag;
  double x, y, z, xaux, yaux, zaux;
  FILE *output;

  // Set number of threads in OpenMP
  int ncores = omp_get_num_procs();
  omp_set_num_threads(ncores);
  omp_set_nested(1);

  /* Open output PDB file (<PDB>.KVFinder.output.pdb) */
  output = fopen(output_pdb, "w");
  fprintf(output, "MODEL     %4.d\n", 1);

  for (count = 1, tag = 2; tag <= ncav + 2; tag++)
#pragma omp parallel default(none)                                             \
    shared(A, S, M, HP, sina, sinb, cosa, cosb, h, ncav, tag, count, m, n, o,  \
           output, kvp_mode, X1, Y1, Z1),                                      \
    private(i, j, k, x, y, z, xaux, yaux, zaux)
  {
#pragma omp for schedule(static) collapse(3) ordered nowait
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        for (k = 0; k < o; k++) {
          // Check if cavity point with value tag
          if (A[i][j][k] == tag) {
            // Convert 3D grid coordinates to real coordinates
            x = i * h;
            y = j * h;
            z = k * h;
            xaux = x * cosb + y * sina * sinb - z * cosa * sinb;
            yaux = y * cosa + z * sina;
            zaux = x * sinb - y * sina * cosb + z * cosa * cosb;
            xaux += X1;
            yaux += Y1;
            zaux += Z1;

            /* Save cavity point coordinates */
#pragma omp critical
            if (S[i][j][k] == tag) {

              /* Write each cavity point */
              fprintf(output,
                      "ATOM  %5.d  HA  K%c%c   259    %8.3lf%8.3lf%8.3lf"
                      "%6.2lf%6.2lf\n",
                      count % 100000, 65 + (((S[i][j][k] - 2) / 26) % 26),
                      65 + ((S[i][j][k] - 2) % 26), xaux, yaux, zaux,
                      HP[i][j][k], M[i][j][k]);

            } else {
              if (kvp_mode)
                fprintf(output,
                        "ATOM  %5.d  H   K%c%c   259    %8.3lf%8.3lf%8.3lf"
                        "%6.2lf%6.2lf\n",
                        count % 100000,
                        65 + (((abs(A[i][j][k]) - 2) / 26) % 26),
                        65 + ((abs(A[i][j][k]) - 2) % 26), xaux, yaux, zaux,
                        HP[i][j][k], M[i][j][k]);
              else if (_filter_cavity(A, m, n, o, i, j, k) != 0)
                fprintf(output,
                        "ATOM  %5.d  H   K%c%c   259    %8.3lf%8.3lf%8.3lf"
                        "%6.2lf%6.2lf\n",
                        count % 100000,
                        65 + (((abs(A[i][j][k]) - 2) / 26) % 26),
                        65 + ((abs(A[i][j][k]) - 2) % 26), xaux, yaux, zaux,
                        HP[i][j][k], M[i][j][k]);
            }
            count++;
          }
        }
  }

  fprintf(output, "END\n");
  fprintf(output, "ENDMDL\n");

  // Close file
  fclose(output);
}

/* Clean memory */

/*
 * Function: _free_igrid
 * ---------------------
 *
 * Free integer 3D grid.
 *
 * A: integer 3D grid
 * m: x grid units
 * n: y grid units
 * o: z grid units
 *
 */
void free_igrid(int ***A, int m, int n, int o) {
  int i, j;

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
      free(A[i][j]);
    free(A[i]);
  }

  free(A);
}

/*
 * Function: _free_dgrid
 * ---------------------
 *
 * Free integer 3D grid.
 *
 * M: double 3D grid
 * m: x grid units
 * n: y grid units
 * o: z grid units
 *
 */
void free_dgrid(double ***M, int m, int n, int o) {
  int i, j;

  for (i = 0; i < m; i++) {

    for (j = 0; j < n; j++)
      free(M[i][j]);
    free(M[i]);
  }

  free(M);
}

/*
 * Function: _free_dgrid
 * ---------------------
 *
 * Free linked list with nodes.
 *
 */
void free_node() {
  /* Declare variables */
  node *p;

  while (V != NULL) {
    p = V;
    V = V->next;
    free(p);
  }
}

/* Cavity hydropathy */
char *resn[] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
double scale[20] = {-0.64, 2.6,  0.8,   0.92,  -0.3,  0.87,  0.76,
                    -0.49, 0.41, -1.42, -1.09, 1.54,  -0.66, -1.22,
                    -0.12, 0.18, 0.05,  -0.83, -0.27, -1.11};

/*
 * Function: get_hydrophobicity_value
 * ----------------------------------
 *
 * Get hydrophobicity scale value for a target residue name
 *
 * resname: target residue name
 * resn: 1D-array hydrophobicity scale residues names
 * scale: 1D-array of hydrophocity scale values
 *
 */
double get_hydrophobicity_value(char *resname, char *resn[], double *scale) {
  int i;

  // Get hydrophobicity value
  for (i = 0; i < 20; i++)
    if (strcmp(resname, resn[i]) == 0)
      return scale[i];

  return 0.0;
}

/*
 * Function: project_hydropathy
 * ---------------------------
 *
 * Map a hydrophobicity scale per surface point of detected cavities.
 *
 * HP: hydrophobicity scale 3D grid
 * S: surface points 3D grid
 * m: x grid units
 * n: y grid units
 * o: z grid units
 * h: 3D grid spacing (A)
 * probe: Probe In size (A)
 * X1: x coordinate of P1
 * Y1: y coordinate of P1
 * Z1: z coordinate of P1
 *
 */
void project_hydropathy(double ***HP, int ***S, int m, int n, int o, double h,
                        double probe, double X1, double Y1, double Z1) {
  int i, j, k;
  double x, y, z, xaux, yaux, zaux, distance, H, ***ref;
  atom *p;

  // Initiliaze 3D grid for residues distances
  ref = dgrid(m, n, o);

  /* Loop around PDB linked list */
  for (p = v; p != NULL; p = p->next) {

    /* Standardize each position */
    x = (p->x - X1) / h;
    y = (p->y - Y1) / h;
    z = (p->z - Z1) / h;
    xaux = x * cosb + z * sinb;
    yaux = y;
    zaux = (-x) * sinb + z * cosb;
    x = xaux;
    y = yaux * cosa - zaux * sina;
    z = yaux * sina + zaux * cosa;

    /* Create a variable for space occupied by probe and radius of atom */
    H = (probe + p->radius) / h;

    for (i = floor(x - H); i <= ceil(x + H); i++)
      for (j = floor(y - H); j <= ceil(y + H); j++)
        for (k = floor(z - H); k <= ceil(z + H); k++) {
          if (i < m && i >= 0 && j < n && j >= 0 && k < o && k >= 0)
            // Found a surface point
            if (S[i][j][k] > 1) {
              // Calculate distance bewteen atom and surface point
              distance = sqrt(pow(i - x, 2) + pow(j - y, 2) + pow(k - z, 2));
              // Check if surface point was not checked before
              if (ref[i][j][k] == 0.0) {
                ref[i][j][k] = distance;
                HP[i][j][k] = get_hydrophobicity_value(
                    _code2residue(p->resname), resn, scale);
              }
              // Check if this atom is closer to the previous one assigned
              else if (ref[i][j][k] > distance) {
                ref[i][j][k] = distance;
                HP[i][j][k] = get_hydrophobicity_value(
                    _code2residue(p->resname), resn, scale);
              }
            }
        }
  }

  // Free 3D grid for residues distances
  free_dgrid(ref, m, n, o);
}

/*
 * Function: estimate_average_hydropathy
 * -------------------------------------
 *
 * Calculate average hydropathy of detected cavities.
 *
 * avgh: empty array of average hydropathy
 * ncav: number of cavities
 * hydropathy: hydrophobicity scale 3D grid
 * surface: surface points 3D grid
 * nx: x grid units
 * ny: y grid units
 * nz: z grid units
 * nthreads: number of threads for OpenMP
 *
 */
void estimate_average_hydropathy(double ***HP, int ***S, int m, int n, int o,
                                 int ncav) {
  int i, j, k, *pts;
  double *avgh;

  /* Set number of processes in OpenMP */
  int ncores = omp_get_num_procs() - 1;
  omp_set_num_threads(ncores);
  omp_set_nested(1);

  /* Initialize area object */
  avgh = (double *)calloc(ncav, sizeof(double));
  pts = (int *)calloc(ncav, sizeof(int));
  for (i = 0; i < ncav; i++) {
    pts[i] = 0;
    avgh[i] = 0.0;
  }

#pragma omp parallel default(none), shared(avgh, HP, S, pts, m, n, o),         \
    private(i, j, k)
  {
#pragma omp for collapse(3) ordered
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        for (k = 0; k < o; k++) {
#pragma omp critical
          if (S[i][j][k] > 1) {
            pts[S[i][j][k] - 2]++;
            avgh[S[i][j][k] - 2] += HP[i][j][k];
          }
        }
  }

  for (i = 0; i < ncav; i++)
    avgh[i] /= pts[i];

  /* Save area in KVFinder results struct */
  for (i = 0; i < ncav; i++) {
    KVFinder_results[i].avg_hydropathy = avgh[i];
  }

  // Free array with number of points per cavity
  free(pts);
  free(avgh);
}
