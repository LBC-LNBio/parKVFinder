/* This file contains all the grid functions utilized by KVFinder. The KVFinder grid data is represented by a matrix.
The most complex functions are briefly explained in the source code.*/

/* Import native modules */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <omp.h>

/* Import custom modules */
#include "dictionaryprocessing.h"
#include "pdbprocessing.h"
#include "resultsprocessing.h"
#include "matrixprocessing.h"

double
max (double a,
     double b)
{

	if (a > b)
	    return a;
	else
	    return b;

}

double
min (double a,
     double b)
{

	if (a < b)
	    return a;
	else
	    return b;

}

/* Check if a given cavity point on the grid is next to a protein point */
int
check_pos (int ***A,
           int i,
           int j,
           int k,
           int m,
           int n,
           int o)
{
    /* Declare variables */
	int a, b, c;

	/* Loop around one position in each direction from starting point */
	for (a = i-1; a <= i+1; a++)
        for (b = j-1; b <= j+1; b++)
			for (c = k-1; c <= k+1; c++) {

			    /* If point inside analysis box, do ... */
				if (a < 0 || b < 0 || c < 0 || a > m-1 || b > n-1 || c > o-1);
				else
				    /* If point next to a protein point, return True */
				    if (A[a][b][c] == 0 || A[a][b][c] == -2)
				        return 1;

			}

	return 0;

}

/* Check if a given cavity point on the grid is next to a protein point */
int
check_pos2 (int ***A,
            int i,
            int j,
            int k,
            int m,
            int n,
            int o)
{
    /* Declare variables */
	int a, b, c;

	/* Loop around one position in each direction from starting point */
	for (a = i-1; a <= i+1; a++)
		for (b = j-1; b <= j+1; b++)
			for (c = k-1; c <= k+1; c++) {

			    /*If point inside analysis box, do...*/
				if (a < 0 || b < 0 || c < 0 || a > m-1 || b > n-1 || c > o-1);
				else
				    /*If point next to a protein point, return tag number*/
				    if (A[a][b][c] == 0)
				        return A[i][j][k];

			}

	return 0;

}

/* When the cavities are grouped, each cavity is identified by an integer tag. This function checks, for a given
identified cavity point, if there is a undefined cavity point around it */
int
check_pos3 (int ***A,
            int i,
            int j,
            int k,
            int m,
            int n,
            int o)
{
	/* Declare variables */
	int a, b, c;

	/* Loop around one position in each direction from starting point */
	for (a = i-1; a <= i+1; a++)
		for (b = j-1; b <= j+1; b++)
			for (c = k-1; c <= k+1; c++) {

			    /* If point inside analysis box, do ... */
				if (a < 0 || b < 0 || c < 0 || a > m-1 || b > n-1 || c > o-1);
				else
				    /*If point next to a protein point, return tag number*/
				    if (A[a][b][c] > 1)
				        return A[a][b][c];

			}

	return 0;

}

/* Check if a given cavity point on the grid is next to a protein point */
int
check_pos4 (int ***A,
            int i,
            int j,
            int k,
            int m,
            int n,
            int o)
{

	/* Check if a protein point(0) is next to a cavity point(>=1) */
	if (i-1 >= 0)
	    if (A[i-1][j][k] == 0)
	        return A[i][j][k];
	if (i+1 < m)
	    if (A[i+1][j][k] == 0)
	        return A[i][j][k];
	if (j-1 >= 0)
	    if (A[i][j-1][k] == 0)
	        return A[i][j][k];
	if (j+1 < n)
	    if (A[i][j+1][k] == 0)
	        return A[i][j][k];
	if (k-1 >= 0)
	    if (A[i][j][k-1] == 0)
	        return A[i][j][k];
	if (k+1 < o)
	    if (A[i][j][k+1] == 0)
	        return A[i][j][k];

	return -1;

}

/* Check if a given cavity point on the grid is next to a medium point */
int
contact_protein_kv (int ***A,
                    int i,
                    int j,
                    int k,
                    int m,
                    int n,
                    int o)
{

	/* Check if a medium point(-1) is next to a cavity point(>=1) */
	if (i-1 >= 0)
	    if (A[i-1][j][k] == -1)
	        return -A[i][j][k];
	if (i+1 < m)
	    if (A[i+1][j][k] == -1)
	        return -A[i][j][k];
	if (j-1 >= 0)
	    if (A[i][j-1][k] == -1)
	        return -A[i][j][k];
	if (j+1 < n)
	    if (A[i][j+1][k] == -1)
	        return -A[i][j][k];
	if (k-1 >= 0)
	    if (A[i][j][k-1] == -1)
	        return -A[i][j][k];
	if (k+1 < o)
	    if (A[i][j][k+1] == -1)
	        return -A[i][j][k];

	return A[i][j][k];

}

/* Mark the grid leaving a probe size around the protein */
void
Matrix_fill (int ***A,
              int m,
              int n,
              int o,
              double h,
              double probe,
              double X1,
              double Y1,
              double Z1)
{

    /* Declare variables */
	int i, j, k;
	double distance, x, y, z, xaux, yaux, zaux, H, x1, y1, z1;
	atom *p;

	/* Loop around PDB linked list */
	for (p = v; p != NULL; p = p->next) {

		/* Standardize each position */
		x1 = (p->x-X1)/h; y1 = (p->y-Y1)/h; z1 = (p->z-Z1)/h;
		xaux = x1*cosb + z1*sinb; yaux = y1; zaux = -x1*sinb + z1*cosb;
		x1 = xaux; y1 = yaux*cosa - zaux*sina; z1 = yaux*sina + zaux*cosa;

        /* Create a variable for space occupied by probe and radius of atom */
		H = (probe+p->radius)/h;

		/* Loop around space occupied by probe and radius of atom from atom position */
		for (i = floor (x1-H); i <= ceil (x1+H); i++)
			for (j = floor (y1-H); j <= ceil (y1+H); j++)
				for (k = floor (z1-H)+1; k <= ceil (z1+H)+1; k++) {

					/* Get absolute distance between protein and grid point inside box */
					distance = sqrt ((i-x1)*(i-x1) + (j-y1)*(j-y1) + (k-z1)*(k-z1));
					/* Mark the grid with 0, leaving a probe size around the protein */
					if (distance < H && i >= 0 && i < m && j >= 0 && j < n && k >= 0 && k < o)
					    A[i][j][k] = 0;

				}

	}

}

/* Create a list of residues contacting each cavity */
void
Matrix_search (int ***A,
               int ***S,
               int m,
               int n,
               int o,
               double h,
               double probe,
               double X1,
               double Y1,
               double Z1,
               int ncav)
{
    /* Declare variables */
	int i, j, k, oldnum = 0, cont = 2;
	double xaux, yaux, zaux, H, x1, y1, z1, dist;
	atom *p;

	while (cont <= ncav+1) {

		/* Loop around PDB linked list */
		for(p = v; p != NULL; p = p->next){

			/* Standardize each position */
			x1 = (p->x-X1)/h; y1 = (p->y-Y1)/h; z1 = (p->z-Z1)/h;
			xaux = x1*cosb + z1*sinb; yaux = y1; zaux = -x1*sinb + z1*cosb;
			x1 = xaux; y1 = yaux*cosa - zaux*sina; z1 = yaux*sina + zaux*cosa;

			/* Create a variable for space occupied by probe and radius of atom */
			H = (probe+p->radius)/h;

			/* Loop around space occupied by probe and radius of atom from atom position */
			for (i = floor (x1-H); i <= ceil (x1+H); i++)
				for (j = floor (y1-H); j <= ceil (y1+H); j++)
					for (k = floor (z1-H); k <= ceil (z1+H); k++)
					    /* If inside box, do ... */
						if (i < m-1 && i > 0 && j < n-1 && j > 0 && k < o-1 && k > 0) {

                            /* If point contain tag number, do ... */
							if (A[i][j][k] == cont) {
								dist = sqrt (pow ((x1 - i), 2) + pow ((y1 - j), 2) + pow ((z1 - k), 2));
								if (dist <= H) {
									/*Residue number range from 1 to 9999 in a PDB file*/
									/*Letters range from 65 to 90 in ASCII table*/
									/*Residue number should be unique for specified cavity index (cont)*/
									if (oldnum != p->resnum && p->resnum <= 9999 && p->chain >= 65 && p->chain <= 90) {
										/* Insert residue information inside KVFinder_results */
										insert_res (p->resnum,
										           p->chain,
										           p->res_name,
										           cont-2);
									}
									oldnum = p->resnum;
								}
							}

						}

			}

        /* Increment cavity tag */
		cont++;

	}

}

/* Pass the probe around the surface of the protein */
void
Matrix_surf (int ***A,
             int m,
             int n,
             int o,
             double h,
             double radius)
{

    /* Declare variables */
	int i, j, k, i2, j2, k2, aux = ((int) radius/h) + 1;
	double norm;

    /* Set number of processes in OpenMP */
	int ncores = omp_get_num_procs () - 1;
	omp_set_num_threads (ncores);

/* Create a parallel region */
#pragma omp parallel default(none), shared(A,m,n,o,aux,radius,h), private(i,j,k,i2,j2,k2,norm)
{
/* Create a parallel loop, collapsing 3 loops inside 1, which will send values to next loop before it ends */
#pragma omp for schedule(dynamic) collapse(3) nowait
	/* Loop around the search box */
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < o; k++) {

                /* If a given cavity point on the grid is next to a protein point, do ... */
				if (A[i][j][k] == 1 && check_pos (A, i, j, k, m, n, o)) {

					/* Loop around space occupied by radius of atom from atom position */
					for (i2 = i-aux; i2 <= i+aux; i2++)
						for (j2 = j-aux; j2 <= j+aux; j2++)
							for (k2 = k-aux; k2 <= k+aux; k2++) {

		                        /* If point inside analysis box, do ... */
								if (i2 > 0 && j2 > 0 && k2 > 0 && i2 < m-1 && j2 < n-1 && k2 < o-1) {

								    /* Get absolute distance between point inside radius from atom and atom point */
									norm = sqrt ((i - i2)*(i - i2) + (j - j2)*(j - j2) + (k - k2)*(k - k2));
									/* If distance inside radius and point is a cavity, do ... */
									if (norm < (radius/h) && A[i2][j2][k2] == 0)
									    /* Mark space occupied by a big probe size from protein surface */
										A[i2][j2][k2] = -2;

								}

							}

				}

			}

/* Create a parallel loop, collapsing 3 loops inside 1, which will receive values from the above loop before it ends */
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

/* Mark the surface points on the grid */
void
Matrix_surface (int ***A,
                int ***S,
                int m,
                int n,
                int o,
                double h,
                double X1,
                double Y1,
                double Z1)
{
    /* Declare variables */
	int i, j, k, iterator;

	/* Loop around the search box */
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < o; k++) {

			    /* A[i][j][k] has an integer tag (cavity tag) */
				if (A[i][j][k] > 1) {

					/* Apply filters */
					/* Mark surface points */
					S[i][j][k] = check_pos4 (A, i, j, k, m, n, o);
					/* Mark frontier points between cavity and medium */
					A[i][j][k] = contact_protein_kv (A, i, j, k, m, n, o);

				}
				else {

				    /* A[i][j][k] is a protein point(0), pass a protein point to S[i][j][k] */
					if (A[i][j][k] == 0)
					    S[i][j][k] = 0;
					/* A[i][j][k] is not a protein or a cavity point, pass medium point to S[i][j][k] */
					else
					    S[i][j][k] = -1;

				}

			}

}

/* Adjust the search space for a given ligand */
void
Matrix_adjust (int ***A,
               int m,
               int n,
               int o,
               double h,
               double limit,
               double X1,
               double Y1,
               double Z1)
{
    /* Declare variables */
	int i, j, k, flag, aux;
	double distance, x, y, z, xaux, yaux, zaux;
	atom *p;

	/* Loop around analysis box */
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < o; k++) {

			flag = 0;

			/* Loop around ligand information linked list */
			for (p = v; p != NULL && !flag; p = p->next) {

				/* Get ligand atom coordinates */
				x = i*h; y = j*h; z = k*h;
				xaux = x*cosb + y*sina*sinb - z*cosa*sinb;
				yaux =          y*cosa      + z*sina;
				zaux = x*sinb - y*sina*cosb + z*cosa*cosb;
				xaux += X1; yaux += Y1; zaux += Z1;

				/* Get distance between protein and grid point inside box */
				distance = sqrt ((xaux-p->x)*(xaux-p->x) + (yaux-p->y)*(yaux-p->y) + (zaux-p->z)*(zaux-p->z));

                /* Mark Point (i,j,k) is inside limited region */
				if (distance < fabs (limit))
				    flag = 1;

			}

            /* Cavity point is not inside ligand search space */
			if (flag == 0 && A[i][j][k] && limit > 0.0)
				A[i][j][k] = -1;

            /*  Cavity point is not inside ligand search space */
			if (flag == 1 && A[i][j][k] && limit < 0.0)
				A[i][j][k] = -1;

			}

}

/* Filter search space based on box adjusment mode.
Analyze if points are outside the user defined search space and points outside it are excluded */
void
Matrix_filter (int ***A,
               int ***S,
               int m,
               int n,
               int o,
               double h,
               double bX1,
               double bY1,
               double bZ1,
               double bX2,
               double bY2,
               double bZ2,
               double norm1)
{
    /* Declare variables */
	int i, j, k;
	double aux, normB;

    /* Set number of processes in OpenMP */
	int ncores = omp_get_num_procs () - 1;
	omp_set_num_threads (ncores);

	normB = sqrt ((bX2-bX1) * (bX2-bX1) + (bY2-bY1) * (bY2-bY1) + (bZ2-bZ1) * (bZ2-bZ1));
	aux = (int)(norm1 - normB) / (2 * h);

/* Create a parallel region */
#pragma omp parallel default(shared), private(i,j,k)
{

	for (i = 0; i <= aux; i++)
/* Create a parallel loop, collapsing 2 loops inside 1, which will send values to next loop before it ends */
#pragma omp for collapse(2) nowait
		for (j = 0; j < n; j++)
			for (k = 0; k < o; k++)
			    A[i][j][k] = -1;
			    S[i][j][k] = -1;

	for (i = m-1; i >= m-aux-1; i--)
/* Create a parallel loop, collapsing 2 loops inside 1, which will send values to next loop before it ends */
#pragma omp for collapse(2) nowait
		for (j = 0; j < n; j++)
			for (k = 0; k < o; k++)
				A[i][j][k] = -1;
				S[i][j][k] = -1;

	for (j = 0; j <= aux; j++)
/* Create a parallel loop, collapsing 2 loops inside 1, which will send values to next loop before it ends */
#pragma omp for collapse(2) nowait
		for (i = 0; i < m; i++)
			for (k = 0; k < o; k++)
				A[i][j][k] = -1;
				S[i][j][k] = -1;

	for (j = n-1;j >= n-aux-1; j--)
/* Create a parallel loop, collapsing 2 loops inside 1, which will send values to next loop before it ends */
#pragma omp for collapse(2) nowait
		for (i = 0; i < m; i++)
			for (k = 0; k < o; k++)
				A[i][j][k] = -1;
				S[i][j][k] = -1;

	for (k = 0; k <= aux; k++)
/* Create a parallel loop, collapsing 2 loops inside 1, which will send values to next loop before it ends */
#pragma omp for collapse(2) nowait
		for (j = 0; j < n; j++)
			for (i = 0; i < m; i++)
				A[i][j][k] = -1;
				S[i][j][k] = -1;

	for (k = o-1; k >= o-aux-1; k--)
/* Create a parallel loop, collapsing 2 loops inside 1, which will send values to next loop before it ends */
#pragma omp for collapse(2) nowait
		for (j = 0; j < n; j++)
			for (i = 0; i < m; i++)
				A[i][j][k] = -1;
				S[i][j][k] = -1;

}

}

/* Export the cavity points to PDB file */
void
Matrix_export (int ***A,
               int ***S,
               int kvp_mode,
               int m,
               int n,
               int o,
               double h,
               int ncav,
               char *output_base_name,
               char *output_pdb,
               double X1,
               double Y1,
               double Z1)
{
    /* Declare variables */
	int i, j, k, iterator, count = 1, control = 1, tag = 2;
	double x, y, z, xaux, yaux, zaux, MAX = 0.0;
	FILE *output, *output_scale;

	/* Open output PDB file (<PDB>.KVFinder.output.pdb) */
	output = fopen (output_pdb, "w");

	while (control == 1) {
	    control = 0;

		/* Loop around grid */
		for (i = 0; i < m; i++)
			for (j = 0; j < n; j++)
				for (k = 0; k < o; k++) {

                    /* If grid point is a tag cavity point, do ... */
					if (abs (A[i][j][k]) == tag) {

						control = 1;

						/* Convert grid coordinates to real coordinates */
						x = i*h; y = j*h; z = k*h;
						xaux = x*cosb + y*sina*sinb - z*cosa*sinb;
						yaux =          y*cosa      + z*sina;
						zaux = x*sinb - y*sina*cosb + z*cosa*cosb;
						xaux += X1; yaux += Y1; zaux += Z1;

						/* Save cavity point coordinates */
						if (abs (S[i][j][k]) == tag) {

							/* Write each cavity point */
							fprintf (output,
							         "ATOM  %5.d  HS  K%c%c   259    %8.3lf%8.3lf%8.3lf  1.00%6.2lf\n",
							         count,
							         65+(((S[i][j][k]-2)/26)%26),
							         65+((S[i][j][k]-2)%26),
							         xaux,
							         yaux,
							         zaux,
							         0.0);

						}
						else {
							if (kvp_mode)
							    fprintf(output,
							            "ATOM  %5.d  H   K%c%c   259    %8.3lf%8.3lf%8.3lf  1.00%6.2lf\n",
							            count,
							            65 + (((abs(A[i][j][k]) - 2) / 26) % 26),
							            65 + ((abs(A[i][j][k]) - 2) % 26),
							            xaux,
							            yaux,
							            zaux,
							            0.0);
							else
							    if (check_pos2(A, i, j, k, m, n, o) != 0)
							        fprintf(output,
							                "ATOM  %5.d  H   K%c%c   259    %8.3lf%8.3lf%8.3lf  1.00%6.2lf\n",
							                count,
							                65 + (((abs(A[i][j][k]) - 2) / 26) % 26),
							                65 + ((abs(A[i][j][k]) - 2) % 26),
							                xaux,
							                yaux,
							                zaux,
							                0.0);
						}
						count++;

						/* If count equal to 100,000, restart count */
						if (count == 100000)
						    count = 1;

					}
				}

				/* Next cavity tag */
				tag++;

	}

    /* Close output PDB file */
	fclose (output);

}

/* Mark as invalid all the cavities with volume smaller than a user defined threshold */
void
remove_cavity (int ***A,
               int m,
               int n,
               int o,
               int tag)
{
    /* Declare variables */
	int i, j, k;

    /* Set number of processes in OpenMP */
	int ncores = omp_get_num_procs () - 1;
	omp_set_num_threads (ncores);


/*Create a parallel region*/
#pragma omp parallel default(shared)
/* Create a parallel loop, collapsing 3 loops inside 1 */
#pragma omp for collapse(3)
    /* Loop around the search box */
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < o; k++)
			    /*Remove tag*/
				if (A[i][j][k] == tag)
				    A[i][j][k] = 0;

}

/* Depth-First Search algorithm */
void
DFS (int ***A,
     int a,
     int b,
     int c,
     int m,
     int n,
     int o,
     int tag,
     double h)
{
    /* Declare variables */
	int i, j, k;

	/* Ignore points in border */
	if (a == 0 || a == m-1 || b == 0 || b == n-1 || c == 0 || c == o-1)
		return;

    /* If point is a cavity point, do ... */
	if (A[a][b][c] == 1 && !flagr) {

		A[a][b][c] = tag;
		volume += 1.0;

		/* Split big cavities */
		if (volume == 20000)
		    flagr = 1;

		if (!flagr) {

			/* Loop around one position in each direction from starting point */
			for (i = a-1 ; i <= a+1 ; i++)
				for (j = b-1 ; j<= b+1 ; j++)
					for (k = c-1; k <= c+1; k++)
					    /* Recursive call */
						DFS (A, i, j, k, m, n, o, tag, h);

		}

	}

}

/* Apply Depth-First Search in accessible nodes and cavity untagged points (1). Calculate cavity volume based on cavity
points with the same numeric tag.
Due to memory restrictions, the recursion is divided for big cavities */
int
DFS_search (int ***A,
            int m,
            int n,
            int o,
            double h,
            double filter)
{
	/* Declare variables */
	int i, j, k, r, t, y, tag = 1;
	double vol_aux = 0.0;
	node *p;
	flagr = 0;

	/* NULL pointer to last item in volume linked list */
	V = NULL;

	/*Loops around analysis box*/
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < o; k++)
			        /* If point is a cavity point, do ... */
					if (A[i][j][k] == 1) {

					volume = 0.0;
					tag++;
					/* Call DFS algorithm */
					DFS (A, i, j, k, m, n, o, tag, h);
					vol_aux = volume;

					/* Loop for big cavities */
					while (flagr) {

						vol_aux = 0;
						/* Loop around search box */
						for (r = 0; r < m; r++)
							for (t = 0; t < n; t++)
								for (y = 0; y < o; y++) {
									flagr = 0;
									/* Save volume in vol_aux */
									vol_aux += volume;
									volume = 0.0;
									/* For a given identified cavity point, check if there is a unidentified cavity
									point around it */
									if(A[r][t][y] == 1 && check_pos3 (A, r, t, y, m, n, o) == tag)
									    /* Call DFS algorithm */
									    DFS (A, r, t, y, m, n, o, tag, h);

								}
					}
					/* Cavity volume */
					volume = vol_aux;

                    /* If volume is less than threshold, remove cavity */
					if (volume*h*h*h < filter) {
						remove_cavity (A, m, n, o, tag);
						tag--;
					}
					else {
						/* Append item to volume linked list */
						p = malloc (sizeof (node));
						p->volume = volume*h*h*h;
						p->pos = tag;

						if (V != NULL)
						    p->next = V;
						else
						    p->next = NULL;
						V = p;
					}
				}

    /* Return number of cavities */
	return tag;

}

/* Check number of faces in contact with a protein point and return pondered area value */
void
check_faces (int ***S,
             int i,
             int j,
             int k,
             double h,
             double *area)
{
    /* Declare variables */
	int tag = S[i][j][k], faces = 0;
	double weight = 0.0;

	/* If face accessible to protein point, increment faces */
	if (S[i - 1][j][k] == 0)
	    faces++;
	if (S[i + 1][j][k] == 0)
	    faces++;
	if (S[i][j - 1][k] == 0)
	    faces++;
	if (S[i][j + 1][k] == 0)
	    faces++;
	if (S[i][j][k - 1] == 0)
	    faces++;
	if (S[i][j][k + 1] == 0)
	    faces++;

	/* Attribute weight based on voxel class */
	switch (faces) {
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
			    (S[i][j][k + 1] == 0 && S[i][j][k - 1] == 0))
			        weight = 2;
			else
			    weight = 1.5879;
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

	/* Save value in area vector inside its respective cavity */
	area[tag - 2] += h * h * weight;

}

/* Area calculation: KVFinder method */
void
Area_search (int ***S,
             int m,
             int n,
             int o,
             double h,
             int tag)
{
    /* Declare variables */
	int *area, i, j, k;

    /* Set number of processes in OpenMP */
	int ncores = omp_get_num_procs () - 1;
	omp_set_num_threads (ncores);

	/* Allocate memory to area object */
	area = (int*) calloc (tag, sizeof (int));

	for (i = 0; i < tag; i++)
	    area[i] = 0;

/* Create a parallel loop, collapsing 3 loops inside 1 and schedule dynamic allocation of threads */
#pragma omp parallel for shared (S, i, j, k, m, n, o) collapse (3) reduction (+: area[:tag])
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < o; k++){
				if (S[i][j][k] > 1)
					area[S[i][j][k] - 2]++;
			}

    /* Save area in KVFinder results struct */
	for (i = 0; i < tag; i++)
	    KVFinder_results[i].area = area[i] * h * h;

	 /* Free area object from memory */
	 free(area);

}

/* Compare two matrices, selecting the points where the smaller probe passed and the bigger one does not */
void
Matrix_subtract (int ***S,
                 int ***A,
                 int m,
                 int n,
                 int o,
                 double h,
                 double removal_distance)
{
    /* Declare variables */
	int i, j, k, i2, j2, k2, grid_spaces = ceil (removal_distance / h);

    /* Set number of processes in OpenMP */
	int ncores = omp_get_num_procs () - 1;
	omp_set_num_threads (ncores);

/* Create a parallel region */
#pragma omp parallel default(none), shared(A,S,grid_spaces,m,n,o,i,j,k), private(j2,i2,k2)
{
/* Create a parallel loop, collapsing 3 loops inside 1 and schedule dynamic allocation of threads */
#pragma omp for schedule(dynamic) collapse(3)
	/* Loop around the search box */
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < o; k++) {

                /* If point is a cavity, do ... */
				if (S[i][j][k]) {

					/*Loops around space occupied by probe from atom position*/
					for(i2 = i-grid_spaces; i2 <= i+grid_spaces; i2++)
						for(j2 = j-grid_spaces; j2 <= j+grid_spaces; j2++)
							for(k2 = k-grid_spaces; k2 <= k+grid_spaces; k2++)
								/*If inside box, do... */
								if(i2 >= 0 && i2 < m && j2 >= 0 && j2 < n && k2 >= 0 && k2 < o)
								    /* Mark points where big probe passed in cavities in A */
								    A[i2][j2][k2] = -1;

				}

			}

}

}

void
filter_outliers (int ***A,
                 int m,
                 int n,
                 int o)
{

	/* Declare variables */
	int i, j, k, faces2protein;

    /* Set number of processes in OpenMP */
	int ncores = omp_get_num_procs () - 1;
	omp_set_num_threads (ncores);

	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
		    for (k = 0; k < o; k++) {

                if (A[i][j][k] == 1) {

                    /* Initialize counter */
                    faces2protein = 0;

                    /* Check if a protein point (0) or a medium point (-1) is next to a cavity point (>=1) */
                    if (i-1 >= 0)
                        if (A[i-1][j][k] == 0 || A[i-1][j][k] == -1)
                            faces2protein++;
                    if (i+1 < m)
                        if (A[i+1][j][k] == 0 || A[i+1][j][k] == -1)
                            faces2protein++;
                    if (j-1 >= 0)
                        if (A[i][j-1][k] == 0 || A[i][j-1][k] == -1)
                            faces2protein++;
                    if (j+1 < n)
                        if (A[i][j+1][k] == 0 || A[i][j+1][k] == -1)
                            faces2protein++;
                    if (k-1 >= 0)
                        if (A[i][j][k-1] == 0 || A[i][j][k-1] == -1)
                            faces2protein++;
                    if (k+1 < o)
                        if (A[i][j][k+1] == 0 || A[i][j][k+1] == -1)
                            faces2protein++;

                    /* Cavity point is a medium point */
                    if (faces2protein == 6)
                        A[i][j][k] = -1;

                }

		    }

}

char*
combine (const char *s1,
         const char *s2)
{
    /* Declare variables */
	const size_t len1 = strlen (s1);
	const size_t len2 = strlen (s2);
	char *result = malloc (len1 + len2 + 1); /* +1 for the zero-terminator */

	memcpy (result, s1, len1);
	memcpy (result + len1, s2, len2 + 1); /* +1 to copy the null-terminator */

	return result;

}

/* Free int*** matrix from memory */
void
free_matrix (int ***A,
             int m,
             int n,
             int o)
{
    /* Declare variables */
	int i, j;

	for (i = 0; i < m; i++) {

		for (j = 0; j < n; j++)
		    free (A[i][j]);
		free (A[i]);

	}

	free (A);

}

/* Free node (volume and tag linked list) from memory */
void
free_node ()
{
    /* Declare variables */
	node *p;

	while (V != NULL) {
	   p = V;
	   V = V->next;
	   free (p);
	}

}

