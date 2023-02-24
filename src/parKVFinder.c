/* This file contains the main script of KVFinder.
It is responsible for controlling the execution flux of KVFinder software,
according to the user definitions. A brief explanation of the most important
parts may be found in the source code */

/* Import builtin modules */
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

/* Import custom modules */
/* WARNING: keep this importing order */
#include "utils.h"
#include "fileprocessing.h"
#include "gridprocessing.h"
#include "argparser.h"

/* Main function */
int main(int argc, char **argv) {

  /* Evaluate elapsed time */
  struct timeval toc, tic;
  gettimeofday(&tic, NULL);

  /* Get date and time */
  time_t t = time(NULL);
  struct tm *tm = localtime(&t);

  /* Command to allow nested thread creation */
  omp_set_nested(1);

  /* Setting number of cores */
  int ncores = omp_get_num_procs() - 1;
  omp_set_num_threads(ncores);

  /* Global variables */
  double h, probe_in, probe_out, volume_cutoff, ligand_cutoff, removal_distance,
      norm1, norm2, norm3, Vvoxel, multiple;
  double X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4;
  double bX1, bY1, bZ1, bX2, bY2, bZ2, bX3, bY3, bZ3, bX4, bY4, bZ4;
  int ligand_mode, surface_mode, whole_protein_mode, resolution_mode, box_mode,
      kvp_mode;
  static int verbose_flag = 0;
  int m, n, o, i, j, k, ncav, tablesize;
  char TABLE[500][4], PDB_NAME[500], LIGAND_NAME[500], dictionary_name[500],
      OUTPUT[500], BASE_NAME[500];
  char boxmode_flag[6], resolution_flag[7], whole_protein_flag[6], mode_flag[6],
      surface_flag[6], step_flag[6], kvpmode_flag[6];
  char log_buffer[4096], *output, *output_folder, *output_pdb, *output_results,
      *pdb_name;
  vdw *DIC[500];
  FILE *parameters_file, *log_file;
  atom *p;
  int ***A, ***S;
  double ***M;

  if (argc == 1) {
    /* Check if parameters.toml exists */
    if (access("parameters.toml", F_OK)) {

      /* Print header */
      print_header();
      /* Print usage */
      print_usage();
      /* File does not exist */
      exit(-1);
    }
    /* Read parameters TOML files */
    parameters *parameters =
        readTOML("parameters.toml"); /* Read TOML file inside struct TOML */

    /* Save TOML parameters from struct TOML parameters to KVFinder variables */
    X1 = parameters->X1;
    Y1 = parameters->Y1;
    Z1 = parameters->Z1;
    X2 = parameters->X2;
    Y2 = parameters->Y2;
    Z2 = parameters->Z2;
    X3 = parameters->X3;
    Y3 = parameters->Y3;
    Z3 = parameters->Z3;
    X4 = parameters->X4;
    Y4 = parameters->Y4;
    Z4 = parameters->Z4;
    bX1 = parameters->bX1;
    bY1 = parameters->bY1;
    bZ1 = parameters->bZ1;
    bX2 = parameters->bX2;
    bY2 = parameters->bY2;
    bZ2 = parameters->bZ2;
    bX3 = parameters->bX3;
    bY3 = parameters->bY3;
    bZ3 = parameters->bZ3;
    bX4 = parameters->bX4;
    bY4 = parameters->bY4;
    bZ4 = parameters->bZ4;
    strcpy(PDB_NAME, parameters->PDB_NAME);
    strcpy(OUTPUT, parameters->OUTPUT);
    strcpy(BASE_NAME, parameters->BASE_NAME);
    strcpy(LIGAND_NAME, parameters->LIGAND_NAME);
    strcpy(dictionary_name, parameters->dictionary_name);
    strcpy(resolution_flag, parameters->resolution_flag);
    volume_cutoff = parameters->volume_cutoff;
    ligand_cutoff = parameters->ligand_cutoff;
    removal_distance = parameters->removal_distance;
    probe_in = parameters->probe_in;
    probe_out = parameters->probe_out;
    h = parameters->h;
    whole_protein_mode = parameters->whole_protein_mode;
    box_mode = parameters->box_mode;
    surface_mode = parameters->surface_mode;
    kvp_mode = parameters->kvp_mode;
    ligand_mode = parameters->ligand_mode;

    /* Free struct TOML */
    free(parameters);
  } else {
    /* Save command line arguments inside KVFinder variables */
    verbose_flag =
        argparser(argc, argv, &box_mode, &kvp_mode, &ligand_mode, &surface_mode,
                  &whole_protein_mode, PDB_NAME, LIGAND_NAME, dictionary_name,
                  OUTPUT, BASE_NAME, resolution_flag, &h, &probe_in, &probe_out,
                  &volume_cutoff, &ligand_cutoff, &removal_distance, &X1, &Y1,
                  &Z1, &X2, &Y2, &Z2, &X3, &Y3, &Z3, &X4, &Y4, &Z4, &bX1, &bY1,
                  &bZ1, &bX2, &bY2, &bZ2, &bX3, &bY3, &bZ3, &bX4, &bY4, &bZ4);
  }
  /* Set step size (h) and resolution_mode */
  if (!strcmp(resolution_flag, "Off"))
    resolution_mode = 0;
  else
    h = _resolution2step(resolution_flag);
  resolution_mode = 1;

  /* Load vdW dictionary */
  tablesize = _get_residues_information(
      dictionary_name,
      TABLE); /*Read dictionary file and return number of residues*/
  if (verbose_flag)
    fprintf(stdout, "> Loading atomic dictionary file\n");
  read_vdw(dictionary_name, DIC, tablesize); /*Dictionary Loaded*/

  /* Preparing files paths */
  if (OUTPUT[strlen(OUTPUT) - 1] == '/') {
    if (OUTPUT[strlen(OUTPUT) - 2] == '/')
      OUTPUT[strlen(OUTPUT) - 1] = '\0';
    output = _combine(OUTPUT, "KV_Files/");
  } else {
    if (OUTPUT[0] == '\0')
      output = _combine(OUTPUT, "KV_Files/");
    else
      output = _combine(OUTPUT, "/KV_Files/");
  }
  pdb_name = realpath(PDB_NAME, NULL);

  /* Create KV_Files folder */
  mkdir(output, S_IRWXU);

  /* Create log_file */
  log_file = fopen(_combine(output, "KVFinder.log"),
                   "a+"); /* Open log file and append information */
  memset(log_buffer, '\0', sizeof(log_buffer)); /* Create buffer */
  setvbuf(log_file, log_buffer, _IOFBF,
          4096); /* Define buffer as writing buffer of size 4096 */

  /* Create BASE_NAME folder for the running analysis */
  output =
      _combine(output, BASE_NAME); /* Include BASE_NAME folder in KV_Files */
  mkdir(output, S_IRWXU);          /* Create BASE_NAME folder in KV_Files */
  output_folder =
      _combine(output, "/"); /* Include bar after appending BASE_NAME */
  output = _combine(
      output_folder,
      BASE_NAME); /* Include BASE_NAME to output path in BASE_NAME folder */

  /* Create output PDB and results file names */
  output_results = _combine(
      output,
      ".KVFinder.results.toml"); /* Create a output path to results file */
  output_pdb = _combine(
      output,
      ".KVFinder.output.pdb"); /* Create a output path to output PDB file */

  /* Print in shell the PDB path that is running in KVFinder */
  if (!verbose_flag)
    fprintf(stdout, "[PID %u] Running parKVFinder for: %s\n", getpid(),
            pdb_name);

  fprintf(log_file,
          "==========\tSTART\tRUN\t=========\n\n"); /* Print the start tag for
                                                       the run inside log */
  fprintf(log_file, "Date and time: %s\n", asctime(tm));
  fprintf(log_file, "Running parKVFinder for: %s\n",
          pdb_name); /* Print path to PDB file */
  fprintf(log_file, "Dictionary: %s\n",
          dictionary_name); /* Print path to dictionary */
  /* if Resolution is Off, print Resolution and Step size */
  if (strcmp(resolution_flag, "Off") == 0)
    fprintf(log_file, "Resolution: %s\nInput step size: %.2lf\n",
            resolution_flag, h);
  /* Else, print resolution */
  else
    fprintf(log_file, "Resolution: %s\n", resolution_flag);

  if (verbose_flag)
    fprintf(stdout, "> Calculating grid dimensions\n");

  if (whole_protein_mode) {

    X1 = 999999;
    Y1 = 999999;
    Z1 = 999999;
    X2 = -999999;
    Y3 = -999999;
    Z4 = -999999;

    /* Create a linked list (dictionary) for PDB file | saves only position
     * (x,y,z) and chain */
    soft_read_pdb(PDB_NAME, 0, 0);

    /*Reduces box to protein size*/
    for (p = v; p != NULL; p = p->next) {

      if (p->x < X1)
        X1 = (p->x);
      if (p->y < Y1)
        Y1 = (p->y);
      if (p->z < Z1)
        Z1 = (p->z);
      if (p->x > X2)
        X2 = (p->x);
      if (p->y > Y3)
        Y3 = (p->y);
      if (p->z > Z4)
        Z4 = (p->z);
    }

    /* Free van der Waals radius dictionary from memory */
    _free_atom();

  } /* If it will be used a user defined search space, without the Probe Out
       Adjustment define its limits */
  else if (box_mode) {
  };

  /* Predefined resolution:
  - Low: h = 0.6A
  - Medium: h = 0.5A
  - High h = 0.25A */
  if (resolution_mode) {
    fprintf(log_file, "Chosen step size = %.2lf\n", h);
  }
  /* Convert step size to string */
  sprintf(step_flag, "%.2lf", h);

  if (whole_protein_mode) {

    /* Prepare vertices */
    X1 = X1 - probe_out - h;
    Y1 = Y1 - probe_out - h;
    Z1 = Z1 - probe_out - h;
    X2 = X2 + probe_out + h;
    Y3 = Y3 + probe_out + h;
    Z4 = Z4 + probe_out + h;
    X3 = X1;
    X4 = X1;
    Y2 = Y1;
    Y4 = Y1;
    Z2 = Z1;
    Z3 = Z1;
  }

  /* Defines the grid inside the search space | Point 1 is the reference of our
   * grid */
  /* Calculate distances between points in box and reference */
  norm1 = sqrt((X2 - X1) * (X2 - X1) + (Y2 - Y1) * (Y2 - Y1) +
               (Z2 - Z1) * (Z2 - Z1));
  norm2 = sqrt((X3 - X1) * (X3 - X1) + (Y3 - Y1) * (Y3 - Y1) +
               (Z3 - Z1) * (Z3 - Z1));
  norm3 = sqrt((X4 - X1) * (X4 - X1) + (Y4 - Y1) * (Y4 - Y1) +
               (Z4 - Z1) * (Z4 - Z1));

  /* Grid steps based on distances */
  /*X axis*/
  if (fmod(norm1, h) != 0) {
    m = (int)(norm1 / h) + 1;
  } else {
    m = (int)(norm1 / h);
  }
  /*Y axis*/
  if (fmod(norm2, h) != 0) {
    n = (int)(norm2 / h) + 1;
  } else {
    n = (int)(norm2 / h);
  }
  /*Z axis*/
  if (fmod(norm3, h) != 0) {
    o = (int)(norm3 / h) + 1;
  } else {
    o = (int)(norm3 / h);
  }

  /* Calculate data used on the spatial manipulation of the protein */
  sina = (Y4 - Y1) / norm3;
  cosa = (Y3 - Y1) / norm2;
  sinb = (Z2 - Z1) / norm1;
  cosb = (X2 - X1) / norm1;

  /* Calculate Voxel volume */
  Vvoxel = h * h * h;

  /* Write voxel volume inside log file */
  fprintf(log_file, "Unitary box volume: %.4lf Angstrons^3\n", Vvoxel);
  /* Write m, n and o axis inside log file */
  fprintf(log_file, "m: %d\tn: %d\to: %d\n", m, n, o);
  /* Write sina, sinb, cosa and cosb inside log file */
  fprintf(log_file, "sina: %.2lf\tsinb: %.2lf\t\ncosa: %.2lf\tcosb: %.2lf\n",
          sina, sinb, cosa, cosb);

  if (verbose_flag)
    fprintf(stdout, "> Reading PDB coordinates\n");

  /* Protein Coordinates Extraction */
  /* Create a linked list for PDB information */
  /* Save coordinates (x,y,z), atom radius, residue number and chain */
  if (read_pdb(PDB_NAME, DIC, tablesize, TABLE, probe_in, m, n, o, h, X1, Y1,
               Z1, &log_file)) {

    if (verbose_flag)
      fprintf(stdout, "> Creating grid\n");

    /* Matrix Allocation and Initialization */
    /* int ***A: Grid representing empty spaces and surface points along marked
    by small probe int ***S: Grid representing empty spaces and surface points
    along marked by big probe double ***M: Grid representing depth in each
    cavity point */
    A = (int ***)calloc(m, sizeof(int **));
    S = (int ***)calloc(m, sizeof(int **));
    M = (double ***)calloc(m, sizeof(double **));
    for (i = 0; i < m; i++) {
      A[i] = (int **)calloc(n, sizeof(int *));
      S[i] = (int **)calloc(n, sizeof(int *));
      M[i] = (double **)calloc(n, sizeof(double *));
      for (j = 0; j < n; j++) {
        A[i][j] = (int *)calloc(o, sizeof(int));
        S[i][j] = (int *)calloc(o, sizeof(int));
        M[i][j] = (double *)calloc(o, sizeof(double));
        for (k = 0; k < o; k++) {
          A[i][j][k] = 1;
          S[i][j][k] = 1;
          M[i][j][k] = 0.0;
        }
      }
    }

    if (verbose_flag)
      fprintf(stdout, "> Filling grid with probe in surface\n");
    /* Mark the grid with 0, leaving a small probe size around the protein */
    SAS(A, m, n, o, h, probe_in, X1, Y1, Z1);

    /* Mark space occupied by a small probe size from protein surface */
    if (surface_mode) {
      SES(A, m, n, o, h, probe_in);
    }

    if (verbose_flag)
      fprintf(stdout, "> Filling grid with probe out surface\n");
    /* Mark the grid with 0, leaving a big probe size around the protein */
    SAS(S, m, n, o, h, probe_out, X1, Y1, Z1);
    /* Mark space occupied by a big probe size from protein surface */
    SES(S, m, n, o, h, probe_out);

    if (verbose_flag)
      fprintf(stdout, "> Defining biomolecular cavities\n");
    /* Mark points where small probe passed and big probe did not */
    subtract(A, S, m, n, o, h, removal_distance);

    /* Ligand adjustment mode */
    if (ligand_mode) {
      if (verbose_flag)
        fprintf(stdout, "> Adjusting ligand\n");
      /* Free linked list (dictionary) from memory */
      _free_atom();
      /* Creates a linked list for Ligand information* | saves position (x,y,z),
       * atom radius, resnumber, chain */
      read_pdb(LIGAND_NAME, DIC, tablesize, TABLE, probe_in, m, n, o, h, X1, Y1,
               Z1, &log_file);
      /* Mark regions that do not belong to ligand_cutoff */
      adjust2ligand(A, m, n, o, h, ligand_cutoff, X1, Y1, Z1);
      /* Free linked list (dictionary) from memory */
      _free_atom();
      /* Create a linked list for PDB information* | saves position (x,y,z),
       * atom radius, resnumber, chain */
      read_pdb(PDB_NAME, DIC, tablesize, TABLE, probe_in, m, n, o, h, X1, Y1,
               Z1, &log_file);
    }

    /* The points outside the user defined search space are excluded here */
    if (box_mode) {
      if (verbose_flag)
        fprintf(stdout, "> Filtering grid points\n");

      filter2box(A, m, n, o, h, bX1, bY1, bZ1, bX2, bY2, bZ2, norm1);
    }

    /*Remove outlier points*/
    filter_noise(A, m, n, o);

    /* Grouping Cavities and calculating Volume and */
    if (verbose_flag)
      fprintf(stdout, "> Clustering cavities and calculating volume\n");
    ncav = clustering(A, m, n, o, h, volume_cutoff);

    if (ncav > 0) {
      /* Create KVFinder_results structure */
      KVFinder_results = (KVresults *)calloc(ncav, sizeof(KVresults));
      cavity = (coords *)calloc(ncav, sizeof(coords));
      boundary = (coords *)calloc(ncav, sizeof(coords));

      /* Pass volume to KVFinder_results structure */
      node *p;
      for (p = V; p != NULL; p = p->next)
        KVFinder_results[(p->pos)].volume = p->volume;
      free_node();
      free(p);
      free(V);

      /* Defining surface points and calculating area*/
      if (verbose_flag)
        fprintf(stdout, "> Defining surface points and calculating area\n");
      filter_surface(A, S, m, n, o);
      area(S, m, n, o, h, ncav);

      /* Define interface residues for each cavity */
      if (verbose_flag)
        fprintf(stdout, "> Retrieving interface residues\n");
      interface(A, m, n, o, h, probe_in, ncav, X1, Y1, Z1);

      /* Computing depth */
      if (verbose_flag)
        fprintf(stdout,
                "> Defining cavity-bulk boundary and calculating depth\n");
      filter_boundary(A, m, n, o, ncav);
      depth(A, M, m, n, o, h, ncav);

      /* Turn ON(1) filled cavities option */
      if (verbose_flag)
        fprintf(stdout, "> Writing cavities PDB file\n");
      /* Export Cavities PDB */
      export(output_pdb, A, S, M, kvp_mode, m, n, o, h, ncav, X1, Y1, Z1);

      /* Write results file */
      if (verbose_flag)
        fprintf(stdout, "> Writing results file\n");
      write_results(output_results, pdb_name, output_pdb, LIGAND_NAME, h, ncav);

    } else {
      fprintf(stdout, "> parKVFinder found no cavities!\n");
    }

    /*Free data structures used for depth calculation*/
    _free_atom();
    free(cavity);
    free(boundary);
    free_igrid(A, m, n, o);
    free_igrid(S, m, n, o);
    free_dgrid(M, m, n, o);
  }

  /*Evaluate elapsed time*/
  gettimeofday(&toc, NULL);
  printf("done!\n");
  printf(
      "Elapsed time: %.2lf seconds\n",
      (double)(toc.tv_sec - tic.tv_sec) +
          ((toc.tv_usec - tic.tv_usec) / 1000000.0)); /*Print the elapsed time*/
  fprintf(log_file, "Elapsed time: %.2lf seconds\n",
          (double)(toc.tv_sec - tic.tv_sec) +
              ((toc.tv_usec - tic.tv_usec) /
               1000000.0)); /*Print the elapsed time for the run inside log*/
  fprintf(log_file,
          "\n==========\tEND\tRUN\t=========\n\n"); /*Print the end tag for the
                                                       run inside log*/

  /*Write log.txt*/
  fflush(log_file);
  fclose(log_file);

  return 0;
}
