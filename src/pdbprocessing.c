/* This file contains all the functions used to process PDB files used as input
 * for KVFinder */

/* Import native modules */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Import custom modules  */
#include "utils.h"
#include "dictionaryprocessing.h"
#include "pdbprocessing.h"

/* Convert string to double */
void convert(char S[50], double *coord) {
  /* Declare variables */
  int i;
  char AUX[50] = "";

  for (i = 0; S[i] != '\0' && i < 50; i++)
    ;
  // extract(S, 0, i + 1, AUX);
  _extract(S, strlen(S), AUX, strlen(AUX), 0, i+1);
  *coord = atof(AUX);
}

/* Create a linked list for PDB information */
void insert_atom(double x, double y, double z, double radius, int resnumber,
                 char chain, char res_name) {
  /* Declare variables */
  atom *p;

  if (v == NULL) {
    v = malloc(sizeof(atom));
    v->x = x;
    v->y = y;
    v->z = z;
    v->radius = radius;
    v->resnum = resnumber;
    v->chain = chain;
    v->res_name = res_name;
    v->next = NULL;
  } else {
    p = malloc(sizeof(atom));
    p->x = x;
    p->y = y;
    p->z = z;
    p->radius = radius;
    p->resnum = resnumber;
    p->chain = chain;
    p->res_name = res_name;
    p->next = v->next;
    v->next = p;
  }
}

/* Free linked list for PDB information from memory */
void free_atom() {
  /* Declare variables */
  atom *p;

  while (v != NULL) {
    p = v;
    v = v->next;
    free(p);
  }
}

/* PDB_load function gets DIC object that contain residues linked lists (motif
and radius), tablesize object that contains number of residues, TABLE[i][j]
object that contain all residues symbols, path to PDB file, Probe In size,
grid size (m,n,o), step size and P1 (X1, Y1, Z1) that is the reference of the
grid */
/* Create a linked list for PDB information and save coordinates (x,y,z), atom
 * radius, residue number and chain */
int PDB_load(dict *DIC[TABLE_SIZE], int tablesize,
             char TABLE[TABLE_SIZE][RES_SIZE], char PDB_NAME[NAME_MAX],
             double probe, int m, int n, int o, double h, double X1, double Y1,
             double Z1, FILE **log_file) {
  /* Declare variables */
  int flag = 1, i, j, number;
  char AUX[50] = "", LINE[100] = "", X[50] = "", Y[50] = "",
       Z[50] = "", RES[RES_SIZE], ATOM_TYPE[ATOM_SIZE], ATOM_SYMBOL[2], chain;
  double x, y, z, x1, y1, z1, xaux, yaux, zaux, radius;
  FILE *arqPDB;

  /* Open PDB file */
  arqPDB = fopen(PDB_NAME, "r");

  /* PDB file not found */
  if (arqPDB == NULL) {

    /* Print error and exit*/
    printf("\033[0;31mError:\033[0m PDB file not found!\n");
    exit(-1);

  }
  /* PDB file found */
  else {

    /* Parse PDB */
    /* While PDB file is not over, do ... */
    while (_read_line(arqPDB, LINE, 100)) {

      /* Extract Record Name */
      // extract(LINE, 0, 6, AUX);
      _extract(LINE, strlen(LINE), AUX, strlen(AUX), 0, 6);
      /* If Record Name is equal to ATOM or HETATM, do ... */
      if (!strcmp(AUX, "ATOM  ") || !strcmp(AUX, "HETATM")) {
        for (i = 12, j = 0; i < 17; i++) {
          if (LINE[i] != ' ') {

            /* Get ATOM name */
            ATOM_TYPE[j] = LINE[i];
            j++;
          }
        }

        /* Mark last position of ATOM name */
        ATOM_TYPE[j] = '\0';

        /* Get residue name */
        RES[0] = LINE[17];
        RES[1] = LINE[18];
        RES[2] = LINE[19];
        RES[3] = '\0';

        /* Get residue sequence number */
        number = 0;
        _extract(LINE, strlen(LINE), AUX, strlen(AUX), 22, 26);
        number = atoi(AUX);

        /* Extract atom symbol */
        // extract(LINE, 76, 78, ATOM_SYMBOL);
        _extract(LINE, strlen(LINE), ATOM_SYMBOL, strlen(ATOM_SYMBOL), 76, 78);

        /* Extract chain identifier */
        chain = _extract_chain(LINE, strlen(LINE));

        /* Extract x coordinate */
        _extract(LINE, strlen(LINE), X, strlen(X), 30, 38);
        x = atof(X);
        /* Extract y coordinate */
        _extract(LINE, strlen(LINE), Y, strlen(Y), 38, 46);
        y = atof(Y);
        /* Extract z coordinate */
        _extract(LINE, strlen(LINE), Z, strlen(Z), 46, 54);
        z = atof(Z);

        /* Get radius for an atom in a specific residue based on VdW radius
         * dictionary */
        radius = dictionary_consult_radius(RES, ATOM_TYPE, DIC, tablesize,
                                           TABLE, ATOM_SYMBOL, log_file);

        /* Calculate coordinate (x1, y1, z1) for atom */
        x1 = (x - X1) / h;
        y1 = (y - Y1) / h;
        z1 = (z - Z1) / h;
        xaux = x1 * cosb + z1 * sinb;
        yaux = y1;
        zaux = -x1 * sinb + z1 * cosb;
        x1 = xaux;
        y1 = yaux * cosa - zaux * sina;
        z1 = yaux * sina + zaux * cosa;

        /* Create a linked list only for atoms inside search box */
        if (x1 > 0.0 - (probe + radius) / h &&
            x1 < (double)m + (probe + radius) / h &&
            y1 > 0.0 - (probe + radius) / h &&
            y1 < (double)n + (probe + radius) / h &&
            z1 > 0.0 - (probe + radius) / h &&
            z1 < (double)o + (probe + radius) / h) {

          /* Save coordinates (x,y,z), radius, residue number and chain */
          insert_atom(x, y, z, radius, number, chain, convertRES(RES));
        }
      }
    }
  }

  /* Close PDB file */
  fclose(arqPDB);

  /* Return flag indicating file has been read */
  return flag;
}

/* PDB_load2 function gets DIC object that contain residues linked lists (motif
and radius), tablesize object that contains number of residues, TABLE[i][j]
object that contain all residues symbols, and path to PDB file*/
/* Create a linked list for PDB information and save coordinates (x,y,z)*/
int PDB_load2(char PDB_NAME[NAME_MAX]) {
  /* Declare variables */
  int flag = 1, i;
  double x, y, z;
  char AUX[50] = "", X[50] = "", Y[50] = "", Z[50] = "", LINE[100], chain;
  FILE *arqPDB;

  /* NULL pointer to last item in PDB information linked list */
  v = NULL;

  /* Open PDB file */
  arqPDB = fopen(PDB_NAME, "r");

  /* PDB file not found */
  if (arqPDB == NULL) {

    /* Print error and exit */
    printf("\033[0;31mError:\033[0m PDB file not found!\n");
    exit(-1);

  }
  /* PDB file found */
  else {

    /* Parse PDB */
    /* While PDB file is not over, do... */
    while (_read_line(arqPDB, LINE, 100)) {

      /* Extract Record Name */
      // extract(LINE, 0, 6, AUX);
      _extract(LINE, strlen(LINE), AUX, strlen(AUX), 0, 6);
      /* If Record Name is equal to ATOM or HETATM, do ... */
      if (!strcmp(AUX, "ATOM  ") || !strcmp(AUX, "HETATM")) {

        /* Extract x coordinate */
        // extract(LINE, 30, 38, X);
        _extract(LINE, strlen(LINE), X, strlen(X), 30, 38);
        /* Extract y coordinate */
        // extract(LINE, 38, 46, Y);
        _extract(LINE, strlen(LINE), Y, strlen(Y), 38, 46);
        /* Extract z axis coordinate */
        // extract(LINE, 46, 54, Z);
        _extract(LINE, strlen(LINE), Z, strlen(Z), 46, 54);
        /* Remove whitespaces */
        trim(X, ' ');
        trim(Y, ' ');
        trim(Z, ' ');
        /* Convert coordinates from strings to doubles */
        convert(X, &x);
        convert(Y, &y);
        convert(Z, &z);

        /* Save coordinate (x,y,z) */
        insert_atom(x, y, z, 0.0, 0, 0, 0);
      }
    }
  }

  /* Close PDB file */
  fclose(arqPDB);

  /* Return flag indicating file has been read */
  return flag;
}

/* PDB_load3 gets path to PDB file */
/* Creates a linked list for PDB information */
int PDB_load3(char PDB_NAME[NAME_MAX]) {
  /* Declare variables */
  int flag = 1, i, resnumber;
  double x, y, z;
  char AUX[50] = "", X[50] = "", Y[50] = "", Z[50] = "",
       LINE[100], chain;
  FILE *arqPDB;

  /* NULL pointer to last item in PDB information linked list */
  v = NULL;

  /* Open PDB file */
  arqPDB = fopen(PDB_NAME, "r");

  /* PDB file not found */
  if (arqPDB == NULL) {

    /* Print error and exit */
    printf("\033[0;31mError:\033[0m PDB file not found!\n");
    exit(-1);

  } else {

    /* Parse PDB */
    /* While PDB is not over, do ... */
    while (_read_line(arqPDB, LINE, 100)) {

      /* Extract Record Name */
      // extract(LINE, 0, 6, AUX);
      _extract(LINE, strlen(LINE), AUX, strlen(AUX), 0, 6);
      /* If Record Name is equal to ATOM or HETATM, do ... */
      if (!strcmp(AUX, "ATOM  ") || !strcmp(AUX, "HETATM")) {

        /*Get residue sequence number*/
        resnumber = 0;
        _extract(LINE, strlen(LINE), AUX, strlen(AUX), 22, 26);
        resnumber = atoi(AUX);

        /* Extract x coordinate */
        // extract(LINE, 30, 38, X);
        _extract(LINE, 100, X, strlen(X), 30, 38);
        /* Extract y coordinate */
        // extract(LINE, 38, 46, Y);
        _extract(LINE, 100, Y, strlen(Y), 38, 46);
        /* Extract z coordinate */
        // extract(LINE, 46, 54, Z);
        _extract(LINE, 100, Z, strlen(Z), 46, 54);

        /* Remove whitespaces */
        trim(X, ' ');
        trim(Y, ' ');
        trim(Z, ' ');
        /* Convert coordinates from strings to doubles */
        convert(X, &x);
        convert(Y, &y);
        convert(Z, &z);

        /* Extract chain identifier */
        chain = _extract_chain(LINE, strlen(LINE));

        /* Save coordinate (x,y,z), residue number and chain */
        insert_atom(x, y, z, 0.0, resnumber, chain, 0);
      }
    }
  }

  /* Close PDB file */
  fclose(arqPDB);

  /* Return flag indicating file has been read */
  return flag;
}
