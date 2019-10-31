/* This file contains all the functions used to process PDB files used as input for KVFinder */

/* Import native modules */
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Import custom modules  */
#include "dictionary.h"
#include "PDBprocessing.h"

/*====================================================================================================================*/
/*==============================================  Auxiliary Functions  ===============================================*/

/* Read PDB file inside LINE */
int
get_line (FILE *arq,
          char LINE[100])
{
    /* Declare variables */
	int i;

	/* Read PDB file inside LINE[100] */
	for (i = 0, LINE[i] = getc (arq); LINE[i] != EOF && LINE[i] != '\n' && i < 100 ; i++, LINE[i] = getc (arq));

	/* If LINE[100] is not \n, read PDB file until LINE[100] assume \n or EOF value */
	if (i == 100 && LINE[i] != '\n')
		for (; LINE[i] != '\n' && LINE[i] != EOF; LINE[i] = getc (arq));

    /* If PDB file is over, return False */
	if (LINE[i] == EOF)
	    return 0;
	/* If PDB file is not over, return True */
	else
	    return 1;

}

/* Extract from LINE[a] until LINE[b] to a char array S[] of size b-a */
void
extract (char LINE[100],
         int a,
         int b,
         char S[50])
{
    /* Declare variables */
	int i;

	/* Test input */
	if (b-a > 50 || b > 100 || a > b)
	    return;

    /* Extract process */
	for (i = a; i < b; i++)
	    S[i-a] = LINE[i];

	/* Mark last position in S[] */
	if (i < 50)
	    S[i] = '\0';

}

/* Extract chain identifier in LINE[21] */
char
extractChain (char LINE[100])
{
	return LINE[21];
}

/* Convert string to double */
void
convert (char S[50],
         double* coord)
{
    /* Declare variables */
	int i;
	char AUX[50] = "";

	for (i = 0; S[i] != '\0' && i < 50; i++);
	extract (S, 0, i+1, AUX);
	*coord = atof (AUX);
}

/* Create a linked list for PDB information */
void
insert_atom (double x,
             double y,
             double z,
             double radius,
             int resnumber,
             char chain,
             char res_name)
{
    /* Declare variables */
	atom *p;

	if (v == NULL) {
		v = malloc (sizeof (atom));
		v->x = x;
		v->y = y;
		v->z = z;
		v->radius = radius;
		v->resnum = resnumber;
		v->chain = chain;
		v->res_name = res_name;
		v->next = NULL;
	}
	else {
		p = malloc (sizeof (atom));
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
void
free_atom ()
{
    /* Declare variables */
	atom *p;

	while (v != NULL) {
		p = v;
		v = v->next;
		free(p);
	}

}

/*==============================================  Auxiliary Functions  ===============================================*/
/*====================================================================================================================*/

/*====================================================================================================================*/
/*================================================  Main Functions  ==================================================*/

/* PDB_load function gets DIC object that contain residues linked lists (motif and radius), tablesize object that
contains number of residues, TABLE[i][j] object that contain all residues symbols, path to PDB file, Probe In size,
grid size (m,n,o), step size and P1 (X1, Y1, Z1) that is the reference of the grid */
/* Create a linked list for PDB information and save coordinates (x,y,z), atom radius, residue number and chain */
int
PDB_load (dict *DIC[TABLE_SIZE],
          int tablesize,
          char TABLE[TABLE_SIZE][RES_SIZE],
          char PDB_NAME[NAME_MAX],
          double probe,
          int m,
          int n,
          int o,
          double h,
          double X1,
          double Y1,
          double Z1,
          FILE **log_file)
{
    /* Declare variables */
	int flag = 1, i, j, number;
	char AUX[50] = "",
	     AUX2[50] = "",
	     LINE[100] = "",
	     X[50] = "",
	     Y[50] = "",
	     Z[50] = "",
	     RES[RES_SIZE],
	     ATOM_TYPE[ATOM_SIZE],
	     ATOM_SYMBOL[2],
	     chain;
	double x, y, z, x1, y1, z1, xaux, yaux, zaux, radius;
	FILE *arqPDB;

    /* Open PDB file */
	arqPDB = fopen (PDB_NAME, "r");

	/* PDB file not found */
	if (arqPDB == NULL) {

	    /* Print error */
		printf ("Reading Error: PDB file not found!\n");
		/* Close PDB file */
		fclose (arqPDB);
		/* Return flag indicating file has not been found */
		return !flag;

	}
	/* PDB file found */
	else {

		/* Parse PDB */
		/* While PDB file is not over, do ... */
		while (get_line (arqPDB, LINE)) {

		    /* Extract Record Name */
		    extract(LINE, 0, 6, AUX);
		    /* If Record Name is equal to ATOM or HETATM, do ... */
			if (!strcmp (AUX, "ATOM  ") || !strcmp (AUX, "HETATM")) {
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
				RES[0] = LINE[17]; RES[1] = LINE[18]; RES[2] = LINE[19]; RES[3] = '\0';

				/* Get residue sequence number */
				number = 0;
				extract (LINE, 22, 26, AUX2);
				number = atoi (AUX2);

                /* Extract atom symbol */
				extract (LINE, 76, 78, ATOM_SYMBOL);

				/* Extract chain identifier */
				chain = extractChain (LINE);

				/* Extract x coordinate */
				extract (LINE, 30, 38, X);
				/* Extract y coordinate */
				extract (LINE, 38, 46, Y);
				/* Extract z coordinate */
				extract (LINE, 46, 54, Z);
                /* Remove whitespaces */
				trim (X, ' ');  trim (Y, ' ');  trim (Z, ' ');
				/* Convert coordinates from strings to doubles */
				convert (X, &x); convert (Y, &y); convert (Z, &z);

				/* Get radius for an atom in a specific residue based on VdW radius dictionary */
				radius = dictionary_consult_radius (RES, ATOM_TYPE, DIC, tablesize, TABLE, ATOM_SYMBOL, log_file);

				/* Calculate coordinate (x1, y1, z1) for atom */
				x1 = (x-X1)/h; y1 = (y-Y1)/h; z1 = (z-Z1)/h;
				xaux = x1*cosb + z1*sinb; yaux = y1; zaux = -x1*sinb + z1*cosb;
				x1 = xaux; y1 = yaux*cosa - zaux*sina; z1 = yaux*sina + zaux*cosa;

				/* Create a linked list only for atoms inside search box */
				if (x1 > 0.0-(probe+radius)/h  &&
				    x1 < (double)m+(probe+radius)/h &&
				    y1 > 0.0-(probe+radius)/h &&
				    y1 < (double)n+(probe+radius)/h &&
				    z1 > 0.0-(probe+radius)/h &&
				    z1 < (double)o+(probe+radius)/h) {

				        /* Save coordinates (x,y,z), radius, residue number and chain */
				        insert_atom (x, y, z, radius, number, chain, convertRES (RES));
				    }

			}

		}

	}

	/* Close PDB file */
	fclose (arqPDB);

    /* Return flag indicating file has been read */
	return flag;

}

/* PDB_load2 function gets DIC object that contain residues linked lists (motif and radius), tablesize object that
contains number of residues, TABLE[i][j] object that contain all residues symbols, and path to PDB file*/
/* Create a linked list for PDB information and save coordinates (x,y,z)*/
int
PDB_load2(char PDB_NAME[NAME_MAX])
{
    /* Declare variables */
	int flag = 1, i;
	double x, y, z;
	char AUX[50] = "", X[50] = "", Y[50] = "", Z[50] = "", LINE[100], chain;
	FILE *arqPDB;

    /* NULL pointer to last item in PDB information linked list */
	v = NULL;

     /* Open PDB file */
	arqPDB = fopen (PDB_NAME,"r");

	/* PDB file not found */
	if (arqPDB == NULL) {

	    /* Print error */
		printf ("Reading Error: PDB file not found!\n");
		/* Close PDB file */
		fclose (arqPDB);
		/* Return flag indicating file has not been found */
		return !flag;

	}
	/* PDB file found */
	else {

		/* Parse PDB */
		/* While PDB file is not over, do... */
		while (get_line (arqPDB, LINE)) {

            /* Extract Record Name */
			extract (LINE, 0, 6, AUX);
			/* If Record Name is equal to ATOM or HETATM, do ... */
			if (!strcmp (AUX,"ATOM  ") || !strcmp (AUX,"HETATM")) {

                /* Extract x coordinate */
				extract (LINE, 30, 38, X);
				/* Extract y coordinate */
				extract (LINE, 38, 46, Y);
				/* Extract z axis coordinate */
				extract (LINE, 46, 54, Z);
				/* Remove whitespaces */
				trim(X, ' ');  trim(Y, ' ');  trim(Z, ' ');
				/* Convert coordinates from strings to doubles */
				convert (X, &x); convert (Y, &y); convert (Z, &z);

                /* Save coordinate (x,y,z) */
				insert_atom (x, y, z, 0.0, 0, 0, 0);

			}

		}

	}

	/* Close PDB file */
	fclose (arqPDB);

    /* Return flag indicating file has been read */
	return flag;

}

/* PDB_load3 gets path to PDB file */
/* Creates a linked list for PDB information */
int
PDB_load3(char PDB_NAME[NAME_MAX])
{
    /* Declare variables */
	int flag = 1, i, resnumber;
	double x, y, z;
	char AUX[50] = "", AUX2[50] = "", X[50] = "", Y[50] = "", Z[50] = "", LINE[100], chain;
	FILE *arqPDB;

    /* NULL pointer to last item in PDB information linked list */
	v = NULL;

    /* Open PDB file */
	arqPDB = fopen (PDB_NAME,"r");


	/* PDB file not found */
	if (arqPDB == NULL) {

	    /* Print error */
		printf ("Reading Error: PDB file not found!\n");
		/* Close PDB file */
		fclose (arqPDB);
		/* Return flag indicating file has not been found */
		return !flag;

	}
	else {

		/* Parse PDB */
		/* While PDB is not over, do ... */
		while (get_line (arqPDB, LINE)) {

            /* Extract Record Name */
			extract (LINE, 0, 6, AUX);
			/* If Record Name is equal to ATOM or HETATM, do ... */
			if (!strcmp (AUX,"ATOM  ") || !strcmp (AUX,"HETATM")) {

				/*Get residue sequence number*/
				resnumber = 0;
				extract (LINE, 22, 26, AUX2);
				resnumber = atoi (AUX2);

				/* Extract x coordinate */
				extract (LINE, 30, 38, X);
				/* Extract y coordinate */
				extract (LINE, 38, 46, Y);
				/* Extract z coordinate */
				extract (LINE, 46, 54, Z);

				/* Remove whitespaces */
				trim (X, ' ');  trim (Y, ' ');  trim (Z, ' ');
				/* Convert coordinates from strings to doubles */
				convert (X, &x); convert (Y, &y); convert (Z, &z);

                /* Extract chain identifier */
				chain = extractChain (LINE);

				/* Save coordinate (x,y,z), residue number and chain */
				insert_atom (x, y, z, 0.0, resnumber, chain, 0);

			}

		}

	}

    /* Close PDB file */
	fclose (arqPDB);

    /* Return flag indicating file has been read */
	return flag;

}
