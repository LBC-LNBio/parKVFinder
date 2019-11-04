/* This file contains the functions used by KVFinder to process the user defined dictionaries, fundamental data to the
software */

/* Import native modules */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Import custom modules */
#include "dictionaryprocessing.h"
#include "resultsprocessing.h"

/* Remove a char c from a string S[50] */
void
trim (char S[50],
      char c)
{
    /* Declare variables */
	int i,j;

	for (i = 0; S[i] != '\0'; i++)
		if (S[i] == c)
		    for (j = i; S[j] != '\0'; j++)
		        S[j] = S[j+1];

}

/* Convert 3-letter residue system to 1-letter residue system */
char
convertRES (char RES[RES_SIZE])
{

	if (!strcmp (RES, "ALA"))
	    return 'A';
	if (!strcmp (RES, "ARG"))
	    return 'R';
	if (!strcmp (RES, "ASN"))
	    return 'N';
	if (!strcmp (RES, "ASP"))
	    return 'D';
	if (!strcmp (RES, "CYS"))
	    return 'C';
	if (!strcmp (RES, "GLN"))
	    return 'Q';
	if (!strcmp (RES, "GLU"))
	    return 'E';
	if (!strcmp (RES, "GLY"))
	    return 'G';
	if (!strcmp (RES, "HIS"))
	    return 'H';
	if (!strcmp (RES, "ILE"))
	    return 'I';
	if (!strcmp (RES, "LEU"))
	    return 'L';
	if (!strcmp (RES, "LYS"))
    	return 'K';
	if (!strcmp (RES, "MET"))
	    return 'M';
	if (!strcmp (RES, "PHE"))
	    return 'F';
	if (!strcmp (RES, "PRO"))
	    return 'P';
	if (!strcmp (RES, "SER"))
	    return 'S';
	if (!strcmp (RES, "THR"))
	    return 'T';
	if (!strcmp (RES, "TRP"))
	    return 'W';
	if (!strcmp (RES, "TYR"))
	    return 'Y';
	if (!strcmp (RES, "VAL"))
	    return 'V';

	return 'X';

}

/* Read dictionary file and saves residue names inside TABLE matrix and return number of residues */
int
define_table (char TABLE[TABLE_SIZE][RES_SIZE],
              char dictionary_name[DIC_NAME_MAX])
{
    /* Declare variables */
	int i = 0, j = 0;
	char AUX[50];
	FILE *arq;

    /* Open dictionary file */
	arq = fopen (dictionary_name, "r");

	/* File has not been found */
	if (arq == NULL)
	    printf ("Reading Error: Residues dictionary file not found!\n");
	/* File has been found */
	else
	    /* While EOF not read in AUX object, do ... */
	    while (fscanf (arq, "%s", AUX) != EOF) {
			if (AUX[0] == '>') {
                /* Copy each letter of residue symbol inside TABLE */
				for (j = 1; j < strlen(AUX); j++)
				    TABLE[i][j - 1] = AUX[j];
				/* Put each symbol at the end of residue symbol */
				TABLE[i][j - 1] = '\0';
				/* Increment line iterator */
				i++;
			}
	}

	/* Close dictionary file */
	fclose (arq);

    /* Return number of residues in TABLE */
	return i;

}

/* Load a VdW radius linked. Each item contain the radius of an atom of a specific residue. */
int
dictionary_load (dict *DIC[TABLE_SIZE],
                 int tablesize,
                 char dictionary_name[DIC_NAME_MAX])
{
    /* Declare variables */
	int i, flag = 1;
	char AUX[50];
	FILE *dictionary_file;
	/* Dictionary structure: {symbol, radius, *next} */
	dict *p;

     /* Open dictionary file */
	dictionary_file = fopen (dictionary_name, "r");
	if (dictionary_file == NULL) {

	    /* Print error */
		printf ("Dictionary reading error! Please select a valid dictionary filename and try again.\n");
		 /* Close dictionary file */
		fclose (dictionary_file);

		/* Return flag indicating file has not been found */
		return !flag;

	}

    /* from zero to tablesize object, attribute NULL to each point of DIC */
	for (i = 0; i < tablesize; i++)
	    DIC[i] = NULL;

    /* Read lines in dictionary file until reaches EOF and i < tablesize */
	for (i = -1; fscanf (dictionary_file, "%s", AUX) != EOF && i < tablesize; ) {

		if (AUX[0] == '>')
		    /* Count a residue */
		    i++;
		else {

		    /* Allocate a dict space for p in memory */
			p = malloc (sizeof (dict));
			/* Eliminate whitespaces between symbol and its radius value */
			trim (AUX, ' ');
			/* Read radius value and save inside struct p inside radius */
			fscanf (dictionary_file, "%lf", &p->radius);
			/* Copy AUX string and paste inside struct p inside symbol */
			strcpy(p->symbol, AUX);
			/* NULL pointer to last item inserted in linked list */
			p->next = NULL;

			/* Create a linked list for accessing each residue information i indicates the residue. So, DIC[i] is the
			linked list of the residue i */
			if (DIC[i] == NULL)
			    /* Allocate a dict space for DIC[i] in memory */
			    DIC[i] = malloc (sizeof (dict));
			else
			    p->next = DIC[i]->next;
			DIC[i]->next = p;

		}

	}

	/* Close dictionary file */
	fclose (dictionary_file);

    /* Return flag indicating file has been found */
	return flag;

}

/* dictionary_consult_radius gets RES object that contains residues names from PDB file, ATOM_TYPE contains atom name,
tablesize that contains number of residues in dictionary file, TABLE[i][j] object that contains all residues symbols,
ATOM_SYMBOL contains atom symbol */
/* Return the radius of an atom of a residue */
double
dictionary_consult_radius (char RES[RES_SIZE],
                           char ATOM_TYPE[ATOM_SIZE],
                           dict *DIC[TABLE_SIZE],
                           int tablesize,
                           char TABLE[TABLE_SIZE][RES_SIZE],
                           char ATOM_SYMBOL[2],
                           FILE **log_file)
{
    /* Declare variables */
	int i;
	double value = 0.0;
	char ATOM_TYPE_AUX[ATOM_SIZE], RES_AUX[RES_SIZE];
	dict *p;

	for (i = 0; i < tablesize && value == 0.0; i++) {

	    /* If TABLE[i] residue is equal to RES residue, do ... */
		if (!strcmp (RES,TABLE[i])) {

			/* Loop inside residue linked list and look for ATOM name, then retrieve radius to value */
			for (p = DIC[i]->next; p != NULL && value == 0.0; p = p->next)
				if (!strcmp (ATOM_TYPE, p->symbol)) {

				    /* Retrieve radius of an atom of a residue */
					value = p->radius;
					/* Return radius value */
					return value;

				}

		}
	}

    /* If radius not found inside linked list, do ... */
	if (!value) {

		fprintf (*log_file, "Warning: Atom radius not found in dictionary:%s %s!\n", ATOM_TYPE, RES);
		strcpy (RES_AUX, "GEN");

		if (ATOM_SYMBOL[0] > 90 || ATOM_SYMBOL[0] < 65) {

			ATOM_TYPE_AUX[0] = ATOM_SYMBOL[1];
			ATOM_TYPE_AUX[1] = '\0';

		}
		else {

			ATOM_TYPE_AUX[0] = ATOM_SYMBOL[0];
			ATOM_TYPE_AUX[1] = ATOM_SYMBOL[1];
			ATOM_TYPE_AUX[2] = '\0';

		}

		/* Loop inside generic residue linked list and look for an atom, then retrieve its radius to value */
		for (i = 0; i < tablesize && value == 0.0; i++) {

			if (!strcmp (RES_AUX, TABLE[i])) {

				for (p = DIC[i]->next; p != NULL && value == 0.0; p = p->next) {

					if (!strcmp (ATOM_TYPE_AUX, p->symbol)) {

					    /* If radius value is found, print generic atom and its radius value */
						value = p->radius;
						/* Print warning in log file */
						fprintf(*log_file,
						        "Warning: Using generic atom %s radius value %.2lf\n",
						        ATOM_TYPE_AUX,
						        p->radius);
						/* Return radius value */
						return value;

					}

				}

			}

		}

	/* Print warning in log file */
	fprintf (*log_file,
	         "Warning: Radius data not found for atom %s. This atom will be excluded from analysis.\n",
	         ATOM_TYPE_AUX);

	}

    /* Atom excluded from analysis and return radius 0.0 */
	return value;

}


/* Set Vvoxel, step size and resolution_mode based on resolution_flag (Low, Medium, High, Off) */
void
resolution_input (char flag_in[7],
                  double* Vvoxel,
                  int* resolution_mode,
                  double* h)
{
	/* Declare variables */
	*resolution_mode = 1;

	switch (flag_in[0]){

		case 'L':
		    *Vvoxel = 0.2;
		    *h = 0.1;
		    break;

		case 'M':
		    *Vvoxel = 0.1;
		    *h = 0.1;
		    break;

		case 'H':
		    *Vvoxel = 0.01;
		    *h = 0.1;
		    break;

		default:
		    *Vvoxel = (*h)*(*h)*(*h);
		    *resolution_mode = 0;
		    break;

	}
}
