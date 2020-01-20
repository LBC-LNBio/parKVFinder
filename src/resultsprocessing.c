/* This file contains the functions used by KVFinder to process the user defined dictionaries, fundamental data to the
software*/

/* Import native modules */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Import custom modules */
#include "dictionaryprocessing.h"
#include "resultsprocessing.h"

/* Insert residue information (resnum, chain, resname, resclass) in linked list (res_info) inside KVresults structure */
void
insert_res (int resnum,
            char chain,
            char resname,
            int kvnum)
{
	/* Declare variables */
	residues_info *res_info;

	/* Special case I: Create data structure */
	if (KVFinder_results[kvnum].res_info == NULL) {

		/* Allocate memory */
		KVFinder_results[kvnum].res_info = malloc (sizeof (residues_info));

		/* Save residue information */
		KVFinder_results[kvnum].res_info->resnum = resnum;
		KVFinder_results[kvnum].res_info->chain = chain;
		KVFinder_results[kvnum].res_info->resname = resname;

		/* Identify head */
		KVFinder_results[kvnum].res_info->next = NULL;

	}
	else {

		/* Allocate memory for generic structure */
		res_info = malloc (sizeof (residues_info));

		/* Save residue information in generic structure*/
		res_info->resnum = resnum;
		res_info->chain = chain;
		res_info->resname = resname;

		/* Special case II: Insertion in the head end */
		if (KVFinder_results[kvnum].res_info->resnum >= resnum) {

			/* Pass generic structure to linked list */
			res_info->next = KVFinder_results[kvnum].res_info;
			KVFinder_results[kvnum].res_info = res_info;

		}
		else {

			/* Locate node before the point of insertion */
			while (KVFinder_results[kvnum].res_info->next != NULL && KVFinder_results[kvnum].res_info->resnum < resnum)
				KVFinder_results[kvnum].res_info = KVFinder_results[kvnum].res_info->next;

			/* Pass generic structure to linked list */
			res_info->next = KVFinder_results[kvnum].res_info->next;
			KVFinder_results[kvnum].res_info->next = res_info;

		}
	}
}

void
write_results (char *output_results,
               char *pdb_name,
               char *output_pdb,
               char LIGAND_NAME[500],
               char resolution_flag[7],
               char step_flag[6],
               int ncav)
{

	/* Declare variables */
	FILE *results_file;
	char results[1024];
	int kvnum, iterator;

    /* Open KVFinder.results.toml */#pragma omp for schedule(dynamic)
	results_file = fopen (output_results, "w");
	/* Save memory for buffer */
	memset (results, '\0', sizeof (results));
	/* Create buffer */
	setvbuf (results_file, results, _IOFBF, 1024);

	/* Write results file */
	/* File header */
	fprintf (results_file, "# TOML results file for parKVFinder software\n\ntitle = \"parKVFinder results file\"\n\n");

	/* Files paths */
	fprintf (results_file,
	         "[FILES_PATH]\nINPUT = \"%s\"\nOUTPUT = \"%s\"\nLIGAND = \"%s\"\n\n",
	         pdb_name,
	         output_pdb,
	         LIGAND_NAME);

	/* Parameters */
	fprintf (results_file, "[PARAMETERS]\nRESOLUTION = \"%s\"\nSTEP_SIZE = %s\n\n", resolution_flag, step_flag);

	/* Results header */
	fprintf (results_file, "[RESULTS]\n# Volume, area and interface residues information for each cavity\n");

	/* Volume */
	fprintf (results_file, "\n\t[RESULTS.VOLUME]\n\t# Volume unit is cubic angstrom\n");
	for (kvnum = 0; kvnum < ncav; kvnum++) {

		fprintf (results_file,
		         "\tK%c%c = %.2lf\n",
		         65+(((kvnum)/26)%26),
		         65+((kvnum)%26),
		         KVFinder_results[kvnum].volume);

	}

	/* Surface Area */
	fprintf (results_file, "\n\t[RESULTS.AREA]\n\t# Area unit is square angstrom\n");
	for (kvnum = 0; kvnum < ncav; kvnum++) {

		fprintf (results_file,
		         "\tK%c%c = %.2lf\n",
		         65+(((kvnum)/26)%26),
		         65+((kvnum)%26),
		         KVFinder_results[kvnum].area);

	}

	/* Interface Residues */
	fprintf (results_file, "\n\t[RESULTS.RESIDUES]\n\t# Interface residues for each cavity\n");
	fprintf (results_file, "\t# [\"residue number\",\"chain identifier\",\"residue name\"]\n");
	for (kvnum = 0; kvnum < ncav; kvnum++) {

		fprintf (results_file,
		         "\tK%c%c = [",
		         65+(((kvnum)/26)%26),
		         65+((kvnum)%26));

		for (t = KVFinder_results[kvnum].res_info; t != NULL; t = t->next) {

			if (t->next != NULL)
			    fprintf (results_file,
			             "[\"%d\",\"%c\",\"%c\"],",
			             t->resnum,
			             t->chain,
			             t->resname);

			else
			    fprintf (results_file,
			             "[\"%d\",\"%c\",\"%c\"]",
			             t->resnum,
			             t->chain,
			             t->resname);

		}

		fprintf (results_file, "]\n");

	}

    /* Print KVFinder.results.toml */
	fflush(results_file);

	/* Close KVFinder.results.toml */
	fclose(results_file);

}