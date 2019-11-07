/* This file contains all the functions used to process TOML parameters file used as input for KVFinder */

/* Import native modules */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Import custom modules */
#include "tomlprocessing.h"

/*Reads TOML file inside LINE*/
int
get_toml_line (FILE *arq,
           char LINE[500])
{
    /* Declare variables */
	int i;

	/* Read PDB file inside LINE[500] */
	for (i = 0, LINE[i] = getc (arq); LINE[i] != EOF && LINE[i] != '\n' && i < 500 ; i++, LINE[i] = getc(arq));

	/* If LINE[500] is not \n, read PDB file until LINE[100] assume \n or EOF value */
	if (i == 500 && LINE[i] != '\n')
	    for (; LINE[i] != '\n' && LINE[i] != EOF; LINE[i] = getc (arq));

    /* If PDB file is over, return False */
	if (LINE[i] == EOF)
	    return 0;
	/* If PDB file is not over, return True */
	else
	    return 1;

}

/* Extract from LINE[a] until LINE[b] to a char array S[100] of size b-a */
void
extract_toml_line (char LINE[500],
          int a,
          int b,
          char S[500])
{
    /* Declare variables */
	int i;

	/* Tests inputs */
	if (b > 500 || a > b)
	    return;

    /* Extract process */
	for (i = a; i < b; i++)
	    S[i - a] = LINE[i];

	/* Mark last position in S[] */
	if (i < 500)
	    S[i] = '\0';

}

/* Remove a char c from a string S[100] */
void
trim2 (char S[500],
       char c)
{
    /* Declare variables */
	int i,j;

	for (i = 0; S[i] != '\0'; i++) {
		if (S[i] == c)
		    for (j = i; S[j] != '\0'; j++)
		        S[j] = S[j+1];
	}

}

/* Fill string S[100] with 0s */
void
init500 (char S[500])
{
    /* Declare variables */
	int i;

	for (i = 0; i < 500; i++)
	    S[i] = '\0';

}

/* Fill vector S[2] with 0s */
void
init2 (int S[2])
{
    /* Declare variables */
	int i;

	for (i = 0; i < 2; i++)
	    S[i] = 0;

}

/* Fill vector S[4] with 0s */
void
init4 (int S[4])
{
    /* Declare variables */
	int i;

	for (i = 0; i < 4; i++)
	    S[i] = 0;

}

/* Convert true or false input to int input */
int
TF_input (char flag_in[6])
{

    if (strcmp (flag_in, "false") == 0 || strcmp (flag_in, "0") == 0)
        return 0;
	else
	    return 1;

}

/* Read TOML file and pass to struct TOML */
toml*
readTOML (toml *p,
          char *path)
{
    /* Declare variables */
	FILE *parameters_file;
	int i = 0, j = 0, m = 0, n = 0, max = 0;
	int equal[4], virg[2], chave[2], end, config_bar[2];
	char LINE[500], config[500] = "", value[500] = "", key[500] = "";
	char keys[3][500], values[3][500];

	/* Allocate memory to struct TOML */
	p = (toml*) malloc (sizeof (toml));

    /* Open parameters TOML file */
	parameters_file = fopen (path, "r");

    /* TOML file not found */
	if (parameters_file == NULL) {

	    /* Print error and exit */
		printf ("\033[0;31mError:\033[0m Invalid parameters file! Please select a valid parameters and try again.\n");
		exit (-1);

	}
	/* TOML file found */
	else {
	    /* While TOML file is not over, do... */
        while (get_toml_line (parameters_file, LINE)) {

			HANDLE_LAST_LINE:
			/* Remove tabs */
            trim2 (LINE, '\t');
            /* Remove whitespaces */
             trim2 (LINE, ' ');

            /* If a '#' is found, ignore line */
            if (LINE[0] == '#' || LINE[0] == ' ');
            else {

				for (max = 0; LINE[max] != '\n'; max++);

                /* Loop through all positions in LINE */
                for (j = 0; j < max; j++) {

                    /* Get positions of equal symbols in LINE */
                    if (LINE[j] == '=') {

                        equal[n] = j;
                        n++;

                    }

                    /* Get positions of comma symbols in LINE */
                    if (LINE[j] == ',') {

                        virg[m] = j;
                        m++;

                    }

                    /* Get position of newline symbol in LINE */
                    if (LINE[j] == '\n')
                        end = j;

                    /* Get position of '[' symbol in LINE */
                    if (LINE[j] == '[')
                        config_bar[0] = j;

                    /* Get position of ']' symbol in LINE */
                    if (LINE[j] == ']')
                        config_bar[1] = j;

                    /* Get position of '{' symbol in LINE */
                    if (LINE[j] == '{')
                        chave[0] = j;

                    /* Get position of '}' symbol in LINE */
                    if (LINE[j] == '}')
                        chave[1] = j;

                }

                /* If vector config_bar are filled, do ... */
                if (config_bar[1] != 0) {

                    /* Restart string config[100] */
                    init500 (config);
                    /* Extract name of configuration */
                    extract_toml_line (LINE, config_bar[0] + 1, config_bar[1], config);

                }

                /* If equal symbol is found, do ... */
                if (equal[0] != 0) {

                    /* If bracket is found, do ... */
                    if (chave[0] != 0 && chave[1] != 0) {

                        /* Extract first key inside keys[0] */
                        extract_toml_line (LINE, chave[0] + 1, equal[1], keys[0]);

                        /* Extract first value inside values[0] */
                        extract_toml_line (LINE, equal[1] + 1, virg[0], values[0]);

                        /* Extract second key inside keys[1] */
                        extract_toml_line (LINE, virg[0] + 1, equal[2], keys[1]);

                        /* Extract second value inside values[1] */
                        extract_toml_line (LINE, equal[2] + 1, virg[1], values[1]);

                        /* Extract third key inside keys[2] */
                        extract_toml_line (LINE, virg[1] + 1, equal[3], keys[2]);

                        /* Extract third value inside values[2] */
                        extract_toml_line (LINE, equal[3] + 1, chave[1], values[2]);

                        /* Save values inside struct TOML */
                        if (strcmp (keys[0], "bX1") == 0) {
                            (p->bX1) = atof (values[0]);
                            (p->bY1) = atof (values[1]);
                            (p->bZ1) = atof (values[2]);
                        }
                        if (strcmp (keys[0], "bX2") == 0) {
                            (p->bX2) = atof (values[0]);
                            (p->bY2) = atof (values[1]);
                            (p->bZ2) = atof (values[2]);
                        }
                        if (strcmp (keys[0], "bX3") == 0) {
                            (p->bX3) = atof (values[0]);
                            (p->bY3) = atof (values[1]);
                            (p->bZ3) = atof (values[2]);
                        }
                        if (strcmp (keys[0], "bX4") == 0) {
                            (p->bX4) = atof (values[0]);
                            (p->bY4) = atof (values[1]);
                            (p->bZ4) = atof (values[2]);
                        }
                        if (strcmp (keys[0], "X1") == 0) {
                            (p->X1) = atof (values[0]);
                            (p->Y1) = atof (values[1]);
                            (p->Z1) = atof (values[2]);
                        }
                        if (strcmp (keys[0], "X2") == 0) {
                            (p->X2) = atof (values[0]);
                            (p->Y2) = atof (values[1]);
                            (p->Z2) = atof (values[2]);
                        }
                        if (strcmp (keys[0], "X3") == 0) {
                            (p->X3) = atof (values[0]);
                            (p->Y3) = atof (values[1]);
                            (p->Z3) = atof (values[2]);
                        }
                        if (strcmp (keys[0], "X4") == 0) {
                            (p->X4) = atof (values[0]);
                            (p->Y4) = atof (values[1]);
                            (p->Z4) = atof (values[2]);
                        }

                    }
                    else {

                        /* Extract key */
                        extract_toml_line (LINE, 0, equal[0], key);

                        /* Extract value */
                        if (end < max) end = max;
                        extract_toml_line (LINE, equal[0] + 1, end, value);

                        /* Remove quotation marks */
                        trim2 (value, '\"');
                        trim2 (value, '\"');

                        /* Save values inside struct TOML */
                        if (strcmp (key, "output") == 0)
                            strcpy ((p->OUTPUT), value);
                        if (strcmp (key, "base_name") == 0)
                            strcpy( (p->BASE_NAME), value);
                        if (strcmp (key, "pdb") == 0)
                            strcpy ((p->PDB_NAME), value);
                        if (strcmp (key, "ligand") == 0)
                            strcpy ((p->LIGAND_NAME), value);
                        if (strcmp (key, "dictionary") == 0)
                            strcpy ((p->dictionary_name), value);
                        if (strcmp (key, "whole_protein_mode") == 0)
                            (p->whole_protein_mode) = TF_input (value);
                        if (strcmp (key, "resolution_mode") == 0)
                            strcpy( (p->resolution_flag), value);
                        if (strcmp (key, "box_mode") == 0)
                            (p->box_mode) = TF_input (value);
                        if (strcmp (key, "surface_mode") == 0)
                            (p->surface_mode) = TF_input (value);
                        if (strcmp (key, "kvp_mode") == 0)
                            (p->kvp_mode) = TF_input (value);
                        if (strcmp (key, "ligand_mode") == 0)
                            (p->ligand_mode) = TF_input (value);
                        if (strcmp (key, "step_size") == 0)
                            (p->h) = atof (value);
                        if (strcmp (key, "probe_in") == 0)
                            (p->probe_in) = atof (value);
                        if (strcmp (key, "probe_out") == 0)
                            (p->probe_out) = atof (value);
                        if (strcmp (key, "volume_cutoff") == 0)
                            (p->volume_cutoff) = atof (value);
                        if (strcmp (key, "ligand_cutoff") == 0)
                            (p->ligand_cutoff) = atof (value);
                        if (strcmp (key, "removal_distance") == 0)
                            (p->removal_distance) = atof (value);

                    }

                }

            }

            /* Restart all varibles */
            init500 (keys[0]);
            init500 (keys[1]);
            init500 (keys[2]);
            init500 (values[0]);
            init500 (values[1]);
            init500 (values[2]);
            init500 (LINE);
            init500 (key);
            init500 (value);
            init2 (chave);
            init2 (virg);
            init4 (equal);
            config_bar[0] = 0;
            config_bar[1] = 0;
            end = 0;
            n = 0;
            m = 0;

        }

	}

	/* Check if last line read found EOF in the middle of the line */
	for (i = 0; LINE[i] != EOF; i++);
	if (i > 0){
		LINE[i] = '\n';
		goto HANDLE_LAST_LINE;
	}

    /* Return struct TOML */
	return p;

}