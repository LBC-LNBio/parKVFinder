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

	    /* Print error */
		printf ("Parameters reading error! Please select a valid parameters and try again.\n");
		/* Close TOML file */
		fclose (parameters_file);
		/* Exit code and return -1 */
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

//int
//main(
//int argc,
//char *argv[])
//{
//
//	/* Declare variables */
//
//	double h, probe_in, probe_out, volume_cutoff, ligand_cutoff, removal_distance;
//	double X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4;
//	double bX1, bY1, bZ1, bX2, bY2, bZ2, bX3, bY3, bZ3, bX4, bY4, bZ4;
//	int ligand_mode, surface_mode, whole_protein_mode, box_mode, kvp_mode;
//	int m, n, o;
//	char PDB_NAME[500], LIGAND_NAME[500], dictionary_name[500], OUTPUT[500],
//	BASE_NAME[500], resolution_flag[7];
//
//	/* Read parameters TOML files */
//	toml *parameters = readTOML (param, "parameters.toml");
//
//	/* Save TOML parameters from struct TOML parameters to KVFinder variables */
//	X1 = parameters->X1; Y1 = parameters->Y1; Z1 = parameters->Z1;
//	X2 = parameters->X2; Y2 = parameters->Y2; Z2 = parameters->Z2;
//	X3 = parameters->X3; Y3 = parameters->Y3; Z3 = parameters->Z3;
//	X4 = parameters->X4; Y4 = parameters->Y4; Z4 = parameters->Z4;
//	bX1 = parameters->bX1; bY1 = parameters->bY1; bZ1 = parameters->bZ1;
//	bX2 = parameters->bX2; bY2 = parameters->bY2; bZ2 = parameters->bZ2;
//	bX3 = parameters->bX3; bY3 = parameters->bY3; bZ3 = parameters->bZ3;
//	bX4 = parameters->bX4; bY4 = parameters->bY4; bZ4 = parameters->bZ4;
//	strcpy (PDB_NAME, parameters->PDB_NAME);
//	strcpy (OUTPUT, parameters->OUTPUT);
//	strcpy (BASE_NAME, parameters->BASE_NAME);
//	strcpy (LIGAND_NAME, parameters->LIGAND_NAME);
//	strcpy (dictionary_name, parameters->dictionary_name);
//	strcpy (resolution_flag, parameters->resolution_flag);
//	volume_cutoff = parameters->volume_cutoff;
//	ligand_cutoff = parameters->ligand_cutoff;
//	removal_distance = parameters->removal_distance;
//	probe_in = parameters->probe_in;
//	probe_out = parameters->probe_out;
//	h = parameters->h;
//	whole_protein_mode = parameters->whole_protein_mode;
//	box_mode = parameters->box_mode;
//	surface_mode = parameters->surface_mode;
//	kvp_mode = parameters->kvp_mode;
//	ligand_mode = parameters->ligand_mode;
//
//	if (argc >= 2){
//		int flag = atoi (argv[1]);
//		if (flag) {
//			char* arr[2] = {"false", "true"};
//			printf("\n\n\n[FILES_PATH]\n");
//			printf("dictionary = \"%s\"\n", dictionary_name);
//			printf("pdb = \"%s\"\n", PDB_NAME);
//			printf("output = \"%s\"\n", OUTPUT);
//			printf("base_name = \"%s\"\n", BASE_NAME);
//			printf("ligand = \"%s\"\n", LIGAND_NAME);
//			printf("[SETTINGS]\n");
//			printf("\t[SETTINGS.modes]]\n");
//			printf("\twhole_protein_mode = %s\n", arr[whole_protein_mode]);
//			printf("\tbox_mode = %s\n", arr[box_mode]);
//			printf("\tsurface_mode = %s\n", arr[surface_mode]);
//			printf("\tkvp_mode = %s\n", arr[kvp_mode]);
//			printf("\tligand_mode = %s\n", arr[ligand_mode]);
//			printf("\t[SETTINGS.step_size]\n");
//			printf("\tstep_size = %lf\n", h);
//			printf("\t[SETTINGS.probes]\n");
//			printf("\tprobe_in = %lf\n", probe_in);
//			printf("\tprobe_out = %lf\n", probe_out);
//			printf("\t[SETTINGS.cutoffs]\n");
//			printf("\tvolume_cutoff = %lf\n", volume_cutoff);
//			printf("\tligand_cutoff = %lf\n", ligand_cutoff);
//			printf("\tremoval_distance = %lf\n", removal_distance);
//			printf("\t[SETTINGS.visiblebox]\n");
//			printf("\tbP1 = {bX1 = %.2lf, bY1 = %.2lf, bZ1 = %.2lf}\n", bX1, bY1, bZ1);
//			printf("\tbP2 = {bX2 = %.2lf, bY2 = %.2lf, bZ2 = %.2lf}\n", bX2, bY2, bZ2);
//			printf("\tbP3 = {bX3 = %.2lf, bY3 = %.2lf, bZ3 = %.2lf}\n", bX3, bY3, bZ3);
//			printf("\tbP4 = {bX4 = %.2lf, bY4 = %.2lf, bZ4 = %.2lf}\n", bX4, bY4, bZ4);
//			printf("\t[SETTINGS.internalbox]\n");
//			printf("\tP1 = {X1 = %.2lf, Y1 = %.2lf, Z1 = %.2lf}\n", X1, Y1, Z1);
//			printf("\tP2 = {X2 = %.2lf, Y2 = %.2lf, Z2 = %.2lf}\n", X2, Y2, Z2);
//			printf("\tP3 = {X3 = %.2lf, Y3 = %.2lf, Z3 = %.2lf}\n", X3, Y3, Z3);
//			printf("\tP4 = {X4 = %.2lf, Y4 = %.2lf, Z4 = %.2lf}\n", X4, Y4, Z4);
//	}
//
//
//
//	}
//
//return 0;
//
//}