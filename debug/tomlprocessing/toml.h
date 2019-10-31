/* This file contains all the functions used to process TOML parameters file used as input for KVFinder */

#ifndef TOML_H
#define TOML_H

typedef struct TOML {
	double X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4;
	double bX1, bY1, bZ1, bX2, bY2, bZ2, bX3, bY3, bZ3, bX4, bY4, bZ4;
	double h, probe_in, probe_out, volume_cutoff, ligand_cutoff, removal_distance;
	char PDB_NAME[500], LIGAND_NAME[500], dictionary_name[500], OUTPUT[500], BASE_NAME[500], resolution_flag[7];
	int whole_protein_mode, box_mode, surface_mode, kvp_mode, ligand_mode;

} toml;

toml *param;

toml* readTOML(toml *p, char *path);
int TF_input(char flag_in[6]);
void init2(int S[2]);
void init4(int S[4]);
void init100(char S[100]);
void trim2(char S[100], char c);
void extract2(char LINE[200], int a, int b, char S[100]);
int get_line2(FILE *arq, char LINE[200]);

#endif
