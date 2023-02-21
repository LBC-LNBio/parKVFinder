/* This file contains all the functions used to process TOML parameters file
 * used as input for KVFinder */

#ifndef TOMLPROCESSING_H
#define TOMLPROCESSING_H

#define DIC_NAME_MAX 500
#define NAME_MAX 500

typedef struct TOML {
  double X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4;
  double bX1, bY1, bZ1, bX2, bY2, bZ2, bX3, bY3, bZ3, bX4, bY4, bZ4;
  double h, probe_in, probe_out, volume_cutoff, ligand_cutoff, removal_distance;
  char PDB_NAME[NAME_MAX], LIGAND_NAME[NAME_MAX], dictionary_name[DIC_NAME_MAX],
      OUTPUT[NAME_MAX], BASE_NAME[NAME_MAX], resolution_flag[7];
  int whole_protein_mode, box_mode, surface_mode, kvp_mode, ligand_mode;

} toml;

toml *param;

toml *readTOML(toml *p, char *path);
int TF_input(char flag_in[6]);
void init500(char S[500]);
void extract_toml_line(char LINE[500], int a, int b, char S[500]);
int get_toml_line(FILE *arq, char LINE[500]);

#endif
