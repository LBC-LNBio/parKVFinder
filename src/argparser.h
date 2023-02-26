#ifndef ARGPARSER_H
#define ARGPARSER_H

void create_custom_box(char *box_name, double *Xmin, double *Xmax, double *Ymin,
                       double *Ymax, double *Zmin, double *Zmax);
void create_residues_box(char *box_name, double *Xmin, double *Xmax,
                         double *Ymin, double *Ymax, double *Zmin, double *Zmax,
                         double padding, char PDB_NAME[500]);
void print_version();
void print_header();
void print_usage();
void print_options();
void print_help();
int argparser(int argc, char **argv, int *box_mode, int *kvp_mode,
              int *ligand_mode, int *surface_mode, int *whole_protein_mode,
              char PDB_NAME[500], char LIGAND_NAME[500],
              char dictionary_name[500], char OUTPUT[500], char BASE_NAME[500],
              char resolution_flag[7], double *h, double *probe_in,
              double *probe_out, double *volume_cutoff, double *ligand_cutoff,
              double *removal_distance, double *X1, double *Y1, double *Z1,
              double *X2, double *Y2, double *Z2, double *X3, double *Y3,
              double *Z3, double *X4, double *Y4, double *Z4, double *bX1,
              double *bY1, double *bZ1, double *bX2, double *bY2, double *bZ2,
              double *bX3, double *bY3, double *bZ3, double *bX4, double *bY4,
              double *bZ4);
int check_input(char *optarg, char *error);

#endif
