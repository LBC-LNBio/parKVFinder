/* This file contains all the command line functions utilized by KVFinder */

#ifndef ARGPARSER_H
#define ARGPARSER_H

/* Define custom functions */
char TF (int mode,
         char **tomlmode);
char *get_file_extension (char *filename);
void create_custom_box (char *box_name,
                        double *Xmin,
                        double *Xmax,
                        double *Ymin,
                        double *Ymax,
                        double *Zmin,
                        double *Zmax);
void create_residues_box (char *box_name,
                          double *Xmin,
                          double *Xmax,
                          double *Ymin,
                          double *Ymax,
                          double *Zmin,
                          double *Zmax,
                          double padding,
                          char PDB_NAME[NAME_MAX]);
void init25 (char S[25]);
void init75 (char S[75]);
void init80 (char S[80]);
void format_variable (char* variable,
                      int a,
                      int b,
                      char AUX[25]);
void format_text (char* description,
                  int a,
                  int b,
                  char AUX[75]);
void format_title (char *title,
                   int a,
                   int b,
                   char AUX[100]);
void print_line (char *title);
void print_arguments (char *variable,
                      char *description,
                      char *extra);
void print_version ();
void print_header ();
void print_usage();
void print_options();
void print_help();
void print_toml (char *toml_name,
                 char OUTPUT[500],
                 char BASE_NAME[500],
                 char dictionary_name[500],
                 char PDB_NAME[500],
                 char LIGAND_NAME[500],
                 int whole_protein_mode,
                 char resolution_mode[7],
                 int box_mode,
                 int surface_mode,
                 int kvp_mode,
                 int ligand_mode,
                 double h,
                 double probe_in,
                 double probe_out,
                 double volume_cutoff,
                 double ligand_cutoff,
                 double removal_distance,
                 double X1,
                 double Y1,
                 double Z1,
                 double X2,
                 double Y2,
                 double Z2,
                 double X3,
                 double Y3,
                 double Z3,
                 double X4,
                 double Y4,
                 double Z4,
                 double bX1,
                 double bY1,
                 double bZ1,
                 double bX2,
                 double bY2,
                 double bZ2,
                 double bX3,
                 double bY3,
                 double bZ3,
                 double bX4,
                 double bY4,
                 double bZ4);
int argparser (int argc,
               char **argv,
               int *box_mode,
               int *kvp_mode,
               int *ligand_mode,
               int *surface_mode,
               int *whole_protein_mode,
               char PDB_NAME[500],
               char LIGAND_NAME[500],
               char dictionary_name[500],
               char OUTPUT[500],
               char BASE_NAME[500],
               char resolution_flag[7],
               double *h,
               double *probe_in,
               double *probe_out,
               double *volume_cutoff,
               double *ligand_cutoff,
               double *removal_distance,
               double *X1,
               double *Y1,
               double *Z1,
               double *X2,
               double *Y2,
               double *Z2,
               double *X3,
               double *Y3,
               double *Z3,
               double *X4,
               double *Y4,
               double *Z4,
               double *bX1,
               double *bY1,
               double *bZ1,
               double *bX2,
               double *bY2,
               double *bZ2,
               double *bX3,
               double *bY3,
               double *bZ3,
               double *bX4,
               double *bY4,
               double *bZ4);
int check_input (char *optarg,
                 char *error);

#endif
