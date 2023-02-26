#ifndef FILEPROCESSING_H
#define FILEPROCESSING_H

#include "utils.h"

/* parKVFinder parameter file processing */
void extract_toml_line(char LINE[500], int a, int b, char S[500]);
int get_toml_line(FILE *arq, char LINE[500]);
parameters *readTOML(char *path);
void write_parameters(char *toml_name, char OUTPUT[500], char BASE_NAME[500],
                      char dictionary_name[500], char PDB_NAME[500],
                      char LIGAND_NAME[500], int whole_protein_mode,
                      char resolution_mode[7], int box_mode, int surface_mode,
                      int kvp_mode, int ligand_mode, double h, double probe_in,
                      double probe_out, double volume_cutoff,
                      double ligand_cutoff, double removal_distance, double X1,
                      double Y1, double Z1, double X2, double Y2, double Z2,
                      double X3, double Y3, double Z3, double X4, double Y4,
                      double Z4, double bX1, double bY1, double bZ1, double bX2,
                      double bY2, double bZ2, double bX3, double bY3,
                      double bZ3, double bX4, double bY4, double bZ4);

/* van der Waals file processing */
double _get_vdw_radius(char RESIDUE[4], char ATOM_TYPE[4], vdw *DIC[500],
                       int tablesize, char TABLE[500][4], char ATOM_SYMBOL[2],
                       FILE **log_file);
int _get_residues_information(char dictionary_name[500], char TABLE[500][4]);
int read_vdw(char dictionary_name[500], vdw *DIC[500], int tablesize);

/* Protein DataBank (PDB) file processing */
void _insert_atom(atom **head, atom *new);
atom *_create_atom(double x, double y, double z, double radius, int resnumber,
                   char resname, char chain);
int soft_read_pdb(char PDB_NAME[500], int has_resnum, int has_chain);
int read_pdb(char PDB_NAME[500], vdw *DIC[500], int tablesize,
             char TABLE[500][4], double probe, int m, int n, int o, double h,
             double X1, double Y1, double Z1, FILE **log_file);
void _free_atom();

/* parKVFinder results file processing */
void write_results(char *output_results, char *pdb_name, char *output_pdb,
                   char LIGAND_NAME[500], double h, int ncav);

#endif