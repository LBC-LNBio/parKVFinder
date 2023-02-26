#ifndef GRIDPROCESSING_H
#define GRIDPROCESSING_H

/* Grid initialization */
int ***igrid(int m, int n, int o);
double ***dgrid(int m, int n, int o);

/* Molecular representation */
int check_protein_neighbours(int ***A, int i, int j, int k, int m, int n,
                             int o);
void SAS(int ***A, int m, int n, int o, double h, double probe, double X1,
         double Y1, double Z1);
void SES(int ***A, int m, int n, int o, double h, double probe);

/* Cavity detection (Probe In - Probe Out) */
void subtract(int ***A, int ***S, int m, int n, int o, double h,
              double removal_distance);
void filter_noise(int ***A, int m, int n, int o);

/* Ligand adjustment */
void adjust2ligand(int ***A, int m, int n, int o, double h, double limit,
                   double X1, double Y1, double Z1);

/* Box adjustment */
void filter2box(int ***A, int m, int n, int o, double h, double bX1, double bY1,
                double bZ1, double bX2, double bY2, double bZ2, double norm1);

/* Cavity clustering and volume estimation */
int check_unclustered_neighbours(int ***A, int m, int n, int o, int i, int j,
                                 int k);
void remove_cavity(int ***A, int m, int n, int o, int tag);
void DFS(int ***A, int m, int n, int o, int i, int j, int k, int tag);
int clustering(int ***A, int m, int n, int o, double h, double volume_cutoff);

/* Cavity surface and area estimation */
int define_surface_points(int ***A, int m, int n, int o, int i, int j, int k);
void filter_surface(int ***A, int ***S, int m, int n, int o);
double check_voxel_class(int ***S, int i, int j, int k);
void area(int ***S, int m, int n, int o, double h, int ncav);

/* Constitutional characterization */
residues_info *_create_residue(int resnumber, char resname, char chain);
void _insert_residue(residues_info **head, residues_info *new);
void _remove_duplicate_residue(residues_info *head);
void interface(int ***A, int m, int n, int o, double h, double probe, int ncav,
               double X1, double Y1, double Z1);

/* Cavity boundary and depth estimation */
int define_boundary_points(int ***A, int m, int n, int o, int i, int j, int k);
void filter_boundary(int ***A, int m, int n, int o, int ncav);
void remove_boundary(int ***A, int m, int n, int o, int ncav);
void depth(int ***A, double ***M, int m, int n, int o, double h, int ncav);

/* Cavity hydropathy */
double get_hydrophobicity_value(char *resname, char *resn[], double *scale);
void project_hydropathy(double ***HP, int ***S, int m, int n, int o, double h,
                        double probe, double X1, double Y1, double Z1);
void estimate_average_hydropathy(double ***HP, int ***S, int m, int n, int o,
                                 int ncav);

/* Export cavity PDB file */
int _filter_cavity(int ***A, int m, int n, int o, int i, int j, int k);
void export(char *output_pdb, int ***A, int ***S, double ***M, double ***HP,
            int kvp_mode, int m, int n, int o, double h, int ncav, double X1,
            double Y1, double Z1);

/* Clean memory */
void free_igrid(int ***A, int m, int n, int o);
void free_dgrid(double ***M, int m, int n, int o);
void free_node();

#endif
