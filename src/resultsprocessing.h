/* This file contains the functions used by KVFinder to process the user defined
dictionaries, fundamental data to the software*/

#ifndef RESULTSPROCESSING_H
#define RESULTSPROCESSING_H

typedef struct COORDINATES{
	double Xmin;
	double Xmax;
	double Ymin;
	double Ymax;
	double Zmin;
	double Zmax;
} coords;

typedef struct RESIDUES_INFORMATION {
  int resnum;
  char resname;
  char chain;
  struct RESIDUES_INFORMATION *next;
} residues_info;

typedef struct KVFINDER_RESULTS {
  /* Geometric properties */
  double volume;
  double area;
  double max_depth;
  double avg_depth;
  /* Residues */
  residues_info *res_info;
} KVresults;

/* Declare structs */
residues_info *t;
KVresults *KVFinder_results;
coords* kvcoords;
coords* frontiercoords;

/* Define custom functions */
void insert_res(int resnum, char chain, char resname, int kvnum);
void write_results(char *output_results, char *pdb_name, char *output_pdb,
                   char LIGAND_NAME[500], char resolution_flag[7],
                   char step_flag[6], int ncav);

#endif
