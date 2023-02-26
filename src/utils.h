#ifndef UTILS_H
#define UTILS_H

/* Structs */

/*
 * Struct: PARAMETERS
 * ------------------
 *
 * A struct containing parKVFinder parameters
 *
 */
typedef struct PARAMETERS {
  double X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4,
      Z4; /* Internal box vertices */
  double bX1, bY1, bZ1, bX2, bY2, bZ2, bX3, bY3, bZ3, bX4, bY4,
      bZ4; /* Visible box vertices */
  double h, probe_in, probe_out, volume_cutoff, ligand_cutoff,
      removal_distance; /* Detection parameters */
  char PDB_NAME[500], LIGAND_NAME[500], dictionary_name[500], OUTPUT[500],
      BASE_NAME[500], resolution_flag[7]; /* File paths */
  int whole_protein_mode, box_mode, surface_mode, kvp_mode,
      ligand_mode; /* Detection modes */

} parameters;

/*
 * Struct: VDW
 * -----------
 *
 * A struct containing atom name and radius
 *
 * radius: atom radius
 * symbol: atom name
 * next: pointer to the next VDW struct
 *
 */
typedef struct VDW {
  double radius;
  char symbol[6];
  struct VDW *next;
} vdw;

/*
 * Struct: ATOM
 * ------------
 *
 * A struct containing atomic information
 *
 * x: X-axis coordinate
 * y: Y-axis coordinate
 * z: Z-axis coordinate
 * radius: atom radius
 * resnumber: residue number
 * resname: residue name
 * chain: chain identifier
 * next: pointer to next ATOM struct
 *
 */
typedef struct ATOM {
  double x;
  double y;
  double z;
  double radius;
  int resnumber;
  char resname;
  char chain;
  struct ATOM *next;
} atom;

/*
 * Struct: COORDINATES
 * -------------------
 *
 * A struct contatining box coordinates
 *
 * Xmin: X-axis origin coordinate
 * Ymin: Y-axis origin coordinate
 * Zmin: Z-axis origin coordinate
 * Xmax: X-axis maximum coordinate
 * Ymax: Y-axis maximum coordinate
 * Zmax: Z-axis maximum coordinate
 *
 */
typedef struct COORDINATES {
  double Xmin;
  double Xmax;
  double Ymin;
  double Ymax;
  double Zmin;
  double Zmax;
} coords;

/*
 * Struct: RESIDUES_INFORMATION
 * ----------------------------
 *
 * A struct containing residue information (residue number, residue name and
 * chain identifier)
 *
 * resnum: residue number
 * resname: residue name
 * chain: chain identifier
 * next: pointer to the next RESIDUES_INFORMATION struct
 *
 */
typedef struct RESIDUES_INFORMATION {
  int resnumber;
  char resname;
  char chain;
  struct RESIDUES_INFORMATION *next;
} residues_info;

/*
 * Struct: KVFINDER_RESULTS
 * ------------------------
 *
 * A struct containing parKVFinder results
 *
 * volume: cavities volume
 * area: cavities area
 * max_depth: cavities maximum depth
 * avg_depth: cavities average depth
 * res_info: interface residues information (residue number, residue name and
 * chain identifier)
 *
 */
typedef struct KVFINDER_RESULTS {
  double volume;
  double area;
  double max_depth;
  double avg_depth;
  double avg_hydropathy;
  residues_info *res_info;
} KVresults;

/*
 * Struct: node
 * ------------
 *
 * A struct containing volume information
 * 
 * volume: cavity volume
 * pos: cavity index
 * struct node* next: pointer to next linked list node
 *
 */
typedef struct NODE {
  double volume;
  int pos;
  struct NODE *next;
} node;

/* Global variables */
double sina, sinb, cosa, cosb;
int big, volume;
atom *v;
node *V;
residues_info *t;
KVresults *KVFinder_results;
coords *cavity, *boundary;

/* Functions */
double max(double a, double b);
double min(double a, double b);
double _resolution2step(char flag[]);
char *_combine(const char *s1, const char *s2);
char _residue2code(char RESIDUE[]);
char *_code2residue(char RESIDUE);
char *_get_file_extension(char *fn);
int _int2toml(int boolean, char **flag);
int _read_line(FILE *arq, char LINE[], int size);
int _toml2int(char flag[6]);
void _extract(char FROM[], int nF, char TO[], int nT, int start, int end);
void _initialize_string(char FROM[], int nF);
void _remove_char(char FROM[], int nF, char c);

#endif
