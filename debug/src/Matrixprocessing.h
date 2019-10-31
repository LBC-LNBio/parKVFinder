#ifndef MATRIXPROCESSING_H
#define MATRIXPROCESSING_H

#define NAME_MAX 500

/* Define type node for a NODE structure */
typedef struct NODE {
	double volume;
	int pos;
	struct NODE *next;
} node;

/* Declare structs */
node *V;
double volume;
int flagr;

/* Define custom functions */
double min (double a,
            double b);
double max (double a,
            double b);
int check_pos (int ***A,
               int i,
               int j,
               int k,
               int m,
               int n,
               int o);
int check_pos2 (int ***A,
                int i,
                int j,
                int k,
                int m,
                int n,
                int o);
int check_pos3 (int ***A,
                int i,
                int j,
                int k,
                int m,
                int n,
                int o);
int check_pos4 (int ***A,
                int i,
                int j,
                int k,
                int m,
                int n,
                int o);
int contact_protein_kv (int ***A,
                        int i,
                        int j,
                        int k,
                        int m,
                        int n,
                        int o);
int DFS_search (int ***A,
                int m,
                int n,
                int o,
                double h,
                double filter);
char* combine (const char *s1,
               const char *s2);
void Matrix_fill (int ***A,
                   int m,
                   int n,
                   int o,
                   double h,
                   double probe,
                   double X1,
                   double Y1,
                   double Z1);
void Matrix_export (int ***A,
                    int ***S,
                    int kvp_mode,
                    int m,
                    int n,
                    int o,
                    double h,
                    int ncav,
                    char *output_base_name,
                    char *output_pdb,
                    double X1,
                    double Y1,
                    double Z1);
void Matrix_surf (int ***A,
                  int m,
                  int n,
                  int o,
                  double h,
                  double radius);
void Matrix_surface (int ***A,
                     int ***S,
                     int m,
                     int n,
                     int o,
                     double h,
                     double X1,
                     double Y1,
                     double Z1);
void Matrix_search (int ***A,
                    int ***S,
                    int m,
                    int n,
                    int o,
                    double h,
                    double probe,
                    double X1,
                    double Y1,
                    double Z1,
                    int ncav);
void Matrix_adjust (int ***A,
                    int m,
                    int n,
                    int o,
                    double h,
                    double limit,
                    double X1,
                    double Y1,
                    double Z1);
void Matrix_filter (int ***A,
                    int ***S,
                    int m,
                    int n,
                    int o,
                    double h,
                    double bX1,
                    double bY1,
                    double bZ1,
                    double bX2,
                    double bY2,
                    double bZ2,
                    double norm1);
void DFS (int ***A,
          int a,
          int b,
          int c,
          int m,
          int n,
          int o,
          int tag,
          double h);
void Matrix_subtract (int ***S,
                      int ***A,
                      int m,
                      int n,
                      int o,
                      double h,
                      double removal_distance);
void remove_cavity (int ***A,
                    int m,
                    int n,
                    int o,
                    int tag);
void filter_outliers (int ***A,
                      int m,
                      int n,
                      int o);
void Area_search (int ***S,
                 int m,
                 int n,
                 int o,
                 double h,
                 int tag);
void check_faces (int ***S,
                  int i,
                  int j,
                  int k,
                  double h,
                  double *area);
void free_matrix (int ***A,
                  int m,
                  int n,
                  int o);
void free_node ();

#endif
