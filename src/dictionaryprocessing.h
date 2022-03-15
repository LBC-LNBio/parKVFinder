#ifndef DICTIONARYPROCESSING_H
#define DICTIONARYPROCESSING_H

#define DIC_NAME_MAX 500
#define TABLE_SIZE 500
#define RES_SIZE 4
#define ATOM_SIZE 4

/* Define type dict for a DICTIONARY structure */
typedef struct DICTIONARY {
  double radius;
  char symbol[6];
  struct DICTIONARY *next;
} dict;

/* Define custom functions */
void trim(char S[50], char c);
void resolution_input(char flag_in[7], double *Vvoxel, int *resolution_mode,
                      double *h);
char convertRES(char RES[RES_SIZE]);
int dictionary_load(dict *DIC[TABLE_SIZE], int tablesize,
                    char dictionary_name[DIC_NAME_MAX]);
int define_table(char TABLE[TABLE_SIZE][RES_SIZE],
                 char dictionary_name[DIC_NAME_MAX]);
double dictionary_consult_radius(char RES[RES_SIZE], char ATOM_TYPE[ATOM_SIZE],
                                 dict *DIC[TABLE_SIZE], int tablesize,
                                 char TABLE[TABLE_SIZE][RES_SIZE],
                                 char ATOM_SYMBOL[2], FILE **log_file);

#endif
