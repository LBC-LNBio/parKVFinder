#ifndef UTILS_H
#define UTILS_H

char *combine(const char *s1, const char *s2);
double max(double a, double b);
double min(double a, double b);
int _read_line(FILE *arq, char LINE[], int size);
void _extract(char FROM[], int nF, char TO[], int nT, int start, int end);
void _remove_char(char FROM[], int nF, char c);
char _convert_residue_code(char RESIDUE[]);
void _initialize_string(char FROM[], int nF);
int _toml2int(char flag[6]);
int _int2toml(int boolean, char **flag);
char *_get_file_extension(char *fn);

#endif
