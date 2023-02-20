#ifndef UTILS_H
#define UTILS_H

char *combine(const char *s1, const char *s2);
double max(double a, double b);
double min(double a, double b);
int _read_line(FILE *arq, char LINE[], int size);
void _extract(char FROM[], int nF, char TO[], int nT, int start, int end);
char _extract_chain(char FROM[], int nF);

#endif
