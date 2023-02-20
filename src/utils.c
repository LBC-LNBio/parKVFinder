#include <stdio.h>  //FILE, EOF
#include <stdlib.h> //malloc, calloc
#include <string.h> //memcpy, size_t, strlen

/*
 * Function: combine
 * -----------------
 *
 * Combine two strings into one
 *
 * s1: string 1
 * s2: string 2
 *
 * returns: combined string
 *
 */
char *combine(const char *s1, const char *s2)
{
    /* Declare variables */
    const size_t len1 = strlen(s1);
    const size_t len2 = strlen(s2);
    char *result = malloc(len1 + len2 + 1); /* +1 for the zero-terminator */

    memcpy(result, s1, len1);
    memcpy(result + len1, s2, len2 + 1); /* +1 to copy the null-terminator */

    return result;
}

/*
 * Function: max
 * -------------
 *
 * Get the maximum value of two doubles
 *
 * a: double
 * b: double
 *
 * returns: maximum between a and b
 *
 */
double max(double a, double b)
{

    if (a > b)
        return a;
    else
        return b;
}

/*
 * Function: min
 * -------------
 *
 * Get the minimum value of two doubles
 *
 * a: double
 * b: double
 *
 * returns: minimum between a and b
 *
 */
double min(double a, double b)
{

    if (a < b)
        return a;
    else
        return b;
}

/*
 * Function: _read_line
 * --------------------
 *
 * Read a line from a file to a string (vector).
 *
 * file: file opened by fopen
 * LINE: flexible array
 * size: size of flexible array
 *
 * returns: file is over (1) or not over (0)
 *
 */
int _read_line(FILE *file, char LINE[], int size)
{
    /* Declare variables */
    int i;

    /* Read PDB file inside LINE[100] */
    for (i = 0, LINE[i] = getc(file); LINE[i] != EOF && LINE[i] != '\n' && i < size;
         i++, LINE[i] = getc(file))
        ;

    /* If LINE[100] is not \n, read PDB file until LINE[100] assume \n or EOF
     * value */
    if (i == size && LINE[i] != '\n')
        for (; LINE[i] != '\n' && LINE[i] != EOF; LINE[i] = getc(file))
            ;

    /* If PDB file is over, return False */
    if (LINE[i] == EOF)
        return 0;
    /* If PDB file is not over, return True */
    else
        return 1;
}

/*
 * Function: _extract
 * ------------------
 *
 * Extract from FROM[start] until FROM[end] to a string TO
 *
 * FROM: flexible array to extract a snippet
 * nF: size of FROM
 * TO: flexible array to move the snippet
 * nT: size of TO
 * start: start of snippet in FROM
 * end: end of snippet in TO
 *
 */
void _extract(char FROM[], int nF, char TO[], int nT, int start, int end)
{
    int i;

    if (end > start || nF > end || (end - start) <= nT)
        for (i = start; i < end; i++)
            TO[i - start] = FROM[i];

    if (i < nT)
        TO[i] = '\0';
}

/*
 * Function: _remove_char
 * ----------------------
 *
 * Remove char from string
 *
 * FROM: flexible array with a string
 * c: char to remove from string
 *
 */
void _remove_char(char FROM[], int nF, char c)
{
    int i, j;

    for (i = 0; FROM[i] != '\0'; i++) {
        if (FROM[i] == c)
            for (j = i; FROM[j] != '\0'; j++)
                FROM[j] = FROM[j + 1];
  }
}
