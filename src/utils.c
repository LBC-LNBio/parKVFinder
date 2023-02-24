#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Functions */

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
double max(double a, double b) {

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
double min(double a, double b) {

  if (a < b)
    return a;
  else
    return b;
}

/*
 * Function: _resolution2step
 * --------------------------
 *
 * Set grid spacing base on resolution mode
 *
 * flag: resolution mode (Off, Low, Medium, High)
 *
 * returns: grid spacing (step size)
 *
 */
double _resolution2step(char flag[]) {

  switch (flag[0]) {

  case 'L':
    return 0.6;
    break;

  case 'M':
    return 0.5;
    break;

  case 'H':
    return 0.25;
    break;

  default:
    break;
  }
}

/*
 * Function: _combine
 * ------------------
 *
 * Combine two strings into one
 *
 * s1: string 1
 * s2: string 2
 *
 * returns: combined string
 *
 */
char *_combine(const char *s1, const char *s2) {
  /* Declare variables */
  const size_t len1 = strlen(s1);
  const size_t len2 = strlen(s2);
  char *result = malloc(len1 + len2 + 1); /* +1 for the zero-terminator */

  memcpy(result, s1, len1);
  memcpy(result + len1, s2, len2 + 1); /* +1 to copy the null-terminator */

  return result;
}

/*
 * Function: _convert_residue_code
 * -------------------------------
 *
 * Covert 3-letter code to 1-letter code
 *
 * RESIDUE: residue in 3-letter code
 *
 * returns: residue in 1-letter code
 *
 */
char _convert_residue_code(char RESIDUE[]) {

  if (!strcmp(RESIDUE, "ALA"))
    return 'A';
  if (!strcmp(RESIDUE, "ARG"))
    return 'R';
  if (!strcmp(RESIDUE, "ASN"))
    return 'N';
  if (!strcmp(RESIDUE, "ASP"))
    return 'D';
  if (!strcmp(RESIDUE, "CYS"))
    return 'C';
  if (!strcmp(RESIDUE, "GLN"))
    return 'Q';
  if (!strcmp(RESIDUE, "GLU"))
    return 'E';
  if (!strcmp(RESIDUE, "GLY"))
    return 'G';
  if (!strcmp(RESIDUE, "HIS"))
    return 'H';
  if (!strcmp(RESIDUE, "ILE"))
    return 'I';
  if (!strcmp(RESIDUE, "LEU"))
    return 'L';
  if (!strcmp(RESIDUE, "LYS"))
    return 'K';
  if (!strcmp(RESIDUE, "MET"))
    return 'M';
  if (!strcmp(RESIDUE, "PHE"))
    return 'F';
  if (!strcmp(RESIDUE, "PRO"))
    return 'P';
  if (!strcmp(RESIDUE, "SER"))
    return 'S';
  if (!strcmp(RESIDUE, "THR"))
    return 'T';
  if (!strcmp(RESIDUE, "TRP"))
    return 'W';
  if (!strcmp(RESIDUE, "TYR"))
    return 'Y';
  if (!strcmp(RESIDUE, "VAL"))
    return 'V';

  return 'X';
}

/*
 * Function: _get_file_extension
 * -----------------------------
 *
 * Get file extension
 *
 * fn: filename
 *
 * returns: file extension
 *
 */
char *_get_file_extension(char *fn) {
  char *dot = strrchr(fn, '.');

  if (!dot || dot == fn)
    return "";
  return dot + 1;
}

/*
 * Function: _int2toml
 * -------------------
 *
 * Convert 1 or 0 to 'true' or 'false'
 *
 * boolean: 1 or 0
 * flag: TOML boolean ('true' or 'false')
 *
 */
int _int2toml(int boolean, char **flag) {

  if (boolean)
    *flag = "true";
  else
    *flag = "false";
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
int _read_line(FILE *file, char LINE[], int size) {
  /* Declare variables */
  int i;

  /* Read PDB file inside LINE[100] */
  for (i = 0, LINE[i] = getc(file);
       LINE[i] != EOF && LINE[i] != '\n' && i < size; i++, LINE[i] = getc(file))
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
 * Function: _toml2int
 * -------------------
 *
 * Convert 'true' or 'false' to 1 or 0
 *
 * flag: TOML boolean
 *
 */
int _toml2int(char flag[6]) {

  if (strcmp(flag, "false") == 0 || strcmp(flag, "0") == 0)
    return 0;
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
void _extract(char FROM[], int nF, char TO[], int nT, int start, int end) {
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
void _remove_char(char FROM[], int nF, char c) {
  int i, j;

  for (i = 0; FROM[i] != '\0'; i++) {
    if (FROM[i] == c)
      for (j = i; FROM[j] != '\0'; j++)
        FROM[j] = FROM[j + 1];
  }
}

/*
 * Function: _initialize_string
 * ----------------------------
 *
 * Fill string FROM[] with '\0'
 *
 * FROM: flexible array with a string
 * nF: size of FROM
 *
 */
void _initialize_string(char FROM[], int nF) {
  int i;

  for (i = 0; i < nF; i++)
    FROM[i] = '\0';
}
