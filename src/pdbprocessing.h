/*    This file is part of KVFinder.

    KVFinder is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    KVFinder is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with KVFinder.  If not, see <http://www.gnu.org/licenses/>.

    The KVFinder software was developed by:
    Joao Victor da Silva Guerra
    Saulo Henrique Pires de Oliveira
    Felipe Augusto Nunes Ferraz
    Rodrigo Vargas Honorato
    Jose Xavier Neto
    Tiago Jose Paschoal Sobreira
    Paulo Sergio Lopes de Oliveira

    National Center of Energy and Material Research - CNPEM
    National Laboratory of Biosciences - LNBio
    Campinas, Brazil - P.O. Box 6192 - CEP 13083-970, Campinas - SP

    Contact: paulo.oliveira@lnbio.cnpem.br
    KVFinder website:
   http://lnbio.cnpem.br/bioinformatics/main/software/KVFinder*/

#ifndef PDBPROCESSING_H
#define PDBPROCESSING_H

/* Maximum length of PDB filename */
#define NAME_MAX 500

/* Define type atom for a ATOM structure */
typedef struct ATOM {
  double x;
  double y;
  double z;
  double radius;
  int resnum;
  char chain;
  char res_name;
  struct ATOM *next;
} atom;

/* Declare variables */
atom *v;
double sina, sinb, cosa, cosb;

int PDB_load(dict *DIC[TABLE_SIZE], int tablesize,
             char TABLE[TABLE_SIZE][RES_SIZE], char PDB_NAME[NAME_MAX],
             double probe, int m, int n, int o, double h, double X1, double Y1,
             double Z1, FILE **log_file);
int PDB_load2(char PDB_NAME[NAME_MAX]);
int PDB_load3(char PDB_NAME[NAME_MAX]);
char extractChain(char LINE[100]);
void convert(char S[50], double *coord);
void extract(char LINE[100], int a, int b, char S[]);
void insert_atom(double x, double y, double z, double radius, int resnumber,
                 char chain, char res_name);
void free_atom();

#endif