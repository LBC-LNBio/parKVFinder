#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "utils.h"

#include "fileprocessing.h"

/* parKVFinder parameter file processing */

/*Reads TOML file inside LINE*/
int get_toml_line(FILE *arq, char LINE[500]) {
  /* Declare variables */
  int i;

  /* Read PDB file inside LINE[500] */
  for (i = 0, LINE[i] = getc(arq); LINE[i] != EOF && LINE[i] != '\n' && i < 500;
       i++, LINE[i] = getc(arq))
    ;

  /* If LINE[500] is not \n, read PDB file until LINE[100] assume \n or EOF
   * value */
  if (i == 500 && LINE[i] != '\n')
    for (; LINE[i] != '\n' && LINE[i] != EOF; LINE[i] = getc(arq))
      ;

  /* If PDB file is over, return False */
  if (LINE[i] == EOF)
    return 0;
  /* If PDB file is not over, return True */
  else
    return 1;
}

/* Extract from LINE[a] until LINE[b] to a char array S[100] of size b-a */
void extract_toml_line(char LINE[500], int a, int b, char S[500]) {
  /* Declare variables */
  int i;

  /* Tests inputs */
  if (b > 500 || a > b)
    return;

  /* Extract process */
  for (i = a; i < b; i++)
    S[i - a] = LINE[i];

  /* Mark last position in S[] */
  if (i < 500)
    S[i] = '\0';
}

/* Read TOML file and pass to struct TOML */
parameters *readTOML(char *path) {
  /* Declare variables */
  FILE *parameters_file;
  char LINE[500], config[500] = "", value[500] = "", key[500] = "";
  char *vb, *ib, *p1, *p2, *p3, *p4;
  int i = 0, equal = 0, flag_visiblebox = 0, flag_internalbox = 0, count = 0;
  int config_bar[2];

  /* Allocate memory to struct TOML */
  parameters *p = (parameters *)malloc(sizeof(parameters));

  /* Open parameters TOML file */
  parameters_file = fopen(path, "r");

  /* TOML file not found */
  if (parameters_file == NULL) {

    /* Print error and exit */
    printf("\033[0;31mError:\033[0m Invalid parameters file! Please select a "
           "valid parameters and try again.\n");
    exit(-1);

  } else {
    /* While TOML file is not over, do... */
    while (get_toml_line(parameters_file, LINE)) {

    HANDLE_LAST_LINE:
      /* Remove tabs */
      _remove_char(LINE, strlen(LINE), '\t');
      /* Remove whitespaces */
      _remove_char(LINE, strlen(LINE), ' ');

      /* If a '#' is found, ignore line */
      if (LINE[0] == '#' || LINE[0] == '\n' || LINE[0] == ' ')
        ;
      else {

        for (i = 0; i < strlen(LINE); i++) {

          /* Get position of newline symbol in LINE */
          if (LINE[i] == '=')
            equal = i;

          /* Get position of '[' symbol in LINE */
          if (LINE[i] == '[')
            config_bar[0] = i;

          /* Get position of ']' symbol in LINE */
          if (LINE[i] == ']')
            config_bar[1] = i;
        }

        /* If vector config_bar are filled, do ... */
        if (config_bar[1] != 0) {

          /* Restart string config[100] */
          _initialize_string(config, strlen(config));
          /* Extract name of configuration */
          extract_toml_line(LINE, config_bar[0] + 1, config_bar[1], config);

          vb = strstr(config, ".visiblebox.");
          ib = strstr(config, ".internalbox.");

          /* If configuration is visiblebox, do ... */
          if (vb != NULL) {
            flag_visiblebox = 1;
            count = 0;
          }

          /* If configuration is internalbox, do ... */
          if (ib != NULL) {
            flag_internalbox = 1;
            count = 0;
          }
        }

        /* If equal symbol is found, do ... */
        if (equal != 0) {

          /* Extract key */
          extract_toml_line(LINE, 0, equal, key);

          /* Extract value */
          extract_toml_line(LINE, equal + 1, strlen(LINE) - 1, value);

          /* If bracket is found, do ... */
          if (flag_visiblebox || flag_internalbox) {

            p1 = strstr(config, ".p1");
            p2 = strstr(config, ".p2");
            p3 = strstr(config, ".p3");
            p4 = strstr(config, ".p4");

            /* Visible box */
            if (flag_visiblebox) {

              /* Count number of coordinates read */
              count++;

              /* Visible point 1 : bP1 */
              if (p1 != NULL) {

                /* Get x, y or z coordinates */
                if (strcmp(key, "x") == 0)
                  (p->bX1) = atof(value);
                if (strcmp(key, "y") == 0)
                  (p->bY1) = atof(value);
                if (strcmp(key, "z") == 0)
                  (p->bZ1) = atof(value);
              }

              /* Visible point 2 : bP2 */
              if (p2 != NULL) {

                if (strcmp(key, "x") == 0)
                  (p->bX2) = atof(value);
                if (strcmp(key, "y") == 0)
                  (p->bY2) = atof(value);
                if (strcmp(key, "z") == 0)
                  (p->bZ2) = atof(value);
              }

              /* Visible point 3 : bP3 */
              if (p3 != NULL) {

                if (strcmp(key, "x") == 0)
                  (p->bX3) = atof(value);
                if (strcmp(key, "y") == 0)
                  (p->bY3) = atof(value);
                if (strcmp(key, "z") == 0)
                  (p->bZ3) = atof(value);
              }

              /* Visible point 4 : bP4 */
              if (p4 != NULL) {

                if (strcmp(key, "x") == 0)
                  (p->bX4) = atof(value);
                if (strcmp(key, "y") == 0)
                  (p->bY4) = atof(value);
                if (strcmp(key, "z") == 0)
                  (p->bZ4) = atof(value);
              }

              if (count == 3)
                flag_visiblebox = 0;
            }

            /* Visible box */
            if (flag_internalbox) {

              count++;

              /* Internal point 1 : P1 */
              if (p1 != NULL) {

                /* Get x, y or z coordinates */
                if (strcmp(key, "x") == 0)
                  (p->X1) = atof(value);
                if (strcmp(key, "y") == 0)
                  (p->Y1) = atof(value);
                if (strcmp(key, "z") == 0)
                  (p->Z1) = atof(value);
              }

              /* Internal point 2 : P2 */
              if (p2 != NULL) {

                if (strcmp(key, "x") == 0)
                  (p->X2) = atof(value);
                if (strcmp(key, "y") == 0)
                  (p->Y2) = atof(value);
                if (strcmp(key, "z") == 0)
                  (p->Z2) = atof(value);
              }

              /* Internal point 3 : P3 */
              if (p3 != NULL) {

                if (strcmp(key, "x") == 0)
                  (p->X3) = atof(value);
                if (strcmp(key, "y") == 0)
                  (p->Y3) = atof(value);
                if (strcmp(key, "z") == 0)
                  (p->Z3) = atof(value);
              }

              /* Internal point 4 : P4 */
              if (p4 != NULL) {

                if (strcmp(key, "x") == 0)
                  (p->X4) = atof(value);
                if (strcmp(key, "y") == 0)
                  (p->Y4) = atof(value);
                if (strcmp(key, "z") == 0)
                  (p->Z4) = atof(value);
              }

              if (count == 3)
                flag_internalbox = 0;
            }

          } else {

            /* Remove quotation marks */
            _remove_char(value, strlen(value), '\"');
            _remove_char(value, strlen(value), '\"');

            /* Save values inside struct TOML */
            if (strcmp(key, "output") == 0)
              strcpy((p->OUTPUT), value);
            if (strcmp(key, "base_name") == 0)
              strcpy((p->BASE_NAME), value);
            if (strcmp(key, "pdb") == 0)
              strcpy((p->PDB_NAME), value);
            if (strcmp(key, "ligand") == 0)
              strcpy((p->LIGAND_NAME), value);
            if (strcmp(key, "dictionary") == 0)
              strcpy((p->dictionary_name), value);
            if (strcmp(key, "whole_protein_mode") == 0)
              (p->whole_protein_mode) = _toml2int(value);
            if (strcmp(key, "resolution_mode") == 0)
              strcpy((p->resolution_flag), value);
            if (strcmp(key, "box_mode") == 0)
              (p->box_mode) = _toml2int(value);
            if (strcmp(key, "surface_mode") == 0)
              (p->surface_mode) = _toml2int(value);
            if (strcmp(key, "kvp_mode") == 0)
              (p->kvp_mode) = _toml2int(value);
            if (strcmp(key, "ligand_mode") == 0)
              (p->ligand_mode) = _toml2int(value);
            if (strcmp(key, "step_size") == 0)
              (p->h) = atof(value);
            if (strcmp(key, "probe_in") == 0)
              (p->probe_in) = atof(value);
            if (strcmp(key, "probe_out") == 0)
              (p->probe_out) = atof(value);
            if (strcmp(key, "volume_cutoff") == 0)
              (p->volume_cutoff) = atof(value);
            if (strcmp(key, "ligand_cutoff") == 0)
              (p->ligand_cutoff) = atof(value);
            if (strcmp(key, "removal_distance") == 0)
              (p->removal_distance) = atof(value);
          }
        }
      }

      /* Restart variables */
      _initialize_string(LINE, strlen(LINE));
      _initialize_string(key, strlen(key));
      _initialize_string(value, strlen(value));
      equal = 0;
      config_bar[0] = 0;
      config_bar[1] = 0;
    }
  }

  /* Check if last line read found EOF in the middle of the line */
  for (i = 0; LINE[i] != EOF; i++)
    ;
  if (i > 0) {
    LINE[i] = '\n';
    goto HANDLE_LAST_LINE;
  }

  /* Return struct TOML */
  return p;
}

/*
 * Function: write_parameters
 * --------------------------
 *
 * Write TOML-formatted parameters file of parKVFinder.
 *
 * toml_name: a path to a TOML-formatted file.
 * OUTPUT: output directory.
 * BASE_NAME: base name.
 * dictionary_name: path of van der Waals dictionary file
 * (KVFinder_PATH/dictionary).
 * PDB_NAME: a path to a target PDB file.
 * LIGAND_NAME: a path to a target ligand PDB file.
 * whole_protein_mode: Whether to use the whole protein to build the 3D grid.
 * resolution_mode: A flag to a pre-defined grid spacing: Low (0.6), Medium
 * (0.5), High (0.25).
 * box_mode: Whether a custom grid is applied.
 * surface_mode: Surface representation. SES (0) or SAS (1).
 * kvp_mode: Filled (1) or filtered (0) cavities.
 * h: Grid spacing (A) (step).
 * probe_in: Probe In size (A).
 * probe_out: Probe Out size (A).
 * volume_cutoff: Volume filter for detected cavities (A3).
 * ligand_cutoff: A radius to limit a space around a ligand (A).
 * removal_distance: A length to be removed from the cavity-bulk frontier (A).
 * Xi, Yi, Zi: Internal box coordinates.
 * bXi, bYi, bZi: Visible box coordinates.
 *
 */
void write_parameters(char *toml_name, char OUTPUT[500], char BASE_NAME[500],
                      char dictionary_name[500], char PDB_NAME[500],
                      char LIGAND_NAME[500], int whole_protein_mode,
                      char resolution_mode[7], int box_mode, int surface_mode,
                      int kvp_mode, int ligand_mode, double h, double probe_in,
                      double probe_out, double volume_cutoff,
                      double ligand_cutoff, double removal_distance, double X1,
                      double Y1, double Z1, double X2, double Y2, double Z2,
                      double X3, double Y3, double Z3, double X4, double Y4,
                      double Z4, double bX1, double bY1, double bZ1, double bX2,
                      double bY2, double bZ2, double bX3, double bY3,
                      double bZ3, double bX4, double bY4, double bZ4) {
  /* Declare variables */
  FILE *toml_file;
  char buffer[1024], *wpmode, *bmode, *smode, *kmode, *lmode, *hmode, *emode;

  /* Create KV_Files directory */
  mkdir(_combine(OUTPUT, "KV_Files/"), S_IRWXU);

  /* Open TOML file */
  toml_file = fopen(toml_name, "w");
  /* Create a buffer */
  memset(buffer, '\0', sizeof(buffer));
  /* Define buffer as writing buffer of size 1024 */
  setvbuf(toml_file, buffer, _IOFBF, 1024);

  /* Print TOML file */
  fprintf(toml_file, "# TOML configuration file for parKVFinder software.\n");
  fprintf(toml_file, "\ntitle = \"parKVFinder parameters file\"\n");

  fprintf(toml_file, "\n[FILES_PATH]\n");
  fprintf(toml_file,
          "# The path of van der Waals radii dictionary for parKVFinder.\n");
  fprintf(toml_file, "dictionary = \"%s\"\n", dictionary_name);
  fprintf(toml_file, "# The path of the input PDB file.\n");
  fprintf(toml_file, "pdb = \"%s\"\n", PDB_NAME);
  fprintf(toml_file, "# The path of the output directory.\n");
  fprintf(toml_file, "output = \"%s\"\n", OUTPUT);
  fprintf(toml_file, "# Base name for output files.\n");
  fprintf(toml_file, "base_name = \"%s\"\n", BASE_NAME);
  fprintf(toml_file, "# Path for the ligand's PDB file.\n");
  fprintf(toml_file, "ligand = \"%s\"\n", LIGAND_NAME);

  fprintf(toml_file, "\n[SETTINGS]\n");
  fprintf(toml_file, "# Settings for parKVFinder software\n");
  fprintf(toml_file, "\n\t[SETTINGS.modes]\n");
  fprintf(toml_file, "\t# Whole Protein mode defines the search space as the "
                     "whole protein.\n");
  _int2toml(whole_protein_mode, &wpmode);
  fprintf(toml_file, "\twhole_protein_mode = %s\n", wpmode);
  fprintf(toml_file, "\t# Box Adjustment mode defines the search space as a "
                     "box that includes a specific region.\n");
  _int2toml(box_mode, &bmode);
  fprintf(toml_file, "\tbox_mode = %s\n", bmode);
  fprintf(toml_file, "\t# Resolution mode implicitly sets the step size (grid "
                     "spacing) of the 3D grid.\n");
  fprintf(toml_file,
          "\t# If set to High, sets a voxel volume of 0.2. If set to Medium, "
          "sets a voxel volume of 0.1. If set to Low, sets a voxel volume of "
          "0.01. If set to Off, the step size must be set explicitly.\n");
  fprintf(toml_file, "\tresolution_mode = \"%s\"\n", resolution_mode);
  fprintf(toml_file, "\t# Surface mode defines the type of surface "
                     "representation to be applied, van der Waals molecular "
                     "surface (true) or solvent accessible surface (false).\n");
  _int2toml(surface_mode, &smode);
  fprintf(toml_file, "\tsurface_mode = %s\n", smode);
  fprintf(toml_file, "\t# Cavity output mode defines whether cavities are "
                     "exported to the output PDB file as filled cavities "
                     "(true) or filtered cavities (false).\n");
  _int2toml(kvp_mode, &kmode);
  fprintf(toml_file, "\tkvp_mode = %s\n", kmode);
  fprintf(toml_file, "\t# Ligand adjustment mode defines the search space "
                     "around the ligand.\n");
  _int2toml(ligand_mode, &lmode);
  fprintf(toml_file, "\tligand_mode = %s\n", lmode);

  fprintf(toml_file, "\n\t[SETTINGS.step_size]\n");
  fprintf(toml_file, "\t# Sets the 3D grid spacing. It directly affects "
                     "accuracy and runtime.\n");
  fprintf(toml_file, "\tstep_size = %.2lf\n", h);

  fprintf(toml_file, "\n\t[SETTINGS.probes]\n");
  fprintf(toml_file, "\t# parKVFinder works with a two sized probe system. A "
                     "smaller probe, called Probe In, and a bigger one, called "
                     "Probe Out, rolls around the protein.\n");
  fprintf(toml_file, "\t# Points reached by the Probe In, but not the Probe "
                     "Out are considered cavity points.\n");
  fprintf(toml_file, "\t# Sets Probe In diameter. Default: 1.4 angstroms.\n");
  fprintf(toml_file, "\tprobe_in = %.2lf\n", probe_in);
  fprintf(toml_file, "\t# Sets Probe Out diameter. Default: 4.0 angstroms.\n");
  fprintf(toml_file, "\tprobe_out = %.2lf\n", probe_out);

  fprintf(toml_file, "\n\t[SETTINGS.cutoffs]\n");
  fprintf(toml_file, "\t# Sets a volume cutoff for the detected cavities. "
                     "Default: 5.0 angstroms.\n");
  fprintf(toml_file, "\tvolume_cutoff = %.2lf\n", volume_cutoff);
  fprintf(toml_file,
          "\t# Sets a distance cutoff for a search space around the ligand in "
          "ligand adjustment mode. Default: 5.0 angstroms.\n");
  fprintf(toml_file, "\tligand_cutoff = %.2lf\n", ligand_cutoff);
  fprintf(toml_file, "\t# Sets a removal distance for the cavity frontier, "
                     "which is defined by comparing Probe In and Probe Out "
                     "surfaces. Default: 2.4 angstroms.\n");
  fprintf(toml_file, "\tremoval_distance = %.2lf\n", removal_distance);

  fprintf(toml_file, "\n\t[SETTINGS.visiblebox]\n");
  fprintf(toml_file,
          "\t# Coordinates of the vertices that define the visible 3D grid. "
          "Only four points are required to define the search space.\n");
  fprintf(toml_file, "\n\t[SETTINGS.visiblebox.p1]\n");
  fprintf(toml_file, "\tx = %.2f\n", bX1);
  fprintf(toml_file, "\ty = %.2f\n", bY1);
  fprintf(toml_file, "\tz = %.2f\n", bZ1);
  fprintf(toml_file, "\n\t[SETTINGS.visiblebox.p2]\n");
  fprintf(toml_file, "\tx = %.2f\n", bX2);
  fprintf(toml_file, "\ty = %.2f\n", bY2);
  fprintf(toml_file, "\tz = %.2f\n", bZ2);
  fprintf(toml_file, "\n\t[SETTINGS.visiblebox.p3]\n");
  fprintf(toml_file, "\tx = %.2f\n", bX3);
  fprintf(toml_file, "\ty = %.2f\n", bY3);
  fprintf(toml_file, "\tz = %.2f\n", bZ3);
  fprintf(toml_file, "\n\t[SETTINGS.visiblebox.p4]\n");
  fprintf(toml_file, "\tx = %.2f\n", bX4);
  fprintf(toml_file, "\ty = %.2f\n", bY4);
  fprintf(toml_file, "\tz = %.2f\n", bZ4);

  fprintf(toml_file, "\n\t[SETTINGS.internalbox]\n");
  fprintf(toml_file,
          "\t# Coordinates of the internal 3D grid. Used for calculations.\n");
  fprintf(toml_file, "\n\t[SETTINGS.internalbox.p1]\n");
  fprintf(toml_file, "\tx = %.2f\n", X1);
  fprintf(toml_file, "\ty = %.2f\n", Y1);
  fprintf(toml_file, "\tz = %.2f\n", Z1);
  fprintf(toml_file, "\n\t[SETTINGS.internalbox.p2]\n");
  fprintf(toml_file, "\tx = %.2f\n", X2);
  fprintf(toml_file, "\ty = %.2f\n", Y2);
  fprintf(toml_file, "\tz = %.2f\n", Z2);
  fprintf(toml_file, "\n\t[SETTINGS.internalbox.p3]\n");
  fprintf(toml_file, "\tx = %.2f\n", X3);
  fprintf(toml_file, "\ty = %.2f\n", Y3);
  fprintf(toml_file, "\tz = %.2f\n", Z3);
  fprintf(toml_file, "\n\t[SETTINGS.internalbox.p4]\n");
  fprintf(toml_file, "\tx = %.2f\n", X4);
  fprintf(toml_file, "\ty = %.2f\n", Y4);
  fprintf(toml_file, "\tz = %.2f\n", Z4);

  fflush(toml_file);
  fclose(toml_file);
}

/* van der Waals file processing */

/*
 * Function: _get_vdw_radius
 * -------------------------
 *
 * Get van der Waals radius for an atom of a residue
 *
 * RESIDUE: residue in 3-letter code
 * ATOM_TYPE: atom name
 * DIC: a vector containing atom name and atom radius per residue type
 * tablesize: number of residue types
 * TABLE: a vector containing residues information (index and residue name)
 * ATOM_SYMBOL: atom symbol
 * log_file: path to log file (KVFinder.log)
 *
 * returns: success (1) or fail (0) to read van der Waals radii file
 *
 */
double _get_vdw_radius(char RESIDUE[4], char ATOM_TYPE[4], vdw *DIC[500],
                       int tablesize, char TABLE[500][4], char ATOM_SYMBOL[2],
                       FILE **log_file) {
  /* Declare variables */
  int i;
  double value = 0.0;
  char ATOM_TYPE_AUX[4], RES_AUX[4];
  vdw *p;

  for (i = 0; i < tablesize && value == 0.0; i++) {

    /* If TABLE[i] residue is equal to RES residue, do ... */
    if (!strcmp(RESIDUE, TABLE[i])) {

      /* Loop inside residue linked list and look for ATOM name, then retrieve
       * radius to value */
      for (p = DIC[i]->next; p != NULL && value == 0.0; p = p->next)
        if (!strcmp(ATOM_TYPE, p->symbol)) {

          /* Retrieve radius of an atom of a residue */
          value = p->radius;
          /* Return radius value */
          return value;
        }
    }
  }

  /* If radius not found inside linked list, do ... */
  if (!value) {

    fprintf(*log_file, "Warning: Atom radius not found in dictionary:%s %s!\n",
            ATOM_TYPE, RESIDUE);
    strcpy(RES_AUX, "GEN");

    if (ATOM_SYMBOL[0] > 90 || ATOM_SYMBOL[0] < 65) {

      ATOM_TYPE_AUX[0] = ATOM_SYMBOL[1];
      ATOM_TYPE_AUX[1] = '\0';

    } else {

      ATOM_TYPE_AUX[0] = ATOM_SYMBOL[0];
      ATOM_TYPE_AUX[1] = ATOM_SYMBOL[1];
      ATOM_TYPE_AUX[2] = '\0';
    }

    /* Loop inside generic residue linked list and look for an atom, then
     * retrieve its radius to value */
    for (i = 0; i < tablesize && value == 0.0; i++) {

      if (!strcmp(RES_AUX, TABLE[i])) {

        for (p = DIC[i]->next; p != NULL && value == 0.0; p = p->next) {

          if (!strcmp(ATOM_TYPE_AUX, p->symbol)) {

            /* If radius value is found, print generic atom and its radius value
             */
            value = p->radius;
            /* Print warning in log file */
            fprintf(*log_file,
                    "Warning: Using generic atom %s radius value %.2lf\n",
                    ATOM_TYPE_AUX, p->radius);
            /* Return radius value */
            return value;
          }
        }
      }
    }

    /* Print warning in log file */
    fprintf(*log_file,
            "Warning: Radius data not found for atom %s. This atom will be "
            "excluded from analysis.\n",
            ATOM_TYPE_AUX);
  }

  /* Atom excluded from analysis and return radius 0.0 */
  return value;
}

/*
 * Function: _get_residues_information
 * ----------------------------------
 *
 * Read a van der Waals radii dictionary file ($KVFinder_PATH/dictionary) and
 * saves residues information inside TABLE vector and return number of residues.
 *
 * TABLE: a vector containing residues information (index and residue name)
 * tablesize: number of residue types
 * dictionary_name: path of van der Waals dictionary file
 * ($KVFinder_PATH/dictionary)
 *
 * returns: number of residues
 *
 */
int _get_residues_information(char dictionary_name[500], char TABLE[500][4]) {
  /* Declare variables */
  int i = 0, j = 0;
  char AUX[50];
  FILE *arq;

  /* Open dictionary file */
  arq = fopen(dictionary_name, "r");

  /* File has not been found */
  if (arq == NULL) {

    /* Print error and exit */
    fprintf(stderr,
            "\033[0;31mError:\033[0m Residues dictionary file not found!\n");
    exit(-1);

  }
  /* File has been found */
  else
    /* While EOF not read in AUX object, do ... */
    while (fscanf(arq, "%s", AUX) != EOF) {
      if (AUX[0] == '>') {
        /* Copy each letter of residue symbol inside TABLE */
        for (j = 1; j < strlen(AUX); j++)
          TABLE[i][j - 1] = AUX[j];
        /* Put each symbol at the end of residue symbol */
        TABLE[i][j - 1] = '\0';
        /* Increment line iterator */
        i++;
      }
    }

  /* Close dictionary file */
  fclose(arq);

  /* Return number of residues in TABLE */
  return i;
}

/*
 * Function: read_vdw
 * ------------------
 *
 * Read a van der Waals radii dictionary file ($KVFinder_PATH/dictionary) to a
 * linked list.
 *
 * DIC: a vector containing atom name and atom radius per residue type
 * tablesize: number of residue types
 * dictionary_name: path of van der Waals dictionary file
 * ($KVFinder_PATH/dictionary)
 *
 * returns: success (1) or fail (0) to read van der Waals radii file
 *
 */
int read_vdw(char dictionary_name[500], vdw *DIC[500], int tablesize) {
  /* Declare variables */
  int i, flag = 1;
  char AUX[50];
  FILE *dictionary_file;
  /* Dictionary structure: {symbol, radius, *next} */
  vdw *p;

  /* Open dictionary file */
  dictionary_file = fopen(dictionary_name, "r");
  if (dictionary_file == NULL) {

    /* Print error and exit */
    fprintf(stderr, "\033[0;31mError:\033[0m Invalid dictionary file. Please "
                    "select a valid dictionary filename and try again.\n");
    exit(-1);
  }

  /* from zero to tablesize object, attribute NULL to each point of DIC */
  for (i = 0; i < tablesize; i++)
    DIC[i] = NULL;

  /* Read lines in dictionary file until reaches EOF and i < tablesize */
  for (i = -1; fscanf(dictionary_file, "%s", AUX) != EOF && i < tablesize;) {

    if (AUX[0] == '>')
      /* Count a residue */
      i++;
    else {

      /* Allocate a vdw space for p in memory */
      p = malloc(sizeof(vdw));
      /* Read radius value and save inside struct p inside radius */
      fscanf(dictionary_file, "%lf", &p->radius);
      /* Copy AUX string and paste inside struct p inside symbol */
      strcpy(p->symbol, AUX);
      /* NULL pointer to last item inserted in linked list */
      p->next = NULL;

      /* Create a linked list for accessing each residue information i indicates
      the residue. So, DIC[i] is the linked list of the residue i */
      if (DIC[i] == NULL)
        /* Allocate a vdw space for DIC[i] in memory */
        DIC[i] = malloc(sizeof(vdw));
      else
        p->next = DIC[i]->next;
      DIC[i]->next = p;
    }
  }

  /* Close dictionary file */
  fclose(dictionary_file);

  /* Return flag indicating file has been found */
  return flag;
}

/* Protein DataBank (PDB) file processing */

/*
 * Function: _create_atom
 * ----------------------
 *
 * Create a atom node
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
 * returns: atom node with atomic information
 */
atom *_create_atom(double x, double y, double z, double radius, int resnumber,
                   char resname, char chain) {
  atom *new = (atom *)malloc(sizeof(atom));

  new->x = x;
  new->y = y;
  new->z = z;
  new->radius = radius;
  new->resnumber = resnumber;
  new->chain = chain;
  new->resname = resname;
  new->next = NULL;

  return new;
}

/*
 * Function: _insert_atom
 * ----------------------
 *
 * Insert atom node in linked list
 *
 * head: pointer to linked list head
 * new: atom node
 *
 */
void _insert_atom(atom **head, atom *new) {
  atom *current;

  if (*head == NULL || (*head)->resnumber >= new->resnumber) {
    new->next = *head;
    *head = new;
  } else {
    current = *head;
    while (current->next != NULL && current->next->resnumber < new->resnumber) {
      current = current->next;
    }
    new->next = current->next;
    current->next = new;
  }
}

/*
 * Function: soft_read_pdb
 * -----------------------
 *
 * Soft read atomic information of a target PDB file
 *
 * PDB_NAME: path to a target PDB file
 * has_resnumber: whether to save residue number
 * has_chain: whether to save chain identifier
 *
 */
int soft_read_pdb(char PDB_NAME[500], int has_resnumber, int has_chain) {
  /* Declare variables */
  int flag = 1, i, resnumber;
  double x, y, z;
  char AUX[10] = "", X[10] = "", Y[10] = "", Z[10] = "", LINE[100],
       CHAIN[10] = "";
  atom *new;
  FILE *arqPDB;

  /* NULL pointer to last item in PDB information linked list */
  v = NULL;

  /* Open PDB file */
  arqPDB = fopen(PDB_NAME, "r");

  /* PDB file not found */
  if (arqPDB == NULL) {

    /* Print error and exit */
    printf("\033[0;31mError:\033[0m PDB file not found!\n");
    exit(-1);

  } else {

    /* Parse PDB */
    /* While PDB is not over, do ... */
    while (_read_line(arqPDB, LINE, 100)) {

      /* Extract Record Name */
      _extract(LINE, strlen(LINE), AUX, strlen(AUX), 0, 6);
      /* If Record Name is equal to ATOM or HETATM, do ... */
      if (!strcmp(AUX, "ATOM  ") || !strcmp(AUX, "HETATM")) {

        /*Get residue sequence number*/
        if (has_resnumber) {
          _extract(LINE, strlen(LINE), AUX, strlen(AUX), 22, 26);
          resnumber = atoi(AUX);
        }

        /* Extract x coordinate */
        _extract(LINE, strlen(LINE), X, strlen(X), 30, 38);
        x = atof(X);
        /* Extract y coordinate */
        _extract(LINE, strlen(LINE), Y, strlen(Y), 38, 46);
        y = atof(Y);
        /* Extract z coordinate */
        _extract(LINE, strlen(LINE), Z, strlen(Z), 46, 54);
        z = atof(Z);

        /* Extract chain identifier */
        if (has_chain)
          _extract(LINE, strlen(LINE), CHAIN, strlen(CHAIN), 21, 22);

        /* Save coordinate (x,y,z), residue number and chain */
        new = _create_atom(x, y, z, 0.0, resnumber, 0, CHAIN[0]);
        if (has_resnumber && has_chain)
          _insert_atom(&v, new);
        else
          _insert_atom(&v, new);
      }
    }
  }

  /* Close PDB file */
  fclose(arqPDB);

  /* Return flag indicating file has been read */
  return flag;
}

/*
 * Function: read_pdb
 * ------------------
 *
 * Read atomic information of a target PDB file
 *
 * PDB_NAME: path to a target PDB file
 * DIC: a vector containing atom name and atom radius per residue type
 * tablesize: number of residue types
 * TABLE: a vector containing residues information (index and residue name)
 * probe: probe size
 * m: x grid units
 * n: y grid units
 * o: z grid units
 * h: grid spacing (step)
 * X1: x coordinate of P1
 * Y1: y coordinate of P1
 * Z1: z coordinate of P1
 * log_file: path to log file (KVFinder.log)
 *
 */
int read_pdb(char PDB_NAME[500], vdw *DIC[500], int tablesize,
             char TABLE[500][4], double probe, int m, int n, int o, double h,
             double X1, double Y1, double Z1, FILE **log_file) {
  /* Declare variables */
  int flag = 1, i, j, number;
  char AUX[10] = "", LINE[100] = "", X[10] = "", Y[10] = "", Z[10] = "",
       RESIDUE[10] = "", ATOM_TYPE[10] = "", ATOM_SYMBOL[10] = "",
       CHAIN[10] = "";
  double x, y, z, x1, y1, z1, xaux, yaux, zaux, radius;
  atom *new;
  FILE *arqPDB;

  /* Open PDB file */
  arqPDB = fopen(PDB_NAME, "r");

  /* PDB file not found */
  if (arqPDB == NULL) {

    /* Print error and exit*/
    printf("\033[0;31mError:\033[0m PDB file not found!\n");
    exit(-1);

  }
  /* PDB file found */
  else {

    /* Parse PDB */
    /* While PDB file is not over, do ... */
    while (_read_line(arqPDB, LINE, 100)) {

      /* Extract Record Name */
      _extract(LINE, strlen(LINE), AUX, strlen(AUX), 0, 6);
      /* If Record Name is equal to ATOM or HETATM, do ... */
      if (!strcmp(AUX, "ATOM  ") || !strcmp(AUX, "HETATM")) {
        /* Get atom type */
        _extract(LINE, strlen(LINE), ATOM_TYPE, strlen(ATOM_TYPE), 12, 17);
        _remove_char(ATOM_TYPE, strlen(ATOM_TYPE), ' ');
        _remove_char(ATOM_TYPE, strlen(ATOM_TYPE), ' ');

        /* Get residue name */
        _extract(LINE, strlen(LINE), RESIDUE, strlen(RESIDUE), 17, 20);
        _remove_char(RESIDUE, strlen(RESIDUE), ' ');
        _remove_char(ATOM_TYPE, strlen(ATOM_TYPE), ' ');

        /* Get residue sequence number */
        _extract(LINE, strlen(LINE), AUX, strlen(AUX), 22, 26);
        number = atoi(AUX);

        /* Extract atom symbol */
        _extract(LINE, strlen(LINE), ATOM_SYMBOL, strlen(ATOM_SYMBOL), 76, 78);
        _remove_char(ATOM_SYMBOL, strlen(ATOM_SYMBOL), ' ');

        /* Extract chain identifier */
        _extract(LINE, strlen(LINE), CHAIN, strlen(CHAIN), 21, 22);

        /* Extract x coordinate */
        _extract(LINE, strlen(LINE), X, strlen(X), 30, 38);
        x = atof(X);
        /* Extract y coordinate */
        _extract(LINE, strlen(LINE), Y, strlen(Y), 38, 46);
        y = atof(Y);
        /* Extract z coordinate */
        _extract(LINE, strlen(LINE), Z, strlen(Z), 46, 54);
        z = atof(Z);

        /* Get radius for an atom in a specific residue based on VdW radius
         * dictionary */
        radius = _get_vdw_radius(RESIDUE, ATOM_TYPE, DIC, tablesize, TABLE,
                                 ATOM_SYMBOL, log_file);

        /* Calculate coordinate (x1, y1, z1) for atom */
        x1 = (x - X1) / h;
        y1 = (y - Y1) / h;
        z1 = (z - Z1) / h;
        xaux = x1 * cosb + z1 * sinb;
        yaux = y1;
        zaux = -x1 * sinb + z1 * cosb;
        x1 = xaux;
        y1 = yaux * cosa - zaux * sina;
        z1 = yaux * sina + zaux * cosa;

        /* Create a linked list only for atoms inside search box */
        if (x1 > 0.0 - (probe + radius) / h &&
            x1 < (double)m + (probe + radius) / h &&
            y1 > 0.0 - (probe + radius) / h &&
            y1 < (double)n + (probe + radius) / h &&
            z1 > 0.0 - (probe + radius) / h &&
            z1 < (double)o + (probe + radius) / h) {

          /* Save coordinates (x,y,z), radius, residue number and chain */
          new = _create_atom(x, y, z, radius, number,
                             _residue2code(RESIDUE), CHAIN[0]);
          _insert_atom(&v, new);
        }
      }
    }
  }

  /* Close PDB file */
  fclose(arqPDB);

  /* Return flag indicating file has been read */
  return flag;
}

/*
 * Function: _free_atom
 * --------------------
 *
 * Free linked list with PDB atomic information
 *
 */
void _free_atom() {
  atom *p;

  while (v != NULL) {
    p = v;
    v = v->next;
    free(p);
  }
}

/* parKVFinder results file processing */

/*
 * Function: write_results
 * -----------------------
 *
 * Write parKVFinder TOML-formatted results to file
 *
 * output_results: a path to a TOML-formatted results file
 * pdb_name: path to target PDB file
 * output_pdb: path to cavity PDB file
 * LIGAND_NAME: path to target ligand PDB file
 * h: grid spacing (step)
 * ncav: number of cavities
 *
 */
void write_results(char *output_results, char *pdb_name, char *output_pdb,
                   char LIGAND_NAME[500], double h, int ncav) {

  /* Declare variables */
  FILE *results_file;
  char results[1024];
  int kvnum, iterator;

  /* Open KVFinder.results.toml */
  results_file = fopen(output_results, "w");
  /* Save memory for buffer */
  memset(results, '\0', sizeof(results));
  /* Create buffer */
  setvbuf(results_file, results, _IOFBF, 1024);

  /* Write results file */
  /* File header */
  fprintf(results_file,
          "# TOML results file for parKVFinder software\n\ntitle = "
          "\"parKVFinder results file\"\n\n");

  /* Files paths */
  fprintf(results_file,
          "[FILES_PATH]\nINPUT = \"%s\"\nOUTPUT = \"%s\"\nLIGAND = \"%s\"\n\n",
          pdb_name, output_pdb, LIGAND_NAME);

  /* Parameters */
  fprintf(results_file, "[PARAMETERS]\nSTEP = %.2lf\n\n", h);

  /* Results header */
  fprintf(results_file,
          "[RESULTS]\n# Volume, area, depth and interface residues "
          "information for each cavity\n");

  /* Volume */
  fprintf(results_file,
          "\n\t[RESULTS.VOLUME]\n\t# Volume unit is cubic angstrom\n");
  for (kvnum = 0; kvnum < ncav; kvnum++) {

    fprintf(results_file, "\tK%c%c = %.2lf\n", 65 + (((kvnum) / 26) % 26),
            65 + ((kvnum) % 26), KVFinder_results[kvnum].volume);
  }

  /* Surface Area */
  fprintf(results_file,
          "\n\t[RESULTS.AREA]\n\t# Area unit is square angstrom\n");
  for (kvnum = 0; kvnum < ncav; kvnum++) {

    fprintf(results_file, "\tK%c%c = %.2lf\n", 65 + (((kvnum) / 26) % 26),
            65 + ((kvnum) % 26), KVFinder_results[kvnum].area);
  }

  /* Maximum Depth */
  fprintf(results_file,
          "\n\t[RESULTS.MAX_DEPTH]\n\t# Maximum depth unit is angstrom\n");
  for (kvnum = 0; kvnum < ncav; kvnum++) {

    fprintf(results_file, "\tK%c%c = %.2lf\n", 65 + (((kvnum) / 26) % 26),
            65 + ((kvnum) % 26), KVFinder_results[kvnum].max_depth);
  }

  /* Average Depth */
  fprintf(results_file,
          "\n\t[RESULTS.AVG_DEPTH]\n\t# Average depth unit is angstrom\n");
  for (kvnum = 0; kvnum < ncav; kvnum++) {
    fprintf(results_file, "\tK%c%c = %.2lf\n", 65 + (((kvnum) / 26) % 26),
            65 + ((kvnum) % 26), KVFinder_results[kvnum].avg_depth);
  }

  /* Average Hydropathy */
  fprintf(results_file,
          "\n\t[RESULTS.AVG_HYDROPATHY]\n\t# Average hydropathy has no unit\n");
  for (kvnum = 0; kvnum < ncav; kvnum++) {
    fprintf(results_file, "\tK%c%c = %.2lf\n", 65 + (((kvnum) / 26) % 26),
            65 + ((kvnum) % 26), KVFinder_results[kvnum].avg_hydropathy);
  }
      fprintf(results_file,
          "\tEisenbergWeiss = [ %.2lf, %.2lf,]\n", -1.42, 2.6);

  /* Interface Residues */
  fprintf(results_file,
          "\n\t[RESULTS.RESIDUES]\n\t# Interface residues for each cavity\n");
  fprintf(results_file,
          "\t# [\"residue number\",\"chain identifier\",\"residue name\"]\n");
  for (kvnum = 0; kvnum < ncav; kvnum++) {

    fprintf(results_file, "\tK%c%c = [", 65 + (((kvnum) / 26) % 26),
            65 + ((kvnum) % 26));

    for (t = KVFinder_results[kvnum].res_info; t != NULL; t = t->next) {

      if (t->next != NULL)
        fprintf(results_file, "[\"%d\",\"%c\",\"%c\"],", t->resnumber, t->chain,
                t->resname);

      else
        fprintf(results_file, "[\"%d\",\"%c\",\"%c\"]", t->resnumber, t->chain,
                t->resname);
    }

    fprintf(results_file, "]\n");
  }

  /* Print KVFinder.results.toml */
  fflush(results_file);

  /* Close KVFinder.results.toml */
  fclose(results_file);
}
