/* This file contains all the functions used to process TOML parameters file
 * used as input for KVFinder */

/* Import native modules */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Import custom modules */
#include "utils.h"
#include "tomlprocessing.h"

/* Read TOML file and pass to struct TOML */
toml *readTOML(toml *p, char *path)
{
  /* Declare variables */
  FILE *parameters_file;
  char LINE[500], config[100] = "", value[500] = "", key[100] = "";
  char *vb, *ib, *p1, *p2, *p3, *p4;
  int i = 0, equal = 0, flag_visiblebox = 0, flag_internalbox = 0, count = 0;
  int config_bar[2];

  /* Allocate memory to struct TOML */
  p = (toml *)malloc(sizeof(toml));

  /* Open parameters TOML file */
  parameters_file = fopen(path, "r");

  /* TOML file not found */
  if (parameters_file == NULL)
  {

    /* Print error and exit */
    printf("\033[0;31mError:\033[0m Invalid parameters file! Please select a "
           "valid parameters and try again.\n");
    exit(-1);
  }
  else
  {
    /* While TOML file is not over, do... */
    while (_read_line(parameters_file, LINE, 500))
    {

    HANDLE_LAST_LINE:
      /* Remove tabs */
      _remove_char(LINE, strlen(LINE), '\t');
      /* Remove whitespaces */
      _remove_char(LINE, strlen(LINE), ' ');

      /* If a '#' is found, ignore line */
      if (LINE[0] == '#' || LINE[0] == '\n' || LINE[0] == ' ')
        ;
      else
      {

        for (i = 0; i < strlen(LINE); i++)
        {

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
        if (config_bar[1] != 0)
        {

          /* Restart string config[500] */
          _initialize_string(config, strlen(config));
          /* Extract name of configuration */
          _extract(LINE, strlen(LINE), config, strlen(config), config_bar[0] + 1, config_bar[1]);

          vb = strstr(config, ".visiblebox.");
          ib = strstr(config, ".internalbox.");

          /* If configuration is visiblebox, do ... */
          if (vb != NULL)
          {
            flag_visiblebox = 1;
            count = 0;
          }

          /* If configuration is internalbox, do ... */
          if (ib != NULL)
          {
            flag_internalbox = 1;
            count = 0;
          }
        }

        /* If equal symbol is found, do ... */
        if (equal != 0)
        {

          /* Extract key */
          _extract(LINE, strlen(LINE), key, strlen(key), 0, equal);

          /* Extract value */
          _extract(LINE, strlen(LINE), value, strlen(value), equal + 1, strlen(LINE) - 1);

          /* If bracket is found, do ... */
          if (flag_visiblebox || flag_internalbox)
          {

            p1 = strstr(config, ".p1");
            p2 = strstr(config, ".p2");
            p3 = strstr(config, ".p3");
            p4 = strstr(config, ".p4");

            /* Visible box */
            if (flag_visiblebox)
            {

              /* Count number of coordinates read */
              count++;

              /* Visible point 1 : bP1 */
              if (p1 != NULL)
              {

                /* Get x, y or z coordinates */
                if (strcmp(key, "x") == 0)
                  (p->bX1) = atof(value);
                if (strcmp(key, "y") == 0)
                  (p->bY1) = atof(value);
                if (strcmp(key, "z") == 0)
                  (p->bZ1) = atof(value);
              }

              /* Visible point 2 : bP2 */
              if (p2 != NULL)
              {

                if (strcmp(key, "x") == 0)
                  (p->bX2) = atof(value);
                if (strcmp(key, "y") == 0)
                  (p->bY2) = atof(value);
                if (strcmp(key, "z") == 0)
                  (p->bZ2) = atof(value);
              }

              /* Visible point 3 : bP3 */
              if (p3 != NULL)
              {

                if (strcmp(key, "x") == 0)
                  (p->bX3) = atof(value);
                if (strcmp(key, "y") == 0)
                  (p->bY3) = atof(value);
                if (strcmp(key, "z") == 0)
                  (p->bZ3) = atof(value);
              }

              /* Visible point 4 : bP4 */
              if (p4 != NULL)
              {

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
            if (flag_internalbox)
            {

              count++;

              /* Internal point 1 : P1 */
              if (p1 != NULL)
              {

                /* Get x, y or z coordinates */
                if (strcmp(key, "x") == 0)
                  (p->X1) = atof(value);
                if (strcmp(key, "y") == 0)
                  (p->Y1) = atof(value);
                if (strcmp(key, "z") == 0)
                  (p->Z1) = atof(value);
              }

              /* Internal point 2 : P2 */
              if (p2 != NULL)
              {

                if (strcmp(key, "x") == 0)
                  (p->X2) = atof(value);
                if (strcmp(key, "y") == 0)
                  (p->Y2) = atof(value);
                if (strcmp(key, "z") == 0)
                  (p->Z2) = atof(value);
              }

              /* Internal point 3 : P3 */
              if (p3 != NULL)
              {

                if (strcmp(key, "x") == 0)
                  (p->X3) = atof(value);
                if (strcmp(key, "y") == 0)
                  (p->Y3) = atof(value);
                if (strcmp(key, "z") == 0)
                  (p->Z3) = atof(value);
              }

              /* Internal point 4 : P4 */
              if (p4 != NULL)
              {

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
          }
          else
          {

            /* Remove quotation marks */
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
  if (i > 0)
  {
    LINE[i] = '\n';
    goto HANDLE_LAST_LINE;
  }

  /* Return struct TOML */
  return p;
}