/* This file contains all the command line functions utilized by KVFinder */

/* Import native modules */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <sys/stat.h>

/* Import customs modules */
#include "dictionaryprocessing.h"
#include "matrixprocessing.h"
#include "tomlprocessing.h"
#include "pdbprocessing.h"

/* Get residues box coordinates (Xmin, Xmax, Ymin, Ymax, Zmin, Zmax) from residues box file */
void
create_residues_box (char *box_name,
                     double *Xmin,
                     double *Xmax,
                     double *Ymin,
                     double *Ymax,
                     double *Zmin,
                     double *Zmax,
                     double padding,
                     char PDB_NAME[NAME_MAX])
{
	/* Declare variables */
	int resnum;
	char chain[2];
	FILE *box_file;
	atom *p;

	/* Prepare coordinate values */
	*Xmin = 999999; *Ymin = 999999; *Zmin = 999999;
	*Xmax = -999999; *Ymax = -999999; *Zmax = -999999;

	/* Load PDB coordinates */
	PDB_load3 (PDB_NAME);

	/* Open box file */
	box_file = fopen (box_name, "r");

    /* Read file by each element (resnum_chain) */
    while (!feof (box_file)) {

        fscanf (box_file, "%d_%s[^\t][^\n]", &resnum, chain);

	    /* Check coordinates */
        for (p = v; p != NULL; p = p->next) {

            /* Find RESNUM and CHAIN */
            if (p->resnum == resnum && (p->chain == chain[0] || p->chain == chain[1])) {

                /* Update min and max coordinates */
                if (p->x < *Xmin)
                    *Xmin = p->x;
                if (p->x > *Xmax)
                    *Xmax = p->x;
                if (p->y < *Ymin)
                    *Ymin = p->y;
                if (p->y > *Ymax)
                    *Ymax = p->y;
                if (p->z < *Zmin)
                    *Zmin = p->z;
                if (p->z > *Zmax)
                    *Zmax = p->z;

            }

        }

    }

    if (*Xmax < *Xmin || *Zmax < *Zmin || *Zmax < *Zmin) {
        fprintf(stderr, "Error: Residues provided in residues box file has not been found in the PDB file!\n");
        exit(-1);
    }

	/* Free PDB coordinates */
	free_atom();

	/* Prepare coordinates and add padding in each direction */
	*Xmin = *Xmin - padding;
	*Xmax = *Xmax + padding;
	*Ymin = *Ymin - padding;
	*Ymax = *Ymax + padding;
	*Zmin = *Zmin - padding;
	*Zmax = *Zmax + padding;

}

/* Get custom box coordinates (Xmin, Xmax, Ymin, Ymax, Zmin, Zmax) from custom box file */
void
create_custom_box (char *box_name,
                   double *Xmin,
                   double *Xmax,
                   double *Ymin,
                   double *Ymax,
                   double *Zmin,
                   double *Zmax)
{
	/* Declare variables */
	int i = 0;
	double coord[6];
	FILE *box_file;

	/* Open box file */
	box_file = fopen (box_name, "r");

    /* Read file until EOF */
    while (!feof(box_file) && i < 6){
        fscanf(box_file, "%lf[^\t][^\n]", &coord[i]);
	    i++;
    }

    /* Pass values to coordinates variables */
    *Xmin = coord[0];
    *Xmax = coord[1];
    *Ymin = coord[2];
    *Ymax = coord[3];
    *Zmin = coord[4];
    *Zmax = coord[5];

}

/* Get file extension */
char
*get_file_extension (char *filename)
{
    /* Declare variables */
	char *dot = strrchr (filename, '.');

	if (!dot || dot == filename)
	    return "";
	return dot+1;

}

/* Convert int to true or false input */
char
TF (int mode,
    char **tomlmode)
{
	if (mode)
	    *tomlmode = "true";
	else
	    *tomlmode = "false";

}

/* Print TOML file inside KV_Files folder */
void
print_toml (char *toml_name,
            char OUTPUT[500],
            char BASE_NAME[500],
            char dictionary_name[500],
            char PDB_NAME[500],
            char LIGAND_NAME[500],
            int whole_protein_mode,
            char resolution_mode[7],
            int box_mode,
            int surface_mode,
            int kvp_mode,
            int ligand_mode,
            double h,
            double probe_in,
            double probe_out,
            double volume_cutoff,
            double ligand_cutoff,
            double removal_distance,
            double X1,
            double Y1,
            double Z1,
            double X2,
            double Y2,
            double Z2,
            double X3,
            double Y3,
            double Z3,
            double X4,
            double Y4,
            double Z4,
            double bX1,
            double bY1,
            double bZ1,
            double bX2,
            double bY2,
            double bZ2,
            double bX3,
            double bY3,
            double bZ3,
            double bX4,
            double bY4,
            double bZ4)
{
    /* Declare variables */
	FILE *toml_file;
	char buffer[1024], *wpmode, *bmode, *smode, *kmode, *lmode, *hmode, *emode;

    /* Create KV_Files directory */
	mkdir(combine(OUTPUT, "KV_Files/"), S_IRWXU);

    /* Open TOML file */
	toml_file = fopen (toml_name, "w");
	/* Create a buffer */
	memset (buffer, '\0', sizeof (buffer));
	/* Define buffer as writing buffer of size 1024 */
	setvbuf (toml_file, buffer, _IOFBF, 1024);

	/* Print TOML file */
	fprintf (toml_file, "# TOML configuration file for parKVFinder software.\n");
	fprintf (toml_file,"\ntitle = \"parKVFinder parameters file\"\n");

	fprintf (toml_file, "\n[FILES_PATH]\n");
	fprintf (toml_file, "# The path of van der Waals radii dictionary for parKVFinder.\n");
	fprintf (toml_file, "dictionary = \"%s\"\n", dictionary_name);
	fprintf (toml_file, "# The path of the input PDB file.\n");
	fprintf (toml_file, "pdb = \"%s\"\n", PDB_NAME);
	fprintf (toml_file, "# The path of the output directory.\n");
	fprintf (toml_file, "output = \"%s\"\n", OUTPUT);
	fprintf (toml_file, "# Base name for output files.\n");
	fprintf (toml_file, "base_name = \"%s\"\n", BASE_NAME);
	fprintf (toml_file, "# Path for the ligand's PDB file.\n");
	fprintf (toml_file, "ligand = \"%s\"\n", LIGAND_NAME);

	fprintf (toml_file, "\n[SETTINGS]\n");
	fprintf (toml_file, "# Settings for parKVFinder software\n");
	fprintf (toml_file, "\n\t[SETTINGS.modes]\n");
	fprintf (toml_file, "\t# Whole Protein mode defines the search space as the whole protein.\n");
	TF (whole_protein_mode, &wpmode);
	fprintf (toml_file, "\twhole_protein_mode = %s\n", wpmode);
	fprintf (toml_file, "\t# Box Adjustment mode defines the search space as a box that includes a specific region.\n");
	TF (box_mode, &bmode);
	fprintf (toml_file, "\tbox_mode = %s\n", bmode);
	fprintf (toml_file, "\t# Resolution mode implicitly sets the step size (grid spacing) of the 3D grid.\n");
	fprintf (toml_file, "\t# If set to High, sets a voxel volume of 0.2. If set to Medium, sets a voxel volume of 0.1. If set to Low, sets a voxel volume of 0.01. If set to Off, the step size must be set explicitly.\n");
	fprintf (toml_file, "\tresolution_mode = \"%s\"\n", resolution_mode);
	fprintf (toml_file, "\t# Surface mode defines the type of surface representation to be applied, van der Waals molecular surface (true) or solvent accessible surface (false).\n");
	TF (surface_mode, &smode);
	fprintf (toml_file, "\tsurface_mode = %s\n", smode);
	fprintf (toml_file, "\t# Cavity output mode defines whether cavities are exported to the output PDB file as filled cavities (true) or filtered cavities (false).\n");
	TF (kvp_mode, &kmode);
	fprintf (toml_file, "\tkvp_mode = %s\n", kmode);
	fprintf (toml_file, "\t# Ligand adjustment mode defines the search space around the ligand.\n");
	TF (ligand_mode, &lmode);
	fprintf (toml_file, "\tligand_mode = %s\n", lmode);

	fprintf (toml_file, "\n\t[SETTINGS.step_size]\n");
	fprintf (toml_file, "\t# Sets the 3D grid spacing. It directly affects accuracy and runtime.\n");
	fprintf (toml_file, "\tstep_size = %.2lf\n", h);

	fprintf (toml_file, "\n\t[SETTINGS.probes]\n");
	fprintf (toml_file, "\t# parKVFinder works with a two sized probe system. A smaller probe, called Probe In, and a bigger one, called Probe Out, rolls around the protein.\n");
	fprintf (toml_file, "\t# Points reached by the Probe In, but not the Probe Out are considered cavity points.\n");
	fprintf (toml_file, "\t# Sets Probe In diameter. Default: 1.4 angstroms.\n");
	fprintf (toml_file, "\tprobe_in = %.2lf\n", probe_in);
	fprintf (toml_file, "\t# Sets Probe Out diameter. Default: 4.0 angstroms.\n");
	fprintf (toml_file, "\tprobe_out = %.2lf\n", probe_out);

	fprintf (toml_file, "\n\t[SETTINGS.cutoffs]\n");
	fprintf (toml_file, "\t# Sets a volume cutoff for the detected cavities. Default: 5.0 angstroms.\n");
	fprintf (toml_file, "\tvolume_cutoff = %.2lf\n", volume_cutoff);
	fprintf (toml_file, "\t# Sets a distance cutoff for a search space around the ligand in ligand adjustment mode. Default: 5.0 angstroms.\n");
	fprintf (toml_file, "\tligand_cutoff = %.2lf\n", ligand_cutoff);
	fprintf (toml_file, "\t# Sets a removal distance for the cavity frontier, which is defined by comparing Probe In and Probe Out surfaces. Default: 2.4 angstroms.\n");
	fprintf (toml_file, "\tremoval_distance = %.2lf\n", removal_distance);

	fprintf (toml_file, "\n\t[SETTINGS.visiblebox]\n");
	fprintf (toml_file, "\t# Coordinates of the vertices that define the visible 3D grid. Only four points are required to define the search space.\n");
	fprintf (toml_file, "\tbP1 = {bX1 = %.2lf, bY1 = %.2lf, bZ1 = %.2lf}\n", bX1, bY1, bZ1);
	fprintf (toml_file, "\tbP2 = {bX2 = %.2lf, bY2 = %.2lf, bZ2 = %.2lf}\n", bX2, bY2, bZ2);
	fprintf (toml_file, "\tbP3 = {bX3 = %.2lf, bY3 = %.2lf, bZ3 = %.2lf}\n", bX3, bY3, bZ3);
	fprintf (toml_file, "\tbP4 = {bX4 = %.2lf, bY4 = %.2lf, bZ4 = %.2lf}\n", bX4, bY4, bZ4);

	fprintf (toml_file, "\n\t[SETTINGS.internalbox]\n");
	fprintf (toml_file, "\t# Coordinates of the internal 3D grid. Used for calculations.\n");
	fprintf (toml_file, "\tP1 = {X1 = %.2lf, Y1 = %.2lf, Z1 = %.2lf}\n", X1, Y1, Z1);
	fprintf (toml_file, "\tP2 = {X2 = %.2lf, Y2 = %.2lf, Z2 = %.2lf}\n", X2, Y2, Z2);
	fprintf (toml_file, "\tP3 = {X3 = %.2lf, Y3 = %.2lf, Z3 = %.2lf}\n", X3, Y3, Z3);
	fprintf (toml_file, "\tP4 = {X4 = %.2lf, Y4 = %.2lf, Z4 = %.2lf}\n", X4, Y4, Z4);

	fflush (toml_file);
	fclose (toml_file);

}

/* Check if numeric inputs are in correct format */
int
check_input (char *optarg,
             char *error)
{
    /* Declare variables */
	int i, hitDecimal = 0, len;

	/* Get option length */
	len = strlen (optarg);

	/* Loop through string option */
	for (i = 0; i < len; i++) {

		/* Count number of points in numeric option */
		if (optarg[i] == '.' && !hitDecimal)
		    hitDecimal = 1;
		else
		    /* Check if all chars, that are not a point, are not digits */
		    if (!isdigit (optarg[i])) {

		        /* Found a char that is not a digit, print error*/
			    fprintf (stderr, "%s", error);
			    /* Exit code */
			    exit (-1);

		}

	}

	/* Function worked, return 1 */
	return 1;

}

void
print_version ()
{
	printf ("parKVFinder (parallel KVFinder) v1.0\n");
}

void
init25 (char S[25])
{
    /* Declare variables */
	int i;

	for (i = 0; i < 25; i++)
	    S[i] = ' ';

}

void
init55 (char S[55])
{
    /* Declare variables */
	int i;

	for (i = 0; i < 55; i++)
	    S[i] = ' ';

}

void
init80 (char S[80])
{
    /* Declare variables */
	int i;

	for (i = 0; i < 80; i++)
	    S[i] = ' ';

}

void
print_header ()
{

	fprintf (stdout, "parKVFinder (parallel KVFinder) software identifies and describes cavities in\n");
	fprintf (stdout, "target biomolecular structure using a dual probe system.\n");
	fprintf (stdout, "\n");
	fprintf (stdout, "The description includes spatial and constitutional characterization. Spatial \n");
	fprintf (stdout, "description includes shape, volume and area. Constitutional description includes\n");
	fprintf (stdout, "amino acids that form the identified cavities.\n");
	fprintf (stdout, "\n");

}

void
print_usage ()
{

	fprintf (stdout, "Usage: parKVFinder PDB [options],\n");
	fprintf (stdout, "\twhere PDB is a path to a target PDB file.\n");
  	fprintf (stdout, "\n");
	fprintf (stdout, "Options:\n");
	fprintf (stdout, "  -h, --help\n");
	fprintf (stdout, "\t  Display this help message.\n");
	fprintf (stdout, "  -v, --version\n");
	fprintf (stdout, "\t  Display parKVFinder version.\n");
	fprintf (stdout, "  --verbose\n");
	fprintf (stdout, "\t  Print extra information to stdout.\n");
	fprintf (stdout, "\n");

}

void
print_options ()
{

	/* GENERAL KVFINDER PARAMETERS */
  	fprintf (stdout, "General options:\n");
	fprintf (stdout, "  -p, --parameters\t[<.toml>]\n");
	fprintf (stdout, "\t  Define path to parameters file.\n");
	fprintf (stdout, "  -d, --dictionary\t[<dictionary>]\n");
	fprintf (stdout, "\t  Define path to a custom dictionary file.\n");
	fprintf (stdout, "  -r, --resolution\t<enum>\t(Low)\n");
	fprintf (stdout, "\t  Define resolution mode. Options include: Off, Low, Medium and High.\n");
	fprintf (stdout, "  -s, --step\t\t<real>\t\t(0.0)\n");
	fprintf (stdout, "\t  Define step size (grid spacing).\n");
	fprintf (stdout, "  -i, --probe_in\t<real>\t\t(1.4)\n");
	fprintf (stdout, "\t  Define probe in size.\n");
	fprintf (stdout, "  -o, --probe_out\t<real>\t\t(4.0)\n");
	fprintf (stdout, "\t  Define probe out size.\n");
	fprintf (stdout, "  --volume_cutoff\t<real>\t\t(5.0)\n");
	fprintf (stdout, "\t  Define cavities volume filter.\n");
	fprintf (stdout, "  --removal_distance\t<real>\t\t(2.4)\n");
	fprintf (stdout, "\t  Define removal distance when comparing probes surfaces.\n");
	fprintf (stdout, "  -t, --template\t\t\t(paramters.toml)\n");
	fprintf (stdout, "\t  Create a parameter file template with defined parameters in current\n");
	fprintf (stdout, "\t  working directory.\n");
	fprintf (stdout, "\n");
	/* BOX ADJUSTMENT PARAMETERS */
	fprintf (stdout, "Box adjustment options:\n");
	fprintf (stdout, "  -B, --box\n");
	fprintf (stdout, "\t  Define a search box mode where parKVFinder will detect cavities.\n");
	fprintf (stdout, "  --custom_box\t\t[<file>]\n");
	fprintf (stdout, "\t  Define a custom search box based on a file containing the minimum and \n");
	fprintf (stdout, "\t  maximum cartesian values of each axis in angstrom.\n");
	fprintf (stdout, "  --residues_box\t[<file>]\n");
	fprintf (stdout, "\t  Automatically set a search box based a file containing a tab-separated\n");
	fprintf (stdout, "\t  list of residues.\n");
	fprintf (stdout, "  --padding\t\t<real>\t\t(3.5)\n");
	fprintf (stdout, "\t  Define residues box padding. Adds a length in each box direction.\n");
	fprintf (stdout, "\n");
	/* SURFACE MODE */
	fprintf (stdout, "Surface options:\n");
	fprintf (stdout, "  -S, --surface\t\t<enum>\t(VdW)\n");
	fprintf (stdout, "\t  Define a surface representation. Options include: SAS and VdW. SAS\n");
	fprintf (stdout, "\t  specifies solvent accessible surface. VdW specifies van der Waals\n");
	fprintf (stdout, "\t  molecular surface.\n");
	fprintf (stdout, "\n");
	/* LIGAND ADJUSTMENT PARAMETERS */
	fprintf (stdout, "Ligand options:\n");
	fprintf (stdout, "  -L, --ligand\t\t[<.pdb>]\n");
	fprintf (stdout, "\t  Define path to ligand PDB file.\n");
	fprintf (stdout, "  --ligand_cutoff\t<real>\t\t(5.0)\n");
	fprintf (stdout, "\t  Define ligand radius distance cutoff.\n");
	fprintf (stdout, "\n");

}

void
print_help()
{

	fprintf (stdout, "================================================================================\n");
	fprintf (stdout, "============================= parKVFinder help menu ============================\n");
	print_header ();
	print_usage ();
	print_options ();
	fprintf (stdout, "================================================================================\n");
	fprintf (stdout, "================================================================================\n");

}

int
argparser (int argc,
           char **argv,
           int *box_mode,
           int *kvp_mode,
           int *ligand_mode,
           int *surface_mode,
           int *whole_protein_mode,
           char PDB_NAME[500],
           char LIGAND_NAME[500],
           char dictionary_name[500],
           char OUTPUT[500],
           char BASE_NAME[500],
           char resolution_flag[7],
           double *h,
           double *probe_in,
           double *probe_out,
           double *volume_cutoff,
           double *ligand_cutoff,
           double *removal_distance,
           double *X1,
           double *Y1,
           double *Z1,
           double *X2,
           double *Y2,
           double *Z2,
           double *X3,
           double *Y3,
           double *Z3,
           double *X4,
           double *Y4,
           double *Z4,
           double *bX1,
           double *bY1,
           double *bZ1,
           double *bX2,
           double *bY2,
           double *bZ2,
           double *bX3,
           double *bY3,
           double *bZ3,
           double *bX4,
           double *bY4,
           double *bZ4)
{

    /* Declare variables */
	/* Flag set by ‘--verbose’. */
	static int verbose_flag = 0;
	/* Flag set by '--ligand or -L' */
	*ligand_mode = 0;
	/* Flag set by '--filled or -k' */
	*kvp_mode = 0;
	/* Flag set by '--box or -B' */
	*box_mode = 0;
	*whole_protein_mode = 1;

	/* Get current directory */
	char cwd[256];
	getcwd (cwd, sizeof (cwd));

	/* Declare variables */
	char *parameters_name, *template_name, *toml_name, *box_name;
	double padding;

	/* Declare counters */
	int point = 0, bar = -1, i = 0, j = 0, c;

	/* Declare flags */
	/* Flag for paths */
	int p_flag = 0, d_flag = 0, pdb_flag = 0, l_flag = 0, t_flag = 0;
	/* Flag for strings */
	int r_flag = 0,
	    surface_flag = 0;
	/* Flag for doubles */
	int o_flag = 0,
	    i_flag = 0,
	    s_flag = 0,
	    vc_flag = 0,
	    lc_flag = 0,
	    rd_flag = 0,
	    cb_flag = 0,
	    rb_flag = 0,
	    prb_flag = 0;

	while (1) {

		static struct option long_options[] = {

			/* Options */
			{"verbose", no_argument, &verbose_flag, 1},
			{"version", no_argument, NULL, 'v'},
			{"help", no_argument, NULL, 'h'},
			/* File paths */
			{"parameters", required_argument, NULL, 'p'},
			{"dictionary", required_argument, NULL, 'd'},
			{"ligand", required_argument, NULL, 'L'},
			{"template", no_argument, NULL, 't'},
			/* Modes */
			{"filled", no_argument, NULL, 'K'},
			{"surface", no_argument, NULL, 'S'},
			{"box", no_argument, NULL, 'B'},
			/* Settings */
			{"resolution", required_argument, NULL, 'r'},
			{"step", required_argument, NULL, 's'},
			{"probe_in", required_argument, NULL, 'i'},
			{"probe_out", required_argument, NULL, 'o'},
			{"volume_cutoff", required_argument, NULL, 0},
			{"ligand_cutoff", required_argument, NULL, 0},
			{"removal_distance", required_argument, NULL, 0},
			/* Custom box settings */
			{"residues_box", required_argument, NULL, 0},
			{"padding", required_argument, NULL, 0},
			{"custom_box", required_argument, NULL, 0},

			{NULL, no_argument, NULL, 0}

		};

		/* getopt_long stores the option index here */
		int option_index = 0;

		c = getopt_long (argc, argv, "vhp:d:L:t::r:s:i:o:B", long_options, &option_index);

		if (c == -1)
			break;

		switch (c) {

			/* LONG OPTIONS */
			/* Handle of long options without a short arg */
			case 0:
				/* LIGAND CUTOFF */
				if (strcmp ("ligand_cutoff", long_options[option_index].name) == 0) {
					if (check_input (optarg, "Error: Invalid ligand cutoff input!\n")) {
						*ligand_cutoff = atof (optarg);
						lc_flag = 1;
					}
				}
				/* VOLUME CUTOFF */
				if (strcmp ("volume_cutoff", long_options[option_index].name) == 0) {
					if (check_input (optarg, "Error: Invalid volume cutoff input!\n")) {
						*volume_cutoff = atof (optarg);
						vc_flag = 1;
					}
				}
				/* REMOVAL DISTANCE */
				if (strcmp ("removal_distance", long_options[option_index].name) == 0) {
					if (check_input (optarg, "Error: Invalid removal distance input!\n")) {
						*removal_distance = atof (optarg);
						rd_flag = 1;
					}
				}
				/* BOX MODE PARAMETERS */
				/* residues box */
				if (strcmp ("residues_box", long_options[option_index].name) == 0) {
					box_name = optarg;
					rb_flag = 1;
					if (access (box_name, F_OK)) {
						fprintf (stderr, "Error: Residues list file does not exist!\n");
						exit (-1);
					}
				}
				/* residues box padding */
				if (strcmp ("padding", long_options[option_index].name) == 0) {
					if (check_input (optarg, "Error: Invalid residue box padding value input!\n")) {
						padding = atof (optarg);
						prb_flag = 1;
					}
				}
				/* custom box */
				if (strcmp ("custom_box", long_options[option_index].name) == 0) {
					box_name = optarg;
					cb_flag = 1;
					if (access (box_name, F_OK)) {
						fprintf (stderr, "Error: Custom box file does not exist!\n");
						exit (-1);
					}
				}
				break;

			/* SHORT OPTIONS */
            /* SURFACE MODE */
			case 'S':
				/* If input is SAS, do ... */
				if (strcmp (optarg, "SAS")  == 0) {
					*surface_mode = 0;
				}
				/* If input is VdW, do ... */
				else
				    if (strcmp (optarg, "VdW") == 0) {
					*surface_mode = 1;
				}
				/* If input is not SAS or VdW, print error */
				else {
					fprintf (stderr, "Error: Wrong surface representation selected!\nPossible inputs: SAS, VdW.\n");
					exit (-1);
				}
				surface_flag = 1;
			    break;

            /* BOX MODE */
			case 'B':
				*box_mode = 1;
				*whole_protein_mode = 0;
			    break;

            /* FILLED CAVITIES MODE (KVP MODE) */
			case 'K':
				*kvp_mode = 1;
			    break;

            /* TEMPLATE */
			case 't':
				if (optarg != NULL)
				    template_name = optarg;
				else
				    template_name = "parameters.toml";
				t_flag = 1;
			    break;

            /* LIGAND MODE + LIGAND PATH */
			case 'L':
				snprintf (LIGAND_NAME, 500, "%s", realpath (optarg, NULL));
				*ligand_mode = 1;
				l_flag = 1;
				if (access (LIGAND_NAME, F_OK)) {
					fprintf (stderr, "Error: Ligand PDB file does not exist.\n");
					exit (-1);
				}
			    break;

            /* PROBE IN */
			case 'i':
			    /*Check if input is numeric*/
				if (check_input (optarg, "Error: Invalid probe in input!\n")) {
				    /*Save probe in inside probe_in variable*/
					*probe_in = atof (optarg);
					i_flag = 1;
				}
			    break;

            /* PROBE OUT */
			case 'o':
			    /* Check if input is numeric */
				if (check_input (optarg, "Error: Invalid probe out input!\n")) {
				    /* Save probe out inside probe_out variable */
					*probe_out = atof (optarg);
					o_flag = 1;
				}
			    break;

            /* STEP SIZE */
			case 's':
				/*Check if input is numeric*/
				if (check_input (optarg, "Error: Invalid step size input!\n")) {
				    /* Save step size inside h variable */
					*h = atof (optarg);
					s_flag = 1;
				}
			    break;

            /* RESOLUTION */
			case 'r':
			    /* Save resolution inside resolution_flag variable */
				snprintf (resolution_flag, 7, "%s", optarg);
				/* Check if input is Off, Low, Medium, High */
				if (strcmp (resolution_flag, "Off")  == 0);
				else
				    if (strcmp (resolution_flag, "High") == 0);
				else
				    if (strcmp (resolution_flag, "Medium") == 0);
				else
				    if (strcmp (resolution_flag, "Low") == 0);
				/* If input is not Off, Low, Medium, High, print error */
				else {
					fprintf (stderr, "Error: Wrong resolution selected!\nPossible inputs: Off, Low, Medium, High.\n");
					exit(-1);
				}
				r_flag = 1;
			    break;

            /* DICTIONARY */
			case 'd':
				snprintf (dictionary_name, 500, "%s", realpath (optarg, NULL));
				d_flag = 1;
				if (access (dictionary_name, F_OK)) {
					fprintf (stderr, "Error: Dictionary file does not exist!\n");
					exit (-1);
				}
			    break;

            /* PARAMETERS FILE */
			case 'p':
				parameters_name = optarg;
				p_flag = 1;
				/* Check if parameters file exists */
				if (access(parameters_name, F_OK)) {
					fprintf (stderr, "Error: Parameter file does not exist!\n");
					exit (-1);
				}
				/* Check parameters file extension*/
				if (strcmp (get_file_extension(parameters_name), "toml")) {
					fprintf (stderr, "Error: Wrong parameters file extension!\narg: [\'%s\']\n", parameters_name);
					exit (-1);
				}
			    break;

            /* HELP */
			case 'h':
				print_help ();
				exit(0);

            /* VERSION */
			case 'v':
				print_version ();
				exit (0);

			case ':':
				exit (-1);
			    break;

			case '?':
			/* getopt_long already printed an error message. */
				exit (-1);
			    break;

		}

	}

	/* COMMAND-LINE ARGUMENTS PROCESSING */

	/* Parameters file test conditions
	Parameters file must be set alone */
	if (p_flag) {
		if (d_flag || l_flag || t_flag || r_flag || o_flag || i_flag || s_flag || vc_flag || lc_flag || rd_flag) {
			fprintf (stderr, "Error: Just define parameters file argument!\n");
			exit (-1);
		}
		else {
			/* Read parameters TOML files */
			toml *parameters = readTOML (param, parameters_name); /*Read TOML file inside struct TOML*/

			/* Save TOML parameters from struct TOML parameters to KVFinder variables */
            *X1 = parameters->X1; *Y1 = parameters->Y1; *Z1 = parameters->Z1;
		    *X2 = parameters->X2; *Y2 = parameters->Y2; *Z2 = parameters->Z2;
		    *X3 = parameters->X3; *Y3 = parameters->Y3; *Z3 = parameters->Z3;
		    *X4 = parameters->X4; *Y4 = parameters->Y4; *Z4 = parameters->Z4;
		    *bX1 = parameters->bX1; *bY1 = parameters->bY1; *bZ1 = parameters->bZ1;
		    *bX2 = parameters->bX2; *bY2 = parameters->bY2; *bZ2 = parameters->bZ2;
		    *bX3 = parameters->bX3; *bY3 = parameters->bY3; *bZ3 = parameters->bZ3;
		    *bX4 = parameters->bX4; *bY4 = parameters->bY4; *bZ4 = parameters->bZ4;
            strcpy(PDB_NAME, parameters->PDB_NAME);
			strcpy(OUTPUT, parameters->OUTPUT);
			strcpy(BASE_NAME, parameters->BASE_NAME);
			strcpy(LIGAND_NAME, parameters->LIGAND_NAME);
			strcpy(dictionary_name, parameters->dictionary_name);
			strcpy(resolution_flag, parameters->resolution_flag);
			*volume_cutoff = parameters->volume_cutoff;
			*ligand_cutoff = parameters->ligand_cutoff;
			*removal_distance = parameters->removal_distance;
			*probe_in = parameters->probe_in;
			*probe_out = parameters->probe_out;
			*h = parameters->h;
			*whole_protein_mode = parameters->whole_protein_mode;
			*box_mode = parameters->box_mode;
			*surface_mode = parameters->surface_mode;
			*kvp_mode = parameters->kvp_mode;
			*ligand_mode = parameters->ligand_mode;

            /* Print in shell the PDB path that is running in KVFinder */
	        if (verbose_flag)
	            fprintf (stdout, "[PID %u] Running parKVFinder for: %s\n", getpid (), PDB_NAME);

            /* Free struct TOML */
			free (param);
			/* Loaded parameters, return to main script */
			return verbose_flag;
		}
	}

	/* PROCESS EXTRA ARGUMENTS - PDB FILE PROCESSING */
	/* User provided more than one PDB file */
	if (argc-optind > 1) {

		fprintf (stderr, "Error: Incorrect number of PDB files!\nargs: [\'%s\'", argv[optind++]);
		while (optind < argc) {
			fprintf (stderr, ", \'%s\'", argv[optind++]);
		}
		fprintf (stderr, "]\n\n");
		print_usage ();
		exit (-1);

	}
	else
	    /* User provided one PDB file */
	    if (optind < argc) {

            /* Save PDB full path inside PDB_NAME variable */
            snprintf (PDB_NAME, 500, "%s", realpath (argv[optind], NULL));

            /* Loop through string PDB_NAME */
            for (i = 0; PDB_NAME[i] != '\0'; i++) {
                /* Get position of last bar symbol in PDB_NAME */
                if (PDB_NAME[i] == '/')
                    bar = i;
                /* Get position of last point symbol in PDB_NAME */
                if (PDB_NAME[i] == '.')
                    point = i;
            }

            /* Copy PDB name to OUTPUT */
            for (i = 0, j = 0; i < bar+1; i++, j++)
                OUTPUT[j] = PDB_NAME[i];
            /* Put end symbol in OUTPUT */
            OUTPUT[j] = '\0';

            /* Copy PDB name to OUTPUT */
            for (i = bar+1, j = 0; i < point; i++, j++)
                BASE_NAME[j] = PDB_NAME[i];
            /* Put end symbol in BASE_NAME */
            BASE_NAME[j] = '\0';

            pdb_flag = 1;

            /* Print in shell the PDB path that is running in KVFinder */
	        if (verbose_flag)
	            fprintf (stdout, "[PID %u] Running parKVFinder for: %s\n", getpid (), PDB_NAME);

            /* Check if provided PDB file exist */
            if (access (PDB_NAME, F_OK)) {
                fprintf (stderr, "Error: PDB file does not exist.\n");
                exit (-1);
            }
            /* Check PDB file extension*/
            if (strcmp (get_file_extension (PDB_NAME), "pdb")) {
                fprintf (stderr, "Error: Wrong PDB file extension!\narg: [\'%s\']\n", PDB_NAME);
            }

	    }
        /* User do not provide a PDB file */
        else {
            fprintf (stderr, "Error: Missing path to PDB file!\n\n");
            print_usage ();
            exit (-1);
        }

	/* Set default values for non-defined arguments */
	/* Path to dictionary file */
	if (!d_flag){

	    /* Set default path to dictionary */
		strcpy (dictionary_name, getenv ("KVFinder_PATH"));
		strcat (dictionary_name, "/dictionary");
		if (verbose_flag)
		    fprintf (stdout, "> Setting \'dictionary_name\' to default file: %s\n", dictionary_name);

	}
	/* Probe in */
	if (!i_flag) {

	    /* Set default value to probe in */
		*probe_in = 1.4;
		if (verbose_flag)
		    fprintf(stdout, "> Setting \'probe_in\' to default value: %.2lf\n", *probe_in);

	}
	/* Probe out */
	if (!o_flag) {

	    /* Set default value to probe out */
		*probe_out = 4.0;
		if (verbose_flag)
		    fprintf (stdout, "> Setting \'probe_out\' to default value: %.2lf\n", *probe_out);

	}
	/* Volume cutoff */
	if (!vc_flag) {

	    /* Set default value to volume cutoff */
		*volume_cutoff = 5.0;
		if (verbose_flag)
		    fprintf (stdout, "> Setting \'volume_cutoff\' to default value: %.2lf\n", *volume_cutoff);

	}
	/* Removal distance */
	if (!rd_flag) {

	    /* Set default value to removal distance */
		*removal_distance = 2.4;
		if (verbose_flag)
		    fprintf (stdout, "> Setting \'removal_distance\' to default value: %.2lf\n", *removal_distance);

	}
	/* Step size and Resolution */
	if (!s_flag && !r_flag) {

		*h = 0.0;
		snprintf (resolution_flag, 7, "%s", "Low");
		if (verbose_flag)
		    fprintf (stdout,
		             "> Setting \'step\' (grid spacing) to %.2lf.\n> Setting \'resolution\' to flag: %s\n",
		             *h,
		             resolution_flag);

	}
	else
	    if (s_flag || r_flag) {

            if (s_flag) {

                snprintf (resolution_flag, 7, "%s", "Off");
                if (verbose_flag)
                    fprintf (stdout,
                             "> User set \'step\' (grid spacing) to %.2lf. Setting \'resolution\' to flag: %s\n",
                             *h,
                             resolution_flag);

            }
            if (r_flag) {

                if (strcmp (resolution_flag, "Off") == 0) {
                    fprintf (stderr, "Error: Resolution mode is Off! Step size (grid spacing) must be defined!\n");
                    exit (-1);
                }
                *h = 0.0;
                if (verbose_flag)
                    fprintf (stdout,
                             "> User set \'resolution\' to %s. Setting \'resolution_flag\' to default value: %.2lf\n",
                             resolution_flag,
                             *h);

            }

	}
	else {

		fprintf (stderr, "Error: Resolution mode is On! Step size (grid spacing) should not be defined!\n");
		exit (-1);

	}
	/* Ligand mode */
	if (!lc_flag) {

		*ligand_cutoff = 5.0;
		if(verbose_flag) fprintf(stdout, "> Setting \'ligand_cutoff\' to default value: %.2lf\n", *ligand_cutoff);

	}
	else
	    /* Ligand mode is not set */
		if (!*ligand_mode) {

		    fprintf (stderr, "Error: Path to ligand PDB file is not provided! ");
			fprintf (stderr, "Define path through \'-L\' or \'--ligand_mode\' flag.\n");
			exit(-1);

		}
	/* Surface mode */
    /* Surface mode is not defined */
	if (!surface_flag) {

		*surface_mode = 1;
		if (verbose_flag)
		    fprintf (stdout, "> Chosen \'surface\' representation: VdW\n");

	}
	/* Box mode */
	/* Box mode is Off */
	if (!*box_mode) {

		/* visiblebox coordinates */
		*bX1 = 0; *bY1 = 0; *bZ1 = 0;
		*bX2 = 0; *bY2 = 0; *bZ2 = 0;
		*bX3 = 0; *bY3 = 0; *bZ3 = 0;
		*bX4 = 0; *bY4 = 0; *bZ4 = 0;

		/* internalbox coordinates*/
		*X1 = (*bX1) - (*probe_out); *Y1 = (*bY1) - (*probe_out); *Z1 = (*bZ1) - (*probe_out);
		*X2 = (*bX2) + (*probe_out); *Y2 = (*bY2) - (*probe_out); *Z2 = (*bZ2) - (*probe_out);
		*X3 = (*bX3) - (*probe_out); *Y3 = (*bY3) + (*probe_out); *Z3 = (*bZ3) - (*probe_out);
		*X4 = (*bX4) - (*probe_out); *Y4 = (*bY4) - (*probe_out); *Z4 = (*bZ4) + (*probe_out);

		if (rb_flag || cb_flag) {

		    fprintf (stderr, "Error: Whole protein mode and Box adjustment options are chosen together! ");
			fprintf (stderr, "Define \'box\' flag or remove Box adjustment option (\'custom_box\' or \'residues_box\').\n");
			exit (-1);

		}
		if (prb_flag) {

			fprintf (stderr, "Error: Residue box mode is not chosen! Residue box padding should not be defined.\n");
			exit (-1);

		}
		if (verbose_flag)
		    fprintf (stdout, "> Running parKVFinder for whole biomolecular structure\n");

	}
	/* Box mode is On */
	else {

		if (verbose_flag)
		    fprintf (stdout, "> Running parKVFinder for custom search box.\n");

		if (rb_flag && cb_flag) {

            fprintf (stderr, "Error: Just choose one box adjustment options! ");
			fprintf (stderr, "Define \'residues_box\' or \'custom_box\' options.\n");
			exit (-1);

		}
		else
		    /* Only one box adjustment option is chosen */
		    if (rb_flag || cb_flag) {

		        /* Residues box */
                if (rb_flag) {

                     /* Padding not defined */
                    if (!prb_flag)
                        padding = 3.5;
                    create_residues_box (box_name, bX1, bX2, bY1, bY3, bZ1, bZ4, padding, PDB_NAME);

                }
                /* Custom box */
                if (cb_flag) {

                     /* Padding incorrectly defined */
                    if (prb_flag) {

                        fprintf (stderr, "Error: Residue box padding should not be defined in custom box option!\n");
                        exit (-1);

                    }
                    /* Save Xmin, Xmax, Ymin, Ymax, Zmin, Zmax */
                    create_custom_box (box_name, bX1, bX2, bY1, bY3, bZ1, bZ4);

                }
                /* Create visible box */
                *bX1 = *bX1; *bY1 = *bY1; *bZ1 = *bZ1;
                *bX2 = *bX2; *bY2 = *bY1; *bZ2 = *bZ1;
                *bX3 = *bX1; *bY3 = *bY3; *bZ3 = *bZ1;
                *bX4 = *bX1; *bY4 = *bY1; *bZ4 = *bZ4;
                /* Create internal box */
                *X1 = (*bX1) - (*probe_out); *Y1 = (*bY1) - (*probe_out); *Z1 = (*bZ1) - (*probe_out);
                *X2 = (*bX2) + (*probe_out); *Y2 = (*bY2) - (*probe_out); *Z2 = (*bZ2) - (*probe_out);
                *X3 = (*bX3) - (*probe_out); *Y3 = (*bY3) + (*probe_out); *Z3 = (*bZ3) - (*probe_out);
                *X4 = (*bX4) - (*probe_out); *Y4 = (*bY4) - (*probe_out); *Z4 = (*bZ4) + (*probe_out);

		    }
		    else {

                fprintf (stderr, "Error: Choose a box adjustment option! ");
                fprintf (stderr, "Define \'residues_box\' or \'custom_box\' options.\n");
                exit (-1);

            }

	}

	/* Print template parameters file */
	if (t_flag) {
		print_toml (template_name,
		            OUTPUT,
		            BASE_NAME,
		            dictionary_name,
		            PDB_NAME,
		            LIGAND_NAME,
		            *whole_protein_mode,
		            resolution_flag,
		            *box_mode,
		            *surface_mode,
		            *kvp_mode,
		            *ligand_mode,
		            *h,
		            *probe_in,
		            *probe_out,
		            *volume_cutoff,
		            *ligand_cutoff,
		            *removal_distance,
		            *X1, *Y1, *Z1, *X2, *Y2, *Z2, *X3, *Y3, *Z3, *X4, *Y4, *Z4,
		            *bX1, *bY1, *bZ1, *bX2, *bY2, *bZ2, *bX3, *bY3, *bZ3, *bX4, *bY4, *bZ4);
		exit (0);
	}

	toml_name = combine (combine (combine (OUTPUT, "KV_Files/parameters_"), BASE_NAME), ".toml");
	print_toml (toml_name,
	            OUTPUT,
	            BASE_NAME,
	            dictionary_name,
	            PDB_NAME,
	            LIGAND_NAME,
	            *whole_protein_mode,
	            resolution_flag,
	            *box_mode,
	            *surface_mode,
	            *kvp_mode,
	            *ligand_mode,
	            *h,
	            *probe_in,
	            *probe_out,
	            *volume_cutoff,
	            *ligand_cutoff,
	            *removal_distance,
	            *X1, *Y1, *Z1, *X2, *Y2, *Z2, *X3, *Y3, *Z3, *X4, *Y4, *Z4,
	            *bX1, *bY1, *bZ1, *bX2, *bY2, *bZ2, *bX3, *bY3, *bZ3, *bX4, *bY4, *bZ4);

	return verbose_flag;

}
