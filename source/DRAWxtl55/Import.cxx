// $Id: Import.cxx 1107 2011-01-19 23:53:52Z martin $
//
// Import.cxx - Source module for DRAWxtl V5.5
// Coded using the FLTK 1.1.7 widget set
//
//     Larry W. Finger, Martin Kroeker and Brian Toby
//
//
// This module contains the following routines:
//
//  import_cif - get structural information from a CIF (Crystallographic Information File)
//  import_fdat - get structural information from a CSD- or FDAT-format file
//  import_gsas - get structural information from a GSAS-format file
//  import_schakal - get structural information from a SCHAKAL-format file
//  import_shelx - get structural information from a SHELX-format file
//  import_wien - get structural information from a WIEN2k struct file
//  import_discus - get structural information from a DISCUS file
//  import_pcr - get structural information from a FULLPROF-format file
//  position_cif - position the CIF import file at a specific token
//  dissect_symbol - split space group name for parsing by symop()

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#if defined(WIN32)
#define snprintf _snprintf
#endif
#include <ctype.h>
#include "drawxtl.h"
#include "draw_ext.h"
#include "DRAWxtlViewUI.h"

#include "DRAWxtl_proto.h"

extern int g_Quick;

int position_cif (int numblocks, long startpos, FILE * impin, const char *search_string, char *string);

char *dissect_symbol (const char *line);

void
import_cif (char input[], int in_line, int Quick, int *Block, int frame_no)

/* routine to extract structural information from a CIF */
{
    FILE *impin;

    char string[256], filename[256], tstring[256];

    int U_pos[6], X_pos[3], Id_pos = 0, O_pos = 0, Label_pos = 0;

    int i, j, sg;

    int n_items, items;

    int numblocks;

    long skip_id;

    long startpos;

    int BIJ_flag;

    static char spgps[230][11]={
	"P 1","P -1","P 2","P 21","C 2","P m","P c","C m","C c","P 2/m",
	"P 21/m","C 2/m","P 2/c","P 21/c","C 2/c","P 2 2 2","P 2 2 21","P 21 21 2","P 21 21 21","C 2 2 21",
	"C 2 2 2","F 2 2 2","I 2 2 2","I 21 21 21","P m m 2","P m c 21","P c c 2","P m a 2","P c a 21","P n c 2",
	"P m n 21","P b a 2","P n a 21","P n n 2","C m m 2","C m c 21","C c c 2","A m m 2","A e m 2","A m a 2",
	"A e a 2","F m m 2","F d d 2","I m m 2","I b a 2","I m a 2","P m m m","P n n n","P c c m","P b a n",
	"P m m a","P n n a","P m n a","P c c a","P b a m","P c c n","P b c m","P n n m","P m m n","P b c n",
	"P b c a","P n m a","C m c m","C m c e","C m m m","C c c m","C m m e","C c c e","F m m m","F d d d",
	"I m m m","I b a m","I b c a","I m m a","P 4","P 41","P 42","P 43","I 4","I 41",
	"P -4","I -4","P 4/m","P 42/m","P 4/n","P 42/n","I 4/m","I 41/a","P 4 2 2","P 4 21 2",
	"P 41 2 2","P 41 21 2","P 42 2 2","P 42 21 2","P 43 2 2","P 43 21 2","I 4 2 2","I 41 2 2","P 4 m m","P 4 b m",
	"P 42 c m","P 42 n m","P 4 c c","P 4 n c","P 42 m c","P 42 b c","I 4 m m","I 4 c m","I 41 m d","I 41 c d",
	"P -4 2 m","P -4 2 c","P -4 21 m","P -4 21 c","P -4 m 2","P -4 c 2","P -4 b 2","P -4 n 2","I -4 m 2","I -4 c 2",
	"I -4 2 m","I -4 2 d","P 4/m m m","P 4/m c c","P 4/n b m","P 4/n n c","P 4/m b m","P 4/m n c","P 4/n m m","P 4/n c c",
	"P 42/m m c","P 42/m c m","P 42/n b c","P 42/n n m","P 42/m b c","P 42/m n m","P 42/n m c","P 42/n c m","I 4/m m m","I 4/m c m",
	"I 41/a m d","I 41/a c d","P 3","P 31","P 32","R 3","P -3","R -3","P 3 1 2","P 3 2 1",
	"P 31 1 2","P 31 2 1","P 32 1 2","P 32 2 1","R 3 2","P 3 m 1","P 3 1 m","P 3 c 1","P 3 1 c","R 3 m",
	"R 3 c","P -3 1 m","P -3 1 c","P -3 m 1","P -3 c 1","R -3 m","R -3 c","P 6","P 61","P 65",
	"P 62","P 64","P 63","P -6","P 6/m","P 63/m","P 6 2 2","P 61 2 2","P 65 2 2","P 62 2 2",
	"P 64 2 2","P 63 2 2","P 6 m m","P 6 c c","P 63 c m","P 63 m c","P -6 m 2","P -6 c 2","P -6 2 m","P -6 2 c",
	"P 6/m m m","P 6/m m c","P 63/m c m","P 63/m m c","P 2 3","F 2 3","I 2 3","P 21 3","I 21 3","P m -3",
	"P n -3","F m -3","F d -3","I m -3","P a -3","I a -3","P 4 3 2","P 42 3 2","F 4 3 2","F 41 3 2",
	"I 4 3 2","P 43 3 2","P 41 3 2","I 41 3 2","P -4 3 m","F -4 3 m","I -4 3 m","P -4 3 n","F -4 3 c","I -4 3 d",
	"P m -3 m","P n -3 n","P m -3 n","P n -3 m","F m -3 m","F m -3 c","F d -3 m","F d -3 c","I m -3 m","I a -3 d"};


    g_Quick = Quick;
    memset (U_pos, 0, sizeof (U_pos));
    memset (X_pos, 0, sizeof (X_pos));
    if (in_line) {
	impin = drvui->fpin;
	startpos = ftell(impin);
//	Error_Box ("Cannot read inline CIF data. Please use the import instruction");
//	return;
    } else {
    strcpy (filename, input);
    startpos = 0;
    if (!(impin = fopen (filename, "r"))) {
	if (!Quick)
	    fprintf (drvui->flout, "Unable to open CIF file %s\n", filename);
	Error_Box ("Cannot Open CIF Import File, Run aborted.");
	return;
    }
    }
    memset (string, 0, 255);
    memset (tstring, 0, 255);

    numblocks = *Block;

    if (!numblocks) {
	char dataname[2040] = "";

	for (;;) {
	    if (!get_next_token (string, 256, impin))
		break;		/* out on EOF */
	    if (!strncmp (string, "data_", 5)) {
		char tmp[40], tmp2[10];

		sscanf (string, "%s", tmp);
		sprintf (tmp2, " %d. ", ++numblocks);
		strcat (dataname, tmp2);
		strcat (dataname, tmp);
		strcat (dataname, "\n");
	    }
	}
	if (numblocks > 1) {
	    sprintf (tstring,
		     "This CIF has %d data blocks with labels\n%s.\nPlease enter the number"
		     " of the one to use:", numblocks, dataname);
	    const char *which = fl_input ("%s", "1", tstring);

	    if (!which) {
		if (!in_line) fclose (impin);
		return;		// Cancel button was pressed 
	    }
	    strcpy (string, which);
	    (void) sscanf (string, "%d", &numblocks);
	    *Block = numblocks;	/* return to caller */
	} else {
	    *Block = 1;
	}
    }
    if (position_cif (numblocks, startpos, impin, "_cell_length_a", string)) {
	Error_Box ("Error finding cell (a) in CIF Import File, Run aborted.");
	if (!in_line) fclose (impin);
	return;
    }
    get_next_token (string, 256, impin);
    (void) sscanf (string, "%f", &drvui->lat_con[0]);
    if (position_cif (numblocks, startpos, impin, "_cell_length_b", string)) {
	Error_Box ("Error finding cell (b) in CIF Import File, Run aborted.");
	if (!in_line) fclose (impin);
	return;
    }
    get_next_token (string, 256, impin);
    (void) sscanf (string, "%f", &drvui->lat_con[1]);
    if (position_cif (numblocks, startpos, impin, "_cell_length_c", string)) {
	Error_Box ("Error finding cell (c) in CIF Import File, Run aborted.");
	if (!in_line) fclose (impin);
	return;
    }
    get_next_token (string, 256, impin);
    (void) sscanf (string, "%f", &drvui->lat_con[2]);
    if (position_cif (numblocks, startpos, impin, "_cell_angle_al", string)) {
	Error_Box ("Error finding cell (alpha) in CIF Import File, Run aborted.");
	if (!in_line) fclose (impin);
	return;
    }
    get_next_token (string, 256, impin);
    (void) sscanf (string, "%f", &drvui->lat_con[3]);
    if (position_cif (numblocks, startpos, impin, "_cell_angle_be", string)) {
	Error_Box ("Error finding cell (beta) in CIF Import File, Run aborted.");
	if (!in_line) fclose (impin);
	return;
    }
    get_next_token (string, 256, impin);
    (void) sscanf (string, "%f", &drvui->lat_con[4]);
    if (position_cif (numblocks, startpos, impin, "_cell_angle_ga", string)) {
	Error_Box ("Error finding cell (gamma) in CIF Import File, Run aborted.");
	if (!in_line) fclose (impin);
	return;
    }
    get_next_token (string, 256, impin);
    (void) sscanf (string, "%f", &drvui->lat_con[5]);
    if (!Quick) {
	fprintf (drvui->flout,
		 "******      Cell     %8.5f %8.5f %8.5f %8.3f %8.3f %8.3f\n",
		 drvui->lat_con[0], drvui->lat_con[1], drvui->lat_con[2],
		 drvui->lat_con[3], drvui->lat_con[4], drvui->lat_con[5]);
	if (drvui->lat_con[0] == drvui->lat_con[1]
	    && drvui->lat_con[1] == drvui->lat_con[2]
	    && drvui->lat_con[3] == drvui->lat_con[4]
	    && drvui->lat_con[4] == drvui->lat_con[5]
	    && drvui->lat_con[3] != 0.0 && fabs (drvui->lat_con[3] - 90.0f) > 0.01f) {
	    fprintf (drvui->flout,
		     "***** Warning - rhombohedral settings are not supported.\n");
	    fprintf (drvui->flout,
		     "***** Please transform to standard hexagonal setting.\n");
	    Error_Box ("Rhombohedral cells are not supported\n"
		       "Please transform to standard hexagonal setting first\n");
	}
    }
    i = -1;
    if (position_cif (numblocks, startpos, impin, "_space_group_ssg_name", string)) {
	i = 0;			/* indicate not found */
	if (!position_cif (numblocks, startpos, impin, "_space_group_symop_ssg_id", string))
	    i = -1;
    }
    j = 1;
    sg = 0;

    if (i == 0) {
	if (position_cif (numblocks, startpos, impin, "_symmetry_space_group_name_H-M", string))
  	    if (position_cif (numblocks, startpos, impin, "_symmetry_Int_Tables_number", string))
		if (position_cif (numblocks, startpos, impin, "_space_group_IT_number", string))
		    j = 0;		/* indicate not found */
    }

    if (i == -1 && drvui->modulated != -1) {
	if (!Quick) {
	    fprintf (drvui->flout,
		     "******      this CIF describes a modulated structure -- \n");
	    fprintf (drvui->flout,
		     "******      use the 'average' keyword to display average structure only\n");
	}
	drvui->modulated = 1;
	j = 1;
    }

    if (j != 0) {
	fgets (string, 255, impin);
	for (i = 0; i < 60; i++)
	    if (string[i] == (char) 34)
		string[i] = (char) 39;
	tstring[0] = 's';
	tstring[1] = 'p';
	tstring[2] = 'g';
	tstring[3] = 'r';
	tstring[4] = ' ';
        j=atoi(string);
	if (j>0&& j<231) { 
	    strcat(tstring,spgps[j-1]);
	    symop(tstring);
	    sg=1;
	} else {
	    if (string[0] != '?') {
		if ( strchr (string, (char) 39)) {
		    for (i = 0; i < 60; i++)
			if (string[i] == (char) 39)
			    break;
		    if (strchr (string, ':')) {	// WJJ notation has reference sg between colons
			for (i++; i < 60; i++)
			    if (string[i] == ':')
				break;
		    }
		    for (i++, j = 5; i < 60; i++)
			if (string[i] == (char) 39 || string[i] == (char) 40 || string[i] == ':')
			    break;
			else
			    tstring[j++] = string[i];
		    tstring[j] = '\0';
		} else {
		    Blank_Strip(string);
		    strcat(tstring,string);
		}
		for (i = 0; i < 60; i++) {
		    if (tstring[i] == 'H')
			tstring[i] = ' ';	/* remove any H from rhomb. spgrps */
		    if (tstring[i] == 'S')
			tstring[i] = ' ';	/* remove any S from monoclinic spgrps */
		}
		if (tstring[6] != ' ') {	/* no space after lattice symbol, need to dissect line first */
		    char *spstring = dissect_symbol (tstring + 5);
		    if (!Quick)
			fprintf (drvui->flout, "******  %s\n", spstring);
		    symop (spstring);
		    free  (spstring);
		} else {
		    if (!Quick)
			fprintf (drvui->flout, "******  %s\n", tstring);
		    symop (tstring);
		}
		sg = 1;
	    } else {
		if (!Quick)
		    fprintf (drvui->flout,
			 "******      this CIF does not contain a space group name\n");
	    }
	}
    } else {
	if (!Quick)
	    fprintf (drvui->flout,
		     "******      this CIF does not contain a space group record\n");
    }

    if (sg == 0 && drvui->sys == 0) {
	j = 1;
	if (position_cif (numblocks, startpos, impin, "_symmetry_equiv_pos_site_id", string))
	    skip_id = 0;	// no id token in the symop table
	else
	    skip_id = ftell (impin);	// save position of id token

	if (skip_id == 0) {	// try modern spelling
	    if (position_cif (numblocks, startpos, impin, "_space_group_symop_id", string))
		skip_id = 0;	// no id token in the symop table
	    else
		skip_id = ftell (impin);	// save position of id token
	}
	if (position_cif (numblocks, startpos, impin, "_symmetry_equiv_pos_as_xyz", string))
	    j = 0;		/* indicate not found, but no error */
	if (j == 0) {		/* try modern spelling */
	    if (position_cif
		(numblocks, startpos, impin, "_space_group_symop_operation_xyz", string))
		j = 0;
	}			/* indicate not found, but no error */
	if (j == 0) {
	    if (!Quick)
		fprintf (drvui->flout,
			 "******      this CIF contains neither a space group nor a list of symmetry operations\n");
	} else {
	    char *p, *pp;

	    int k;

	    drvui->acentric = 1;
	    /* initialize symmetry arrays */
	    for (i = 0; i <= 23; ++i) {
		for (j = 0; j <= 2; ++j) {
		    drvui->ts[i][j] = 0.0;
		    for (k = 0; k <= 2; ++k)
			drvui->ss[i][j][k] = 0;
		}
	    }
	    for (i = 0; i <= 2; ++i) {
		drvui->lat_pos[0][i] = 0.0;
		drvui->ss[0][i][i] = 1;
		for (j = 0; j <= 3; ++j)
		    drvui->spg[i][j] = ' ';
	    }
	    drvui->ng = 0;

/* find out our current position (i.e. that of the _pos_as_xyz token) 
   in the CIF file. if it is smaller than that of the _id token, 
   we need not strip id numbers from the symmetry operation table */

	    if (ftell (impin) < skip_id)
		skip_id = 0;

	    fgets (string, 255, impin);

	    while (strlen (string) > 2) {
		fgets (string, 255, impin);
		for (j = 0; j < (int) strlen (string); j++)
		    if (string[j] == (unsigned char) 39
			|| string[j] == (unsigned char) 34)
			string[j] = ' ';	/* single or double quote */
		if (skip_id > 0)
		    Token_Strip (string, 2);
		Blank_Strip (string);
		if (strlen (string) < 3)
		    break;
		if (strstr (string, "-x, -y, -z")) {
		    drvui->acentric = 0;
		    break;
		}
		drvui->ng++;
		if (drvui->ng > 23) {
		    Error_Box
			("Please supply a spacegroup name - either in the CIF file itself\n"
			 "or in a 'spgr' command preceding the 'import' line in the str file");
		    break;
		}
		p = strstr (string, ",");
		*p = '\0';
		getsym (string, drvui->ng - 1, 0);
		*p = ' ';
		p++;
		pp = strstr (p, ",");
		*pp = '\0';
		getsym (p, drvui->ng - 1, 1);
		pp++;
		getsym (pp, drvui->ng - 1, 2);
	    }
	    findsys ();
	    find_lattice_type ();
	}
    }

    if (drvui->modulated >= 1) {	// parse all symmetry lines to preserve superspace symmetry
	j = 1;
	if (position_cif (numblocks, startpos, impin, "_space_group_symop_ssg_id", string))
	    skip_id = 0;	// no id token in the symop table
	else
	    skip_id = ftell (impin);	// save position of id token
	if (position_cif
	    (numblocks, startpos, impin, "_space_group_symop_ssg_operation_algebraic", string))
	    j = 0;		/* indicate not found, but no error */
	if (j == 0) {
	    if (!Quick)
		fprintf (drvui->flout,
			 "******      this CIF does not contain superspace symmetry operations\n");
	} else {
	    char *p, *pp;

	    int k;

	    drvui->acentric = 1;
	    /* initialize symmetry arrays */
	    for (i = 0; i <= 23; ++i) {
		for (j = 0; j <= 2; ++j) {
		    drvui->ts[i][j] = 0.0;
		    drvui->ts_m[i][j] = 0.0;
		    for (k = 0; k <= 2; ++k) {
			drvui->ss[i][j][k] = 0;
			drvui->ss_m[i][j][k] = 0;
		    }
		}
	    }
	    for (i = 0; i <= 2; ++i) {
		drvui->lat_pos[0][i] = 0.0;
		drvui->ss[0][i][i] = 1;
		drvui->ss_m[0][i][i] = 1;
		for (j = 0; j <= 3; ++j)
		    drvui->spg[i][j] = ' ';
	    }
	    for (i = 0; i < 3; ++i) {
		drvui->cell_vec[0][i] = 0.0;
		drvui->cell_vec[1][i] = 0.0;
		drvui->cell_vec[2][i] = 0.0;
	    }
	    drvui->ng = 0;

/* find out our current position (i.e. that of the _symop_ token) 
   in the CIF file. if it is smaller than that of the _id token, 
   we need not strip id numbers from the symmetry operation table */

	    if (ftell (impin) < skip_id)
		skip_id = 0;

	    fgets (string, 255, impin);

	    while (strlen (string) > 2) {
		int kk;

		char temp_string[256];

		fgets (string, 255, impin);
		if (skip_id > 0)
		    Token_Strip (string, 2);
		for (j = 0; j < (int) strlen (string); j++)
		    if (string[j] == (unsigned char) 39
			|| string[j] == (unsigned char) 34)
			string[j] = ' ';
		Blank_Strip (string);
		if (strlen (string) < 3 || string[0] == '_')
		    break;
		if (!drvui->ng) {
		    if (strstr (string, "x5"))	/* get the order of the modulation */
			drvui->modulated = 2;
		    if (strstr (string, "x6"))
			drvui->modulated = 3;
		}
		if (drvui->modulated < 3)	/* fill in parts of ss_m not included in input line */
		    drvui->ss_m[drvui->ng][2][2] = 1;
		if (drvui->modulated < 2)
		    drvui->ss_m[drvui->ng][1][1] = 1;
		drvui->ng++;
		strcpy (temp_string, string);
		strcat (temp_string, ",");
		p = temp_string;
		for (kk = 0; kk < 6; kk++) {
		    pp = strstr (p, ",");
		    *pp = '\0';
		    getsym (p, drvui->ng - 1, kk);
		    p = ++pp;
		    if (strlen (p) < 2)
			kk = 6;
		}
	    }
	    findsys ();
	    find_lattice_type ();
	}
    }

/* skip through file to atom_site */
    if (position_cif (numblocks, startpos, impin, "_atom_site_", string)) {
	Error_Box ("Error finding _atom_site_ in CIF Import File, Run aborted.");
	if (!in_line) fclose (impin);
	return;
    }

    /* first _atom_site_ line in string need to find position of label, x, y, and z */
    j = 0;
    while (!strncmp (string, "_atom_site_", 11)) {
	if (!strncmp (string, "_atom_site_label", 16))
	    Label_pos = j;
	if (!strncmp (string, "_atom_site_fract_x", 18))
	    X_pos[0] = j;
	if (!strncmp (string, "_atom_site_fract_y", 18))
	    X_pos[1] = j;
	if (!strncmp (string, "_atom_site_fract_z", 18))
	    X_pos[2] = j;
	if (!strncmp (string, "_atom_site_occupancy", 20))
	    O_pos = j;
	j++;
	if (!get_next_token (string, 256, impin)) {	/* search for next item */
	    Error_Box ("Error reading CIF Import File, Run aborted.");
	    if (!in_line) (void) fclose (impin);
	    return;
	}
    }
    n_items = j;
    for (;;) {
	if (!strncmp (string, "loop_", 5) || !strncmp (string, "data_", 5)
	    || !strncmp (string, "_", 1))
	    break;
	for (items = 0; items < n_items; items++) {
	    if (items == Label_pos) {
		j = 0;
		drvui->atoms[natom].atom_n = 0;
		strcpy (drvui->atoms[natom].atom_l, "    ");
		for (i = 0; i < (int) strlen (string); i++) {
		    if (string[i] >= '0' && string[i] <= '9') {
			drvui->atoms[natom].atom_n =
			    10 * drvui->atoms[natom].atom_n + (int) string[i] - 48;
		    } else {
			drvui->atoms[natom].atom_l[j++] = string[i];
			if (j > 3)
			    j = 3;
		    }
		}
	    } else if (items == X_pos[0]) {
		(void) sscanf (string, "%f", &drvui->atoms[natom].atom_xyz[0]);
	    } else if (items == X_pos[1]) {
		(void) sscanf (string, "%f", &drvui->atoms[natom].atom_xyz[1]);
	    } else if (items == X_pos[2]) {
		(void) sscanf (string, "%f", &drvui->atoms[natom].atom_xyz[2]);
	    } else if (items == O_pos) {
		(void) sscanf (string, "%f", &drvui->atoms[natom].occupancy);
	    }
	    if (!get_next_token (string, 256, impin) || !strlen (string)) {
//        Error_Box("End of File while reading atoms in CIF Import File. loop_ added.");
		strcpy (string, "loop_");
	    }
	}
	if (!Quick) {
	    fprintf (drvui->flout, "******      Atom %c%c%c%c%3d %8.5f %8.5f %8.5f\n",
		     drvui->atoms[natom].atom_l[0], drvui->atoms[natom].atom_l[1],
		     drvui->atoms[natom].atom_l[2], drvui->atoms[natom].atom_l[3],
		     drvui->atoms[natom].atom_n, drvui->atoms[natom].atom_xyz[0],
		     drvui->atoms[natom].atom_xyz[1], drvui->atoms[natom].atom_xyz[2]);
	    fflush (drvui->flout);
	}
	drvui->atoms[natom].atom_fn = frame_no;
	drvui->atoms[natom].sv_atom_n = drvui->atoms[natom].atom_n;
	drvui->atoms[natom].atom_ismod = 0;
	drvui->atoms[natom].occ_ismod = 0;
	drvui->atoms[natom++].min_occ = 0.;
	check_dynamic_storage ();
    }
    j = 1;
/* skip through file to atom_site_aniso */
    if (position_cif (numblocks, startpos, impin, "_atom_site_aniso", string)) {
	j = 0;			// indicate not found, but no error
    }

    if (j == 1) {		// only do if there are anisotropic temperature factors
/* first _atom_site_aniso line in string
   need to find position of label, U11, etc */
	BIJ_flag = 0;
	j = 0;
	while (!strncmp (string, "_atom_site_aniso", 16)) {
	    if (!strncmp (string, "_atom_site_aniso_label", 22))
		Label_pos = j;
	    if (!strncmp (string, "_atom_site_aniso_U_11", 21))
		U_pos[0] = j;
	    if (!strncmp (string, "_atom_site_aniso_U_22", 21))
		U_pos[1] = j;
	    if (!strncmp (string, "_atom_site_aniso_U_33", 21))
		U_pos[2] = j;
	    if (!strncmp (string, "_atom_site_aniso_U_12", 21))
		U_pos[3] = j;
	    if (!strncmp (string, "_atom_site_aniso_U_13", 21))
		U_pos[4] = j;
	    if (!strncmp (string, "_atom_site_aniso_U_23", 21))
		U_pos[5] = j;
	    if (!strncmp (string, "_atom_site_aniso_B_11", 21)) {
		U_pos[0] = j;
		BIJ_flag = 1;
	    }
	    if (!strncmp (string, "_atom_site_aniso_B_22", 21))
		U_pos[1] = j;
	    if (!strncmp (string, "_atom_site_aniso_B_33", 21))
		U_pos[2] = j;
	    if (!strncmp (string, "_atom_site_aniso_B_12", 21))
		U_pos[3] = j;
	    if (!strncmp (string, "_atom_site_aniso_B_13", 21))
		U_pos[4] = j;
	    if (!strncmp (string, "_atom_site_aniso_B_23", 21))
		U_pos[5] = j;
	    j++;
	    if (!get_next_token (string, 256, impin)) {	/* search for next item */
//        Error_Box("End of File while reading Uij in CIF Import File, Run aborted.");
		if (!in_line) (void) fclose (impin);
		return;
	    }
	}
	n_items = j;
	for (;;) {
	    if (!strncmp (string, "_", 1) || strstr (string, "loop_"))
		break;
	    for (items = 0; items <= n_items - 1; items++) {
		if (Label_pos == items) {
		    drvui->ellips[drvui->n_ellips].ellips_n = 0;
		    strcpy (drvui->ellips[drvui->n_ellips].ellips_l, "    ");	/* initialize ellipsoid name */
		    j = 0;
		    for (i = 0; i < (int) strlen (string); i++) {
			if (string[i] >= '0' && string[i] <= '9') {
			    drvui->ellips[drvui->n_ellips].ellips_n =
				10 * drvui->ellips[drvui->n_ellips].ellips_n +
				(int) string[i] - 48;
			} else {
			    drvui->ellips[drvui->n_ellips].ellips_l[j++] = string[i];
			    if (j > 3)
				j = 3;
			}
		    }
		}
		drvui->ellips[drvui->n_ellips].save_el_number =
		    drvui->ellips[drvui->n_ellips].ellips_n;
		drvui->do_ellipsoids = 1;
		if (U_pos[0] == items)
		    (void) sscanf (string, "%f",
				   &drvui->ellips[drvui->n_ellips].ellips[0]);
		if (U_pos[1] == items)
		    (void) sscanf (string, "%f",
				   &drvui->ellips[drvui->n_ellips].ellips[1]);
		if (U_pos[2] == items)
		    (void) sscanf (string, "%f",
				   &drvui->ellips[drvui->n_ellips].ellips[2]);
		if (U_pos[3] == items)
		    (void) sscanf (string, "%f",
				   &drvui->ellips[drvui->n_ellips].ellips[3]);
		if (U_pos[4] == items)
		    (void) sscanf (string, "%f",
				   &drvui->ellips[drvui->n_ellips].ellips[4]);
		if (U_pos[5] == items)
		    (void) sscanf (string, "%f",
				   &drvui->ellips[drvui->n_ellips].ellips[5]);
		if (drvui->auto_ellipse == 1) {
		    if (BIJ_flag)
			drvui->ellips[drvui->n_ellips].ell_type = 1002;	// Bij with ellipcolor
		    else
			drvui->ellips[drvui->n_ellips].ell_type = 1001;	// Uij with ellipcolor
		    memset (drvui->ellips[drvui->n_ellips].ellips_col, 0, 40);
		    strcpy (drvui->ellips[drvui->n_ellips].ellips_col, "Gray30");
		} else {
		    if (BIJ_flag)
			drvui->ellips[drvui->n_ellips].ell_type = 2;	// Bij
		    else
			drvui->ellips[drvui->n_ellips].ell_type = 1;	// Uij
		}
		if (items < n_items - 1)
		    get_next_token (string, 255, impin);
	    }
	    if (!Quick)
		fprintf (drvui->flout,
			 "******      Uij  %c%c%c%c%3d %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
			 drvui->ellips[drvui->n_ellips].ellips_l[0],
			 drvui->ellips[drvui->n_ellips].ellips_l[1],
			 drvui->ellips[drvui->n_ellips].ellips_l[2],
			 drvui->ellips[drvui->n_ellips].ellips_l[3],
			 drvui->ellips[drvui->n_ellips].ellips_n,
			 drvui->ellips[drvui->n_ellips].ellips[0],
			 drvui->ellips[drvui->n_ellips].ellips[1],
			 drvui->ellips[drvui->n_ellips].ellips[2],
			 drvui->ellips[drvui->n_ellips].ellips[3],
			 drvui->ellips[drvui->n_ellips].ellips[4],
			 drvui->ellips[drvui->n_ellips].ellips[5]);
	    drvui->ellips[drvui->n_ellips++].ellips_ismod = 0;
	    check_dynamic_storage ();
	    if (!get_next_token (string, 256, impin)) {
//        Error_Box("Error reading CIF Import File, Run aborted.");
		if (!in_line) (void) fclose (impin);
		return;
	    }
	}
	drvui->auto_ellipse = 0;
    }

    if (drvui->modulated >= 1) {
	int axis = 0;

	int elem = 0;

	int theatom = 0;

	drvui->no_mod_vectors = 0;

/* unit cell wave vectors */

/* skip through file to cell_wave_vector_seq_id line */

	if (position_cif (numblocks, startpos, impin, "_cell_wave_vector_seq_id", string)) {
	    Error_Box
		("Error finding _cell_wave_vector in CIF Import File, Run aborted.");
	    if (!in_line) fclose (impin);
	    return;
	}
	j = 0;
	if (!Quick)
	    fprintf (drvui->flout,
		     "******\n******  Cell Modulation Vector(s)\n****** No.             "
		     "Components\n******\n");
	while (!strncmp (string, "_cell_wave_vector", 17)) {
	    if (!strncmp (string, "_cell_wave_vector_x", 19))
		X_pos[0] = j;
	    if (!strncmp (string, "_cell_wave_vector_y", 19))
		X_pos[1] = j;
	    if (!strncmp (string, "_cell_wave_vector_z", 19))
		X_pos[2] = j;
	    j++;
	    if (!get_next_token (string, 256, impin)) {	/* search for next item */
		Error_Box ("Error reading CIF Import File, Run aborted.");
		if (!in_line) (void) fclose (impin);
		return;
	    }
	}
	n_items = j;
	drvui->no_cell_vec = 0;
	for (;;) {
	    if (!strncmp (string, "loop_", 5) || !strncmp (string, "data_", 5)
		|| !strncmp (string, "_", 1))
		break;
	    for (items = 0; items < n_items; items++) {
		if (items == X_pos[0]) {
		    (void) sscanf (string, "%f", &drvui->cell_vec[drvui->no_cell_vec][0]);
		} else if (items == X_pos[1]) {
		    (void) sscanf (string, "%f", &drvui->cell_vec[drvui->no_cell_vec][1]);
		} else if (items == X_pos[2]) {
		    (void) sscanf (string, "%f", &drvui->cell_vec[drvui->no_cell_vec][2]);
		}
		if (!get_next_token (string, 256, impin) || !strlen (string)) {
		    strcpy (string, "loop_");
		}
	    }
	    if (!Quick) {
		fprintf (drvui->flout,
			 "****** %2d     %8.5f %8.5f %8.5f\n", drvui->no_cell_vec + 1,
			 drvui->cell_vec[drvui->no_cell_vec][0],
			 drvui->cell_vec[drvui->no_cell_vec][1],
			 drvui->cell_vec[drvui->no_cell_vec][2]);
		fflush (drvui->flout);
	    }
	    drvui->no_cell_vec++;
	}
	if (drvui->no_cell_vec && !Quick)
	    fprintf (drvui->flout, "******\n");

	/* modulation wave vector(s) */

	/* skip through file to atom_site_Fourier */
	if (position_cif (numblocks, startpos, impin, "_atom_site_Fourier_", string)) {
	    Error_Box
		("Error reading modulation wave vector in CIF Import File, Run aborted.");
	    if (!in_line) fclose (impin);
	    return;
	}
	/* first _atom_site_ line in string need to find position of label, x, y, and z */
	j = 0;
	X_pos[0] = X_pos[1] = X_pos[2] = Label_pos = -1;
	if (!Quick)
	    fprintf (drvui->flout,
		     "******\n******  Atom Site Modulation Vector(s)\n******  No.     Components         Cell Wave Components\n******\n");
	while (!strncmp (string, "_atom_site_Fourier_", 19)) {
	    if (!strncmp (string, "_atom_site_Fourier_wave_vector_seq_id", 37))
		Label_pos = j;
	    if (!strncmp (string, "_atom_site_Fourier_wave_vector_x", 32))
		X_pos[0] = j;
	    if (!strncmp (string, "_atom_site_Fourier_wave_vector_y", 32))
		X_pos[1] = j;
	    if (!strncmp (string, "_atom_site_Fourier_wave_vector_z", 32))
		X_pos[2] = j;
	    j++;
	    if (!get_next_token (string, 256, impin)) {	/* search for next item */
		Error_Box ("Error reading CIF Import File, Run aborted.");
		if (!in_line) (void) fclose (impin);
		return;
	    }
	}
	n_items = j;
	for (;;) {
	    int no_z, max_z = 5;

	    int no_y, max_y = 5;

	    int no_x, max_x = 5;

	    if (drvui->no_cell_vec < 3)
		max_z = 0;
	    if (drvui->no_cell_vec < 2)
		max_y = 0;
	    if (!strncmp (string, "loop_", 5) || !strncmp (string, "data_", 5))
		break;
	    j = drvui->no_mod_vectors;
	    drvui->modulate_gbl[j].modvector[0] = 0.;
	    drvui->modulate_gbl[j].modvector[1] = 0.;
	    drvui->modulate_gbl[j].modvector[2] = 0.;
	    for (items = 0; items <= n_items - 1; items++) {
		if (X_pos[0] == items)
		    (void) sscanf (string, "%f", &drvui->modulate_gbl[j].modvector[0]);
		if (X_pos[1] == items)
		    (void) sscanf (string, "%f", &drvui->modulate_gbl[j].modvector[1]);
		if (X_pos[2] == items)
		    (void) sscanf (string, "%f", &drvui->modulate_gbl[j].modvector[2]);
		if (items < n_items)
		    get_next_token (string, 256, impin);
	    }
	    for (no_z = -max_z; no_z <= max_z; no_z++) {
		for (no_y = -max_y; no_y <= max_y; no_y++) {
		    for (no_x = -max_x; no_x <= max_x; no_x++) {
			if (vec_dif
			    (no_x, drvui->cell_vec[0], no_y, drvui->cell_vec[1], no_z,
			     drvui->cell_vec[2], drvui->modulate_gbl[j].modvector)) {
			    drvui->modulate_gbl[j].vector_mult[0] = no_x;	/* have components of modvector */
			    drvui->modulate_gbl[j].vector_mult[1] = no_y;
			    drvui->modulate_gbl[j].vector_mult[2] = no_z;
			    goto done;
			}
		    }
		}
	    }
	    Error_Box
		("Unable to decombine modulation wave vector into cell vector parts.");
	  done:
	    if (!Quick) {
		fprintf (drvui->flout, "****** %2d %8.5f %8.5f %8.5f  %3d %3d %3d\n",
			 j + 1, drvui->modulate_gbl[j].modvector[0],
			 drvui->modulate_gbl[j].modvector[1],
			 drvui->modulate_gbl[j].modvector[2],
			 drvui->modulate_gbl[j].vector_mult[0],
			 drvui->modulate_gbl[j].vector_mult[1],
			 drvui->modulate_gbl[j].vector_mult[2]);
	    }
	    drvui->no_mod_vectors++;
	    check_dynamic_storage ();
	}

	/* displacement modulation - fourier series */

	/* skip through file to atom_site_displace_Fourier */
	(void) position_cif (numblocks, startpos, impin, "_atom_site_displace_Fourier_", string);
	/* first _atom_site_ line in string need to find position of label and parameters */
	j = 0;
	while (!strncmp (string, "_atom_site_displace_Fourier_", 28)) {
	    if (!strncmp (string, "_atom_site_displace_Fourier_atom_site_label", 43))
		Label_pos = j;
	    if (!strncmp (string, "_atom_site_displace_Fourier_wave_vector_seq_id", 46))
		Id_pos = j;
	    if (!strncmp (string, "_atom_site_displace_Fourier_axis", 32))
		X_pos[0] = j;
	    if (!strncmp (string, "_atom_site_displace_Fourier_param_cos", 37))
		X_pos[1] = j;
	    if (!strncmp (string, "_atom_site_displace_Fourier_param_sin", 37))
		X_pos[2] = j;
	    j++;
	    if (!get_next_token (string, 256, impin)) {	/* search for next item */
		Error_Box ("Error reading CIF Import File, Run aborted.");
		if (!in_line) (void) fclose (impin);	// we saw the header, so there should be data
		return;
	    }
	}
	if (j > 0) {
	    if (!Quick)
		fprintf (drvui->flout,
			 "******\n****** Atom Site Fourier Displacement Items\n****** \n"
			 "****** Atom  axis ID      cos      sin\n******\n");
	    n_items = j;
	    drvui->no_site_displace = 0;
	    for (;;) {
		if (!strncmp (string, "_", 1) || strstr (string, "loop_")
		    || !strlen (string))
		    break;
		for (items = 0; items <= n_items - 1; items++) {
		    if (Label_pos == items) {
			int modnum = 0;

			char modl[5];

			strcpy (modl, "    ");	/* initialize name */
			j = 0;
			for (i = 0; i < (int) strlen (string); i++) {
			    if (string[i] >= '0' && string[i] <= '9') {
				modnum = 10 * modnum + (int) string[i] - 48;
			    } else {
				modl[j++] = string[i];
				if (j > 3)
				    j = 3;
			    }
			}
			for (j = 0; j < natom; j++)
			    if (check_atom_name (modl, drvui->atoms[j].atom_l)
				&& modnum == drvui->atoms[j].atom_n) {
				theatom = j;
				drvui->atoms[j].atom_ismod |= 1;
				drvui->atoms[j].occ_ismod = 0;
				break;
			    }
		    }
		    if (X_pos[0] == items) {
			axis = 0;
			if (!strncmp (string, "y", 1))
			    axis = 1;
			else if (!strncmp (string, "z", 1))
			    axis = 2;
		    }
		    drvui->modulate_3x[drvui->no_site_displace].atom_modpar_axis = axis;
		    drvui->modulate_3x[drvui->no_site_displace].atom_modpar_atom =
			theatom;
		    if (Id_pos == items)
			(void) sscanf (string, "%d",
				       &drvui->modulate_3x[drvui->no_site_displace].
				       atom_modpar_id);
		    if (X_pos[1] == items)
			(void) sscanf (string, "%f",
				       &drvui->modulate_3x[drvui->no_site_displace].
				       atom_modpar[0]);
		    if (X_pos[2] == items)
			(void) sscanf (string, "%f",
				       &drvui->modulate_3x[drvui->no_site_displace].
				       atom_modpar[1]);

		    if (items < n_items)
			if (!get_next_token (string, 256, impin))
			    strcpy (string, "_loop");

		}
		if (!Quick) {
		    char caxis[3][2] = { "X", "Y", "Z" };
		    char axis2[2];

		    strcpy (axis2,
			    caxis[drvui->modulate_3x[drvui->no_site_displace].
				  atom_modpar_axis]);
		    fprintf (drvui->flout, "****** %4s%2d   %s%3d   %8.5f %8.5f\n",
			     drvui->atoms[theatom].atom_l, drvui->atoms[theatom].atom_n,
			     axis2,
			     drvui->modulate_3x[drvui->no_site_displace].atom_modpar_id,
			     drvui->modulate_3x[drvui->no_site_displace].atom_modpar[0],
			     drvui->modulate_3x[drvui->no_site_displace].atom_modpar[1]);
		}
		drvui->no_site_displace++;
		check_dynamic_storage ();
	    }
	}

	/* displacement modulation - sawtooth function */

	/* skip through file to atom_site_displace_special */
	(void) position_cif (numblocks, startpos, impin, "_atom_site_displace_special_", string);
	/* first _atom_site_ line in string need to find position of label and parameters */
	j = 0;
	Label_pos = U_pos[0] = U_pos[1] = U_pos[2] = U_pos[3] = U_pos[4] = 0;
	while (!strncmp (string, "_atom_site_displace_special_", 28)) {
	    if (!strncmp (string, "_atom_site_displace_special_func_atom_site_label", 48))
		Label_pos = j;
	    if (!strncmp (string, "_atom_site_displace_special_func_sawtooth_ax", 44))
		U_pos[0] = j;
	    if (!strncmp (string, "_atom_site_displace_special_func_sawtooth_ay", 44))
		U_pos[1] = j;
	    if (!strncmp (string, "_atom_site_displace_special_func_sawtooth_az", 44))
		U_pos[2] = j;
	    if (!strncmp (string, "_atom_site_displace_special_func_sawtooth_c", 43))
		U_pos[3] = j;
	    if (!strncmp (string, "_atom_site_displace_special_func_sawtooth_w", 43))
		U_pos[4] = j;
	    j++;
	    if (!get_next_token (string, 256, impin)) {	/* search for next item */
		Error_Box ("Error reading CIF Import File, Run aborted.");
		if (!in_line) (void) fclose (impin);	// we saw the header, so there should be data
		return;
	    }
	}
	if (j > 0) {
	    if (!Quick)
		fprintf (drvui->flout,
			 "******\n****** Atom Site Sawtooth Displacement Items\n****** \n"
			 "****** Atom         ax       ay       az     center   width\n******\n");
	    n_items = j;
	    for (;;) {
		if (!strncmp (string, "_", 1) || strstr (string, "loop_")
		    || !strlen (string))
		    break;
		for (items = 0; items <= n_items - 1; items++) {
		    if (Label_pos == items) {
			int modnum = 0;

			char modl[5];

			strcpy (modl, "    ");	/* initialize name */
			j = 0;
			for (i = 0; i < (int) strlen (string); i++) {
			    if (string[i] >= '0' && string[i] <= '9') {
				modnum = 10 * modnum + (int) string[i] - 48;
			    } else {
				modl[j++] = string[i];
				if (j > 3)
				    j = 3;
			    }
			}
			for (j = 0; j < natom; j++)
			    if (check_atom_name (modl, drvui->atoms[j].atom_l)
				&& modnum == drvui->atoms[j].atom_n) {
				theatom = j;
				drvui->atoms[j].atom_ismod |= 2;
				drvui->atoms[j].occ_ismod = 0;
				break;
			    }
		    }
		    if (U_pos[0] == items)
			(void) sscanf (string, "%f",
				       &drvui->modulate_x[theatom].atom_mod_sawtooth[0]);
		    if (U_pos[1] == items)
			(void) sscanf (string, "%f",
				       &drvui->modulate_x[theatom].atom_mod_sawtooth[1]);
		    if (U_pos[2] == items)
			(void) sscanf (string, "%f",
				       &drvui->modulate_x[theatom].atom_mod_sawtooth[2]);
		    if (U_pos[3] == items)
			(void) sscanf (string, "%f",
				       &drvui->modulate_x[theatom].atom_mod_sawtooth[3]);
		    if (U_pos[4] == items)
			(void) sscanf (string, "%f",
				       &drvui->modulate_x[theatom].atom_mod_sawtooth[4]);

		    if (items < n_items)
			if (!get_next_token (string, 256, impin))
			    strcpy (string, "_loop");

		}
		if (!Quick) {
		    fprintf (drvui->flout,
			     "****** %4s%2d   %8.5f %8.5f %8.5f %8.5f %8.5f\n",
			     drvui->atoms[theatom].atom_l, drvui->atoms[theatom].atom_n,
			     drvui->modulate_x[theatom].atom_mod_sawtooth[0],
			     drvui->modulate_x[theatom].atom_mod_sawtooth[1],
			     drvui->modulate_x[theatom].atom_mod_sawtooth[2],
			     drvui->modulate_x[theatom].atom_mod_sawtooth[3],
			     drvui->modulate_x[theatom].atom_mod_sawtooth[4]);
		}
	    }
	}

	/* occupancy modulation - fourier series */

	/* skip through file to atom_site_occ_Fourier */
	(void) position_cif (numblocks, startpos, impin, "_atom_site_occ_Fourier_", string);
	/* first _atom_site_ line in string need to find position of label, x, y, and z */
	j = 0;
	while (!strncmp (string, "_atom_site_occ_Fourier_", 23)) {
	    if (!strncmp (string, "_atom_site_occ_Fourier_atom_site_label", 38))
		Label_pos = j;
	    if (!strncmp (string, "_atom_site_occ_Fourier_wave_vector_seq_id", 41))
		Id_pos = j;
	    if (!strncmp (string, "_atom_site_occ_Fourier_param_cos", 32))
		X_pos[1] = j;
	    if (!strncmp (string, "_atom_site_occ_Fourier_param_sin", 32))
		X_pos[2] = j;
	    j++;
	    if (!get_next_token (string, 256, impin)) {	/* search for next item */
		Error_Box ("Error reading CIF Import File, Run aborted.");
		if (!in_line) (void) fclose (impin);
		return;
	    }
	}
	if (j > 0) {
	    if (!Quick)
		fprintf (drvui->flout,
			 "******\n****** Atom Site Fourier Occupancy Items\n******\n"
			 "******  Atom  axis    cos      sin\n******\n");
	    n_items = j;
	    drvui->no_site_occ = 0;
	    for (;;) {
		if (!strncmp (string, "_", 1) || strstr (string, "loop_")
		    || !strlen (string))
		    break;
		for (items = 0; items <= n_items - 1; items++) {
		    if (Label_pos == items) {
			int modnum = 0;

			char modl[5];

			strcpy (modl, "    ");	/* initialize name */
			j = 0;
			for (i = 0; i < (int) strlen (string); i++) {
			    if (string[i] >= '0' && string[i] <= '9') {
				modnum = 10 * modnum + (int) string[i] - 48;
			    } else {
				modl[j++] = string[i];
				if (j > 3)
				    j = 3;
			    }
			}
			for (j = 0; j < natom; j++)
			    if (check_atom_name (modl, drvui->atoms[j].atom_l)
				&& modnum == drvui->atoms[j].atom_n) {
				theatom = j;
				drvui->atoms[j].occ_ismod |= 1;
				break;
			    }
		    }
		    drvui->modulate_x[drvui->no_site_occ].atom_occpar_atom = theatom;

		    if (Id_pos == items)
			(void) sscanf (string, "%d",
				       &drvui->modulate_x[drvui->no_site_occ].
				       atom_occpar_id);

		    if (X_pos[1] == items)
			(void) sscanf (string, "%f",
				       &drvui->modulate_x[drvui->no_site_occ].
				       atom_occpar[0]);

		    if (X_pos[2] == items)
			(void) sscanf (string, "%f",
				       &drvui->modulate_x[drvui->no_site_occ].
				       atom_occpar[1]);

		    if (items < n_items)
			if (!get_next_token (string, 256, impin))
			    strcpy (string, "_loop");

		}
		if (!Quick) {
		    fprintf (drvui->flout, "****** %4s%d %5d %8.5f %8.5f\n",
			     drvui->atoms[theatom].atom_l, drvui->atoms[theatom].atom_n,
			     drvui->modulate_x[drvui->no_site_occ].atom_occpar_id,
			     drvui->modulate_x[drvui->no_site_occ].atom_occpar[0],
			     drvui->modulate_x[drvui->no_site_occ].atom_occpar[1]);
		}
		drvui->no_site_occ++;
		check_dynamic_storage ();
	    }
	}

	/* occupancy modulation : crenel function */

	/* skip through file to atom_site_occ_special */
	(void) position_cif (numblocks, startpos, impin, "_atom_site_occ_special_", string);
	/* first _atom_site_ line in string need to find position of label, x, y, and z */
	j = 0;
	while (!strncmp (string, "_atom_site_occ_special_", 23)) {
	    if (!strncmp (string, "_atom_site_occ_special_func_atom_site_label", 43))
		Label_pos = j;
	    if (!strncmp (string, "_atom_site_occ_special_func_crenel_c", 36))
		X_pos[1] = j;
	    if (!strncmp (string, "_atom_site_occ_special_func_crenel_w", 36))
		X_pos[2] = j;
	    j++;
	    if (!get_next_token (string, 256, impin)) {	/* search for next item */
		Error_Box ("Error reading CIF Import File, Run aborted.");
		if (!in_line) (void) fclose (impin);
		return;
	    }
	}
	if (j > 0) {
	    n_items = j;
	    for (;;) {
		if (!strncmp (string, "_", 1) || !strncmp (string, "loop_", 5)
		    || !strlen (string))
		    break;
		for (items = 0; items <= n_items - 1; items++) {
		    if (Label_pos == items) {
			int modnum = 0;

			char modl[5];

			strcpy (modl, "    ");	/* initialize name */
			j = 0;
			for (i = 0; i < (int) strlen (string); i++) {
			    if (string[i] >= '0' && string[i] <= '9') {
				modnum = 10 * modnum + (int) string[i] - 48;
			    } else {
				modl[j++] = string[i];
				if (j > 3)
				    j = 3;
			    }
			}
			for (j = 0; j < natom; j++)
			    if (check_atom_name (modl, drvui->atoms[j].atom_l)
				&& modnum == drvui->atoms[j].atom_n) {
				theatom = j;
				drvui->atoms[j].occ_ismod |= 2;
				break;
			    }
		    }

		    if (X_pos[1] == items)
			(void) sscanf (string, "%f",
				       &drvui->modulate_x[theatom].atom_occ_crenel[0]);
		    if (X_pos[2] == items)
			(void) sscanf (string, "%f",
				       &drvui->modulate_x[theatom].atom_occ_crenel[1]);
		    if (items < n_items)
			if (!get_next_token (string, 256, impin))
			    strcpy (string, "_loop");

		}
	    }
	    for (j = 0; j < natom; j++)
		if (drvui->atoms[j].occ_ismod == 2) {	// convert crenel parameters to min/max values 
		    if (!Quick) {
			fprintf (drvui->flout,
				 "*** Occupational modulation parameters for %s%d:\n",
				 drvui->atoms[j].atom_l, drvui->atoms[j].atom_n);
			fprintf (drvui->flout, "*** \t crenel offset %7.5f width %7.5f\n",
				 drvui->modulate_x[j].atom_occ_crenel[0],
				 drvui->modulate_x[j].atom_occ_crenel[1]);
		    }
		    float w = drvui->modulate_x[j].atom_occ_crenel[1] / 2.0f;

		    drvui->modulate_x[j].atom_occ_crenel[1] =
			drvui->modulate_x[j].atom_occ_crenel[0] + w;
		    drvui->modulate_x[j].atom_occ_crenel[0] -= w;
		}
	}

	/* modulation of thermal parameters - fourier series */

	/* skip through file to atom_site_U_Fourier */
	(void) position_cif (numblocks, startpos, impin, "_atom_site_U_Fourier_", string);
	/* first _atom_site_ line in string need to find position of label and parameters */
	j = 0;
	while (!strncmp (string, "_atom_site_U_Fourier_", 21)) {
	    if (!strncmp (string, "_atom_site_U_Fourier_atom_site_label", 36))
		Label_pos = j;
	    if (!strncmp (string, "_atom_site_U_Fourier_tens_elem", 30))
		X_pos[0] = j;
	    if (!strncmp (string, "_atom_site_U_Fourier_wave_vector_seq_id", 39))
		Id_pos = j;
	    if (!strncmp (string, "_atom_site_U_Fourier_param_cos", 30))
		X_pos[1] = j;
	    if (!strncmp (string, "_atom_site_U_Fourier_param_sin", 30))
		X_pos[2] = j;
	    j++;
	    if (!get_next_token (string, 256, impin)) {	/* search for next item */
		Error_Box ("Error reading CIF Import File, Run aborted.");
		if (!in_line) (void) fclose (impin);	// we saw the header, so there should be data
		return;
	    }
	}
	if (j > 0) {
	    n_items = j;
	    drvui->no_site_U_terms = 0;
	    for (;;) {
		if (!strncmp (string, "_", 1) || strstr (string, "loop_")
		    || !strlen (string))
		    break;
		for (items = 0; items <= n_items - 1; items++) {
		    if (Label_pos == items) {
			int modnum = 0;

			char modl[5];

			strcpy (modl, "    ");	/* initialize name */
			j = 0;
			for (i = 0; i < (int) strlen (string); i++) {
			    if (string[i] >= '0' && string[i] <= '9') {
				modnum = 10 * modnum + (int) string[i] - 48;
			    } else {
				modl[j++] = string[i];
				if (j > 3)
				    j = 3;
			    }
			}
			for (j = 1; j < drvui->n_ellips; j++)
			    if (check_atom_name (modl, drvui->ellips[j].ellips_l)
				&& modnum == drvui->ellips[j].ellips_n) {
				drvui->ellips[j].ellips_ismod = 1;
				break;
			    }

			for (j = 0; j < natom; j++)
			    if (check_atom_name (modl, drvui->atoms[j].atom_l)
				&& modnum == drvui->atoms[j].atom_n) {
				theatom = j;
				break;
			    }
			drvui->modulate_x[drvui->no_site_U_terms].ellips_modpar_atom =
			    theatom;
		    }
		    if (X_pos[0] == items) {
			elem = 0;	// U11
			if (!strncmp (string, "U22", 3))
			    elem = 1;
			else if (!strncmp (string, "U33", 3))
			    elem = 2;
			else if (!strncmp (string, "U12", 3))
			    elem = 3;
			else if (!strncmp (string, "U13", 3))
			    elem = 4;
			else if (!strncmp (string, "U23", 3))
			    elem = 5;
			drvui->modulate_3t[drvui->no_site_U_terms].ellips_modpar_term =
			    elem;
		    }
		    if (Id_pos == items)
			(void) sscanf (string, "%d",
				       &drvui->modulate_3t[drvui->no_site_U_terms].
				       ellips_modpar_id);
		    if (X_pos[1] == items)
			(void) sscanf (string, "%f",
				       &drvui->modulate_3t[drvui->no_site_U_terms].
				       ellips_modpar[0]);
		    if (X_pos[2] == items)
			(void) sscanf (string, "%f",
				       &drvui->modulate_3t[drvui->no_site_U_terms].
				       ellips_modpar[1]);

		    if (items < n_items)
			if (!get_next_token (string, 256, impin))
			    strcpy (string, "_loop");

		}
		drvui->no_site_U_terms++;
		check_dynamic_storage ();
	    }
	    if (!Quick) {
		char Terms[6][4] = { "U11", "U22", "U33", "U12", "U13", "U23" };
		char out_term[4];

		int id, term;

		fprintf (drvui->flout,
			 "******\n****** Uij Fourier Modulation terms\n******\n"
			 "****** Atom    ID  Term      cos      sin\n******\n");
		for (j = 0; j < drvui->no_site_U_terms; j++) {
		    theatom = drvui->modulate_x[j].ellips_modpar_atom;
		    id = drvui->modulate_3t[j].ellips_modpar_id;
		    term = drvui->modulate_3t[j].ellips_modpar_term;
		    strcpy (out_term, Terms[term]);
		    fprintf (drvui->flout, "****** %4s%d   %2d   %3s %10.5f%9.5f\n",
			     drvui->atoms[theatom].atom_l, drvui->atoms[theatom].atom_n,
			     id, out_term, drvui->modulate_3t[j].ellips_modpar[0],
			     drvui->modulate_3t[j].ellips_modpar[1]);
		}
		fprintf (drvui->flout, "******\n");
	    }
	}
    }

/* get subsystem information for composite crystals */

    make_bmat (drvui->sys, drvui->lat_con, drvui->b_mat, drvui->ginv, drvui->rec_lat_con);	/* create the lattice metric */
    if (!position_cif (numblocks, startpos, impin, "_cell_subsystems_number", string)) {	/* no error if not found */
	float e[6] = { 0, 0, 0, 0, 0, 0 };	/* the (3 + d)D reciprocal basis vectors */
	float mat[6][6];	/* the matrix describing the relations between e and ep */

	int kk;

	float ast = drvui->rec_lat_con[0];

	float bst = drvui->rec_lat_con[1];

	float cst = drvui->rec_lat_con[2];

	float csal = drvui->rec_lat_con[3];

	float csbe = drvui->rec_lat_con[4];

	float csga = drvui->rec_lat_con[5];

	float snalp, snbep, sngap, vol;

	float ap, bp, cp, csalp, csbep, csgap;

	float a, b, c, alpha, beta, gamma;

	int l, m;

	get_next_token (string, 256, impin);	/* this CIF describes a composite crystal */
	sscanf (string, "%d", &drvui->no_subsys);
	if (!Quick)
	    fprintf (drvui->flout, "******\n****** CIF describes a composite crystal "
		     "with %d subsystems\n******\n", drvui->no_subsys);
	fgets (string, 255, impin);
	for (;;) {
	    for (j = 0; j < 6; j++)
		for (kk = 0; kk < 6; kk++)
		    mat[j][kk] = 0.0f;
	    fgets (string, 255, impin);
	    if (!strstr (string, "_cell_subsystem"))
		break;		/* skip through all the W matrix lines - FIXME? */
	}
	for (j = 0; j < 3; j++)
	    e[j] = drvui->rec_lat_con[j];	/* first three elements of e are a*, b* and c* */
	for (j = 0; j < drvui->modulated; j++)	/* remainder are e4*, e5*, and e6* */
	    for (kk = 0; kk < 3; kk++)
		e[j + 3] += drvui->cell_vec[j][kk] * drvui->rec_lat_con[kk];
	for (j = 0; j < drvui->no_subsys; j++) {
	    for (i = 0; i < 3 + drvui->modulated; i++) {
		for (kk = 0; kk < 3 + drvui->modulated; kk++) {
		    get_next_token (string, 256, impin);
		    sscanf (string, "%f", &mat[i][kk]);	/* read the appropriate matrix element */
		}
	    }
	    fgets (string, 255, impin);	/* skip the line with ? */
	    for (l = 0; l < 3; l++)
		for (m = 0; m < 3; m++)
		    drvui->subsys_fact[j][l][m] =
			mat[l][m] + mat[l][3] * drvui->cell_vec[0][m]
			+ mat[l][4] * drvui->cell_vec[1][m] +
			mat[l][5] * drvui->cell_vec[2][m];
/* this is really ugly code, but I don't know how to make it pretty */
	    ap = (float) sqrt (ast * ast * drvui->subsys_fact[j][0][0] *
			       drvui->subsys_fact[j][0][0]
			       +
			       bst * bst * drvui->subsys_fact[j][0][1] *
			       drvui->subsys_fact[j][0][1]
			       +
			       cst * cst * drvui->subsys_fact[j][0][2] *
			       drvui->subsys_fact[j][0][2]
			       +
			       2.0f * ast * bst * csga * drvui->subsys_fact[j][0][0] *
			       drvui->subsys_fact[j][0][1]
			       +
			       2.0f * ast * cst * csbe * drvui->subsys_fact[j][0][0] *
			       drvui->subsys_fact[j][0][2]
			       +
			       2.0f * bst * cst * csal * drvui->subsys_fact[j][0][1] *
			       drvui->subsys_fact[j][0][2]);
	    bp = (float) sqrt (ast * ast * drvui->subsys_fact[j][1][0] *
			       drvui->subsys_fact[j][1][0]
			       +
			       bst * bst * drvui->subsys_fact[j][1][1] *
			       drvui->subsys_fact[j][1][1]
			       +
			       cst * cst * drvui->subsys_fact[j][1][2] *
			       drvui->subsys_fact[j][1][2]
			       +
			       2.0f * ast * bst * csga * drvui->subsys_fact[j][1][0] *
			       drvui->subsys_fact[j][1][1]
			       +
			       2.0f * ast * cst * csbe * drvui->subsys_fact[j][1][0] *
			       drvui->subsys_fact[j][1][2]
			       +
			       2.0f * bst * cst * csal * drvui->subsys_fact[j][1][1] *
			       drvui->subsys_fact[j][1][2]);
	    cp = (float) sqrt (ast * ast * drvui->subsys_fact[j][2][0] *
			       drvui->subsys_fact[j][2][0]
			       +
			       bst * bst * drvui->subsys_fact[j][2][1] *
			       drvui->subsys_fact[j][2][1]
			       +
			       cst * cst * drvui->subsys_fact[j][2][2] *
			       drvui->subsys_fact[j][2][2]
			       +
			       2.0f * ast * bst * csga * drvui->subsys_fact[j][2][0] *
			       drvui->subsys_fact[j][2][1]
			       +
			       2.0f * ast * cst * csbe * drvui->subsys_fact[j][2][0] *
			       drvui->subsys_fact[j][2][2]
			       +
			       2.0f * bst * cst * csal * drvui->subsys_fact[j][2][1] *
			       drvui->subsys_fact[j][2][2]);
	    csgap =
		(drvui->subsys_fact[j][0][0] * drvui->subsys_fact[j][1][0] * ap * ap +
		 drvui->subsys_fact[j][0][1] * drvui->subsys_fact[j][1][1] * bp * bp +
		 drvui->subsys_fact[j][0][2] * drvui->subsys_fact[j][1][2] * cp * cp +
		 ast * bst * csga * (drvui->subsys_fact[j][0][0] *
				     drvui->subsys_fact[j][1][1]
				     +
				     drvui->subsys_fact[j][0][1] *
				     drvui->subsys_fact[j][1][0])
		 +
		 ast * cst * csbe * (drvui->subsys_fact[j][0][0] *
				     drvui->subsys_fact[j][1][2]
				     +
				     drvui->subsys_fact[j][0][2] *
				     drvui->subsys_fact[j][1][0])
		 +
		 bst * cst * csal * (drvui->subsys_fact[j][0][1] *
				     drvui->subsys_fact[j][1][2]
				     +
				     drvui->subsys_fact[j][0][2] *
				     drvui->subsys_fact[j][1][1]))
		/ (ap * bp);
	    sngap = (float) sqrt (1.0 - csgap * csgap);
	    csbep = (drvui->subsys_fact[j][0][0] * drvui->subsys_fact[j][2][0] * ap * ap
		     + drvui->subsys_fact[j][0][1] * drvui->subsys_fact[j][2][1] * bp * bp
		     + drvui->subsys_fact[j][0][2] * drvui->subsys_fact[j][2][2] * cp * cp
		     +
		     ast * bst * csga * (drvui->subsys_fact[j][0][0] *
					 drvui->subsys_fact[j][2][1]
					 +
					 drvui->subsys_fact[j][0][1] *
					 drvui->subsys_fact[j][2][0])
		     +
		     ast * cst * csbe * (drvui->subsys_fact[j][0][0] *
					 drvui->subsys_fact[j][2][2]
					 +
					 drvui->subsys_fact[j][0][2] *
					 drvui->subsys_fact[j][2][0])
		     +
		     bst * cst * csal * (drvui->subsys_fact[j][0][1] *
					 drvui->subsys_fact[j][2][2]
					 +
					 drvui->subsys_fact[j][0][2] *
					 drvui->subsys_fact[j][2][1]))
		/ (ap * cp);
	    snbep = (float) sqrt (1.0 - csbep * csbep);
	    csalp = (drvui->subsys_fact[j][1][0] * drvui->subsys_fact[j][2][0] * ap * ap
		     + drvui->subsys_fact[j][1][1] * drvui->subsys_fact[j][2][1] * bp * bp
		     + drvui->subsys_fact[j][1][2] * drvui->subsys_fact[j][2][2] * cp * cp
		     +
		     ast * bst * csga * (drvui->subsys_fact[j][1][0] *
					 drvui->subsys_fact[j][2][1]
					 +
					 drvui->subsys_fact[j][1][1] *
					 drvui->subsys_fact[j][2][0])
		     +
		     ast * cst * csbe * (drvui->subsys_fact[j][1][0] *
					 drvui->subsys_fact[j][2][2]
					 +
					 drvui->subsys_fact[j][1][2] *
					 drvui->subsys_fact[j][2][0])
		     +
		     bst * cst * csal * (drvui->subsys_fact[j][1][1] *
					 drvui->subsys_fact[j][2][2]
					 +
					 drvui->subsys_fact[j][1][2] *
					 drvui->subsys_fact[j][2][1]))
		/ (bp * cp);
	    snalp = (float) sqrt (1.0 - csalp * csalp);
	    vol =
		ap * bp * cp * (float) sqrt (1.0 - csalp * csalp - csbep * csbep -
					     csgap * csgap - csalp * csbep * csgap);
	    a = bp * cp * snalp / vol;
	    b = ap * cp * snbep / vol;
	    c = ap * bp * sngap / vol;
	    alpha =
		(float) (180.0 / PI * acos ((csbep * csgap - csalp) / (snbep * sngap)));
	    beta =
		(float) (180.0 / PI * acos ((csalp * csgap - csbep) / (snalp * sngap)));
	    gamma =
		(float) (180.0 / PI * acos ((csalp * csbep - csgap) / (snalp * snbep)));
	    drvui->subsys_vol[j] = vol;	/* save reciprocal volume */

	    if (!Quick)
		fprintf (drvui->flout, "****** The cell for subsystem %d is:\n"
			 "******  a = %8.5f, b = %8.5f, c = %8.5f, Vol = %9.3f\n"
			 "******  alpha = %8.3f, beta = %8.3f,"
			 " gamma = %8.3f\n******\n", j + 1, a, b, c, 1.0f / vol, alpha,
			 beta, gamma);
	}
    }
    if (!in_line)
	(void) fclose (impin);
    else   
	position_cif (numblocks, startpos, impin, "# End of data", string);
}

/* ************************************************************** */
/* ************************************************************** */

void
import_fdat (char input[], int in_line, int Quick, int frame_no)

/* routine to extract structural information from a CSD 'FDAT' data file */
{
    FILE *impin;

    char string[256];		/* temporary string */

    char filename[256];		/* string for filename */

    int i, j, k, l;		/* temporary */

    char *spstring;		/* space group symbol */

    char lat_s[6][7];		/*temporary storage for lattice constants */

    int ldec[6];		/*temporary storage for decimal places of lattice constants */

    char number[4];

    char ncons[4], nrfacs[4], nrems[4], ndiss[4], nerrs[4], noprs[4], nrads[4], nats[4];

    char sats[4], nbnds[4];

    int ncards, nrfac = 0, nrem = 0, ndis = 0, nerr = 0, nopr = 0, nrad =
	0, nat = 0, sat = 0, ncon = 0, cellflag;
    int intflag, coordflag, centflag, errflag = 0;

    int nskip;			/* lines to skip between lattice constants and atoms */

    int natomin;

    memset (string, 0, 255);
    memset (nrfacs, 0, 4);
    memset (nrems, 0, 4);
    memset (ndiss, 0, 4);
    memset (ncons, 0, 4);
    memset (number, 0, 4);
    memset (nerrs, 0, 4);
    memset (noprs, 0, 4);
    memset (nrads, 0, 4);
    memset (nats, 0, 4);
    memset (sats, 0, 4);
    memset (nbnds, 0, 4);
    for (i = 0; i < 6; i++)
	memset (lat_s[i], 0, 7);
    if (in_line) {
	impin = drvui->fpin;
    } else {
	(void) sscanf (input, "%*s %s", filename);	/* get file name */
	if (!(impin = fopen (filename, "r"))) {
	    Error_Box ("Cannot Open CSD/FDAT Import File, Run aborted.");
	    return;
	}
    }

    if (!Quick)
	fprintf (drvui->flout, "****** Now interpreting CSD-style data\n");

    if (!fgets (input, 255, impin)) {	/* read a line */
	Error_Box ("Error reading CSD/FDAT file.");
	fclose (impin);
	return;
    }

/* the first line should contain the REFCODE of the database entry */
/* and information about the quality of the data and the length of */
/* the following sections                                          */

    if (strncmp (input, "#", 1) == 0) {

	(void) sscanf (input, "%23c%3d%3c%3c%3c%3c%3c%3c%3c%3c%3c%3c%1d%1d%1d%1d%1d",
		       string, &ncards, nrfacs, nrems, ndiss, nerrs, noprs, nrads,
		       nats, sats, nbnds, ncons, &cellflag, &intflag, &coordflag,
		       &centflag, &errflag);

	nrfac = atoi (nrfacs);
	ncon = atoi (ncons);
	nrem = atoi (nrems);
	ndis = atoi (ndiss);
	nerr = atoi (nerrs);
	nopr = atoi (noprs);
	nrad = atoi (nrads);
	nat = atoi (nats);
	sat = atoi (sats);
	input[0] = ' ';
	input[9] = '\0';	/* truncate title */
	if (!Quick)
	    fprintf (drvui->flout, "******      %s\n", input);
	if (!cellflag) {
	    Error_Box ("Error - this CSD entry has neither CELL nor COORDINATES");
	    fclose (impin);
	    return;
	}

	if (!coordflag) {
	    Error_Box ("Error - this CSD entry has no COORDINATES");
	    fclose (impin);
	    return;
	}

	if (errflag) {
	    Error_Box
		("Warning - this CSD entry has the ERROR flag set -\n at least some coordinates may be wrong");
	}

	if (!fgets (input, 255, impin)) {	/* read a line */
	    Error_Box ("Error reading CSD file.");
	    fclose (impin);
	    return;
	}
    }

/* lattice constants come in a fixed FORTRAN format followed by integers */
/* encoding precision, esd's, etc and finally by the spacegroup symbol */

    for (i = 0; i <= 5; ++i)
	drvui->lat_con[i] = 0.0f;

    (void) sscanf (input, "%6c%6c%6c%6c%6c%6c%1d%1d%1d%1d%1d%1d%21c%s", lat_s[0],
		   lat_s[1], lat_s[2], lat_s[3], lat_s[4], lat_s[5], &ldec[0],
		   &ldec[1], &ldec[2], &ldec[3], &ldec[4], &ldec[5], string, string);

    for (i = 0; i < 6; ++i) {
	drvui->lat_con[i] = (float) (atof (lat_s[i]) / pow (10, (double) ldec[i]));
    }

    if (drvui->sys == 0) {
	spstring = dissect_symbol (string);
	if (!Quick)
	    fprintf (drvui->flout, "******      Space group *%s*\n", spstring);
	symop (spstring);
	free (spstring);
    }

    if (!Quick)
	fprintf (drvui->flout,
		 "******      Cell     %8.5f %8.5f %8.5f %8.3f %8.3f %8.3f\n",
		 drvui->lat_con[0], drvui->lat_con[1], drvui->lat_con[2],
		 drvui->lat_con[3], drvui->lat_con[4], drvui->lat_con[5]);

/*skip R value, any textual information, symmetry codes and radii data */

    nskip = (nrfac + nrem + ndis + nerr) / 80 + 1;	/* skip comments (80chars per line) */
    nskip += nopr / 5 + 1;	/* skip symmetry operators (five per line) */
    nskip += nrad / 16 + 1;	/* skip atomic radii (sixteen per line) */


    for (i = 0; i < nskip; i++) {
	if (!fgets (input, 255, impin)) {	/* read a line */
	    Error_Box ("Error reading file while skipping records.");
	    if (!Quick) {
		fprintf (drvui->flout,
			 "nrfac:%d, nrem:%d, ndis:%d,nerr:%d, nopr:%d, nrad:%d, nskip:%d\n",
			 nrfac, nrem, ndis, nerr, nopr, nrad, nskip);
		fclose (impin);
		return;
	    }
	}
    }

    natomin = natom;		/* number of atoms already in list */

/* read coordinates until empty atom symbol encountered */

    for (i = 0; i < (nat + sat) / 3 + 1; i++) {	/*for all expected atom records */

	if (!fgets (input, 255, impin)) {	/* read a line */
	    Error_Box ("Error reading CSD file while reading atom coordinates.");
	    fclose (impin);
	    return;
	}

/* atoms come in groups of three */
	(void) sscanf (input,
		       "%c%c%c%c %7f %7f %7f %c%c%c%c %7f %7f %7f %c%c%c%c %7f %7f %7f",
		       &drvui->atoms[natom].atom_l[0], &drvui->atoms[natom].atom_l[1],
		       &drvui->atoms[natom].atom_l[2], &drvui->atoms[natom].atom_l[3],
		       &drvui->atoms[natom].atom_xyz[0], &drvui->atoms[natom].atom_xyz[1],
		       &drvui->atoms[natom].atom_xyz[2],
		       &drvui->atoms[natom + 1].atom_l[0],
		       &drvui->atoms[natom + 1].atom_l[1],
		       &drvui->atoms[natom + 1].atom_l[2],
		       &drvui->atoms[natom + 1].atom_l[3],
		       &drvui->atoms[natom + 1].atom_xyz[0],
		       &drvui->atoms[natom + 1].atom_xyz[1],
		       &drvui->atoms[natom + 1].atom_xyz[2],
		       &drvui->atoms[natom + 2].atom_l[0],
		       &drvui->atoms[natom + 2].atom_l[1],
		       &drvui->atoms[natom + 2].atom_l[2],
		       &drvui->atoms[natom + 2].atom_l[3],
		       &drvui->atoms[natom + 2].atom_xyz[0],
		       &drvui->atoms[natom + 2].atom_xyz[1],
		       &drvui->atoms[natom + 2].atom_xyz[2]);

/* and the atom list is terminated by numerical data */
	if (drvui->atoms[natom].atom_l[0] < 'A') {
	    break;
	}

/* scale coordinates */
	for (j = 0; j < 3; j++) {
	    drvui->atoms[natom].atom_xyz[j] /= 100000.0f;
	    drvui->atoms[natom + 1].atom_xyz[j] /= 100000.0f;
	    drvui->atoms[natom + 2].atom_xyz[j] /= 100000.0f;
	}

/* separate numbers from labels and remove trailing identifiers */

	for (l = 0; l < 3; l++) {
	    k = 0;
	    for (j = 1; j < 4; j++) {
		number[k] = ' ';
		if (drvui->atoms[natom + l].atom_l[j] >= '0'
		    && drvui->atoms[natom + l].atom_l[j] <= '9') {
		    number[k] = drvui->atoms[natom + l].atom_l[j];
		    k++;
		}
		if (k)
		    drvui->atoms[natom + l].atom_l[j] = ' ';
	    }
	    number[k] = '\0';
	    drvui->atoms[natom + l].atom_n = atoi (number);
	}

	if (Unique_Atom () == 1) {
/* report converted data of first atom */
	    if (!Quick)
		fprintf (drvui->flout, "******      Atom %c%c%c%c%3d %8.5f %8.5f %8.5f\n",
			 drvui->atoms[natom].atom_l[0], drvui->atoms[natom].atom_l[1],
			 drvui->atoms[natom].atom_l[2], drvui->atoms[natom].atom_l[3],
			 drvui->atoms[natom].atom_n, drvui->atoms[natom].atom_xyz[0],
			 drvui->atoms[natom].atom_xyz[1],
			 drvui->atoms[natom].atom_xyz[2]);
	    drvui->atoms[natom].atom_fn = frame_no;
	    drvui->atoms[natom].sv_atom_n = drvui->atoms[natom].atom_n;
	    natom++;
	    check_dynamic_storage ();
	}
	/* first unique */
	/* return now if no second atom on the line */
	if (natom == nat + natomin) {
	    break;		/* got all expected atoms */
	}
	if (drvui->atoms[natom].atom_l[0] == '\0')
	    break;
	if (Unique_Atom () == 1) {
/* report converted data of second atom */
	    if (!Quick)
		fprintf (drvui->flout, "******      Atom %c%c%c%c%3d %8.5f %8.5f %8.5f\n",
			 drvui->atoms[natom].atom_l[0], drvui->atoms[natom].atom_l[1],
			 drvui->atoms[natom].atom_l[2], drvui->atoms[natom].atom_l[3],
			 drvui->atoms[natom].atom_n, drvui->atoms[natom].atom_xyz[0],
			 drvui->atoms[natom].atom_xyz[1],
			 drvui->atoms[natom].atom_xyz[2]);
	    drvui->atoms[natom].atom_fn = frame_no;
	    drvui->atoms[natom].sv_atom_n = drvui->atoms[natom].atom_n;
	    natom++;
	    check_dynamic_storage ();
	}
	/* second unique */
	/* return now if no third atom on the line */
	if (natom == nat + natomin) {
	    break;		/* got all expected atoms */
	}
	if (drvui->atoms[natom].atom_l[0] == '\0')
	    break;
	if (Unique_Atom () == 1) {	/* third unique */
/* report converted data of third atom */
	    if (!Quick)
		fprintf (drvui->flout, "******      Atom %c%c%c%c%3d %8.5f %8.5f %8.5f\n",
			 drvui->atoms[natom].atom_l[0], drvui->atoms[natom].atom_l[1],
			 drvui->atoms[natom].atom_l[2], drvui->atoms[natom].atom_l[3],
			 drvui->atoms[natom].atom_n, drvui->atoms[natom].atom_xyz[0],
			 drvui->atoms[natom].atom_xyz[1],
			 drvui->atoms[natom].atom_xyz[2]);
	    drvui->atoms[natom].atom_fn = frame_no;
	    drvui->atoms[natom].sv_atom_n = drvui->atoms[natom].atom_n;
	    natom++;
	    check_dynamic_storage ();
	    if (natom == nat + natomin)
		break;		/* got all expected atoms */
	}			/* third unique */
    }				/*for all expected atom records */
    for (i++; i < (nat + sat) / 3; i++) {	/* skip remainder of atom records */

	if (!fgets (input, 255, impin)) {	/* read a line */
	    Error_Box ("Error reading file while skipping extra atom records.");
	    if (!in_line)
		(void) fclose (impin);
	    return;
	}
    }
    if (!errflag) {		/* only error-free entries have connectivity records */
	nskip = ncon / 40 + 1;
	if (nat + sat > 100)
	    nskip = ncon / 26 + 1;
	for (i = 0; i < nskip; i++) {
	    if (!fgets (input, 255, impin)) {	/* read a line */
		if (!in_line)
		    (void) fclose (impin);
		return;
	    }
	}
    }
    if (!in_line)
	(void) fclose (impin);

}				/* end of import_fdat */

/* ************************************************************** */
/* ************************************************************** */

void
import_gsas (char input[], int in_line, int Quick, int frame_no)

/* routine to extract structural information from a GSAS exp file */
{
    char string[256], filename[256], temp[5], tstring[20];

    int phase_no = 0, i, j;

    int ii, jj;

    FILE *impin;

    if (in_line) {
	Error_Box ("Cannot read inline GSAS data. Please use the import instruction");
	return;
    }

    ii = strlen (input) - 1;
    if (isdigit (input[ii])) {	// if number at end 
	for (jj = ii; jj > ii - 6; jj--) {
	    if (isspace (input[jj])) {	// separate number
		input[jj] = '\0';
		phase_no = atoi (&input[jj + 1]);
		break;
	    }
	    if (isalpha (input[jj]))
		break;		// unless it is attached to filename
	}
    }
//  (void)sscanf (input, "%s %s %d", filename, filename, &phase_no);    /* get file name and phase number */
    strcpy (filename, input);
    if ((phase_no < 1) || (phase_no > 9)) {
	Error_Box ("Illegal phase number in GSAS import instruction, Run aborted.");
	return;
    }
    if (!(impin = fopen (filename, "r"))) {
	Error_Box ("Cannot Open GSAS Import File, Run aborted.");
	return;
    }
    sprintf (temp, "CRS%d", phase_no);	/* make search object */
    for (;;) {
	if (!fgets (string, 100, impin)) {	/* search for phase info */
	    Error_Box ("Error reading phase info from GSAS Import File, Run aborted.");
	    fclose (impin);
	    return;
	}
	if (strncmp (string, temp, 4) == 0)
	    break;
    }
    for (;;) {			/* search for lattice constants */
	if (!fgets (string, 100, impin)) {
	    Error_Box ("Error reading lattice constants GSAS Import File, Run aborted.");
	    fclose (impin);
	    return;
	}
	(void) sscanf (string, "%s %s", tstring, tstring);
	if (strncmp (tstring, "ABC", 3) == 0) {
	    (void) sscanf (string, "%s %s %f %f %f", tstring, tstring, &drvui->lat_con[0], &drvui->lat_con[1], &drvui->lat_con[2]);	/* read lattice constants */
	    break;
	}
    }				/* search for lattice constants */
    for (;;) {			/* search for angles */
	if (!fgets (string, 100, impin)) {
	    Error_Box ("Error reading lattice angles GSAS Import File, Run aborted.");
	    fclose (impin);
	    return;
	}
	(void) sscanf (string, "%s %s", tstring, tstring);
	if (strncmp (tstring, "ANGLES", 6) == 0) {
	    (void) sscanf (string, "%s %s %f %f %f", tstring, tstring, &drvui->lat_con[3], &drvui->lat_con[4], &drvui->lat_con[5]);	/* read lattice constants */
	    break;
	}
    }				/* search for angles */
    if (!Quick)
	fprintf (drvui->flout,
		 "******      Cell     %8.5f %8.5f %8.5f %8.3f %8.3f %8.3f\n",
		 drvui->lat_con[0], drvui->lat_con[1], drvui->lat_con[2],
		 drvui->lat_con[3], drvui->lat_con[4], drvui->lat_con[5]);
    for (;;) {			/* position to atomic coordinates */
	if (!fgets (string, 100, impin)) {
	    Error_Box
		("Error reading atomic coordinates from GSAS Import File, Run aborted.");
	    fclose (impin);
	    return;
	}
	(void) sscanf (string, "%s %s", tstring, tstring);
	if (strncmp (tstring, "AT", 2) == 0)
	    break;
    }				/* position to atomic coordinates */

/* at this point, positioned at first atom line */
    for (;;) {			/* process atom parameters */
	(void) sscanf (string, "%s %s", tstring, tstring);
	if (strncmp (tstring, "AT", 2) != 0)
	    break;		/* out when atom stuff done */
	(void) sscanf (string, "%s %s %s %s %f %f %f", tstring, tstring, tstring,
		       tstring, &drvui->atoms[natom].atom_xyz[0],
		       &drvui->atoms[natom].atom_xyz[1],
		       &drvui->atoms[natom].atom_xyz[2]);
	drvui->atoms[natom].atom_n = 0;
	for (i = 0; i < 4; i++)
	    drvui->atoms[natom].atom_l[i] = ' ';	/* initialize atom name */
	j = 0;
	if (string[62] >= '0' && string[62] <= '9')
	    drvui->atoms[natom].atom_n = (int) string[62] - 48;
	else
	    drvui->atoms[natom].atom_l[j++] = string[62];
	if (string[63] >= '0' && string[63] <= '9')
	    drvui->atoms[natom].atom_n =
		10 * drvui->atoms[natom].atom_n + (int) string[63] - 48;
	else
	    drvui->atoms[natom].atom_l[j++] = string[63];
	if (string[64] >= '0' && string[64] <= '9')
	    drvui->atoms[natom].atom_n =
		10 * drvui->atoms[natom].atom_n + (int) string[64] - 48;
	else
	    drvui->atoms[natom].atom_l[j++] = string[64];
	if (string[65] >= '0' && string[65] <= '9')
	    drvui->atoms[natom].atom_n =
		10 * drvui->atoms[natom].atom_n + (int) string[65] - 48;
	else
	    drvui->atoms[natom].atom_l[j++] = string[65];
	if (!Quick)
	    fprintf (drvui->flout, "******      Atom %c%c%c%c%3d %8.5f %8.5f %8.5f\n",
		     drvui->atoms[natom].atom_l[0], drvui->atoms[natom].atom_l[1],
		     drvui->atoms[natom].atom_l[2], drvui->atoms[natom].atom_l[3],
		     drvui->atoms[natom].atom_n, drvui->atoms[natom].atom_xyz[0],
		     drvui->atoms[natom].atom_xyz[1], drvui->atoms[natom].atom_xyz[2]);
	(void) fgets (string, 100, impin);
	if (string[74] == 'A') {
	    drvui->ellips[drvui->n_ellips].ell_type = 1;	// Uij
	    (void) sscanf (string, "%s %s %s %f %f %f %f %f %f", tstring, tstring, tstring, &drvui->ellips[drvui->n_ellips].ellips[0], &drvui->ellips[drvui->n_ellips].ellips[1], &drvui->ellips[drvui->n_ellips].ellips[2], &drvui->ellips[drvui->n_ellips].ellips[3], &drvui->ellips[drvui->n_ellips].ellips[4], &drvui->ellips[drvui->n_ellips].ellips[5]);	/* read unique coefficients */
	    for (i = 0; i < 4; i++)
		drvui->ellips[drvui->n_ellips].ellips_l[i] =
		    drvui->atoms[natom].atom_l[i];
	    drvui->ellips[drvui->n_ellips].ellips_n = drvui->atoms[natom].atom_n;
	    drvui->ellips[drvui->n_ellips].save_el_number =
		drvui->ellips[drvui->n_ellips].ellips_n;
	    if (!Quick)
		fprintf (drvui->flout,
			 "******      Uij  %c%c%c%c%3d %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
			 drvui->atoms[natom].atom_l[0], drvui->atoms[natom].atom_l[1],
			 drvui->atoms[natom].atom_l[2], drvui->atoms[natom].atom_l[3],
			 drvui->atoms[natom].atom_n,
			 drvui->ellips[drvui->n_ellips].ellips[0],
			 drvui->ellips[drvui->n_ellips].ellips[1],
			 drvui->ellips[drvui->n_ellips].ellips[2],
			 drvui->ellips[drvui->n_ellips].ellips[3],
			 drvui->ellips[drvui->n_ellips].ellips[4],
			 drvui->ellips[drvui->n_ellips].ellips[5]);
	    drvui->n_ellips++;
	    check_dynamic_storage ();
	}			/* string[74] == 'A' */
	drvui->atoms[natom].atom_fn = frame_no;
	drvui->atoms[natom].sv_atom_n = drvui->atoms[natom].atom_n;
	natom++;
	check_dynamic_storage ();
	if (!fgets (string, 100, impin)) {
	    Error_Box ("Error reading GSAS Import File, Run aborted.");
	    fclose (impin);
	    return;
	}
    }				/* process atom parameters */
    for (;;) {			/* start loop for Space Group information */
	if (!fgets (string, 100, impin)) {
	    Error_Box ("Error reading GSAS Import File, Run aborted.");
	    fclose (impin);
	    return;
	}
	(void) sscanf (string, "%s %s", tstring, tstring);
	if (strncmp (tstring, "SG", 2) == 0)
	    break;		/* out when space group found */
    }				/* start loop for Space Group information */
    string[0] = 's';
    string[1] = 'p';
    string[2] = 'g';
    string[3] = 'r';
    for (i = 3; i < 60; i++)
	string[i + 1] = string[i + 10];
    string[60] = '\0';
    if (!Quick)
	fprintf (drvui->flout, "******      %s\n", string);
    symop (string);
    (void) fclose (impin);
}

/* ************************************************************** */
/* ************************************************************** */

void
import_schakal (char input[], int in_line, int Quick, int frame_no)

/* routine to extract structural information from a SCHAKAL data file */
{
    FILE *impin;

    char string[256];		/* temporary string */

    char filename[256];		/* string for filename */

    int intype;			/* type of datum read */

    int i, j, k;		/* temporary */

    if (in_line) {
	impin = drvui->fpin;
    } else {
//    (void)sscanf (input, "%s %s", filename, filename); /* get file name */
	strcpy (filename, input);
	if (!(impin = fopen (filename, "r"))) {
	    Error_Box ("Cannot Open SCHAKAL Import File, Run aborted.");
	    return;
	}
    }
    if (!Quick)
	fprintf (drvui->flout, "****** Now interpreting SCHAKAL-style data\n");
    while (!feof (impin)) {
	if (!fgets (input, 255, impin)) {	/* read a line */
	    Error_Box ("Error reading SCHAKAL file.");
	    if (!in_line)
		fclose (impin);
	    return;
	}
	intype = 0;
	for (i = 0; i < 4; i++) {
	    if (input[i] == '\0')
		intype = 3;
	}
	if (strncmp (input, "TITL", 4) == 0)
	    intype = 1;		/* echo and ignore */
	if (strncmp (input, "titl ", 4) == 0)
	    intype = 1;		/* echo and ignore */
	if (strncmp (input, "CELL", 4) == 0)
	    intype = 2;
	if (strncmp (input, "cell", 4) == 0)
	    intype = 2;
	if (strncmp (input, "ATOM", 4) == 0)
	    intype = 3;		/* atom coordinate card */
	if (strncmp (input, "atom", 4) == 0)
	    intype = 3;
	if (strncmp (input, "END", 3) == 0)
	    intype = 4;		/* trigger exit */
	if (strncmp (input, "end", 3) == 0)
	    intype = 4;		/* trigger exit */


	switch (intype) {
	case 1:
	    for (i = 0; i < 51; i++)
		if (input[i] == '\n')
		    input[i] = ' ';
	    input[50] = '\0';	/* truncate title */
	    if (!Quick)
		fprintf (drvui->flout, "******      %s\n", input);
	    break;		/* ignore title except for listing */
	case 2:
	    for (i = 3; i <= 5; ++i)
		drvui->lat_con[i] = 0.0f;
	    (void) sscanf (input, "%s %f %f %f %f %f %f", string, &drvui->lat_con[0],
			   &drvui->lat_con[1]
			   , &drvui->lat_con[2], &drvui->lat_con[3], &drvui->lat_con[4],
			   &drvui->lat_con[5]);
	    if (!Quick)
		fprintf (drvui->flout,
			 "******      Cell     %8.5f %8.5f %8.5f %8.3f %8.3f %8.3f\n",
			 drvui->lat_con[0], drvui->lat_con[1], drvui->lat_con[2],
			 drvui->lat_con[3], drvui->lat_con[4], drvui->lat_con[5]);
	    break;
	case 3:		/* read an atom line */
	    (void) sscanf (input, "%s %c%c%c%c %f %f %f", string,
			   &drvui->atoms[natom].atom_l[0], &drvui->atoms[natom].atom_l[1],
			   &drvui->atoms[natom].atom_l[2], &drvui->atoms[natom].atom_l[3],
			   &drvui->atoms[natom].atom_xyz[0],
			   &drvui->atoms[natom].atom_xyz[1],
			   &drvui->atoms[natom].atom_xyz[2]);
	    drvui->atoms[natom].atom_n = 0;
	    k = 0;
	    for (j = 0; j < 4; j++) {	/* extract atom number from SCHAKAL name */
		if (drvui->atoms[natom].atom_l[j] >= '0'
		    && drvui->atoms[natom].atom_l[j] <= '9') {
		    drvui->atoms[natom].atom_n =
			drvui->atoms[natom].atom_n * 10 +
			(int) drvui->atoms[natom].atom_l[j] - 48;
		} else {
		    drvui->atoms[natom].atom_l[k++] = drvui->atoms[natom].atom_l[j];
		}
	    }
	    for (j = k; j < 4; j++)
		drvui->atoms[natom].atom_l[j] = ' ';
	    if (drvui->atoms[natom].atom_n == 0)
		drvui->atoms[natom].atom_n = 1;
	    if (!Quick)
		fprintf (drvui->flout,
			 "******      Atom %c%c%c%c%3d %8.5f %8.5f %8.5f\n",
			 drvui->atoms[natom].atom_l[0], drvui->atoms[natom].atom_l[1],
			 drvui->atoms[natom].atom_l[2], drvui->atoms[natom].atom_l[3],
			 drvui->atoms[natom].atom_n, drvui->atoms[natom].atom_xyz[0],
			 drvui->atoms[natom].atom_xyz[1],
			 drvui->atoms[natom].atom_xyz[2]);
	    drvui->atoms[natom].atom_fn = frame_no;
	    drvui->atoms[natom].sv_atom_n = drvui->atoms[natom].atom_n;
	    natom++;
	    check_dynamic_storage ();
	    break;
	case 4:
	    if (!Quick)
		fprintf (drvui->flout, "****** END - end of SCHAKAL dataset\n");
	    if (!in_line)
		(void) fclose (impin);
	    return;
	default:
	    break;		/* end of default */
	}			/* end of switch */
    }
}				/* end of import_schakal */

/* ************************************************************** */
/* ************************************************************** */

void
import_shelx (char string[], int in_line, int Quick, int frame_no)

/* routine to extract structural information from a SHELX data file */
{
    FILE *impin;

    char input[256];		/* temporary for continuation cards */

    char filename[256];		/* string for filename */

    int intype;			/* type of datum read */

    int i, j, k;		/* temporary */

    float xray;			/* xray wavelength (read, but not used) */

    int read_start;		/* flag for start of atom list */

    int read_uij;		/* flag for uij card */

    float u11, u22;		/* temporary storage for uij */

    float site_occup;		/* occupancy (read, but not used) */

    static char shcmds[][68] = { "ZERR", "LATT", "SYMM", "SFAC", "DISP", "UNIT", "LAUE",
	"REM ", "MORE", "TIME", "OMIT", "SHEL", "BASF", "TWIN",
	"EXTI", "SWAT", "HOPE", "MERG", "SPEC", "RESI", "MOVE",
	"ANIS", "AFIX", "HFIX", "FRAG", "FEND", "EXYZ", "EADP",
	"EQIV", "CONN", "PART", "BIND", "FREE", "DFIX", "DANG",
	"BUMP", "SAME", "SADI", "CHIV", "FLAT", "DELU", "SIMU",
	"DEFS", "ISOR", "NCSY", "SUMP", "L.S.", "CGLS", "BLOC",
	"DAMP", "STIR", "WGHT", "FVAR", "BOND", "CONF", "MPLA",
	"RTAB", "HTAB", "LIST", "ACTA", "SIZE", "TEMP", "WPDB",
	"FMAP", "GRID", "PLAN", "MOLE", "REM\0"
    };
    static char cmd[6] = "     ";	/*temporary storage for start of line */

    float fvarpar[30];

    int nfvar = 0;

    int fv;

    int latt;

    char *p, *pp;

    memset (input, ' ', (size_t) 199);
    latt = 0;

    if (in_line) {
	impin = drvui->fpin;
    } else {
//    (void) sscanf (string, "%*s %s", filename);     /* get file name */
	strcpy (filename, string);
	if (!(impin = fopen (filename, "r"))) {
	    Error_Box ("Cannot Open SHELX Import File, Run aborted.");
	    return;
	}
    }
    read_start = 0;
    read_uij = 0;

    if (drvui->sys == 0) {	/* only do this if no spgr command given */
	/* initialize symmetry arrays */
	for (i = 0; i <= 23; ++i) {
	    for (j = 0; j <= 2; ++j) {
		drvui->ts[i][j] = 0.0;
		for (k = 0; k <= 2; ++k)
		    drvui->ss[i][j][k] = 0;
	    }
	}
	for (i = 0; i <= 2; ++i) {
	    drvui->lat_pos[0][i] = 0.0;
	    drvui->ss[0][i][i] = 1;
	    for (j = 0; j <= 3; ++j)
		drvui->spg[i][j] = ' ';
	}

	/* add trivial x,y,z symmetry */
	drvui->ts[0][0] = 0;
	drvui->ts[0][1] = 0;
	drvui->ts[0][2] = 0;
	drvui->ss[0][0][0] = 1;
	drvui->ss[0][1][1] = 1;
	drvui->ss[0][2][2] = 1;
	drvui->ng = 1;
    }

    if (!Quick)
	fprintf (drvui->flout, "****** Now interpreting SHELX-style data\n");
    while (!feof (impin)) {
	if (!fgets (string, 256, impin)) {	/* read a line */
	    if (!Quick)
		fprintf (drvui->flout, "Error reading input file.");
	    if (!in_line)
		fclose (impin);
	    return;
	}
	strcpy (input, string);
	if (strchr (input, '=')
	    && (strncmp (input, "TITL", 4) && strncmp (input, "REM", 3))) {
	    (*strchr (input, '=')) = '\0';
	    if (!fgets (string, 256, impin)) {
		Error_Box ("Error reading str file.");
		if (!in_line)
		    fclose (impin);
		return;
	    }
	    strcat (input, string);
	}
	intype = 0;
	for (i = 0; i < 4; i++) {
	    cmd[i] = (char) toupper (input[i]);
	    if (input[i] == '\0')
		intype = 3;
	}
	cmd[4] = '\0';

	if (strcmp (cmd, "TITL") == 0)
	    intype = 1;		/* echo and ignore */
	if (strcmp (cmd, "CELL") == 0)
	    intype = 2;

	if (strcmp (cmd, "HKLF") == 0)
	    intype = 4;		/* trigger exit */
	if (strncmp (cmd, "END", 3) == 0)
	    intype = 4;		/* trigger exit */

	if (strncmp (cmd, "FVAR", 4) == 0)
	    intype = 5;		/* store fixed parameters */

	if (drvui->sys == 0) {	/* ignore LATT/SYMM if spgr already seen */
	    if (strncmp (cmd, "LATT", 4) == 0)
		intype = 6;

	    if (strncmp (cmd, "SYMM", 4) == 0)
		intype = 7;
	}

	if (intype == 0) {	/* check against list of 'other' commands */
	    for (i = 0; i < 68; i++) {
		if (strlen (cmd) < 4 || strcmp (cmd, shcmds[i]) == 0) {
		    intype = 3;
		    break;
		}
	    }
	}
	if ((strncmp (cmd, "REM", 3) == 0) || (strncmp (cmd, "    ", 4) == 0))
	    intype = 3;
	switch (intype) {
	case 1:
	    trim_string (input, 51);
	    if (!Quick)
		fprintf (drvui->flout, "******      %s\n", input);
	    break;		/* ignore title except for listing */
	case 2:
	    for (i = 3; i <= 5; ++i)
		drvui->lat_con[i] = 0.0f;
	    (void) sscanf (input, "%*s %f %f %f %f %f %f %f", &xray, &drvui->lat_con[0]
			   , &drvui->lat_con[1], &drvui->lat_con[2], &drvui->lat_con[3],
			   &drvui->lat_con[4]
			   , &drvui->lat_con[5]);
	    if (!Quick)
		fprintf (drvui->flout,
			 "******      Cell     %8.5f %8.5f %8.5f %8.3f %8.3f %8.3f\n",
			 drvui->lat_con[0], drvui->lat_con[1], drvui->lat_con[2],
			 drvui->lat_con[3], drvui->lat_con[4], drvui->lat_con[5]);
	    break;
	case 3:		/* just ignore most SHELX commands */
	    break;
	case 4:
	    read_start = 0;
	    if (!Quick)
		fprintf (drvui->flout, "****** HKLF or END - end of SHELX dataset\n");
	    if (latt == 0 && drvui->sys == 0) {
		drvui->acentric = 0;	/* No LATT line given, spgp = P -1 */
		latt = 1;
		drvui->nlat = 1;
		drvui->nbr = 1;
	    }
	    if (latt == 3)
		drvui->sys = 5;
	    if (latt == 4) {	// cubic only if all axes equal, else orthorhombic
		if (fabs (drvui->lat_con[0] - drvui->lat_con[1]) < 1.e-4
		    && fabs (drvui->lat_con[0] - drvui->lat_con[2]) < 1.e-4)
		    drvui->sys = 6;
		else
		    drvui->sys = 3;
	    }
	    if (drvui->sys == 0)
		findsys ();
	    for (i = 1; i < drvui->n_ellips; i++) {
		if (drvui->ellips[i].ellips[0] <= 0.0f)
		    drvui->ellips[i].ell_type = -100;
	    }
	    if (drvui->autolabel == 2) {
		for (i = 0; i < natom; i++) {
		    drvui->labels[drvui->nlabel].label_fn = drvui->frame_no;
		    drvui->labels[drvui->nlabel].label_x[0] = drvui->atoms[i].atom_xyz[0];
		    drvui->labels[drvui->nlabel].label_x[1] = drvui->atoms[i].atom_xyz[1];
		    drvui->labels[drvui->nlabel].label_x[2] = drvui->atoms[i].atom_xyz[2];
		    int nn = 0;

		    for (j = 0; j < 4; j++) {
			if (drvui->atoms[i].atom_l[j] != ' ')
			    drvui->labels[drvui->nlabel].label_label[nn++] =
				drvui->atoms[i].atom_l[j];
		    }
		    drvui->labels[drvui->nlabel].label_label[nn] = '\0';
		    sprintf (cmd, "%d", drvui->atoms[i].atom_n);
		    strcat (drvui->labels[drvui->nlabel].label_label, cmd);

		    drvui->nlabel++;
		}
		drvui->autolabel = 1;
	    }

	    drvui->auto_ellipse = 0;
	    if (!in_line)
		fclose (impin);
	    return;
	case 5:
	    i = sscanf (input, "%*s %f %f %f %f %f %f %f",
			&fvarpar[nfvar], &fvarpar[nfvar + 1], &fvarpar[nfvar + 2],
			&fvarpar[nfvar + 3], &fvarpar[nfvar + 4], &fvarpar[nfvar + 5],
			&fvarpar[nfvar + 6]);
	    nfvar += 7;
	    if (i < 7)
		nfvar -= (7 - i);
	    break;
	case 6:
	    (void) sscanf (input, "%*s %d", &latt);
	    if (latt < 0) {
		drvui->acentric = 1;
		latt *= -1;
	    } else
		drvui->acentric = 0;
	    switch (latt) {
	    case 1:
		drvui->nlat = 1;
		drvui->nbr = 1;
		break;
	    case 2:
		drvui->nlat = 2;
		drvui->nbr = 6;
		drvui->lat_pos[1][0] = 0.5f;
		drvui->lat_pos[1][1] = 0.5f;
		drvui->lat_pos[1][2] = 0.5f;
		break;
	    case 3:
		drvui->nlat = 3;
		drvui->nbr = 7;
		drvui->lat_pos[1][0] = (1.f / 3.f);
		drvui->lat_pos[1][1] = (2.f / 3.f);
		drvui->lat_pos[1][2] = (2.f / 3.f);
		drvui->lat_pos[2][0] = (2.f / 3.f);
		drvui->lat_pos[2][1] = (1.f / 3.f);
		drvui->lat_pos[2][2] = (1.f / 3.f);
		break;
	    case 4:
		drvui->nlat = 4;
		drvui->nbr = 5;
		drvui->lat_pos[1][0] = 0.0;
		drvui->lat_pos[1][1] = 0.5f;
		drvui->lat_pos[1][2] = 0.5f;
		drvui->lat_pos[2][0] = 0.5f;
		drvui->lat_pos[2][1] = 0.0;
		drvui->lat_pos[2][2] = 0.5f;
		drvui->lat_pos[3][0] = 0.5f;
		drvui->lat_pos[3][1] = 0.5f;
		drvui->lat_pos[3][2] = 0.0;
		break;
	    case 5:
		drvui->nlat = 2;
		drvui->nbr = 2;
		drvui->lat_pos[1][0] = 0.0;
		drvui->lat_pos[1][1] = 0.5f;
		drvui->lat_pos[1][2] = 0.5f;
		break;
	    case 6:
		drvui->nlat = 2;
		drvui->nbr = 3;
		drvui->lat_pos[1][0] = 0.5f;
		drvui->lat_pos[1][1] = 0.0;
		drvui->lat_pos[1][2] = 0.5f;
		break;
	    case 7:
		drvui->nlat = 2;
		drvui->nbr = 4;
		drvui->lat_pos[1][0] = 0.5f;
		drvui->lat_pos[1][1] = 0.5f;
		drvui->lat_pos[1][2] = 0.0;
	    }
	    break;

	case 7:		/* SYMM */
	    drvui->ng++;
	    input[0] = ' ';
	    input[1] = ' ';
	    input[2] = ' ';
	    input[3] = ' ';
	    p = strstr (input, ",");
	    *p = '\0';
	    getsym (input, drvui->ng - 1, 0);
	    p++;
	    pp = strstr (p, ",");
	    *pp = '\0';
	    getsym (p, drvui->ng - 1, 1);
	    pp++;
	    getsym (pp, drvui->ng - 1, 2);
	    break;

	default:
	    site_occup = 1.;
	    u11 = u22 = 0.05f;
	    i = sscanf (input, "%c%c%c%c %d %f %f %f %f %f %f %f %f %f %f",
			&drvui->atoms[natom].atom_l[0], &drvui->atoms[natom].atom_l[1],
			&drvui->atoms[natom].atom_l[2], &drvui->atoms[natom].atom_l[3],
			&drvui->atoms[natom].atom_n, &drvui->atoms[natom].atom_xyz[0],
			&drvui->atoms[natom].atom_xyz[1],
			&drvui->atoms[natom].atom_xyz[2], &site_occup, &u11, &u22,
			&drvui->ellips[drvui->n_ellips].ellips[2],
			&drvui->ellips[drvui->n_ellips].ellips[5],
			&drvui->ellips[drvui->n_ellips].ellips[4],
			&drvui->ellips[drvui->n_ellips].ellips[3]);
	    for (j = 0; j < 3; j++)
		if (fabs (drvui->atoms[natom].atom_xyz[j]) > 15.) {
		    fv = ((int) fabs (drvui->atoms[natom].atom_xyz[j])) / 10;
		    if (drvui->atoms[natom].atom_xyz[j] < 0.)
			drvui->atoms[natom].atom_xyz[j] =
			    (float) (1. -
				     fvarpar[fv - 1] * (drvui->atoms[natom].atom_xyz[j] -
							fv * 10));
		    else
			drvui->atoms[natom].atom_xyz[j] =
			    (float) fvarpar[fv - 1] * (drvui->atoms[natom].atom_xyz[j] -
						       fv * 10);
		} else if (drvui->atoms[natom].atom_xyz[j] > 9.)
		    drvui->atoms[natom].atom_xyz[j] -= 10.0f;

	    if (fabs (site_occup) > 11.) {
		fv = ((int) fabs (site_occup)) / 10;
		if (site_occup < 0)
		    site_occup = (float) (1.0 - fvarpar[fv - 1]);
		else
		    site_occup = fvarpar[fv - 1];
	    }
	    if (fabs (u11) > 11.) {
		fv = ((int) fabs (u11)) / 10;
		if (u11 < 0.)
		    u11 = (float) (1.0 - fvarpar[fv - 1] * (fabs (u11) - fv * 10));
		else
		    u11 = fvarpar[fv - 1] * (u11 - fv * 10);
	    } else if (fabs (u11) > 1.) {
		if (u11 < 0.)
		    u11 += 10.;
		else
		    u11 -= 10.;
	    }

	    drvui->atoms[natom].atom_n = 0;
	    int k = 0;

	    for (j = 0; j < 4; j++)	/* extract atom number if part of SHELX name */
		if (drvui->atoms[natom].atom_l[j] >= '0'
		    && drvui->atoms[natom].atom_l[j] <= '9') {
		    drvui->atoms[natom].atom_n =
			10 * drvui->atoms[natom].atom_n +
			(int) drvui->atoms[natom].atom_l[j] - 48;
		} else {
		    drvui->atoms[natom].atom_l[k++] = drvui->atoms[natom].atom_l[j];
		}
	    for (j = k; j < 4; j++)
		drvui->atoms[natom].atom_l[j] = ' ';
	    drvui->atoms[natom].atom_l[4] = '\0';
	    if (drvui->atoms[natom].atom_n == 0)
		drvui->atoms[natom].atom_n = 1;

	    if (i > 11) {
		strcpy (drvui->ellips[drvui->n_ellips].ellips_l,
			drvui->atoms[natom].atom_l);
		drvui->ellips[drvui->n_ellips].ellips_n = drvui->atoms[natom].atom_n;
		drvui->ellips[drvui->n_ellips].save_el_number =
		    drvui->ellips[drvui->n_ellips].ellips_n;
		if (fabs (u22) > 11.) {
		    fv = ((int) fabs (u22)) / 10;
		    if (u22 < 0.)
			u22 = (float) (1. - fvarpar[fv - 1] * (fabs (u22) - fv * 10));
		    else
			u22 = fvarpar[fv - 1] * (u22 - fv * 10);
		} else if (fabs (u22) > 1.) {
		    if (u22 < 0.)
			u22 += 10.;
		    else
			u22 -= 10.;
		}

		if (fabs (drvui->ellips[drvui->n_ellips].ellips[2]) > 11.) {
		    fv = ((int) fabs (drvui->ellips[drvui->n_ellips].ellips[2])) / 10;
		    if (drvui->ellips[drvui->n_ellips].ellips[2] < 0.)
			drvui->ellips[drvui->n_ellips].ellips[2] =
			    (float) (1.0 -
				     fvarpar[fv -
					     1] *
				     (fabs (drvui->ellips[drvui->n_ellips].ellips[2]) -
				      fv * 10));
		    else
			drvui->ellips[drvui->n_ellips].ellips[2] = fvarpar[fv - 1] *
			    (drvui->ellips[drvui->n_ellips].ellips[2] - fv * 10);
		} else if (fabs (drvui->ellips[drvui->n_ellips].ellips[2]) > 1.) {
		    if (drvui->ellips[drvui->n_ellips].ellips[2] < 0.)
			drvui->ellips[drvui->n_ellips].ellips[2] += 10.;
		    else
			drvui->ellips[drvui->n_ellips].ellips[2] -= 10.;
		}

		if (fabs (drvui->ellips[drvui->n_ellips].ellips[5]) > 11.) {
		    fv = ((int) fabs (drvui->ellips[drvui->n_ellips].ellips[5])) / 10;
		    if (drvui->ellips[drvui->n_ellips].ellips[5] < 0.)
			drvui->ellips[drvui->n_ellips].ellips[5] =
			    (float) (1. -
				     fvarpar[fv -
					     1] *
				     (fabs (drvui->ellips[drvui->n_ellips].ellips[5]) -
				      fv * 10));
		    else
			drvui->ellips[drvui->n_ellips].ellips[5] = fvarpar[fv - 1] *
			    (drvui->ellips[drvui->n_ellips].ellips[5] - fv * 10);
		} else if (fabs (drvui->ellips[drvui->n_ellips].ellips[5]) > 1.) {
		    if (drvui->ellips[drvui->n_ellips].ellips[5] < 0.)
			drvui->ellips[drvui->n_ellips].ellips[5] += 10.;
		    else
			drvui->ellips[drvui->n_ellips].ellips[5] -= 10.;
		}
		if (drvui->ellips[drvui->n_ellips].ellips[4] > 11.) {
		    fv = ((int) fabs (drvui->ellips[drvui->n_ellips].ellips[4])) / 10;
		    if (drvui->ellips[drvui->n_ellips].ellips[4] < 0.)
			drvui->ellips[drvui->n_ellips].ellips[4] =
			    (float) (1.0 -
				     fvarpar[fv -
					     1] *
				     (fabs (drvui->ellips[drvui->n_ellips].ellips[4]) -
				      fv * 10));
		    else
			drvui->ellips[drvui->n_ellips].ellips[4] = fvarpar[fv - 1] *
			    (drvui->ellips[drvui->n_ellips].ellips[4] - fv * 10);
		} else if (fabs (drvui->ellips[drvui->n_ellips].ellips[4]) > 1.) {
		    if (drvui->ellips[drvui->n_ellips].ellips[4] < 0.)
			drvui->ellips[drvui->n_ellips].ellips[4] += 10.;
		    else
			drvui->ellips[drvui->n_ellips].ellips[4] -= 10.;
		}
		if (drvui->ellips[drvui->n_ellips].ellips[3] > 11.) {
		    fv = ((int) fabs (drvui->ellips[drvui->n_ellips].ellips[3])) / 10;
		    if (drvui->ellips[drvui->n_ellips].ellips[3] < 0.)
			drvui->ellips[drvui->n_ellips].ellips[3] =
			    (float) (1.0 -
				     fvarpar[fv -
					     1] *
				     (fabs (drvui->ellips[drvui->n_ellips].ellips[3]) -
				      fv * 10));
		    else
			drvui->ellips[drvui->n_ellips].ellips[3] = fvarpar[fv - 1] *
			    (drvui->ellips[drvui->n_ellips].ellips[3] - fv * 10);
		} else if (fabs (drvui->ellips[drvui->n_ellips].ellips[3]) > 1.) {
		    if (drvui->ellips[drvui->n_ellips].ellips[3] < 0.)
			drvui->ellips[drvui->n_ellips].ellips[3] += 10.;
		    else
			drvui->ellips[drvui->n_ellips].ellips[3] -= 10.;
		}

		drvui->ellips[drvui->n_ellips].ellips[0] = u11;
		drvui->ellips[drvui->n_ellips].ellips[1] = u22;
		drvui->ellips[drvui->n_ellips].ellips_ismod = 0;
		if (!Quick)
		    fprintf (drvui->flout,
			     "******      Atom %c%c%c%c%3d %8.5f %8.5f %8.5f (%8.5f)\n",
			     drvui->atoms[natom].atom_l[0], drvui->atoms[natom].atom_l[1],
			     drvui->atoms[natom].atom_l[2], drvui->atoms[natom].atom_l[3],
			     drvui->atoms[natom].atom_n, drvui->atoms[natom].atom_xyz[0],
			     drvui->atoms[natom].atom_xyz[1],
			     drvui->atoms[natom].atom_xyz[2], site_occup);
		if (!Quick)
		    fprintf (drvui->flout,
			     "******      Uij  %c%c%c%c%3d %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
			     drvui->atoms[natom].atom_l[0], drvui->atoms[natom].atom_l[1],
			     drvui->atoms[natom].atom_l[2], drvui->atoms[natom].atom_l[3],
			     drvui->atoms[natom].atom_n,
			     drvui->ellips[drvui->n_ellips].ellips[0],
			     drvui->ellips[drvui->n_ellips].ellips[1],
			     drvui->ellips[drvui->n_ellips].ellips[2],
			     drvui->ellips[drvui->n_ellips].ellips[3],
			     drvui->ellips[drvui->n_ellips].ellips[4],
			     drvui->ellips[drvui->n_ellips].ellips[5]);
		drvui->ellips[drvui->n_ellips].save_el_number =
		    drvui->ellips[drvui->n_ellips].ellips_n;
		drvui->do_ellipsoids = 1;
		if (drvui->auto_ellipse == 1) {
		    drvui->ellips[drvui->n_ellips].ell_type = 1001;	// Uij with ellipcolor
		    memset (drvui->ellips[drvui->n_ellips].ellips_col, 0, 40);
		    strcpy (drvui->ellips[drvui->n_ellips].ellips_col, "Gray30");
		} else {
		    drvui->ellips[drvui->n_ellips].ell_type = 1;	// Uij
		}
		drvui->n_ellips++;
		check_dynamic_storage ();
	    } else {
		if (i < 10) {
		    if (!Quick)
			fprintf (drvui->flout,
				 "******WARNING - too few parameters given for following atom:\n");
		    if (site_occup < 0)
			site_occup = 0.;
		    if (u11 < 0)
			u11 = 0.;
		}
		if (!Quick)
		    fprintf (drvui->flout,
			     "******      Atom %c%c%c%c%3d %8.5f %8.5f %8.5f (%8.5f) Uiso %8.5f\n",
			     drvui->atoms[natom].atom_l[0], drvui->atoms[natom].atom_l[1],
			     drvui->atoms[natom].atom_l[2], drvui->atoms[natom].atom_l[3],
			     drvui->atoms[natom].atom_n, drvui->atoms[natom].atom_xyz[0],
			     drvui->atoms[natom].atom_xyz[1],
			     drvui->atoms[natom].atom_xyz[2], site_occup, u11);
	    }
	    drvui->atoms[natom].atom_fn = frame_no;
	    drvui->atoms[natom].sv_atom_n = drvui->atoms[natom].atom_n;
	    natom++;
	    check_dynamic_storage ();
	    break;		/* end of default */
	}			/* end of switch */
    }
    if (!in_line)
	(void) fclose (impin);

    if (latt == 3)
	drvui->sys = 5;
    if (latt == 4)
	drvui->sys = 6;

    if (drvui->sys == 0)
	findsys ();
    for (i = 1; i < drvui->n_ellips; i++) {
	if (drvui->ellips[i].ellips[0] <= 0.0f)
	    drvui->ellips[i].ell_type = -100;
    }
    if (drvui->autolabel == 2) {
	for (i = 0; i < natom; i++) {
	    if (drvui->atoms[i].atom_fn != drvui->frame_no)
		continue;
	    drvui->labels[drvui->nlabel].label_fn = drvui->frame_no;
	    drvui->labels[drvui->nlabel].label_x[0] = drvui->atoms[i].atom_xyz[0];
	    drvui->labels[drvui->nlabel].label_x[1] = drvui->atoms[i].atom_xyz[1];
	    drvui->labels[drvui->nlabel].label_x[2] = drvui->atoms[i].atom_xyz[2];
	    int nn = 0;

	    for (j = 0; j < 4; j++) {
		if (drvui->atoms[i].atom_l[j] != ' ')
		    drvui->labels[drvui->nlabel].label_label[nn++] =
			drvui->atoms[i].atom_l[j];
	    }
	    drvui->labels[drvui->nlabel].label_label[nn] = '\0';
	    sprintf (cmd, "%d", drvui->atoms[i].atom_n);
	    strcat (drvui->labels[drvui->nlabel].label_label, cmd);

	    drvui->nlabel++;
	    check_dynamic_storage ();
	}
	drvui->autolabel = 1;
	drvui->auto_ellipse = 0;
    }
}				/* end of import_shelx */

/* ************************************************************** */
/* ************************************************************** */

void
import_wien (char string[], int in_line, int Quick, int frame_no)

/* routine to extract structural information from a WIEN2k struct file */
{
    FILE *impin;

    char filename[256];		/* string for filename */

    char line[85];

    int i, j, k;		/* temporary */

    char spgp[50], *spgr;

    char xs[15], ys[15], zs[15];

    int nat = 0;

    int skip = 0;

    char lattice[4];

    char atoms[4];

    char *spstring;

    if (in_line) {
	impin = drvui->fpin;
    } else {
//    (void) sscanf (string, "%*s %s", filename);     /* get file name */
	strcpy (filename, string);
	if (!(impin = fopen (filename, "r"))) {
	    Error_Box ("Cannot Open WIEN2k Import File, Run aborted.");
	    return;
	}
    }

    if (!Quick)
	fprintf (drvui->flout, "****** Now interpreting WIEN2k-style data\n");

    memset (string, 0, 85);


    if (!fgets (string, 85, impin)) {	/* skip title */
	if (!Quick)
	    fprintf (drvui->flout, "Error reading title.");
	if (!in_line)
	    fclose (impin);
	return;
    }
    trim_string (string, 85);	// kill any ^M or ^J at end of line
    string[85] = '\0';

    if (!Quick)
	fprintf (drvui->flout, "%s\n", string);

    if (!fgets (line, 85, impin)) {
	if (!Quick)
	    fprintf (drvui->flout, "Error reading header.");
	if (!in_line)
	    fclose (impin);
	return;
    }


    memset (spgp, 0, 40);
    memset (atoms, 0, 4);
    sscanf (line, "%s %*s%3c%s", lattice, atoms, spgp);
    natom = atoi (atoms);
    if (natom == 0)
	sscanf (line, "%s %*s %*s %d %*d %s", lattice, &natom, spgp);

    if (strchr (spgp, '_'))
	spgr = strchr (spgp, '_') + 1;
    else
	spgr = spgp;

    if (!Quick)
	fprintf (drvui->flout, "%s Lattice, %d atoms, space group %s\n", lattice, natom,
		 spgr);

    if (drvui->sys == 0) {
	strcpy (line, spgr);
	spstring = dissect_symbol (line);
	symop (spstring);
	free (spstring);
    }				// if spacegroup not already known

    if (!fgets (line, 85, impin)) {	/* read a line */
	if (!Quick)
	    fprintf (drvui->flout, "Error reading units flag.");
	if (!in_line)
	    fclose (impin);
	return;
    }

    if (!fgets (line, 85, impin)) {	/* read a line */
	if (!Quick)
	    fprintf (drvui->flout, "Error reading lattice constants.");
	if (!in_line)
	    fclose (impin);
	return;
    }

    i = sscanf (line, "%f %f %f %f %f %f",
		&drvui->lat_con[0], &drvui->lat_con[1], &drvui->lat_con[2],
		&drvui->lat_con[3], &drvui->lat_con[4], &drvui->lat_con[5]);

    if (i < 6 || drvui->lat_con[4] < 1. || drvui->lat_con[5] < 1.) {	//fortran fixed format lacking space characters
	line[39] = ' ';
	line[49] = ' ';
	i = sscanf (line, "%f %f %f %f %f %f",
		    &drvui->lat_con[0], &drvui->lat_con[1], &drvui->lat_con[2],
		    &drvui->lat_con[3], &drvui->lat_con[4], &drvui->lat_con[5]);
    }

    drvui->lat_con[0] *= (float) 0.529;	/* conversion from a.u. to angstroms */
    drvui->lat_con[1] *= (float) 0.529;
    drvui->lat_con[2] *= (float) 0.529;

    if (!Quick)
	fprintf (drvui->flout, "Lattice constants %f %f %f %5.2f %5.2f %5.2f\n",
		 drvui->lat_con[0], drvui->lat_con[1], drvui->lat_con[2],
		 drvui->lat_con[3], drvui->lat_con[4], drvui->lat_con[5]);

    nat = 0;

    skip = 0;

    while (!feof (impin) && nat < natom) {
	if (!fgets (line, 85, impin)) {	/* read a line */
	    if (!Quick)
		fprintf (drvui->flout, "Error reading atom record.");
	    if (!in_line)
		fclose (impin);
	    return;
	}
	if (!strncmp (line, "ATOM", 4) && skip == 0) {
	    sscanf (line, "%*s %*s %s %s %s", xs, ys, zs);
	    if (strchr(xs,'=')) 
		memmove(xs,strchr(xs,'=')+1,strlen(strchr(xs,'=')));
 	    if (strchr(ys,'=')) 
		memmove(ys,strchr(ys,'=')+1,strlen(strchr(ys,'=')));
	    if (strchr(zs,'=')) 
		memmove(zs,strchr(zs,'=')+1,strlen(strchr(zs,'=')));

	    if (drvui->sys == 5) {	/* struct files for hexagonal lattice have rhombohedral positions */
		drvui->atoms[nat].atom_xyz[0] = (float) atof (xs) * 2.f / 3.f
		    - (float) atof (ys) * 1.f / 3.f
		    - (float) atof (zs) * 1.f / 3.f;
		drvui->atoms[nat].atom_xyz[1] = (float) atof (xs) * 1.f / 3.f
		    + (float) atof (ys) * 1.f / 3.f
		    - (float) atof (zs) * 2.f / 3.f;
		drvui->atoms[nat].atom_xyz[2] = (float) atof (xs) * 1.f / 3.f
		    + (float) atof (ys) * 1.f / 3.f
		    + (float) atof (zs) * 1.f / 3.f;
	    } else {
		drvui->atoms[nat].atom_xyz[0] = (float) atof (xs);
		drvui->atoms[nat].atom_xyz[1] = (float) atof (ys);
		drvui->atoms[nat].atom_xyz[2] = (float) atof (zs);
	    }
	    drvui->atoms[nat].atom_fn = frame_no;
	    skip = 1;
	} else if (strstr (line, "NPT=")) {
	    sscanf (line, "%c%c%c%c",
		    &drvui->atoms[nat].atom_l[0], &drvui->atoms[nat].atom_l[1],
		    &drvui->atoms[nat].atom_l[2], &drvui->atoms[nat].atom_l[3]);
	    drvui->atoms[nat].atom_l[4] = '\0';

	    drvui->atoms[nat].atom_n = 0;
	    k = 0;
	    for (j = 0; j < 4; j++)	/* extract atom number if part of name */
		if (drvui->atoms[nat].atom_l[j] >= '0'
		    && drvui->atoms[nat].atom_l[j] <= '9') {
		    drvui->atoms[nat].atom_n =
			10 * drvui->atoms[nat].atom_n +
			(int) drvui->atoms[nat].atom_l[j] - 48;
		} else {
		    drvui->atoms[nat].atom_l[k++] = drvui->atoms[nat].atom_l[j];
		}
	    for (j = k; j < 4; j++)
		drvui->atoms[nat].atom_l[j] = ' ';
	    drvui->atoms[nat].atom_l[4] = '\0';
	    if (drvui->atoms[nat].atom_n == 0)
		drvui->atoms[nat].atom_n = 1;
	    drvui->atoms[nat].sv_atom_n = drvui->atoms[nat].atom_n;

	    nat++;
            skip = 0;
	}
    }

    if (!Quick) {
	for (i = 0; i < natom; i++)
	    fprintf (drvui->flout, "Atom %c%c%c%c %d  %f %f %f\n",
		     drvui->atoms[i].atom_l[0], drvui->atoms[i].atom_l[1],
		     drvui->atoms[i].atom_l[2], drvui->atoms[i].atom_l[3],
		     drvui->atoms[i].atom_n, drvui->atoms[i].atom_xyz[0],
		     drvui->atoms[i].atom_xyz[1], drvui->atoms[i].atom_xyz[2]);
    }

    if (!in_line)
	(void) fclose (impin);
}				/* end of import_wien */

/* ************************************************************** */
/* ************************************************************** */
void
import_discus (char input[], int in_line, int Quick, int frame_no)

/* routine to extract structural information from a DISCUS file */
{
    FILE *impin;

    char string[200];		/* temporary string */

    char filename[256];		/* string for filename */

    int i, j;			/* temporary */

    char *spstring = NULL;	/* space group symbol */

    int nx, ny, nz, nc;		/* number of cells if structure */

    int structflag = 0;		/* cell or structure */

    int ng = -1;		/* counter for symmetry generators */

    int gss[10][3][3];		/* temporary storage for */

    float gts[10][3];		/* additional generator symmetry */

    memset (string, 0, 100);

    if (in_line) {
	impin = drvui->fpin;
    } else {
//    (void)sscanf (input, "%s %s", filename, filename); /* get file name */
	strcpy (filename, input);
	if (!(impin = fopen (filename, "r"))) {
	    Error_Box ("Cannot Open DISCUS Import File, Run aborted.");
	    return;
	}
    }

    if (!Quick)
	fprintf (drvui->flout, "****** Now interpreting DISCUS-style data\n");

    if (!fgets (input, 100, impin)) {	/* read a line */
	Error_Box ("Error reading DISCUS file - end of file on first line.");
	if (!in_line)
	    (void) fclose (impin);
	return;
    }

    if (strncmp (input, "title", 5)) {
	Error_Box ("Error reading DISCUS file - file does not start with a title.");
	if (!in_line)
	    (void) fclose (impin);
	return;
    }

    while (!feof (impin)) {
	if (!fgets (input, 100, impin)) {	/* read a line */
	    Error_Box ("Error - unexpected end of file while reading DISCUS file.");
	    if (!in_line)
		(void) fclose (impin);
	    return;
	}

	if (!strncmp (input, "cell", 4)) {
	    sscanf (input, "%*s %f%*[ ,]%f%*[ ,]%f%*[ ,]%f%*[ ,]%f%*[ ,]%f",
		    &drvui->lat_con[0], &drvui->lat_con[1], &drvui->lat_con[2],
		    &drvui->lat_con[3], &drvui->lat_con[4], &drvui->lat_con[5]);
	    if (!Quick)
		fprintf (drvui->flout,
			 "******      Cell %5.3f %5.3f %5.3f %5.2f %5.2f %5.2f\n",
			 drvui->lat_con[0], drvui->lat_con[1], drvui->lat_con[2],
			 drvui->lat_con[3], drvui->lat_con[4], drvui->lat_con[5]);
	} else if (!strncmp (input, "ncell", 5)) {
	    structflag = 1;
	    sscanf (input, "%*s %d,%d,%d,%d", &nx, &ny, &nz, &nc);
	    spstring = dissect_symbol ("P1");
	    drvui->lat_con[0] *= (float) nx;
	    drvui->lat_con[1] *= (float) ny;
	    drvui->lat_con[2] *= (float) nz;
	    if (!Quick)
		fprintf (drvui->flout,
			 "******      %d x %d x %d supercell, transforming to P1\n", nx,
			 ny, nz);
	} else if (!strncmp (input, "spcgr", 5)) {
	    sscanf (input, "%*s %s", string);

	    spstring = dissect_symbol (string);

	    if (!Quick)
		fprintf (drvui->flout, "******      Space group %s\n", string);
/* spacegroup might be overridden by ncell, so do not call symop immediately */
	} else if (!strncmp (input, "generator", 9)) {
	    ng++;
	    sscanf (input,
		    "*s %d%*[ ,]%d%*[ ,]%d%*[ ,]%f%*[ ,]%d%*[ ,]%d%*[ ,]%d%*[, ]%f%*[ ,]%d%*[ ,]%d%*[ ,]%d%*[ ,]%f",
		    &gss[ng][0][0], &gss[ng][0][1], &gss[ng][0][2], &gts[ng][0],
		    &gss[ng][1][0], &gss[ng][1][1], &gss[ng][1][2], &gts[ng][1],
		    &gss[ng][2][0], &gss[ng][2][1], &gss[ng][2][2], &gts[ng][2]);
	    if (!Quick) {
		fprintf (drvui->flout, "****** Generator (%d %d %d) (%5.1f)\n",
			 gss[ng][0][0], gss[ng][0][1], gss[ng][0][2], gts[ng][0]);
		fprintf (drvui->flout, "******           (%d %d %d) (%5.1f)\n",
			 gss[ng][1][0], gss[ng][1][1], gss[ng][1][2], gts[ng][1]);
		fprintf (drvui->flout, "******           (%d %d %d) (%5.1f)\n",
			 gss[ng][2][0], gss[ng][2][1], gss[ng][2][2], gts[ng][2]);
	    }
	} else if (!strncmp (input, "ato", 3)) {
	    natom = 0;
	    while (fgets (input, 100, impin) != NULL) {
		if (input[0] == '#')
		    continue;
		if (strstr (input, "molecule")) {
		    if (!Quick)
			fprintf (drvui->flout,
				 "*** molecule keyword is currently not supported !\n");
		    continue;
		}
		memset (drvui->atoms[natom].atom_l, 0, 4);
		sscanf (input, "%4c %f %f %f %*f", drvui->atoms[natom].atom_l,
			&drvui->atoms[natom].atom_xyz[0],
			&drvui->atoms[natom].atom_xyz[1],
			&drvui->atoms[natom].atom_xyz[2]);
		if (structflag == 1) {
		    drvui->atoms[natom].atom_xyz[0] /= (float) nx;
		    drvui->atoms[natom].atom_xyz[1] /= (float) ny;
		    drvui->atoms[natom].atom_xyz[2] /= (float) nz;
		}
		drvui->atoms[natom].atom_fn = frame_no;
		int k = 0;

		for (j = 0; j < 4; j++) {	/* extract atom number if part of name */
		    if (drvui->atoms[natom].atom_l[j] >= '0'
			&& drvui->atoms[natom].atom_l[j] <= '9') {
			drvui->atoms[natom].atom_n =
			    10 * drvui->atoms[natom].atom_n +
			    (int) drvui->atoms[natom].atom_l[j] - 48;
		    } else {
			drvui->atoms[natom].atom_l[k++] = drvui->atoms[natom].atom_l[j];
		    }
		}
		for (j = k; j < 4; j++)
		    drvui->atoms[natom].atom_l[j] = ' ';

		drvui->atoms[natom].atom_l[4] = '\0';
		if (drvui->atoms[natom].atom_n == 0)
		    drvui->atoms[natom].atom_n = 1;
		drvui->atoms[natom].sv_atom_n = drvui->atoms[natom].atom_n;
		natom++;
		check_dynamic_storage ();
	    }			/* here at end of input, i.e. when the last atom has been read */
	    if (!Quick)
		for (i = 0; i < natom; i++)
		    fprintf (drvui->flout, "******      Atom %s %d %f %f %f\n",
			     drvui->atoms[i].atom_l, drvui->atoms[i].atom_n,
			     drvui->atoms[i].atom_xyz[0], drvui->atoms[i].atom_xyz[1],
			     drvui->atoms[i].atom_xyz[2]);

	    symop (spstring);
	    if (ng > -1 && structflag == 0) {
		for (i = 0; i <= ng; i++) {
		    int ig = ++drvui->ng;

		    drvui->ss[ig][0][0] = gss[i][0][0];
		    drvui->ss[ig][0][1] = gss[i][0][1];
		    drvui->ss[ig][0][2] = gss[i][0][2];
		    drvui->ss[ig][1][0] = gss[i][1][0];
		    drvui->ss[ig][1][1] = gss[i][1][1];
		    drvui->ss[ig][1][2] = gss[i][1][2];
		    drvui->ss[ig][2][0] = gss[i][2][0];
		    drvui->ss[ig][2][1] = gss[i][2][1];
		    drvui->ss[ig][2][2] = gss[i][2][2];
		    drvui->ts[ig][0] = gts[i][0];
		    drvui->ts[ig][1] = gts[i][1];
		    drvui->ts[ig][2] = gts[i][2];
		}
	    }
	    if (!in_line)
		(void) fclose (impin);
	    if (spstring) free(spstring);
	    return;
	}			/* end of atom records */
    }				/* main keyword parsing loop */

/* normal exit is through the EOF condition in the atoms loop */
    Error_Box ("Unexpected end of DISCUS file - no atom records found.");
    if (spstring) 
	free (spstring);
    if (!in_line)
	(void) fclose (impin);
    return;
}

/* ************************************************************** */
/* ************************************************************** */

void
import_pcr (char input[], int in_line, int Quick, int frame_no)

/* routine to extract structural information from a FULLPROF pcr file */
{
    char string[256], filename[256], temp[25], tstring[21];

    int phase_no, i, j;

    int ii, jj;

    int nat, anis;

    FILE *impin;

    if (in_line) {
	Error_Box ("Cannot read inline FULLPROF data. Please use the import instruction");
	return;
    }

    phase_no = 0;

    ii = strlen (input) - 1;
    if (isdigit (input[ii])) {	// if number at end 
	for (jj = ii; jj > ii - 6; jj--) {
	    if (isspace (input[jj])) {	// separate number
		input[jj] = '\0';
		phase_no = atoi (&input[jj + 1]);
		break;
	    }
	    if (isalpha (input[jj]))
		break;		// unless it is attached to filename
	}
    }
    strcpy (filename, input);
    if (phase_no < 1) {
	Error_Box ("Illegal phase number in FULLPROF import instruction, Run aborted.");
	return;
    }
    if (!(impin = fopen (filename, "r"))) {
	Error_Box ("Cannot Open FULLPROF Import File, Run aborted.");
	return;
    }
    sprintf (temp, "for PHASE number: %3d", phase_no);	/* make search object */
    for (;;) {
	if (!fgets (string, 200, impin)) {	/* search for phase info */
	    Error_Box
		("Error reading phase info from FULLPROF Import File, Run aborted.");
	    (void) fclose (impin);
	    return;
	}
	if (strstr (string, temp) != 0)
	    break;
    }

    if (!fgets (string, 200, impin)) {
	Error_Box ("Error reading from FULLPROF Import File, Run aborted.");
	(void) fclose (impin);
	return;
    }

    if (!fgets (string, 200, impin)) {
	Error_Box ("Error reading from FULLPROF Import File, Run aborted.");
	(void) fclose (impin);
	return;
    }
    (void) sscanf (string, "%20c", tstring);
    tstring[20] = '\0';
    if (!Quick)
	fprintf (drvui->flout, "******      Phase %d: %s\n", phase_no, tstring);	//title

    do {
	if (!fgets (string, 200, impin)) {
	    Error_Box ("Error reading from FULLPROF Import File, Run aborted.");
	    (void) fclose (impin);
	    return;
	}
	trim_string (string, 200);
    } while (string[0] == '!' || !strlen (string));	//skip any comment lines and empty lines

    nat = 0;
    sscanf (string, "%d", &nat);	// next record must start with number of atoms

    for (;;) {			/* search for spacegroup */
	if (!fgets (string, 200, impin)) {
	    Error_Box ("Error reading spacegroup from FULLPROF File, Run aborted.");
	    (void) fclose (impin);
	    return;
	}
	if (strstr (string, "Space group")) {
	    (void) sscanf (string, "%20c", tstring);
	    tstring[20] = '\0';
	    strcpy (string, "spgr ");
	    strcat (string, tstring);
	    if (!Quick)
		fprintf (drvui->flout, "******      %s\n", string);
	    symop (string);
	    break;
	}
    }

    for (;;) {			/* search for start of atom list */
	if (!fgets (string, 200, impin)) {
	    Error_Box
		("Error reading atom records from FULLPROF Import File, Run aborted.");
	    (void) fclose (impin);
	    return;
	}
	if (!strncmp (string, "!Atom Typ", 9))
	    break;
    }

    do {
	if (!fgets (string, 200, impin)) {
	    Error_Box ("Error reading atom data from FULLPROF Import File, Run aborted.");
	    (void) fclose (impin);
	    return;
	}
	if (string[0] == '!')
	    continue;		//skip any comment lines

	anis = 0;
	if (isalpha (string[0])) {
	    sscanf (string, "%s %*s %f %f %f %*f %*f %*d %*d %d %*d",
		    tstring, &drvui->atoms[natom].atom_xyz[0],
		    &drvui->atoms[natom].atom_xyz[1], &drvui->atoms[natom].atom_xyz[2],
		    &anis);

	    drvui->atoms[natom].atom_n = 0;
	    for (i = 0; i < 4; i++)
		drvui->atoms[natom].atom_l[i] = ' ';	/* initialize atom name */
	    j = 0;
	    if (string[0] >= '0' && string[0] <= '9')
		drvui->atoms[natom].atom_n = (int) string[0] - 48;
	    else
		drvui->atoms[natom].atom_l[j++] = string[0];
	    if (string[1] >= '0' && string[1] <= '9')
		drvui->atoms[natom].atom_n =
		    10 * drvui->atoms[natom].atom_n + (int) string[1] - 48;
	    else
		drvui->atoms[natom].atom_l[j++] = string[1];
	    if (string[2] >= '0' && string[2] <= '9')
		drvui->atoms[natom].atom_n =
		    10 * drvui->atoms[natom].atom_n + (int) string[2] - 48;
	    else
		drvui->atoms[natom].atom_l[j++] = string[2];
	    if (string[3] >= '0' && string[3] <= '9')
		drvui->atoms[natom].atom_n =
		    10 * drvui->atoms[natom].atom_n + (int) string[3] - 48;
	    else
		drvui->atoms[natom].atom_l[j++] = string[3];
	    if (!Quick)
		fprintf (drvui->flout,
			 "******      Atom %c%c%c%c%3d %8.5f %8.5f %8.5f\n",
			 drvui->atoms[natom].atom_l[0], drvui->atoms[natom].atom_l[1],
			 drvui->atoms[natom].atom_l[2], drvui->atoms[natom].atom_l[3],
			 drvui->atoms[natom].atom_n, drvui->atoms[natom].atom_xyz[0],
			 drvui->atoms[natom].atom_xyz[1],
			 drvui->atoms[natom].atom_xyz[2]);
	    drvui->atoms[natom].atom_fn = frame_no;
	    drvui->atoms[natom].sv_atom_n = drvui->atoms[natom].atom_n;
	    drvui->atoms[natom].atom_fn = frame_no;
	    drvui->atoms[natom].sv_atom_n = drvui->atoms[natom].atom_n;
	    fgets (string, 200, impin);	// skip codewords
	    if (anis == 2) {
		strcpy (drvui->ellips[drvui->n_ellips].ellips_l,
			drvui->atoms[natom].atom_l);
		drvui->ellips[drvui->n_ellips].ellips_n = drvui->atoms[natom].atom_n;
		drvui->ellips[drvui->n_ellips].save_el_number =
		    drvui->ellips[drvui->n_ellips].ellips_n;
		drvui->ellips[drvui->n_ellips].ellips_ismod = 0;
		drvui->ellips[drvui->n_ellips].ell_type = 0;
		fgets (string, 200, impin);	// b factors
		if (string[0] == '!')
		    fgets (string, 100, impin);
		sscanf (string, "%f %f %f %f %f %f",
			&drvui->ellips[drvui->n_ellips].ellips[0],
			&drvui->ellips[drvui->n_ellips].ellips[1],
			&drvui->ellips[drvui->n_ellips].ellips[2],
			&drvui->ellips[drvui->n_ellips].ellips[3],
			&drvui->ellips[drvui->n_ellips].ellips[4],
			&drvui->ellips[drvui->n_ellips].ellips[5]);
		if (!Quick)
		    fprintf (drvui->flout,
			     "******      Beta_ij %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
			     drvui->ellips[drvui->n_ellips].ellips[0],
			     drvui->ellips[drvui->n_ellips].ellips[1],
			     drvui->ellips[drvui->n_ellips].ellips[2],
			     drvui->ellips[drvui->n_ellips].ellips[3],
			     drvui->ellips[drvui->n_ellips].ellips[4],
			     drvui->ellips[drvui->n_ellips].ellips[5]);
		drvui->n_ellips++;
		check_dynamic_storage ();
	    }
	    natom++;
	    check_dynamic_storage ();
	}			// if line starts with (atom) label
    } while (natom < nat);

    for (;;) {			/* search for lattice constants */
	if (!fgets (string, 200, impin)) {
	    Error_Box
		("Error reading lattice constants from FULLPROF Import File, Run aborted.");
	    (void) fclose (impin);
	    return;
	}
	if (strstr (string, "alpha") && strstr (string, "beta")
	    && strstr (string, "gamma")) {
	    if (!fgets (string, 200, impin)) {
		Error_Box
		    ("Error reading lattice constants from FULLPROF Import File, Run aborted.");
		(void) fclose (impin);
		return;
	    }
	    sscanf (string, "%f %f %f %f %f %f", &drvui->lat_con[0], &drvui->lat_con[1], &drvui->lat_con[2], &drvui->lat_con[3], &drvui->lat_con[4], &drvui->lat_con[5]);	/* read lattice constants */
	    if (!Quick)
		fprintf (drvui->flout,
			 "******      Cell     %8.5f %8.5f %8.5f %8.3f %8.3f %8.3f\n",
			 drvui->lat_con[0], drvui->lat_con[1], drvui->lat_con[2],
			 drvui->lat_con[3], drvui->lat_con[4], drvui->lat_con[5]);
	    (void) fclose (impin);
	    return;
	}
    }
}

/* ************************************************************** */
/* ************************************************************** */

void
import_exc (char string[], int in_line, int Quick, int frame_no)

/* routine to extract structural information from an Exciting GEOMETRY.OUT file */
{
    FILE *impin;

    char filename[256];		/* string for filename */

    char line[85];

    int i, j, k;		/* temporary */

    int nat = 0;

    int nsites = 0;

    char *p, *pname;

    float scale = 1.0f;

    float scale1 = 1.0f;

    float scale2 = 1.0f;

    float scale3 = 1.0f;

    float avec[3], bvec[3], cvec[3];

    char spstring[10] = "spgr P 1";

    char *name;

    if (in_line) {
	impin = drvui->fpin;
    } else {
//    (void) sscanf (string, "%*s %s", filename);     /* get file name */
	strcpy (filename, string);
	if (!(impin = fopen (filename, "r"))) {
	    Error_Box ("Cannot Open Exciting Import File, Run aborted.");
	    return;
	}
    }

    if (!Quick)
	fprintf (drvui->flout, "****** Now interpreting GEOMETRY.OUT-style data\n");

    memset (string, 0, 85);

    while (1) {
	if (!fgets (line, 85, impin)) {	/* read a line */
	    if (!Quick)
		fprintf (drvui->flout, "End of input reached.\n");
	    if (!in_line)
		fclose (impin);
	    break;
	}
	Blank_Strip (line);	// remove leading spaces
	if (!strncmp (line, "scale1", 6)) {
	    fgets (line, 85, impin);
	    sscanf (line, "%f", &scale1);
	    continue;
	}
	if (!strncmp (line, "scale2", 6)) {
	    fgets (line, 85, impin);
	    sscanf (line, "%f", &scale2);
	    continue;
	}
	if (!strncmp (line, "scale3", 6)) {
	    fgets (line, 85, impin);
	    sscanf (line, "%f", &scale3);
	    continue;
	}
	if (!strncmp (line, "scale", 5)) {
	    fgets (line, 85, impin);
	    sscanf (line, "%f", &scale);
	    continue;
	}

	if (!strncmp (line, "avec", 4)) {
	    fgets (line, 85, impin);
	    sscanf (line, "%f %f %f", &avec[0], &avec[1], &avec[2]);
	    fgets (line, 85, impin);
	    sscanf (line, "%f %f %f", &bvec[0], &bvec[1], &bvec[2]);
	    fgets (line, 85, impin);
	    sscanf (line, "%f %f %f", &cvec[0], &cvec[1], &cvec[2]);
	    continue;
	}

	if (!strncmp (line, "atoms", 5)) {
	    natom = 0;
	    fgets (line, 85, impin);
	    sscanf (line, "%d", &nat);
	    for (j = 0; j < nat; j++) {
		name = (char *) zalloc (40 * sizeof (char));
		fgets (line, 85, impin);
		sscanf (line, "%s", name);
		pname = name + 1;	// skip first quotation mark
		p = strrchr (pname, '/');
		if (p)
		    pname = ++p;
		p = strrchr (pname, '.');
		if (!p)
		    p = strrchr (pname, 39);
		if (p)
		    *p = '\0';	// strip either the .in or the trailing quotation mark
		for (k = 0; k < 3; k++)
		    if (strlen (pname) < 4)
			strcat (pname, " ");

		fgets (line, 85, impin);
		sscanf (line, "%d", &nsites);
		for (k = 0; k < nsites; k++) {
		    fgets (line, 85, impin);
		    sscanf (line, "%f %f %f",
			    &drvui->atoms[natom].atom_xyz[0],
			    &drvui->atoms[natom].atom_xyz[1],
			    &drvui->atoms[natom].atom_xyz[2]);
		    strcpy (drvui->atoms[natom].atom_l, pname);
		    drvui->atoms[natom].atom_n = k + 1;
		    drvui->atoms[natom].atom_fn = frame_no;
		    drvui->atoms[natom].sv_atom_n = drvui->atoms[natom].atom_n;
		    natom++;
		    check_dynamic_storage ();
		}
		free (name);
	    }
	    continue;
	}			// if atoms block
    }				// while not eof

    if (natom == 0) {
	Error_Box("No atoms could be read.\nThis file is probably not in elk.in or GEOMETRY.OUT format.");
    }
    for (k = 0; k < 3; k++) {	// scale lattice vectors and convert to angstroms
	avec[k] *= scale * scale1;
	bvec[k] *= scale * scale2;
	cvec[k] *= scale * scale3;
    }
    drvui->lat_con[0] =
	(float) sqrt (avec[0] * avec[0] + avec[1] * avec[1] + avec[2] * avec[2]);
    drvui->lat_con[1] =
	(float) sqrt (bvec[0] * bvec[0] + bvec[1] * bvec[1] + bvec[2] * bvec[2]);
    drvui->lat_con[2] =
	(float) sqrt (cvec[0] * cvec[0] + cvec[1] * cvec[1] + cvec[2] * cvec[2]);
    drvui->lat_con[3] =
	(float) (180 / M_PI *
		 acos ((bvec[0] * cvec[0] + bvec[1] * cvec[1] + bvec[2] * cvec[2])
		       / (drvui->lat_con[1] * drvui->lat_con[2])));
    drvui->lat_con[4] =
	(float) (180 / M_PI *
		 acos ((avec[0] * cvec[0] + avec[1] * cvec[1] + avec[2] * cvec[2])
		       / (drvui->lat_con[0] * drvui->lat_con[2])));
    drvui->lat_con[5] =
	(float) (180 / M_PI *
		 acos ((avec[0] * bvec[0] + avec[1] * bvec[1] + avec[2] * bvec[2])
		       / (drvui->lat_con[0] * drvui->lat_con[1])));

    drvui->lat_con[0] *= 0.529f;
    drvui->lat_con[1] *= 0.529f;
    drvui->lat_con[2] *= 0.529f;

    symop (spstring);
    if (!Quick)
	fprintf (drvui->flout, "Lattice constants %f %f %f %5.2f %5.2f %5.2f\n",
		 drvui->lat_con[0], drvui->lat_con[1], drvui->lat_con[2],
		 drvui->lat_con[3], drvui->lat_con[4], drvui->lat_con[5]);

    if (!Quick) {
	for (i = 0; i < natom; i++)
	    fprintf (drvui->flout, "Atom %c%c%c%c %d  %f %f %f\n",
		     drvui->atoms[i].atom_l[0], drvui->atoms[i].atom_l[1],
		     drvui->atoms[i].atom_l[2], drvui->atoms[i].atom_l[3],
		     drvui->atoms[i].atom_n, drvui->atoms[i].atom_xyz[0],
		     drvui->atoms[i].atom_xyz[1], drvui->atoms[i].atom_xyz[2]);
    }

}				/* end of import_exc */

/* ************************************************************** */
/* ************************************************************** */

int
position_cif (int numblocks, long startpos, FILE * impin, const char *search_string, char *string)
{
/* routine to search through a CIF until the token in 'search_string' is found 
   If the string is not found, the routine returns 1. If 'search_string' is found,
   the routine returns 0
 */

//    rewind (impin);
    fseek (impin, startpos, SEEK_SET);
    memset (string, 0, 255);
    if (numblocks > 1)
	skip_blocks (numblocks - 1, impin);
    while (strncmp (string, search_string, strlen (search_string))) {
	if (!get_next_token (string, 256, impin)) {
	    return 1;
	}
    }
    return 0;
}


char *
dissect_symbol (const char *line)
{
/* routine to prepare a given space group name for translating with symop  */
    char *spstring;

    int i, j, k, l;

    int k3, k4;

    spstring = (char *) zalloc (50 * sizeof (char));
    spstring[0] = 's';
    spstring[1] = 'p';
    spstring[2] = 'g';
    spstring[3] = 'r';
    spstring[4] = ' ';
    spstring[5] = line[0];
    spstring[6] = ' ';
    j = 6;
    k = 0;
    k3 = 0;
    k4 = 0;
    l = (int) strlen (line);
    for (i = 1; i <= l; i++) {
	switch (line[i]) {
	case ' ':
	    break;
	case '/':
	    spstring[++j] = line[i];
	    break;
	case '-':
	    if (spstring[j] != ' ')
		spstring[++j] = ' ';
	    spstring[++j] = line[i++];
	    spstring[++j] = line[i];
	    break;
	default:
	    j++;
	    /* do we have to collate this symbol (for screw axes or mirror planes) ? */
	    if (line[i] > '0' && line[i] <= '9') {	/* it is a number */
		if (line[0] == 'R' ||	/* no screw axes in R */
		    l == 4 ||	/* no screw axes in space groups with four-letter names */
		    k != 0 ||	/* have just appended to a number */
		    k3 != 0 ||	/* cannot have further screw axes after 3_x */
		    (k4 != 0 && spstring[j - 1] == '3')) {	/* cannot have 3_x after a 4_x */
		    spstring[j] = ' ';	/* so separate it */
		    k = 0;	/* and reset the collation flag */
		    j++;
		} else if (spstring[j - 1] != '-' && (spstring[j - 1] <= line[i]
						      || spstring[j - 1] >= 'a')) {
		    /* neither a minus sign nor a higher number precedes it */
		    spstring[j] = ' ';	/* so separate it */
		    k = 0;	/* and reset the collation flag */
		    j++;
		} else {
		    k++;	/* set flag for collation */
		    if (spstring[j - 1] == '3')
			k3 = 1;	/* if 3_x screw, prohibit further screw axes */
		    if (spstring[j - 1] == '4')
			k4 = 1;	/* if 4_x screw, prohibit 3_x screw axes */
		}
	    } else if (j > 1) {	/* not a number */
		if (spstring[j - 1] != '/' && spstring[j - 1] != ' ') {	/*but not a mirror plane either */
		    spstring[j] = ' ';	/* separate it */
		    j++;
		}
	    }
	    spstring[j] = line[i];	/* now add the symbol */
	    break;
	}			/*switch */
    }				/*for i */
    spstring[j + 1] = '\0';

/* Space group 90, P 4 21 2, is not correctly identified by this code. */
    if (!strcmp (line, "P4212"))
	strcpy (spstring, "spgp P 4 21 2");

    return (spstring);
}
