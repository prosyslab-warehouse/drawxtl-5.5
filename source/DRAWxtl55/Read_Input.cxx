// $Id: Read_Input.cxx 1107 2011-01-19 23:53:52Z martin $
//
// Read_Input.cxx - Source module for DRAWxtl V5.5
// Coded using the FLTK 1.1.6 widget set
//
//     Larry W. Finger, Martin Kroeker and Brian Toby
//
//
// This module contains the following routines:
//
//  Blank_Strip - strips leading spaces from a string
//  check_atom_name - return 1 if input names match, 0 otherwise
//  expand_atom - generates all atoms of a given type within the unit cell
//  findsys - determines crystal system from symmetry information
//  find_lattice_type - interpret a set of x',y',z' symmetry operators to find lattice type
//  get_label - gets 1-4 character atom label from input string
//  getsym - obtain symmetry information from a SHELX-style SYMM card
//  Init_DRAWxtl - initialize all global variables
//  P_to_C - convert ellipsoid probability to scaling parameter
//  read_inp - read and process input ('.str') file
//  set_tf_status - sets the status information regarding temperature factors
//  skip_blocks - skip over data blocks in a CIF
//  Token_Strip - strips tokens from the beginning of a string
//  Transform_POV_Color - change input color info to POV form (if necessary)
//  Transform_VRML_Color - change input color into rgb triple needed by VRML and openGL
//  trim_string - trim strings - remove trailing spaces and control characters
//  Unique_Atom - checks that atom position is unique
//  vec_dif - find if two vectors are equal
//
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

struct tmp_color_struct
{
    char el_color_tmp[256];
};

int Block_CIF = 0;

int g_Quick;

/* ************************************************************** */
/* ************************************************************** */

void
Blank_Strip (char input[])

/* routine to suppress any leading spaces in the input string */
{
    unsigned int i, j;

    for (i = 0; i < strlen(input); i++)
	if (input[i] != ' ')
	    break;
    for (j = 0; j < strlen(input); j++) {
	input[j] = input[i];
	if (input[i++] == '\0')
	    break;
    }
}

/* ************************************************************** */
/* ************************************************************** */

int
check_atom_name (char *name1, char *name2)
{
// determine if two atom names are the same. Return 1 (TRUE) if they are

    if (name1[0] == name2[0] && name1[1] == name2[1] &&
	name1[2] == name2[2] && name1[3] == name2[3])
	return 1;
    else
	return 0;
}

/* ************************************************************** */
/* ************************************************************** */

float
convert_pos (char *string)
{
// routine to convert an atom string to a float. If only digits, a simple sscanf
//  is done. If the position is of the form n/m, each part is decoded and the
// results of the division is returned.

    float result;

    int i, n, m;

    char *pdest;

    char sub_string1[10], sub_string2[10];

    if (!(pdest = strstr (string, "/"))) {
	sscanf (string, "%f", &result);	// no / in the input
    } else {
	i = pdest - string + 1;	// number of chars to copy
	strncpy (sub_string1, string, i);
	sub_string1[i - 1] = '\0';
	strcpy (sub_string2, pdest + 1);	// copy part after /
	sscanf (sub_string1, "%d", &n);
	sscanf (sub_string2, "%d", &m);
	result = (float) n / (float) m;
    }
    return result;
}

/* ************************************************************** */
/* ************************************************************** */

int
end_flip (int value)
{
// routine to flip the endian nature of the input integer
// This routine was copied from a Web article:
//   http://www.codeproject.com/cpp/endianness.asp by Juan Carlos Cobas

    return (((value & 0x000000FF) << 24) + ((value & 0x0000FF00) << 8) +
	    ((value & 0x00FF0000) >> 8) + ((value & 0xFF000000) >> 24));
}

/* ************************************************************** */
/* ************************************************************** */

float
end_flip_real (float value)
{
// routine to flip the endian nature of the input float 
// This routine was copied from a Web article:
//   http://www.codeproject.com/cpp/endianness.asp by Juan Carlos Cobas

    union u
    {
	float vi;
	unsigned char c[4];
    };

    union v
    {
	float ni;
	unsigned char d[4];
    };

    union u un;

    union v vn;

    un.vi = value;
    vn.d[0] = un.c[3];
    vn.d[1] = un.c[2];
    vn.d[2] = un.c[1];
    vn.d[3] = un.c[0];
    return (vn.ni);

}

/* ************************************************************** */
/* ************************************************************** */

void
expand_atom (int natom)

/* routine to expand an atom position and obtain all equipoints within
   the unit cell */
/* natom - number of atom in list to expand */
{
    float xp[3], xpp[3];	/* used in generating position */

    int i, j, k, l;		/* loop variables */

    ncell = 0;			/* clear list */
    for (i = 0; i < drvui->ng; ++i) {	/* symmetry operations */
	for (j = 0; j <= 2; ++j) {
	    xp[j] = drvui->ts[i][j];
	    for (k = 0; k <= 2; ++k)
		xp[j] = xp[j] + drvui->ss[i][j][k] * drvui->atoms[natom].atom_xyz[k];
	}
	for (k = 0; k < drvui->nlat; ++k) {	/*lattice points */
	    for (l = 0; l <= 2; ++l)
		xpp[l] = xp[l] + drvui->lat_pos[k][l];
	    add_to_list (xpp, i, i);
	    if (!drvui->acentric) {
		for (l = 0; l <= 2; ++l)
		    xpp[l] = -xpp[l];
		if (i == 0)
		    add_to_list (xpp, i, -1000);
		else
		    add_to_list (xpp, i, -i);
	    }
	}
    }
}

/* ************************************************************** */
/* ************************************************************** */

void
findsys (void)
/* routine determine crystal system from symmetry information */
{
    int i, j, k;

    int n, kk = 0, nn, *symlist;

    int yzx = 0, zxy = 0, ymxmz = 0, mxmyz = 0, myxmz = 0;

    if (drvui->ng == 1) {	/* only x,y,z => triclinic */
	drvui->sys = 1;
	return;
    }

    if (drvui->ng == 2 && fabs (drvui->lat_con[4] - 90.) < 0.00001) {	/* only two symops => monoclinic */
	drvui->sys = 2;
	return;
    }

    drvui->sys = 3;		/* orthorhombic unless we find specific operators below */

    for (i = 0; i < drvui->ng; ++i) {
	for (j = 0; j <= 2; ++j) {
	    if (fabs (drvui->ts[i][j] - 0.33333333) < 0.00001)
		drvui->sys = 5;
	    if (fabs (drvui->ts[i][j] - 0.16666667) < 0.00001)
		drvui->sys = 5;
	    if (fabs (drvui->ts[i][j] - 0.66666667) < 0.00001)
		drvui->sys = 5;
	    if (fabs (drvui->ts[i][j] - 0.83333333) < 0.00001)
		drvui->sys = 5;
	}
	if (drvui->sys == 5)
	    return;
    }

    if (fabs (drvui->lat_con[5] - 120.0) < 0.00001
	&& fabs (drvui->lat_con[4] - 90.0) < 0.0001) {
	drvui->sys = 5;
	return;
    }
#if 0
    if (fabs (drvui->lat_con[5] - 54.27) < 0.1
	&& fabs (drvui->lat_con[4] - 54.27) < 0.1
	&& fabs (drvui->lat_con[3] - 54.27) < 0.1) {
	drvui->sys = 5;
	return;
    }
#endif

    n = 0;
    symlist = (int *) malloc (drvui->ng * sizeof (int));
    for (i = 0; i < drvui->ng; ++i) {
	for (j = 0; j <= 2; ++j) {
	    if (drvui->ts[i][j] != 0.0)
		break;
	    kk = 0;
	    for (k = 0; k <= 2; ++k)
		if (drvui->ss[i][j][k] != 0)
		    kk++;
	}
	if (kk > 1) {
	    drvui->sys = 5;
	    free (symlist);
	    return;
	}
	symlist[n++] = i;
    }

    for (nn = 0; nn < n; nn++) {
	i = symlist[nn];
	if (drvui->ss[i][0][1] == 1 && drvui->ss[i][1][2] == 1 && drvui->ss[i][2][0] == 1)
	    yzx = 1;
	if (drvui->ss[i][0][2] == 1 && drvui->ss[i][1][0] == 1 && drvui->ss[i][2][1] == 1)
	    zxy = 1;
	if (drvui->ss[i][0][1] == -1 &&
	    drvui->ss[i][1][2] == -1 && drvui->ss[i][2][0] == -1)
	    yzx = 1;
	if (drvui->ss[i][0][2] == -1 &&
	    drvui->ss[i][1][0] == -1 && drvui->ss[i][2][1] == -1)
	    zxy = 1;
	if (drvui->ss[i][0][1] == 1 &&
	    drvui->ss[i][1][0] == -1 && drvui->ss[i][2][2] == -1)
	    ymxmz = 1;
	if (drvui->ss[i][0][0] == -1 &&
	    drvui->ss[i][1][1] == -1 && drvui->ss[i][2][2] == 1)
	    mxmyz = 1;
	if (drvui->ss[i][0][1] == -1 &&
	    drvui->ss[i][1][0] == 1 && drvui->ss[i][2][2] == -1)
	    myxmz = 1;
    }

    if (ymxmz + mxmyz + myxmz == 3)
	drvui->sys = 4;		/*tetragonal */
    if (yzx + zxy == 2)
	drvui->sys = 6;		/*cubic */

    if (drvui->sys == 3 && fabs (drvui->lat_con[4] - 90.0) > 0.0001)
	drvui->sys = 2;
    free (symlist);
}

/* ************************************************************** */
/* ************************************************************** */

void
get_label (char input[], char *c1, char *c2, char *c3, char *c4, int strip)
// Routine to get 1 to 4 character atom ID from input line, and shift input
//   line to left, if strip true, skip over first token
{

    int i, j;

    i = 0;
    if (strip)
	while (input[i] != ' ')
	    input[i++] = ' ';	// skip over command
    while (input[i] == ' ')
	i++;			// skip whitespace after command
    *c1 = input[i];		// return first atom label character
    input[i++] = ' ';
    *c2 = input[i];		// return 2nd atom label char
    if (input[i] != ' ')
	input[i++] = ' ';	// if 2nd not blank, rub it out
    *c3 = input[i];		// return 3rd character
    if (input[i] != ' ')
	input[i++] = ' ';	// if 3rd not blank, rub it out
    *c4 = input[i];		// return 4th character
    if (input[i] != ' ')
	input[i++] = ' ';	// if 4th not blank, rub it out
    while (input[i] == ' ')
	i++;			// skip white space after atom label
    for (j = 0; i < (int) strlen (input); i++) {
	input[j++] = input[i];	// copy rest of command
    }
    input[j] = 0;		// put NULL at end
}				// end of get_label

/* ************************************************************** */
/* ************************************************************** */

void
getsym (char *text, int num, int kk)

/* helper routine to parse SHELX or CIF symmetry operator strings */
{
    int i, j;

    char *txt;

    char *n;

    char *d;

    txt = (char *) zalloc (24 * sizeof (char));
//    if (!g_Quick) fprintf(drvui->flout, "Entering getsym with text = %s, num = %d\n",text,num);

    j = 0;

    if (drvui->modulated >= 1) {
	for (i = 0; i < (int) strlen (text) - 1; i++) {
	    if (text[i] == 'X' || text[i] == 'x') {
		if (text[i + 1] == '2')
		    text[i] = 'Y';
		if (text[i + 1] == '3')
		    text[i] = 'Z';
		if (text[i + 1] == '4')
		    text[i] = 'Q';
		if (text[i + 1] == '5')
		    text[i] = 'R';
		if (text[i + 1] == '6')
		    text[i] = 'S';
		if (text[i] != 'X' && text[i] != 'x')
		    text[++i] = ' ';
	    }
	}
    }

    for (i = 0; i < (int) strlen (text); i++)
	if (text[i] != ' ' && text[i] != '\n')
	    txt[j++] = toupper (text[i]);

    txt[j] = '\0';

    if (kk < 3)
	drvui->ts[num][kk] = 0.;
    else
	drvui->ts_m[num][kk - 3] = 0.;

    n = strstr (txt, "/");
    if (n != NULL) {
	n++;
	d = n - 3;
	if (d < txt || atof (d) == 0.)
	    d = n - 2;
	if (kk < 3)
	    drvui->ts[num][kk] = (float) (atof (d) / atof (n));
	else
	    drvui->ts_m[num][kk - 3] = (float) (atof (d) / atof (n));
    }

    n = strstr (txt, ".");
    if (n != NULL) {
	if (kk < 3) {
	    drvui->ts[num][kk] = (float) atof (n);
	    if (n != txt) {
		if (txt[&txt - &n + 1] == '-')
		    drvui->ts[num][kk] *= -1;
	    }
	} else {
	    drvui->ts_m[num][kk - 3] = (float) atof (n);
	    if (n != txt) {
		if (txt[&txt - &n + 1] == '-')
		    drvui->ts_m[num][kk - 3] *= -1;
	    }
	}
    }

    if (kk < 3) {
	drvui->ss[num][kk][0] = 0;
	drvui->ss[num][kk][1] = 0;
	drvui->ss[num][kk][2] = 0;
    } else {
	drvui->ss_m[num][kk - 3][0] = 0;
	drvui->ss_m[num][kk - 3][1] = 0;
	drvui->ss_m[num][kk - 3][2] = 0;
    }

    for (i = 0; i < (int) strlen (txt); i++) {
	if (txt[i] == 'X') {
	    drvui->ss[num][kk][0] = 1;
	    if (i > 0) {
		if (txt[i - 1] == '-')
		    drvui->ss[num][kk][0] = -1;
	    }
	}
	if (txt[i] == 'Y') {
	    drvui->ss[num][kk][1] = 1;
	    if (i > 0) {
		if (txt[i - 1] == '-')
		    drvui->ss[num][kk][1] = -1;
	    }
	}
	if (txt[i] == 'Z') {
	    drvui->ss[num][kk][2] = 1;
	    if (i > 0) {
		if (txt[i - 1] == '-')
		    drvui->ss[num][kk][2] = -1;
	    }
	}
	if (txt[i] == 'Q') {
	    drvui->ss_m[num][kk - 3][0] = 1;
	    if (i > 0) {
		if (txt[i - 1] == '-')
		    drvui->ss_m[num][kk - 3][0] = -1;
	    }
	}
	if (txt[i] == 'R') {
	    drvui->ss_m[num][kk - 3][1] = 1;
	    if (i > 0) {
		if (txt[i - 1] == '-')
		    drvui->ss_m[num][kk - 3][1] = -1;
	    }
	}
	if (txt[i] == 'S') {
	    drvui->ss_m[num][kk - 3][2] = 1;
	    if (i > 0) {
		if (txt[i - 1] == '-')
		    drvui->ss_m[num][kk - 3][2] = -1;
	    }
	}
    }
    free (txt);
/*
    if (!g_Quick) {
       fprintf(drvui->flout, "ss: %d %d %d / %d %d %d / %d %d %d\n", drvui->ss[num][0][0],
            drvui->ss[num][0][1], drvui->ss[num][0][2], drvui->ss[num][1][0],
            drvui->ss[num][1][1], drvui->ss[num][1][2], drvui->ss[num][2][0],
            drvui->ss[num][2][1], drvui->ss[num][2][2]);
       fprintf(drvui->flout, "ts: %f %f %f\n", drvui->ts[num][0], drvui->ts[num][1], drvui->ts[num][2]);
       fprintf(drvui->flout, "ss_m: %d %d %d / %d %d %d / %d %d %d\n", drvui->ss_m[num][0][0],
            drvui->ss_m[num][0][1], drvui->ss_m[num][0][2], drvui->ss_m[num][1][0],
            drvui->ss_m[num][1][1], drvui->ss_m[num][1][2], drvui->ss_m[num][2][0],
            drvui->ss_m[num][2][1], drvui->ss_m[num][2][2]);
       fprintf(drvui->flout, "ts_m: %f %f %f\n", drvui->ts_m[num][0], drvui->ts_m[num][1], drvui->ts_m[num][2]);
    }
*/
}

/* ************************************************************** */
/* ************************************************************** */

int
get_next_token (char *p, int max_len, FILE * fpin)
{
/* routine to return the next token in a stream. A token is the string
 * characters to the next white space consisting of a space, a tab, or
 * a newline character. Carriage returns are ignored. The returned value
 * is 1 (true) if a token is found and 0 (false) if the end of file is
 * reached */

    char *pp = p;

    int inchar;

    int started = 0;

    int openquote = -1;

    int len = max_len - 2;

    for (;;) {			/* start an infinite loop */
	if ((inchar = getc (fpin)) == EOF)
	    break;		/* here for end-of-file */
	if (inchar == '\t' || inchar == ' ') {
	    if (started && openquote != 1) {
		*pp = '\0';	/* add null terminator */
		return 1;
	    }
	} else if (inchar == '\n') {
	    if (started) {
		*pp = '\0';
		return 1;
	    }
	} else if (inchar == '#' && !started) {	// comment - skip to end of line
	    char comment[256];
	    char *c = comment;
	    for (;;) {
		if ((inchar = getc (fpin)) == EOF)
		    return 0;	/* here for end-of-file */
		*c++=(char)inchar;
		if (inchar == '\n') {
		    *c='\0';
		    if (!strncmp(comment," End of data for", 16)) 
			return 0;
		    break;
		}
	    }
	} else if (inchar != '\r') {
	    *pp++ = (char) inchar;
	    started = 1;
	    if (inchar == (char) 39)
		openquote *= -1;
	    if (!--len)
		return 0;
	}
    }
    return 0;
}

/* ************************************************************** */
/* ************************************************************** */

void
find_lattice_type (void)
{
/* routine to scan symmetry elements, find the lattice type, and fill in the coorcdinates
   of the lattice points */

    int i;

    drvui->nbr = 1;

    for (i = 1; i <= drvui->ng; ++i) {
	if (drvui->ss[i - 1][0][0] == 1 && drvui->ss[i - 1][1][1] == 1 && drvui->ss[i - 1][2][2] == 1) {	/* rotational part of symmetry is x,y,z - centering? */
	    if (drvui->ts[i - 1][0] == 0. && drvui->ts[i - 1][1] == 0.5 && drvui->ts[i - 1][2] == 0.5) {	/* found x,1/2+y,1/2+z */
		if (drvui->nbr == 3 || drvui->nbr == 4)
		    drvui->nbr = 5;	/* x,1/2+y,1/2+z && B || C ==> F */
		else
		    drvui->nbr = 2;	/* x,1/2+y,1/2+z ==> A */
	    }
	    if (drvui->ts[i - 1][0] == 0.5 && drvui->ts[i - 1][1] == 0. && drvui->ts[i - 1][2] == 0.5) {	/* found 1/2+x,y,1/2+z */
		if (drvui->nbr == 2 || drvui->nbr == 4)
		    drvui->nbr = 5;	/* 1/2+x,y,1/2+z && A || C ==> F */
		else
		    drvui->nbr = 3;	/* 1/2+x,y,1/2+z ==> B */
	    }
	    if (drvui->ts[i - 1][0] == 0.5 && drvui->ts[i - 1][1] == 0.5 && drvui->ts[i - 1][2] == 0.) {	/* found 1/2+x,1/2+y,z */
		if (drvui->nbr == 2 || drvui->nbr == 3)
		    drvui->nbr = 5;	/* 1/2+x,1/2+y,z && A || B ==> F */
		else
		    drvui->nbr = 4;	/* 1/2+x,1/2+y,z ==> C */
	    }
	    if (drvui->nbr == 1 && drvui->ts[i - 1][0] == 0.5
		&& drvui->ts[i - 1][1] == 0.5 && drvui->ts[i - 1][2] == 0.5)
		drvui->nbr = 6;	// 0.5+X,0.5+Y,0.5+Z ==> I

	    if (drvui->ts[i - 1][0] < 0.5 && drvui->ts[i - 1][0] > 0.25)
		drvui->nbr = 7;	// 0.33+X ==> R
	}

    }

    switch (drvui->nbr) {
    case 1:
	drvui->nlat = 1;
	break;
    case 2:
	drvui->nlat = 2;
	drvui->lat_pos[1][0] = 0.0;
	drvui->lat_pos[1][1] = 0.5f;
	drvui->lat_pos[1][2] = 0.5f;
	break;
    case 3:
	drvui->nlat = 2;
	drvui->lat_pos[1][0] = 0.5f;
	drvui->lat_pos[1][1] = 0.0;
	drvui->lat_pos[1][2] = 0.5f;
	break;
    case 4:
	drvui->nlat = 2;
	drvui->lat_pos[1][0] = 0.5f;
	drvui->lat_pos[1][1] = 0.5f;
	drvui->lat_pos[1][2] = 0.0;
	break;
    case 5:
	drvui->nlat = 4;
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
    case 6:
	drvui->nlat = 2;
	drvui->lat_pos[1][0] = 0.5f;
	drvui->lat_pos[1][1] = 0.5f;
	drvui->lat_pos[1][2] = 0.5f;
	break;
    case 7:
	drvui->nlat = 3;
	drvui->lat_pos[1][0] = (1.f / 3.f);
	drvui->lat_pos[1][1] = (2.f / 3.f);
	drvui->lat_pos[1][2] = (2.f / 3.f);
	drvui->lat_pos[2][0] = (2.f / 3.f);
	drvui->lat_pos[2][1] = (1.f / 3.f);
	drvui->lat_pos[2][2] = (1.f / 3.f);
    }
}


/* ************************************************************** */
/* ************************************************************** */

/* ************************************************************** */
/* ************************************************************** */

void
Init_DRAWxtl (void)
{
    int j, l, m;

/* Initialize */
    drvui->X_Origin = drvui->Y_Origin = drvui->Z_Origin = 0.5f;
    drvui->X_Boxlim = drvui->Y_Boxlim = drvui->Z_Boxlim = 20.0f;
    drvui->Trans[0] = drvui->Trans[1] = drvui->Trans[2] = 0.0f;
    drvui->automation =0;
    xrot = yrot = zrot = 0.0f;
    boxflag = packflag = clipflag = 0;
    drvui->do_ellipsoids = edges = 0;
    drvui->Ellipsoid_Prob = 0.5f;
    cur_cen[0] = 0.25f;
    cur_cen[1] = 0.25f;
    cur_cen[2] = 0.25f;		// location of cursor (fractional coordinates)
    drvui->cur_step = 0.5f;	// amount to jump cursor in A
    Unit_Cell = M_cameras = no_comment = Labels = Vrml2 = 1;
    X3D = 0;
    Display_axes = drvui->El_Cutout = 0;
    strcpy (drvui->Cutout_color, "Gray20");
    strcpy (drvui->Ellipaxis_color, "Gray20");
    drvui->Ellipaxis_width = 0.;
    drvui->SpMult = drvui->BndMult = Magnification = 1.0f;
    drvui->Ellipsoid_Scale = P_to_C (0.5f);
    Options = 0;
    FourierMapType = 0;
    ShowMapLegend = 0;
    drvui->Fourier2d = 0;
    Map_Info.info_valid = 0;
    Map_Info.map_type = -1;
    Map_Info.res = 4;
    drvui->sys = 0;
    drvui->polylimit = 0.1f;
    drvui->Phong_Value = 0.2f;
    drvui->Phong_Size = 1.0f;
    DepthCue = 0.0f;
    drvui->noshadow = 0;
    drvui->rad_edge = 0.0f;
/*
  drvui->ambient = 0.0f;
  drvui->diffuse = 0.0f;
  drvui->specular = 0.0f;
  drvui->roughness = 0.0f;
*/
    printdist = 3.5f;
    origin[0] = drvui->X_Origin;
    origin[1] = drvui->Y_Origin;
    origin[2] = drvui->Z_Origin;
    boxlim[0] = drvui->X_Boxlim;
    boxlim[1] = drvui->Y_Boxlim;
    boxlim[2] = drvui->Z_Boxlim;
    docell = Unit_Cell;
    Color_Warning = 0;
    drvui->El_Cutout = 0;
    strcpy (drvui->Cutout_color, "");
    rad_cell = 0.02f;
    drvui->Sphere_Mult = drvui->SpMult;
    drvui->Bond_Mult = drvui->BndMult;
    natom = 0;
    drvui->nmag = 0;
    drvui->nedges = 0;
    drvui->nsphere = drvui->npoly = drvui->nbond = drvui->nplane = 1;
    drvui->n_ellips = drvui->ncone = drvui->nlabel = drvui->nbplane = 1;
    drvui->nsurf = drvui->natprop = 1 ;
    drvui->nlabel = 1;
//  POV_Max[0] = POV_Max[1] = POV_Max[2] = -99999.0f;
//  POV_Min[0] = POV_Min[1] = POV_Min[2] = 99999.0f;
    drvui->glback[0] = drvui->glback[1] = drvui->glback[2] = 1.0f;
    slabmode = 0;
    domolcomp = 0;
    drvui->mol_d = 0.0f;
    drvui->modulated = 0;
    for (j = 0; j < 3; j++)
	drvui->phaseshift[j] = 0.0;
    Omit->nomits = 0;
    drvui->no_subsys = 1;	/* initialize the subsystem variables */
    for (l = 0; l < 3; l++) {
	for (m = 0; m < 3; m++)
	    drvui->subsys_fact[0][l][m] = 0.0f;
	drvui->subsys_fact[0][l][l] = 1.0f;	/* set to identity */
    }
    drvui->subsys_ref_volume = 1.0f;
    drvui->subsys_vol[0] = 1.0;
    memset (Omit->omit1, 0, 1000);
    memset (Omit->omit2, 0, 1000);
    memset (drvui->slab_con, 0, sizeof (drvui->slab_con));
    memset (drvui->slab_rot, 0, sizeof (drvui->slab_rot));
    memset (drvui->slab_off, 0, sizeof (drvui->slab_off));
// clear all the color variables
    memset (drvui->col_cell, 0, sizeof (drvui->col_cell));
    memset (drvui->col_edge, 0, sizeof (drvui->col_edge));
    memset (drvui->col_bg, 0, sizeof (drvui->col_bg));
    memset (drvui->Cutout_color, 0, sizeof (drvui->Cutout_color));
    memset (drvui->col_edge, 0, sizeof (drvui->col_edge));

    strcpy (drvui->col_cell, "Black");	/* preset unit-cell colors to Black */
    strcpy (drvui->col_bg, "White");
    strcpy (Map_Info.title, " Title not given in map file");
    if (drvui->voidflag > 0 && drvui->voidmap) {
	for (int j = 0; j < drvui->voidgrid[1]; j++) {
	    for (int k = 0; k < drvui->voidgrid[2]; k++)
		free (drvui->voidmap[j][k]);
	    free (drvui->voidmap[j]);
	}
	free (drvui->voidmap);
	drvui->voidmap = NULL;
    }
    drvui->voidflag = 0;
}

/* ************************************************************** */
/* ************************************************************** */

float
P_to_C (float prob)

/* Routine to convert probability to ellipsoid scaling, C.  Given probability
   'prob', calculate value of C required for ellipsoid,

   C^2 = (delx/sigx)^2 + (dely/sigy)^2 + (delz/sigz)^2,

   where 'prob' is the probability that a random point in the distribution
   will fall inside the ellipsoid.  Table values taken from D.B. Owen,
   "Handbook of Statistical Tables", Addison-Wesley, Reading, Mass., 1962.

   N.B. The present version of the tables only handles probabilities from
   0.01 to 0.99.
*/
#define N_C 99
{
    static double C_table[N_C + 1] =
	{ 0.3389, 0.4299, 0.4951, 0.5479, 0.5932, 0.6334, 0.6699, 0.7035,
	0.7349, 0.7644,
	0.7924, 0.8192, 0.8447, 0.8694, 0.8932, 0.9162, 0.9386, 0.9605,
	0.9818, 1.0026,
	1.0230, 1.0430, 1.0627, 1.0821, 1.1012, 1.1200, 1.1386, 1.1570,
	1.1751, 1.1932,
	1.2110, 1.2288, 1.2464, 1.2638, 1.2812, 1.2985, 1.3158, 1.3330,
	1.3501, 1.3672,
	1.3842, 1.4013, 1.4183, 1.4354, 1.4524, 1.4695, 1.4866, 1.5037,
	1.5209, 1.5382,
	1.5555, 1.5729, 1.5904, 1.6080, 1.6257, 1.6436, 1.6616, 1.6797,
	1.6980, 1.7164,
	1.7351, 1.7540, 1.7730, 1.7924, 1.8119, 1.8318, 1.8519, 1.8724,
	1.8932, 1.9144,
	1.9360, 1.9580, 1.9804, 2.0034, 2.0269, 2.0510, 2.0757, 2.1012,
	2.1274, 2.1544,
	2.1824, 2.2114, 2.2416, 2.2730, 2.3059, 2.3404, 2.3767, 2.4153,
	2.4563, 2.5003,
	2.5478, 2.5997, 2.6571, 2.7216, 2.7955, 2.8829, 2.9912, 3.1365,
	3.3682, 5.9503
    };

    int i;

    double C;

    i = (int) (100.0 * prob);	/* find correct table entry */
    if (i > N_C)
	i = N_C;
/* interpolate (linear) */
    C = C_table[i - 1] + (prob - 0.01 * i) * (C_table[i] - C_table[i - 1]);
    return ((float) C);
}


/* ************************************************************** */
/* ************************************************************** */

void
read_inp (int Quick)

/* routine to read and process input file */
{
    char input[256];		/* input buffer for data read from file */

    char input2[256];		/* input buffer for data read from file */

    char string[256];		/* temporary string */

    char filename[256];		/* used for map file name */

    int intype;			/* type of datum read */

    int i, j, k;		/* temporary */

    char t_name[5];		/* temporary atom name storage */

    char t_color[40];		/* temporary ellipsoid color storage */

    int in_line = 0;		/* flag for inline import of foreign files */

    float t_rad;

    int n_el_tmp = 0;

    int done;

    int tmp_frame_no = 1;	// number of frame being processed

    FILE *tout = 0;

    int curframe;

    drvui->lat_con[0] = 0.0f;

    char atom_pos[3][30];

    int average = 0;

    char modl[5];

    int modnum;

    float avg_occ, min_occ;

    struct tmp_color_struct *el_color_tmp;

/* ************************************************************** */
/* ************************************************************** */
/*
    Commands Implemented

    arrow
    atom
    average
    axislines
    background
    bestplane
    betaij
    bij,Bij
    bond
    bounds
    box
    cell
    clip
    cutout
    dash
    depthcue
    edges
    ellipcolor
    ellipsoids
    finish
    frame
    import
    inline
    labelscale
    labeltext
    list
    lonepair
    lookat
    magtrans
    magnification
    mapcalclimits
    mapcontour
    mapcontour2d
    maplegend
    mapread
    mapregion
    mapslice
    molcomp
    nolabels
    noshadow
    occupancy
    origin
    orthographic
    pack
    phong
    plane
    polyedge
    polysz
    polytolerance
    polylimit
    polyfudge
    polyvert
    rem, REM
    shell
    slab
    special
    spgp, sgrp, spgr
    sphere
    title, titl
    uij,Uij
    values
    vectors    
    view
    voids
    vrml1
    vrml2
    vrml97
    x3d
    xyzoff
    end, END
                                                                  */
/* ************************************************************** */
/* ************************************************************** */

    el_color_tmp =
	(tmp_color_struct *) zalloc (drvui->ellips_alloc *
				     sizeof (struct tmp_color_struct));
    if (!el_color_tmp) {
	Error_Box ("Unable to allocate temporary space.\n");
	exit (0);
    }

    drvui->xyzoff_read = 0;
    drvui->voidflag = 0;
    Map_Info.xlim[0] = 0.0f;
    Map_Info.xlim[1] = 1.0f;
    Map_Info.ylim[0] = 0.0f;
    Map_Info.ylim[1] = 1.0f;
    Map_Info.zlim[0] = 0.0f;
    Map_Info.zlim[1] = 1.0f;
    Map_Info.x4lim[0] = 0.0f;
    Map_Info.x4lim[1] = 0.0f;
    Map_Info.x5lim[0] = 0.0f;
    Map_Info.x5lim[1] = 0.0f;
    Map_Info.x6lim[0] = 0.0f;
    Map_Info.x6lim[1] = 0.0f;
    drvui->numOfFourierContours = 0;
    drvui->frames[drvui->frame_no].slice = 0;
    curframe = drvui->frame_no;
    if (curframe > 1) {
	if (packflag)
	    for (i = 0; i < 6; i++) {
		drvui->frames[curframe].cryst_lim[i] =
		    drvui->frames[curframe - 1].cryst_lim[i];
	} else
	    for (i = 3; i < 6; i++) {
		drvui->frames[curframe].cryst_lim[i] = 1.;
	    }
    }

    for (i = 0; i < 3; i++)
	drvui->xyzoff[i] = 0.0f;
    while (!feof (drvui->fpin)) {
	strcpy (t_color, drvui->ellips[0].ellips_l);
	memset (t_color, 0, 40);
	memset (string, 0, 256);
	memset (input, 0, sizeof (input));
	if (fgets (input, 256, drvui->fpin)) {	/* read a line */
	    trim_string (input, 256);	// kill any ^M or ^J at end of line
	    Blank_Strip (input);	// suppress any leading blanks
	    strcpy (input2, input);	// make a copy for listing file
	    for (i = strlen (input); i < 255; i++) {	//pad line with blanks
		input[i] = ' ';
	    }
	    input[255] = 0;	// force a terminator
	} else {
	    Error_Box ("Error reading file.");
	    free (el_color_tmp);
	    return;
	}
	intype = 1000;
	in_line = 0;
	if (!strncmp (input, "titl", 4) || !strncmp (input, "TITL", 4)) {
	    if (drvui->frame_no == 1)
		intype = 1;	// only echo title in first frame
	    else
		intype = 0;
	}
	if (strncmp (input, "rem", 3) == 0 || strncmp (input, "REM", 3) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 1;	// only echo rem in its frame
	    else
		intype = 0;
	}
	if (strncmp (input, "cell", 4) == 0)
	    intype = 2;
	if (strncmp (input, "sgrp", 4) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 3;
	    else
		intype = 0;	// ignore sgrp command if not current frame
	}
	if (strncmp (input, "spgp", 4) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 3;
	    else
		intype = 0;	// ignore spgp command if not current frame
	}
	if (strncmp (input, "spgr", 4) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 3;
	    else
		intype = 0;	// ignore spgr command if not current frame
	}
	if (strncmp (input, "ellipsoids", 10) == 0)
	    intype = 4;
	if (strncmp (input, "lookat", 6) == 0)
	    intype = 5;
	if (strncmp (input, "atom", 4) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 6;
	    else
		intype = 0;	// ignore atom command if not current frame
	}
	if (strncmp (input, "polysz", 6) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 7;
	    else
		intype = 0;	// ignore polysz command if not current frame
	}
	if (strncmp (input, "sphere", 6) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 8;
	    else
		intype = 0;	// ignore sphere command if not current frame
	}
	if (strncmp (input, "end", 3) == 0)
	    intype = 9;
	if (strncmp (input, "END", 3) == 0)
	    intype = 9;
	if (strncmp (input, "xyzoff", 6) == 0)
	    intype = 10;
	if (strncmp (input, "bond", 4) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 11;
	    else
		intype = 0;	// ignore bond command if not current frame
	}
	if (strncmp (input, "box", 3) == 0)
	    intype = 12;
	if (strncmp (input, "edges", 5) == 0)
	    intype = 13;
	if (strncmp (input, "view", 4) == 0)
	    intype = 14;
	if (strncmp (input, "plane", 5) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 15;
	    else
		intype = 0;	// ignore plane command if not current frame
	}
	if (strncmp (input, "bij", 3) == 0)
	    intype = 16;
	if (strncmp (input, "Bij", 3) == 0)
	    intype = 16;
	if (strncmp (input, "betaij", 6) == 0)
	    intype = 17;
	if (strncmp (input, "uij", 3) == 0)
	    intype = 18;
	if (strncmp (input, "Uij", 3) == 0)
	    intype = 18;
	if (strncmp (input, "bounds", 6) == 0)
	    intype = 1;		//19;
	if (strncmp (input, "pack", 4) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 20;
	    else
		intype = 0;	// ignore pack command if not current frame
	}
	if (strncmp (input, "magnification", 13) == 0)
	    intype = 21;
	if (strncmp (input, "origin", 6) == 0)
	    intype = 22;
	if (strncmp (input, "ellipcolor", 10) == 0)
	    intype = 23;
	if (strncmp (input, "inline", 6) == 0) {
	    if (tmp_frame_no == drvui->frame_no) {
		intype = 24;
	    } else {		// skip ahead to next frame or end
		int end_seen = 0;	// some formats have their own END card

		while (!feof (drvui->fpin) && intype != 0) {
		    if (fgets (input, 256, drvui->fpin)) {	/* read a line */
			trim_string (input, 256);	// kill any ^M or ^J at end of line
			Blank_Strip (input);	// suppress any leading blanks
			input[255] = 0;	// force a terminator
			if (strncmp (input, "frame", 5) == 0) {
			    intype = 37;
			    break;
			}
			if (strncmp (input, "end", 3) == 0) {
			    if (strstr (input2, "schakal") && end_seen == 0) {
				end_seen = 1;
			    } else {
				intype = 9;
				break;
			    }
			}
		    } else {
			Error_Box ("Error reading file.");
			free (el_color_tmp);
			return;
		    }
		}		// advance to next frame or end command
	    }			// if not in current frame 
	}
	if (strncmp (input, "import", 6) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 25;
	    else
		intype = 0;	// ignore import command if not current frame
	}
	if (strncmp (input, "orthographic", 12) == 0)
	    intype = 26;
	if (strncmp (input, "vectors", 7) == 0) {
	    if (tmp_frame_no == 1) {
		intype = 27;
	    } else {
		intype = 0;
	    }
	}
	if (strncmp (input, "phong", 5) == 0)
	    intype = 28;
	if (strncmp (input, "molcomp", 7) == 0)
	    intype = 29;
	if (strncmp (input, "cutout", 6) == 0)
	    intype = 30;
	if (strncmp (input, "depthcue", 8) == 0)
	    intype = 31;
	if (strncmp (input, "vrml1", 5) == 0)
	    intype = 32;
	if (strncmp (input, "vrml2", 5) == 0)
	    intype = 1;
	if (strncmp (input, "vrml97", 6) == 0)
	    intype = 1;
	if (strncmp (input, "polytol", 7) == 0)
	    intype = 33;
	if (strncmp (input, "polylimit", 9) == 0)
	    intype = 33;
	if (strncmp (input, "polyfudge", 9) == 0)
	    intype = 33;
	if (strncmp (input, "clip", 4) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 34;
	    else
		intype = 0;	// ignore clip command if not current frame
	}
	if (strncmp (input, "list", 4) == 0)
	    intype = 35;
	if (strncmp (input, "polyedge", 8) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 36;
	    else
		intype = 0;
	}
	if (strncmp (input, "frame", 5) == 0)
	    intype = 37;
	if (strncmp (input, "axislines", 9) == 0)
	    intype = 38;
	if (strncmp (input, "shell", 5) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 39;
	    else
		intype = 0;	// ignore shell command if not current frame
	}
	if (strncmp (input, "nolabels", 8) == 0)
	    intype = 40;
	if (strncmp (input, "lonepair", 8) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 41;
	    else
		intype = 0;	// ignore lonepair command if not current frame
	}
	if (strncmp (input, "labeltext", 9) == 0)
	    intype = 42;
	if (strncmp (input, "polyvert", 8) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 43;
	    else
		intype = 0;	// ignore polyvert command if not current frame
	}
	if (strncmp (input, "background", 10) == 0)
	    intype = 44;
	if (strncmp (input, "slab", 4) == 0)
	    intype = 45;
	if (strncmp (input, "mag_trans", 9) == 0)
	    intype = 46;
	if (strncmp (input, "special", 7) == 0)
	    intype = 47;
	if (strncmp (input, "noshadow", 8) == 0)
	    intype = 48;
	if (strncmp (input, "finish", 6) == 0)
	    intype = 49;
	if (strncmp (input, "arrow", 5) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 50;
	    else
		intype = 0;	// ignore arrow command if not current frame
	}
	if (strncmp (input, "dash", 4) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 51;
	    else
		intype = 0;	// ignore bond command if not current frame
	}
	if (strncmp (input, "mapcontour", 10) == 0)
	    intype = 52;
	if (strncmp (input, "mapcontour2d", 12) == 0)
	    intype = 61;
	if (strncmp (input, "mapread", 7) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 53;
	    else
		intype = 0;
	}
	if (strncmp (input, "mapregion", 9) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 54;
	    else
		intype = 0;
	}
	if (strncmp (input, "mapcalclimits", 13) == 0)
	    intype = 55;
	if (strncmp (input, "labelscale", 10) == 0)
	    intype = 56;
	if (strncmp (input, "bestplane", 9) == 0)
	    intype = 57;
	if (strncmp (input, "average", 7) == 0)
	    intype = 58;
	if (strncmp (input, "occupancy", 9) == 0)
	    intype = 59;
	if (strncmp (input, "phaseshift", 10) == 0)
	    intype = 60;
	if (strncmp (input, "aimsurf", 7) == 0)
	    intype = 62;
	if (strncmp (input, "voids", 5) == 0)
	    intype = 63;
	if (strncmp (input, "values", 6) == 0)
	    intype = 64;
	if (strncmp (input, "mapslice", 8) == 0) {
	    if (tmp_frame_no == drvui->frame_no)
		intype = 65;
	    else
		intype = 0;
	}
	if (strncmp (input, "qvector", 7) == 0)
	    intype = 66;
	if (strncmp (input, "maplegend", 9) == 0)
	    intype = 67;
	if (strncmp (input, "x3d", 3) == 0)
	    intype = 68;

	switch (intype) {

	case 0:
	    break;
	case 1:
	    break;		// ignore title except for listing
	case 2:		// cell command
	    for (i = 3; i <= 5; ++i)
		drvui->lat_con[i] = 0.0f;
	    sscanf (input, "%s %f %f %f %f %f %f", string, &drvui->lat_con[0]
		    , &drvui->lat_con[1], &drvui->lat_con[2], &drvui->lat_con[3],
		    &drvui->lat_con[4]
		    , &drvui->lat_con[5]);
	    break;
	case 3:
	    symop (input);	/* process space group information */
	    break;
	case 4:		/* ellipsoids command */
	    drvui->do_ellipsoids = 1;
	    (void) sscanf (input, "%s %f", string, &drvui->Ellipsoid_Prob);
	    if (drvui->Ellipsoid_Prob > 1)
		drvui->Ellipsoid_Prob *= 0.01f;	/* Probability to fraction */
	    break;
	case 5:		// read vectors for lookat command
	    sscanf (input, "%s %f %f %f %f %f %f", string, &drvui->lookat_v1[0],
		    &drvui->lookat_v1[1], &drvui->lookat_v1[2], &drvui->lookat_v2[0],
		    &drvui->lookat_v2[1], &drvui->lookat_v2[2]);
	    Options = Options | L_OPT;
	    break;
	case 6:		/* atom command */
	    char t_atom[5];

	    get_label (input, &t_atom[0], &t_atom[1], &t_atom[2], &t_atom[3], 1);
	    t_atom[4] = 0;
	    strcpy (drvui->atoms[natom].atom_l, t_atom);
	    (void) sscanf (input, "%d %s %s %s", &drvui->atoms[natom].atom_n, atom_pos[0],
			   atom_pos[1], atom_pos[2]);
	    drvui->atoms[natom].sv_atom_n = drvui->atoms[natom].atom_n;
	    for (i = 0; i < 3; i++)
		drvui->atoms[natom].atom_xyz[i] = convert_pos (atom_pos[i]);
	    drvui->atoms[natom].atom_fn = tmp_frame_no;
	    natom++;
	    check_dynamic_storage ();
	    break;
	case 7:		/* polysz command */
	case 39:		/* shell command */
	case 43:		/* polyvert command */
	    get_label (input, &drvui->polyhedra[drvui->npoly].poly_l[0],
		       &drvui->polyhedra[drvui->npoly].poly_l[1]
		       , &drvui->polyhedra[drvui->npoly].poly_l[2],
		       &drvui->polyhedra[drvui->npoly].poly_l[3], 1);
	    if (intype == 7) {
		(void) sscanf (input, "%f %39c",
			       &drvui->polyhedra[drvui->npoly].poly_size, string);
	    } else if (intype == 39) {
		(void) sscanf (input, "%f %f %39c",
			       &drvui->polyhedra[drvui->npoly].poly_min,
			       &drvui->polyhedra[drvui->npoly].poly_size, string);
	    } else if (intype == 43) {
		get_label (input, &drvui->polyhedra[drvui->npoly].poly_t[0],
			   &drvui->polyhedra[drvui->npoly].poly_t[1]
			   , &drvui->polyhedra[drvui->npoly].poly_t[2],
			   &drvui->polyhedra[drvui->npoly].poly_t[3], 0);
		(void) sscanf (input, "%f %39c",
			       &drvui->polyhedra[drvui->npoly].poly_size, string);
	    }
	    trim_string (string, 40);
	    if (!strlen (string))
		strcpy (string, "White");
	    strcpy (drvui->polyhedra[drvui->npoly].poly_col, string);
	    if (intype == 7 || intype == 43)
		drvui->polyhedra[drvui->npoly].poly_min = 0.005f;
	    drvui->polyhedra[drvui->npoly].poly_fn = tmp_frame_no;
	    if (intype == 7 || intype == 39)
		strcpy (drvui->polyhedra[drvui->npoly].poly_t, "");
	    drvui->polyhedra[drvui->npoly].poly_rad_edge = 0.0f;
	    drvui->npoly++;
	    check_dynamic_storage ();
	    break;
	case 8:		/* sphere command */
	    get_label (input, &drvui->spheres[drvui->nsphere].sphere_l[0],
		       &drvui->spheres[drvui->nsphere].sphere_l[1]
		       , &drvui->spheres[drvui->nsphere].sphere_l[2],
		       &drvui->spheres[drvui->nsphere].sphere_l[3], 1);
	    char str[9][40], *pnt;

	    double inpval[9];

	    i = sscanf (input, "%s %s %s %s %s %s %s %s %s", str[0], str[1], str[2],
			str[3], str[4], str[5], str[6], str[7], str[8]);
	    for (j = 0; j < i; j++)
		inpval[j] = strtod (str[j], &pnt);
	    drvui->spheres[drvui->nsphere].sphere_n = -1;
	    if (i == 2 || i == 4 || i == 6) {
		drvui->spheres[drvui->nsphere].sphere_size = (float) inpval[0];
		j = 1;
	    } else if (i == 3 || i == 5 || i == 7) {
		drvui->spheres[drvui->nsphere].sphere_n = (int) (inpval[0] + 0.01);
		drvui->spheres[drvui->nsphere].sphere_size = (float) inpval[1];
		j = 2;
	    } else {
		Error_Box ("Improperly formed sphere command.");
	    }
	    strcpy (string, str[j]);
	    for (k = j + 1; k < i; k++) {
		strcat (string, " ");
		strcat (string, str[k]);
	    }
	    trim_string (string, 40);
	    if (!strlen (string))
		strcpy (string, "White");
	    strcpy (drvui->spheres[drvui->nsphere].sphere_col, string);
	    drvui->spheres[drvui->nsphere].sphere_fn = tmp_frame_no;
	    drvui->nsphere++;
	    check_dynamic_storage ();
	    break;
	case 9:		/* end or END */
	    if (!drvui->lat_con[0]) {
		if (errorbox) {	// make sure that the original error message gets a chance to appear
		    int n;

		    for (n = 0; n < 10; n++)
			Fl::wait (1);
		}
		Error_Box ("There is no 'cell' line in your input file.");
	    }
	    label_cell ();
	    if (Quick == 2) {
		drvui->max_frame = tmp_frame_no;
		drvui->Old_Xrot = xrot;
		drvui->Old_Yrot = yrot;
		drvui->Old_Zrot = zrot;
	    }
	    if (average == 1)
		drvui->modulated = -1;
	    int frame;

	    for (frame = 1; frame <= drvui->max_frame; frame++) {
		if (!Quick) {
		    strcpy (input, drvui->Cur_Root);
		    sprintf (input2, ".frm%d", frame);
		    strcat (input, input2);
		    tout = fopen (input, "w");
		}
		for (i = 0; i < natom; i++) {	/* shift atom coordinates for origin offset */
		    if (drvui->atoms[i].atom_fn != frame)
			continue;
		    drvui->atoms[i].TF_status = -1;	// initialize TF type/status
		    for (j = 0; j < 3; j++)
			drvui->atoms[i].atom_xyz[j] -= drvui->xyzoff[j];
		    strncpy (t_atom, drvui->atoms[i].atom_l, 4);
		    t_atom[4] = 0;
		    if (!Quick)
			fprintf (tout, "%s %d\n", t_atom, drvui->atoms[i].sv_atom_n);
		}
		if (!Quick) {
		    fclose (tout);
		}
	    }
	    if (n_el_tmp > 0) {
		strcpy (input, "");
		for (i = 0; i < n_el_tmp; i++) {
		    get_label (el_color_tmp[i].el_color_tmp, &t_name[0], &t_name[1],
			       &t_name[2], &t_name[3], 0);
		    t_name[4] = 0;
		    if (el_color_tmp[i].el_color_tmp[0] == '*') {
			(void) sscanf (el_color_tmp[i].el_color_tmp, "%s %39c", string,
				       t_color);
			j = -1;
		    } else {
			(void) sscanf (el_color_tmp[i].el_color_tmp, "%d %39c", &j,
				       t_color);
		    }
		    trim_string (t_color, 40);
		    if (!strlen (t_color))
			strcpy (t_color, "White");
		    for (k = 1; k < drvui->n_ellips; k++) {	/* find this atom in ellipsoid list */
			if (check_atom_name (t_name, drvui->ellips[k].ellips_l)) {
			    if ((j == -1) || (j == drvui->ellips[k].ellips_n)) {
				drvui->ellips[k].ell_type += 1000;	// >= 1000 means ellipcolor read
				strcpy (drvui->ellips[k].ellips_col, t_color);	/* copy color */
				drvui->ellips[k].save_el_number = j;
			    }
			}
		    }
		}
		n_el_tmp = 0;
	    }
	    free (el_color_tmp);

	    if (drvui->natprop > 0) {

		for (i = 0; i < drvui->natprop; i++) {
		    j = drvui->atprops[i].atprop_n;
		    for (k = 0; k < natom; k++) {	/* find this atom in list */
			if (check_atom_name (drvui->atprops[i].atprop_l, drvui->atoms[k].atom_l)) {
			    if ((j == -1) || (j == drvui->atoms[k].atom_n)) {
				drvui->atoms[k].radius = drvui->atprops[i].radius;
			    }
			}
		    }
		}
	    }

	    for (i = 1; i < drvui->nsphere; i++) {	// remove any auto-generated ellipsoids that
		for (j = 1; j < drvui->n_ellips; j++) {	//     interfere with spheres
		    if (drvui->ellips[j].ell_type == 1001) {
			if (check_atom_name
			    (drvui->spheres[i].sphere_l, drvui->ellips[j].ellips_l)) {
			    if (drvui->spheres[i].sphere_n == -1
				|| drvui->spheres[i].sphere_n ==
				drvui->ellips[j].save_el_number) {
				drvui->ellips[j].ell_type = 1;
			    }
			}
		    }
		}
	    }
	    if (FourierMapType)
		drvui->mainWindow->cursor (FL_CURSOR_WAIT);
	    switch (FourierMapType) {
	    case 1:
		read_grd (FourierFileName, Quick);
		break;
	    case 2:
		read_stf (FourierFileName, Quick);
		break;
	    case 3:
		read_w2k (FourierFileName, Quick);
		break;
	    case 4:
		read_vasp (FourierFileName, Quick);
		break;
	    case 5:
		read_flp (FourierFileName, Quick);
		break;
	    case 6:
		read_fcf (FourierFileName, Quick);
		break;
	    case 7:
		read_dn6 (FourierFileName, Quick);
		break;
	    case 8:
		read_m80 (FourierFileName, Quick);
		break;
	    case 9:
		read_exc (FourierFileName, Quick);
		break;
	    case 10:
		read_m81 (FourierFileName, Quick);
		break;
	    case 11:
		read_xsf (FourierFileName, Quick);
		break;
	    default:
		break;
	    }
	    if (FourierMapType)
		drvui->mainWindow->cursor (FL_CURSOR_DEFAULT);
	    for (frame = 1; frame <= drvui->max_frame; frame++) {
		if (drvui->frames[frame].map_lim_set == 0) {	// mapregion line NOT read
		    xMin = Map_Info.xlim[0];
		    xMax = Map_Info.xlim[1];
		    yMin = Map_Info.ylim[0];
		    yMax = Map_Info.ylim[1];
		    zMin = Map_Info.zlim[0];
		    zMax = Map_Info.zlim[1];
		    x4Val = x5Val = x6Val = 0.0f;
		    drvui->frames[frame].map_lim[0]=xMin;
		    drvui->frames[frame].map_lim[3]=xMax;
		    drvui->frames[frame].map_lim[1]=yMin;
		    drvui->frames[frame].map_lim[4]=yMax;
		    drvui->frames[frame].map_lim[2]=zMin;
		    drvui->frames[frame].map_lim[5]=zMax;
		    drvui->frames[frame].map_lim[6]=x4Val;
		    drvui->frames[frame].map_lim[7]=x5Val;
		    drvui->frames[frame].map_lim[8]=x6Val;
		}
		if (drvui->frames[frame].slice >0 ) { // calculate 2d plane in cartesian space
		    float p[3];
		    vnormalize(drvui->frames[frame].mapnorm);
		    drvui->frames[frame].planeeq[0] = (float)drvui->frames[frame].mapnorm[0];
		    drvui->frames[frame].planeeq[1] = (float)drvui->frames[frame].mapnorm[1];
		    drvui->frames[frame].planeeq[2] = (float)drvui->frames[frame].mapnorm[2];
		    make_bmat (drvui->sys, drvui->lat_con, drvui->b_mat, drvui->ginv,drvui->rec_lat_con);
		    Convert_Cryst_Cart(drvui->b_mat, drvui->frames[frame].mapslice, p, origin);
		    drvui->frames[frame].planeeq[3] = -(float) (drvui->frames[frame].mapnorm[0] * p[0] 
						    + drvui->frames[frame].mapnorm[1] * p[1] 
						    + drvui->frames[frame].mapnorm[2] * p[2]);
		}
	    }
	    if (drvui->subsys_vol[0] == 1.0f)
		drvui->subsys_vol[0] = drvui->subsys_ref_volume;
	    return;		/* thru with end */
	case 10:		/* xyzoff command */
	    (void) sscanf (input, "%s %s %s %s", string, atom_pos[0], atom_pos[1],
			   atom_pos[2]);
	    for (i = 0; i < 3; i++)
		drvui->xyzoff[i] = convert_pos (atom_pos[i]);
	    drvui->xyzoff_read = 1;
	    break;
	case 51:		/* dash command */
	case 11:		/* bond command */
	    if (intype == 51) {
		k = sscanf (input, "%s %d", string, &j);
		if (k == 2)	/* number of segments if given */
		    Token_Strip (input, 1);
		else
		    j = 5;
	    }
	    get_label (input, &drvui->bonds[drvui->nbond].bond_l1[0],
		       &drvui->bonds[drvui->nbond].bond_l1[1]
		       , &drvui->bonds[drvui->nbond].bond_l1[2],
		       &drvui->bonds[drvui->nbond].bond_l1[3], 1);
	    get_label (input, &drvui->bonds[drvui->nbond].bond_l2[0],
		       &drvui->bonds[drvui->nbond].bond_l2[1]
		       , &drvui->bonds[drvui->nbond].bond_l2[2],
		       &drvui->bonds[drvui->nbond].bond_l2[3], 0);
	    (void) sscanf (input, "%f %f %f %39c", &drvui->bonds[drvui->nbond].bond_size,
			   &drvui->bonds[drvui->nbond].bond_min,
			   &drvui->bonds[drvui->nbond].bond_max, string);
	    trim_string (string, 40);
	    if (strlen (string) == 0)
		strcpy (string, "Grey");
	    strcpy (drvui->bonds[drvui->nbond].col_bond, string);
	    drvui->bonds[drvui->nbond].bond_fn = tmp_frame_no;
	    if (intype == 11)
		drvui->bonds[drvui->nbond].bond_style = 0;
	    else
		drvui->bonds[drvui->nbond].bond_style = j;
	    drvui->nbond++;
	    check_dynamic_storage ();
	    break;
	case 12:		/* box command */
	    (void) sscanf (input, "%s %f %39c", string, &rad_cell, t_color);
	    trim_string (t_color, 40);
	    strcpy (drvui->col_cell, t_color);
	    if (strlen (drvui->col_cell) == 0)
		strcpy (drvui->col_cell, "Black");
	    if (rad_cell == 0.0) {
		docell = 0;
	    }
	    break;
	case 13:		/* edges command */
	    (void) sscanf (input, "%s %f %39c", string, &drvui->rad_edge, t_color);
	    trim_string (t_color, 40);
	    strcpy (drvui->col_edge, t_color);
	    if (drvui->rad_edge > 0.005)
		edges = 2;
	    if (strlen (drvui->col_edge) == 0)
		strcpy (drvui->col_edge, "Black");
	    break;
	case 14:		/* view command */
//      if ((Options & V_OPT) == 0)         // only if no -v switch 
	    (void) sscanf (input, "%s %f %f %f", string, &xrot, &yrot, &zrot);
	    break;
	case 15:		/* plane command */
	    get_label (input, &drvui->planes[drvui->nplane].plane_l[0],
		       &drvui->planes[drvui->nplane].plane_l[1]
		       , &drvui->planes[drvui->nplane].plane_l[2],
		       &drvui->planes[drvui->nplane].plane_l[3], 1);
	    (void) sscanf (input, "%f %39c", &drvui->planes[drvui->nplane].plane_size,
			   string);
	    trim_string (string, 40);
	    if (!strlen (string))
		strcpy (string, "White");
	    strcpy (drvui->planes[drvui->nplane].plane_col, string);
	    drvui->planes[drvui->nplane].plane_fn = tmp_frame_no;
	    drvui->nplane++;
	    check_dynamic_storage ();
	    break;
	case 16:		/* bij or Bij command */
	    drvui->ellips[drvui->n_ellips].ell_type = 2;
	case 17:		/* betaij command */
	    if (intype == 17)
		drvui->ellips[drvui->n_ellips].ell_type = 0;
	case 18:		/* uij or Uij command */
	    if (intype == 18)
		drvui->ellips[drvui->n_ellips].ell_type = 1;
	    get_label (input, &drvui->ellips[drvui->n_ellips].ellips_l[0],
		       &drvui->ellips[drvui->n_ellips].ellips_l[1],
		       &drvui->ellips[drvui->n_ellips].ellips_l[2],
		       &drvui->ellips[drvui->n_ellips].ellips_l[3], 1);
	    drvui->ellips[drvui->n_ellips].ellips_l[4] = 0;
	    (void) sscanf (input, "%d %f %f %f %f %f %f %39c", &drvui->ellips[drvui->n_ellips].ellips_n, &drvui->ellips[drvui->n_ellips].ellips[0], &drvui->ellips[drvui->n_ellips].ellips[1], &drvui->ellips[drvui->n_ellips].ellips[2], &drvui->ellips[drvui->n_ellips].ellips[3], &drvui->ellips[drvui->n_ellips].ellips[4], &drvui->ellips[drvui->n_ellips].ellips[5], string);	// read unique coefficients and color
	    drvui->ellips[drvui->n_ellips].save_el_number =
		drvui->ellips[drvui->n_ellips].ellips_n;
	    trim_string (string, 40);
	    strcpy (drvui->ellips[drvui->n_ellips].ellips_col, string);
	    if (!strlen (drvui->ellips[drvui->n_ellips].ellips_col))
		strcpy (drvui->ellips[drvui->n_ellips].ellips_col, "Gray80");
	    drvui->n_ellips++;
	    check_dynamic_storage ();
	    break;
	case 19:		/* bounds command */
//      if ((Options & B_OPT) == 0) {         // only do this if no -b switch 
	    (void) sscanf (input, "%s %f %f %f", string, &boxlim[0], &boxlim[1],
			   &boxlim[2]);
	    boxflag = 1;
//      }
	    break;
	case 20:		/* pack command */
//      if ((Options & P_OPT) == 0) {         // only if no -p switch
	    (void) sscanf (input, "%s %f %f %f %f %f %f", string,
			   &drvui->frames[tmp_frame_no].cryst_lim[0],
			   &drvui->frames[tmp_frame_no].cryst_lim[3],
			   &drvui->frames[tmp_frame_no].cryst_lim[1],
			   &drvui->frames[tmp_frame_no].cryst_lim[4],
			   &drvui->frames[tmp_frame_no].cryst_lim[2],
			   &drvui->frames[tmp_frame_no].cryst_lim[5]);
	    packflag = 1;
//      }
	    break;
	case 21:		/* magnification command */
//      if ((Options & M_OPT) == 0)         // only if no -m switch
	    (void) sscanf (input, "%s %f", string, &Magnification);
	    break;
	case 22:		/* origin command */
//      if ((Options & O_OPT) == 0) {         // only if no -o switch
	    (void) sscanf (input, "%s %f %f %f", string, &origin[0], &origin[1],
			   &origin[2]);
//      }
	    break;
	case 23:		/* ellipcolor command */
	    Token_Strip (input, 1);
	    strcpy (el_color_tmp[n_el_tmp++].el_color_tmp, input);
	    break;
	case 24:		/* inline command */
	    in_line = 1;
	case 25:		/* import command */
	    char filetype[10];

	    done = 0;
	    if (!Quick)
		fprintf (drvui->flout, "*         %s\n======================= Start of import/inline commands\n", input2);	// echo to listing file
	    intype = 0;
	    Token_Strip (input, 1);	/* strip command */
	    sscanf (input, "%s", filetype);
	    Token_Strip (input, 1);	/* strip filetype */
	    trim_string (input, 256);	/* and remove trailing blanks */

	    if (strncmp (filetype, "gsas", 4) == 0) {
		import_gsas (input, in_line, Quick, tmp_frame_no);
		done = 1;
	    }
	    if (strncmp (filetype, "pcr", 3) == 0) {
		import_pcr (input, in_line, Quick, tmp_frame_no);
		done = 1;
	    }
	    if (strncmp (filetype, "cif", 3) == 0) {
		int ii, jj;

		ii = strlen (input) - 1;
		if (isdigit (input[ii])) {	// if number at end 
		    for (jj = ii; jj > ii - 6; jj--) {
			if (isspace (input[jj])) {	// separate number
			    input[jj] = '\0';
			    Block_CIF = atoi (&input[jj + 1]);
			    break;
			}
			if (isalpha (input[jj]))
			    break;	// unless it is attached to filename
		    }
		} else {
		    Block_CIF = 0;
		}
		import_cif (input, in_line, Quick, &Block_CIF, tmp_frame_no);
		done = 1;
	    }
	    if (strncmp (filetype, "schakal", 7) == 0) {
		import_schakal (input, in_line, Quick, tmp_frame_no);
		done = 1;
	    }
	    if (strncmp (filetype, "shakal", 6) == 0) {
		import_schakal (input, in_line, Quick, tmp_frame_no);
		done = 1;
	    }
	    if (strncmp (filetype, "shelx", 5) == 0) {
		import_shelx (input, in_line, Quick, tmp_frame_no);
		done = 1;
	    }
	    if (strncmp (filetype, "fdat", 4) == 0) {
		import_fdat (input, in_line, Quick, tmp_frame_no);
		done = 1;
	    }
	    if (strncmp (filetype, "csd", 3) == 0) {
		import_fdat (input, in_line, Quick, tmp_frame_no);
		done = 1;
	    }
	    if (strncmp (filetype, "ccdf", 4) == 0) {
		import_fdat (input, in_line, Quick, tmp_frame_no);
		done = 1;
	    }
	    if (strncmp (filetype, "wien2k", 6) == 0) {
		import_wien (input, in_line, Quick, tmp_frame_no);
		done = 1;
	    }
	    if (strncmp (filetype, "discus", 6) == 0) {
		import_discus (input, in_line, Quick, tmp_frame_no);
		done = 1;
	    }
	    if (strncmp (filetype, "exciting", 8) == 0
		|| strncmp (filetype, "elk", 3) == 0) {
		import_exc (input, in_line, Quick, tmp_frame_no);
		done = 1;
	    }
	    if (!done) {
		char string[256];

		sprintf (string, "Import/inline source type '%s' not recognized.", input);
		Error_Box (string);
		free (el_color_tmp);
		return;
	    }
	    if (!Quick)
		fprintf (drvui->flout,
			 "======================= End of import/inline commands\n");
	    break;
	case 26:		/* orthographic */
	    M_cameras = 0;
	    break;
	case 27:{		/* unit cell direction vectors */
		float v[3];

		Display_axes = 1;
		offset[0] = offset[1] = offset[2] = 0.0f;
		i = sscanf (input, "%s %f %f %f", string, &v[0], &v[1], &v[2]);
		if (i == 4)
		    for (i = 0; i < 3; i++)
			offset[i] = v[i];
	    }
	    break;
	case 28:		/* phong */
	    (void) sscanf (input, "%s %f %f", string, &drvui->Phong_Value,
			   &drvui->Phong_Size);
	    break;
	case 29:		/* molcomp */
	    (void) sscanf (input, "%s %f", string, &drvui->mol_d);
	    domolcomp = 1;
	    break;
	case 30:		/* cutout */
	    drvui->El_Cutout = 1;
	    (void) sscanf (input, "%s %39c", string, t_color);
	    trim_string (t_color, 40);
	    if (!strlen (t_color))
		strcpy (t_color, "White");
	    strcpy (drvui->Cutout_color, t_color);
	    break;
	case 31:		/* depthcue - thickness of edges scaled by Z coordinate */
	    (void) sscanf (input, "%s %f", string, &DepthCue);
	    break;
	case 32:		/* vrml1 - turn on VRML1 output */
	    Vrml2 = 0;
	    break;
	case 33:		/* polylimit (formerly polyfudge) */
	    (void) sscanf (input, "%s %f", string, &drvui->polylimit);
	    if (!Quick)
		fprintf (drvui->flout,
			 "Tolerance factor for non-planarity of polyhedra faces set to %f. This option should be used with care\n",
			 drvui->polylimit);
	    break;
	case 34:		/* clip command */
	    (void) sscanf (input, "%s %f %f %f %f %f %f", string,
			   &drvui->frames[tmp_frame_no].clip_lim[0],
			   &drvui->frames[tmp_frame_no].clip_lim[3],
			   &drvui->frames[tmp_frame_no].clip_lim[1],
			   &drvui->frames[tmp_frame_no].clip_lim[4],
			   &drvui->frames[tmp_frame_no].clip_lim[2],
			   &drvui->frames[tmp_frame_no].clip_lim[5]);
	    clipflag = 1;
	    break;
	case 35:		/* list - limit for printed distance table */
	    (void) sscanf (input, "%s %f", string, &printdist);
	    if (!Quick)
		fprintf (drvui->flout,
			 "Interatomic distances up to %6.3f input units will be tabulated in the logfile\n",
			 printdist);
	    break;
	case 36:		/* polyedge - enables individual coloring of polyhedral edges */
	    get_label (input, &t_name[0], &t_name[1], &t_name[2], &t_name[3], 1);
	    (void) sscanf (input, "%f %39c", &t_rad, t_color);
	    t_name[4] = 0;
	    trim_string (t_color, 40);
	    if (!strlen (t_color))
		strcpy (t_color, "White");
	    for (i = 1; i < drvui->npoly; i++) {	/* find this atom in poly list */
		if ((t_name[0] == drvui->polyhedra[i].poly_l[0]) &&
		    (t_name[1] == drvui->polyhedra[i].poly_l[1]) &&
		    (t_name[2] == drvui->polyhedra[i].poly_l[2])
		    && (t_name[3] == drvui->polyhedra[i].poly_l[3])) {
		    edges = 1;
		    drvui->polyhedra[i].poly_rad_edge = t_rad;
		    strcpy (drvui->polyhedra[i].poly_col_edge, t_color);
		}
	    }
	    trim_string (t_name, 4);
	    strcpy (drvui->polyedges[drvui->nedges].name, t_name);
	    strcpy (drvui->polyedges[drvui->nedges].color, t_color);
	    drvui->polyedges[drvui->nedges].radius = t_rad;
	    drvui->nedges++;
	    check_dynamic_storage ();
	    break;
	case 37:		/* frame command */
	    /* if the current frame inherited a clip flag, default clip limit to pack range */
	    if (clipflag == 1
		&& drvui->frames[tmp_frame_no].clip_lim[0] ==
		drvui->frames[tmp_frame_no].clip_lim[3]) {
		for (i = 0; i < 6; i++)
		    drvui->frames[tmp_frame_no].clip_lim[i] =
			drvui->frames[tmp_frame_no].cryst_lim[i];
	    }
	    tmp_frame_no++;
	    i = drvui->max_frame;
	    drvui->max_frame = tmp_frame_no;
	    check_dynamic_storage ();
	    drvui->max_frame = i;
	    break;
	case 38:		/* axislines */
	    (void) sscanf (input, "%s %f %39c", string, &drvui->Ellipaxis_width, t_color);
	    trim_string (t_color, 40);
	    strcpy (drvui->Ellipaxis_color, t_color);
	    if (strlen (drvui->Ellipaxis_color) == 0)
		strcpy (drvui->Ellipaxis_color, "Gray20");
	    break;
	case 40:		/* nolabels command */
	    Labels = 0;		/* set value false */
	    break;
	case 41:		/* lonepair command */
	    get_label (input, &drvui->cones[drvui->ncone].cone_l1[0],
		       &drvui->cones[drvui->ncone].cone_l1[1]
		       , &drvui->cones[drvui->ncone].cone_l1[2],
		       &drvui->cones[drvui->ncone].cone_l1[3], 1);
	    (void) sscanf (input, "%d %f %f %f %39c",
			   &drvui->cones[drvui->ncone].numlonepairs,
			   &drvui->cones[drvui->ncone].cone_height,
			   &drvui->cones[drvui->ncone].cone_min,
			   &drvui->cones[drvui->ncone].cone_max, t_color);
	    trim_string (t_color, 40);
	    strcpy (drvui->cones[drvui->ncone].col_cone, t_color);
	    if (strlen (t_color) == 0)
		strcpy (drvui->cones[drvui->ncone].col_cone, "Grey");
	    drvui->cones[drvui->ncone].cone_fn = tmp_frame_no;
	    drvui->ncone++;
	    check_dynamic_storage ();
	    break;
	case 42:		/* labeltext command */
	    memset (drvui->labels[drvui->nlabel].label_label, 0, 64);
	    Token_Strip (input, 1);
	    (void) sscanf (input, "%f %f %f %63c",
			   &(drvui->labels[drvui->nlabel].label_x[0]),
			   &(drvui->labels[drvui->nlabel].label_x[1]),
			   &(drvui->labels[drvui->nlabel].label_x[2]),
			   drvui->labels[drvui->nlabel].label_label);
	    trim_string (drvui->labels[drvui->nlabel].label_label, 64);
	    drvui->labels[drvui->nlabel].label_fn = tmp_frame_no;
	    intype = 0;
	    if (!strcmp (drvui->labels[drvui->nlabel].label_label, "a") && curframe != 1)
		break;
	    if (!strcmp (drvui->labels[drvui->nlabel].label_label, "b") && curframe != 1)
		break;
	    if (!strcmp (drvui->labels[drvui->nlabel].label_label, "c") && curframe != 1)
		break;
	    if (!strcmp (drvui->labels[drvui->nlabel].label_label, "o") && curframe != 1)
		break;
	    intype = 42;
	    drvui->nlabel++;
	    check_dynamic_storage ();
	    break;
	case 44:		// background (color) command
	    Token_Strip (input, 1);
	    strncpy (drvui->col_bg, input, 30);
	    trim_string (drvui->col_bg, 40);
	    if (strlen (drvui->col_bg) < 2)
		strcpy (drvui->col_bg, "White");
	    strcpy (input, drvui->col_bg);
	    Transform_VRML_Color (input);
	    sscanf (input, "%f %f %f", &drvui->glback[0], &drvui->glback[1],
		    &drvui->glback[2]);
	    break;
	case 45:
	    (void) sscanf (input, "%*s %f %f %f %f %f %f %f %f %f %f %f %f %d",
			   &drvui->slab_con[0], &drvui->slab_con[1], &drvui->slab_con[2],
			   &drvui->slab_con[3], &drvui->slab_con[4], &drvui->slab_con[5],
			   &drvui->slab_off[0], &drvui->slab_off[1], &drvui->slab_off[2],
			   &drvui->slab_rot[0], &drvui->slab_rot[1], &drvui->slab_rot[2],
			   &slabmode);

	    break;
	case 46:		// mag_trans command
	    Token_Strip (input, 1);
	    (void) sscanf (input, "%f %f %f %f %f %f %f %f %f", &drvui->mag_matrix[0][0],
			   &drvui->mag_matrix[0][1], &drvui->mag_matrix[0][2],
			   &drvui->mag_matrix[1][0], &drvui->mag_matrix[1][1],
			   &drvui->mag_matrix[1][2], &drvui->mag_matrix[2][0],
			   &drvui->mag_matrix[2][1], &drvui->mag_matrix[2][2]);
	    break;
	case 47:		// special command
	    if (!Quick)
		break;
	    (void) sscanf (input, "%*7c %d %d", &Omit->omit1[Omit->nomits],
			   &Omit->omit2[Omit->nomits]);
	    Omit->nomits++;
	    break;
	case 48:		// noshadow command
	    drvui->noshadow = 1;
	    break;
	case 49:		// finish command
	    (void) sscanf (input, "%*6c %f %f %f %f", &drvui->ambient,
			   &drvui->diffuse, &drvui->specular, &drvui->roughness);
	    break;
	case 50:		// arrow command
	    Token_Strip (input, 1);
	    (void) sscanf (input, "%f %f %f %f %f %f %f %f %39c",
			   &drvui->arrows[drvui->nmag].mag_xp[0],
			   &drvui->arrows[drvui->nmag].mag_xp[1],
			   &drvui->arrows[drvui->nmag].mag_xp[2],
			   &drvui->arrows[drvui->nmag].mag_xc[0],
			   &drvui->arrows[drvui->nmag].mag_xc[1],
			   &drvui->arrows[drvui->nmag].mag_xc[2],
			   &drvui->arrows[drvui->nmag].arrow_length,
			   &drvui->arrows[drvui->nmag].arrow_diam, t_color);
	    trim_string (t_color, 40);
	    if (!strlen (t_color))
		strcpy (t_color, "White");
	    strcpy (drvui->arrows[drvui->nmag].col_arrow, t_color);
	    drvui->arrows[drvui->nmag].arrow_fn = tmp_frame_no;
	    drvui->nmag++;
	    check_dynamic_storage ();
	    break;
	case 52:		/* mapcontour command */
	    {
	    drvui->numOfFourierContours++;
	    check_dynamic_storage ();
	    Token_Strip (input, 1);
	    (void) sscanf (input, "%f %s %69c",
			   &drvui->fourier[drvui->numOfFourierContours].
			   FourierContourLevel, string, t_color);
	    trim_string (t_color, 70);
	    if (!strlen (t_color))
		strcpy (t_color, "White");
	    char *amp  = strchr(t_color,'&');
	    if (amp) { 
		*amp='\0'; 
		amp++; 
		strcpy(drvui->fourier[drvui->numOfFourierContours].FourierBackColor,amp);
	    }
	    strcpy (drvui->fourier[drvui->numOfFourierContours].FourierContourColor,
		    t_color);
	    if (strncmp (string, "mesh", 4) == 0)
		drvui->fourier[drvui->numOfFourierContours].FourierContourSolid = 0;
	    else
		drvui->fourier[drvui->numOfFourierContours].FourierContourSolid = 1;
	    }
	    break;
	case 53:		/* mapread command */
	    char string2[20], res[20];

	    if (!Quick)
		break;
	    Token_Strip (input, 1);
	    i = sscanf(input, "%s %s %s %s", string, filename, string2, res);
	    strcpy (FourierFileName, filename);
	    if (i < 4)
		Map_Info.res = 4; /* No resolution specified - set 0.25 A */
	    else
		sscanf(res, "%d", &Map_Info.res);
	    if (Map_Info.res <= 0) 
		Map_Info.res = 1;
	    if (i < 3)
		strcpy(string2, "Fo"); /* Map type not specified - default to Fo */
	    if (strncmp (string, "grd", 3) == 0) {
		FourierMapType = 1;
	    } else if (strncmp (string, "stf", 3) == 0) {
		FourierMapType = 2;
	    } else if (strncmp (string, "w2k", 3) == 0) {
		FourierMapType = 3;
	    } else if (strncmp (string, "vsp", 3) == 0) {
		FourierMapType = 4;
	    } else if (strncmp (string, "flp", 3) == 0) {
		FourierMapType = 5;
	    } else if (strncmp (string, "fcf", 3) == 0) {
		FourierMapType = 6;
		Map_Info.map_type = 0;
		if (!strcmp (string2, "Fc"))
		    Map_Info.map_type = 1;
		else if (!strcmp (string2, "Fo-Fc"))
		    Map_Info.map_type = 2;
		else if (!strcmp (string2, "2Fo-Fc"))
		    Map_Info.map_type = 3;
		else if (!strcmp (string2, "Fo2"))
		    Map_Info.map_type = 4;
	    } else if (strncmp (string, "dn6", 3) == 0) {
		FourierMapType = 7;
	    } else if (strncmp (string, "m80", 3) == 0) {
		FourierMapType = 8;
		Map_Info.map_type = 0;
		if (!strcmp (string2, "Fc"))
		    Map_Info.map_type = 1;
		else if (!strcmp (string2, "Fo-Fc"))
		    Map_Info.map_type = 2;
		else if (!strcmp (string2, "2Fo-Fc"))
		    Map_Info.map_type = 3;
		else if (!strcmp (string2, "Fo2"))
		    Map_Info.map_type = 4;
	    } else if (strncmp (string, "m81", 3) == 0) {
		FourierMapType = 10;
		Map_Info.map_type = 0;
		if (!strcmp (string2, "Fc"))
		    Map_Info.map_type = 1;
		else if (!strcmp (string2, "Fo-Fc"))
		    Map_Info.map_type = 2;
		else if (!strcmp (string2, "2Fo-Fc"))
		    Map_Info.map_type = 3;
		else if (!strcmp (string2, "Fo2"))
		    Map_Info.map_type = 4;
	    } else if (strncmp (string, "exc", 3) == 0) {
		FourierMapType = 9;
	    } else if (strncmp (string, "xsf", 3) == 0) {
		FourierMapType = 11;
	    } else {
		Error_Box ("mapread command does not specify a legal file type.");
	    }
	    break;
	case 54:		/* mapregion command */
	    Token_Strip (input, 1);
	    (void) sscanf (input, "%f %f %f %f %f %f %f %f %f", &xMin, &xMax, &yMin,
			   &yMax, &zMin, &zMax, &x4Val, &x5Val, &x6Val);
	    drvui->frames[tmp_frame_no].map_lim_set= 1;
	    drvui->frames[tmp_frame_no].map_lim[0]=xMin;
	    drvui->frames[tmp_frame_no].map_lim[3]=xMax;
	    drvui->frames[tmp_frame_no].map_lim[1]=yMin;
	    drvui->frames[tmp_frame_no].map_lim[4]=yMax;
	    drvui->frames[tmp_frame_no].map_lim[2]=zMin;
	    drvui->frames[tmp_frame_no].map_lim[5]=zMax;
	    drvui->frames[tmp_frame_no].map_lim[6]=x4Val;
	    drvui->frames[tmp_frame_no].map_lim[7]=x5Val;
	    drvui->frames[tmp_frame_no].map_lim[8]=x6Val;
	    i = 0;
	    if (fabs (xMax - xMin) < 0.001f)
		i++;
	    if (fabs (yMax - yMin) < 0.001f)
		i++;
	    if (fabs (zMax - zMin) < 0.001f)
		i++;
	    if (i == 1)
		drvui->Fourier2d = 1;
	    break;
	case 55:		/* mapcalclimits command */
	    Token_Strip (input, 1);
	    (void) sscanf (input, "%f %f %f %f %f %f", &Map_Info.xlim[0],
			   &Map_Info.xlim[1], &Map_Info.ylim[0], &Map_Info.ylim[1],
			   &Map_Info.zlim[0], &Map_Info.zlim[1]);
	    break;
	case 56:		/* labelscale command */
	    Token_Strip (input, 1);
	    (void) sscanf (input, "%f", &drvui->label_scale);
	    break;
	case 57:		/* bestplane command */
	    char bplane[15][5];

	    sscanf (input, "%*s %d", &drvui->bplanes[drvui->nbplane].nbatoms);
	    if (drvui->bplanes[drvui->nbplane].nbatoms > 15) {
		Error_Box ("Cannot handle more than 15 atoms in a best plane.");
		if (el_color_tmp)
		    free (el_color_tmp);
		return;
	    }
	    Token_Strip (input, 2);
	    for (i = 0; i < drvui->bplanes[drvui->nbplane].nbatoms; i++) {
		sscanf (input, "%s", bplane[i]);
		Token_Strip (input, 1);
	    }
	    memset (drvui->bplanes[drvui->nbplane].bplane_col, 0, 40);
	    sscanf (input, "%f %f %39c", &drvui->bplanes[drvui->nbplane].bplane_d1,
		    &drvui->bplanes[drvui->nbplane].bplane_d2, t_color);
	    trim_string (t_color, 40);
	    if (!strlen (t_color))
		strcpy (t_color, "White");
	    strcpy (drvui->bplanes[drvui->nbplane].bplane_col, t_color);
	    for (i = 0; i < drvui->bplanes[drvui->nbplane].nbatoms; i++) {
		memset (drvui->bplanes[drvui->nbplane].bplane_t[i], 0, 4);
		strcpy (drvui->bplanes[drvui->nbplane].bplane_t[i], bplane[i]);
		char *p = strstr (bplane[i], "*");

		if (p != NULL) {
		    *p = '\0';
		    drvui->bplanes[drvui->nbplane].bplane_n[i] = -1;
		} else {
		    for (j = 1; j < 4; j++)
			if (bplane[i][j] <= '9') {
			    char *q = drvui->bplanes[drvui->nbplane].bplane_t[i] + j;

			    drvui->bplanes[drvui->nbplane].bplane_n[i] =
				strtol (q, NULL, 10);
			    for (k = j; k < 4; k++)
				drvui->bplanes[drvui->nbplane].bplane_t[i][k] = ' ';
			    break;
			}
		}
	    }
	    drvui->nbplane++;
	    check_dynamic_storage ();
	    break;
	case 58:
	    average = 1;
	    break;
	case 59:
	    Token_Strip (input, 1);
	    memset (modl, 0, 4);
	    sscanf (input, "%s %d %f %f", modl, &modnum, &avg_occ, &min_occ);
	    if (modl[1] == '\0')
		modl[1] = modl[2] = modl[3] = ' ';
	    if (modl[2] == '\0')
		modl[2] = modl[3] = ' ';
	    if (modl[3] == 0)
		modl[3] = ' ';
	    modl[4] = '\0';
	    for (j = 0; j < natom; j++) {
		if (check_atom_name (modl, drvui->atoms[j].atom_l)
		    && modnum == drvui->atoms[j].atom_n) {
		    drvui->atoms[j].occupancy = avg_occ;
		    drvui->atoms[j].min_occ = min_occ;
		    break;
		}
	    }
	    break;
	case 60:
	    sscanf (input, "%*s %f %f %f", &drvui->phaseshift[0], &drvui->phaseshift[1],
		    &drvui->phaseshift[2]);
	    for (j = 0; j < 3; j++)
		drvui->phaseshift[j] = (float) fmod (drvui->phaseshift[j], 1.);
	    break;
	case 61:		/* mapcontour2d command */
	    drvui->numOfFourierContours++;
	    check_dynamic_storage ();
	    Token_Strip (input, 1);
	    (void) sscanf (input, "%f %f %f %39c",
			   &drvui->fourier[drvui->numOfFourierContours].
			   FourierContourLevel,
			   &drvui->fourier[drvui->numOfFourierContours].
			   FourierContourStep,
			   &drvui->fourier[drvui->numOfFourierContours].FourierContourTop,
			   t_color);
	    trim_string (t_color, 40);
	    if (!strlen (t_color))
		strcpy (t_color, "White");
	    strcpy (drvui->fourier[drvui->numOfFourierContours].FourierContourColor,
		    t_color);
	    break;
	case 62:
	    memset (modl, 0, 4);
	    sscanf (input, "%*s %s %d %s %s %39c", modl, &modnum, filename, string,
		    t_color);
	    if (modl[1] == '\0')
		modl[1] = modl[2] = modl[3] = ' ';
	    if (modl[2] == '\0')
		modl[2] = modl[3] = ' ';
	    if (modl[3] == 0)
		modl[3] = ' ';
	    modl[4] = '\0';
	    if (drvui->nsurf == MAX_SURF) {
		if (!Quick)
		    fprintf (drvui->flout,
			     "Maximum number of surfaces exceeded, skipping\n");
		break;
	    }
	    strcpy (drvui->surfatom[drvui->nsurf], modl);
	    drvui->surfnum[drvui->nsurf] = modnum;
	    if (!strncmp (string, "mesh", 4)) {
		drvui->surftype[drvui->nsurf] = 0;
	    } else if (!strncmp (string, "solid", 5)) {
		drvui->surftype[drvui->nsurf] = 1;
	    } else
		drvui->surftype[drvui->nsurf] = 2;
	    trim_string (t_color, 40);
	    if (!strlen (t_color))
		strcpy (t_color, "White");
	    strcpy (drvui->surfcolor[drvui->nsurf], t_color);
	    strcpy (drvui->surffile[drvui->nsurf], filename);
	    read_aim (filename, Quick);
	    drvui->nsurf++;
	    break;
	case 63:		/* voids command */
	    sscanf (input, "%*s %d %f %d %d %d %39c", &j, &drvui->probesize,
		    &drvui->voidgrid[0], &drvui->voidgrid[1], &drvui->voidgrid[2],
		    t_color);
	    trim_string (t_color, 40);
	    if (!strlen (t_color))
		strcpy (t_color, "White");
	    strcpy (drvui->voidcolor, t_color);
	    if (drvui->voidflag >= 0) {
		drvui->voidflag = j;
		if (j == 1) {
		    drvui->voidmap =
			(char ***) zalloc (drvui->voidgrid[0] * sizeof (*drvui->voidmap));
		    for (j = 0; j < drvui->voidgrid[0]; j++) {
			drvui->voidmap[j] =
			    (char **) zalloc (drvui->voidgrid[1] *
					      sizeof (*drvui->voidmap[j]));
			for (k = 0; k < drvui->voidgrid[1]; k++)
			    drvui->voidmap[j][k] =
				(char *) zalloc (drvui->voidgrid[2] *
						 sizeof (*drvui->voidmap[j][k]));
		    }
		}
	    }
	    break;
	case 64:		/* values command */
	    get_label (input, &drvui->atprops[drvui->natprop].atprop_l[0],
		       &drvui->atprops[drvui->natprop].atprop_l[1],
		       &drvui->atprops[drvui->natprop].atprop_l[2],
		       &drvui->atprops[drvui->natprop].atprop_l[3], 1);

	    i = sscanf (input, "%s %s", str[0], str[1]);
	    if (atoi (str[0]) == 0)
		drvui->atprops[drvui->natprop].atprop_n = -1;
	    else
		drvui->atprops[drvui->natprop].atprop_n = atoi(str[0]);
	    drvui->atprops[drvui->natprop].radius = (float)atof (str[1]);
	    drvui->atprops[drvui->natprop].atprop_fn = tmp_frame_no;
	    drvui->natprop++;
	    check_dynamic_storage ();
	    break;
	case 65: /* mapslice command */
	    Token_Strip (input,1);
	    j = sscanf (input,"%f %f %f %f %f %f %d",&drvui->frames[tmp_frame_no].mapslice[0],
		    &drvui->frames[tmp_frame_no].mapslice[1],&drvui->frames[tmp_frame_no].mapslice[2],
		    &drvui->frames[tmp_frame_no].mapnorm[0],&drvui->frames[tmp_frame_no].mapnorm[1],
		    &drvui->frames[tmp_frame_no].mapnorm[2],&drvui->frames[tmp_frame_no].slice);
	    if (j <7) drvui->frames[tmp_frame_no].slice = 1;
	    break;
	case 66: /* qvector command */
	    Token_Strip (input,1);
	    j = drvui->no_cell_vec;
	    sscanf (input,"%f %f %f",  &drvui->cell_vec[j][0],
		&drvui->cell_vec[j][1], &drvui->cell_vec[j][2]);
	    drvui->no_cell_vec++;
	    drvui->modulated++;
	    break;
	case 67:
	    ShowMapLegend = 1;
	    break;
	case 68:
	    Vrml2 = 1;
	    X3D = 1;
	    break;
	default:
	    if (!Quick)
		fprintf (drvui->flout,
			 "********** Next line does not start with a recognizable keyword - ignored\n");
	    break;
	}
	if (!Quick && intype)
	    fprintf (drvui->flout, "*         %s\n", input2);	// echo to listing file
    }
    if (drvui->n_ellips > 1)
	drvui->do_ellipsoids = 1;
    if (packflag != 0) {
	boxlim[0] = boxlim[1] = boxlim[2] = 500.0f;
    }
}				/* end of read_inp */

/* ************************************************************** */
/* ************************************************************** */

void
set_tf_status (void)
{
// set the TF_status flags for the atoms in the list

    int i, j;

    for (i = 0; i < natom; i++) {
	if (drvui->atoms[i].atom_fn != drvui->frame_no)
	    continue;
	drvui->atoms[i].TF_status = -1;
	for (j = 1; j < drvui->n_ellips; j++) {
	    if (check_atom_name (drvui->atoms[i].atom_l, drvui->ellips[j].ellips_l)) {
		if (drvui->ellips[j].ellips_n == -1
		    || drvui->ellips[j].ellips_n == drvui->atoms[i].sv_atom_n) {
		    if (drvui->ellips[j].ell_type > 1000) {
			drvui->atoms[i].TF_status = 1;
			j = drvui->n_ellips;
			break;
		    } else {
			if (drvui->ellips[j].ell_type > 0) {
			    drvui->atoms[i].TF_status = 0;
			    j = drvui->n_ellips;
			    break;
			}
		    }
		}
	    }
	}
	j = drvui->atoms[i].TF_status;
    }
}


/* ************************************************************** */
/* ************************************************************** */

void
skip_blocks (int i, FILE * in)
{
// routine to skip data blocks in a CIF
    char string[256];

    int j;

    memset (string, 0, 255);
    for (j = 0; j <= i; j++) {
	while (strncmp (string, "data_", 5) != 0) {
	    if (!fgets (string, 255, in)) {	/* search for data_ keyword */
		Error_Box ("Error skipping data blocks in CIF Import File, Run aborted.");
		(void) fclose (in);
		return;
	    }
	}
	memset (string, 0, 255);
    }
}


/* ************************************************************** */
/* ************************************************************** */

void
Token_Strip (char string[], int no)

/*  string[] - input output string */
/*  no - number of tokens to strip */
/* routine to strip tokens from the leading part of a string */
{
    int i, j;
    int l = strlen(string);

    for (i = 0, j = 0; i < no; i++) {	/* loop through tokens to be removed */
	while (j < l-1 && string[j++] != ' ')
	    string[j - 1] = ' ';	/* characters to space */
	while (j < l && string[j] == ' ')
	    j++;		/* skip all spaces */
    }
    Blank_Strip (string);	/* move stuff to left */
}

/* ************************************************************** */
/* ************************************************************** */

void
Transform_POV_Color (char *color)
{
// Transform a POV color is specified in the form f1 f2 f2 {filter f} format,
// where {filter} is optional. If a color name is specified, return with no
// action.
    float red, green, blue, filter;

    char working[120], string[20];

    char *p;

    if (!strtod (color, &p) && p == color)
	return;			/* we have the color name form */

    strcpy (working, color);	/* make a working copy */
    (void) sscanf (working, "%f %f %f", &red, &green, &blue);

    if ((p = strstr (working, "filter"))) {	/* if "filter" given, return rgbf<...> form */
	(void) sscanf (p, "%s %f", string, &filter);
	sprintf (working, "rgbf <%.4f,%.4f,%.4f,%.4f>", red, green, blue, filter);
    } else {			/* return rgb <...> form */
	sprintf (working, "rgb <%.4f,%.4f,%.4f>", red, green, blue);
    }
    strcpy (color, working);	/* reformatted color to caller */
}

/* ************************************************************** */
/* ************************************************************** */

void
Transform_VRML_Color (char *input)
{
    char *transp_start;

    char transp[20];

    float fval;

    int i;

    char *a, b[30], line[81];

    char thecolor[50];

    char keyword[80];

    float thered, thegreen, theblue;

    char POV_incpath[255] = "\0";

    FILE *colinc;

/*
   VRML routine to transform VRML (and openGL) colors into RGB triples
   derived from 'colors.inc' of POV

   If the color information is "float1 float2 float3", then extract the rgb
   values from the input.   LWF 9/30/05
*/

    i = strlen (input);
    if (i == 0)
	return;
    if (input[i - 1] < ' ')
	input[i - 1] = 0;	// get rid of any ^J character
    transp[0] = '\0';
    transp_start = strstr (input, "filter");
    if (transp_start) {
	fval = (float) atof (transp_start + 6);
	sprintf (transp, " transparency %5.2f", fval / 2.);
    }
    if (strtod (input, &a) || a != input) {	// rgb form of color
	sscanf (input, "%f %f %f", &thered, &thegreen, &theblue);
	sprintf (input, "%.6f %.6f %.6f", thered, thegreen, theblue);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Black", 5)) {
	strncpy (input, "0 0 0\0", 6);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Brass", 5)) {
	strncpy (input, ".71 .65 .26\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Bronze2", 7)) {
	strncpy (input, ".65 .49 .24\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Bronze", 6)) {
	strncpy (input, ".55 .47 .14\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "BlueViolet", 10)) {
	strncpy (input, ".62352 .372549 .623529\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Blue", 3)) {
	strncpy (input, "0 0 1\0", 6);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Brown", 5)) {
	strncpy (input, ".647059 .164706 .164706\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Copper", 6)) {
	strncpy (input, ".72 .45 .20\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Cyan", 3)) {
	strncpy (input, "0 1 1\0", 6);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "DimGray", 7)) {
	strncpy (input, ".329412 .329412 .329412\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "DimGrey", 7)) {
	strncpy (input, ".329412 .329412 .329412\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Aquamarine", 10)) {
	strncpy (input, ".439216 .858824 .576471\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "CadetBlue", 9)) {
	strncpy (input, ".372549 .623529 .623529\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Coral", 5)) {
	strncpy (input, "1.0 .498039 .0\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "CornflowerBlue", 14)) {
	strncpy (input, ".258824 .258824 .435294\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "DarkGreen", 9)) {
	strncpy (input, ".184314 .309804 .184314\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "DarkOliveGreen", 14)) {
	strncpy (input, ".309804 .309804 .184314\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "DarkOrchid", 10)) {
	strncpy (input, ".6 .196078 .8\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "DarkSlateBlue", 13)) {
	strncpy (input, ".119608 .137255 .556863\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "DarkSlateGray", 13)) {
	strncpy (input, ".184314 .309804 .309804\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "DarkSlateGrey", 13)) {
	strncpy (input, ".184314 .309804 .309804\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "DarkTurquoise", 13)) {
	strncpy (input, ".439216 .576471 .858824\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Firebrick", 9)) {
	strncpy (input, ".556863 .137255 .137255\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "ForestGreen", 11)) {
	strncpy (input, ".137255 .556863 .137255\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Goldenrod", 9)) {
	strncpy (input, ".858824 .858824 .439216\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gold", 4)) {
	strncpy (input, ".8 .498039 .196078\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray05", 6)) {
	strncpy (input, ".05 .05 .05\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray10", 6)) {
	strncpy (input, ".10 .10 .10\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray15", 6)) {
	strncpy (input, ".15 .15 .15\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray20", 6)) {
	strncpy (input, ".20 .20 .20\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray25", 6)) {
	strncpy (input, ".25 .25 .25\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray30", 6)) {
	strncpy (input, ".30 .30 .30\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray35", 6)) {
	strncpy (input, ".35 .35 .35\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray40", 6)) {
	strncpy (input, ".40 .40 .40\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray45", 6)) {
	strncpy (input, ".45 .45 .45\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray50", 6)) {
	strncpy (input, ".50 .50 .50\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray55", 6)) {
	strncpy (input, ".55 .55 .55\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray60", 6)) {
	strncpy (input, ".60 .60 .60\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray65", 6)) {
	strncpy (input, ".65 .65 .65\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray70", 6)) {
	strncpy (input, ".70 .70 .70\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray75", 6)) {
	strncpy (input, ".75 .75 .75\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray80", 6)) {
	strncpy (input, ".80 .80 .80\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray85", 6)) {
	strncpy (input, ".85 .85 .85\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray90", 6)) {
	strncpy (input, ".90 .90 .90\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray95", 6)) {
	strncpy (input, ".95 .95 .95\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Gray", 4)) {
	strncpy (input, ".752941 .752941 .752941\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Grey", 4)) {
	strncpy (input, ".752941 .752941 .752941\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "GreenYellow", 11)) {
	strncpy (input, ".576471 .858824 .439216\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Green", 5)) {
	strncpy (input, "0 1 0\0", 6);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "IndianRed", 9)) {
	strncpy (input, ".309804 .184314 .184314\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Khaki", 5)) {
	strncpy (input, ".623529 .623529 .372549\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "LightBlue", 9)) {
	strncpy (input, ".74902 .847059 .847059\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "LightGray", 9)) {
	strncpy (input, ".658824 .658824 .658824\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "LightGrey", 9)) {
	strncpy (input, ".658824 .658824 .658824\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Light_Purple", 12)) {
	strncpy (input, ".87 .58 .98\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "LightSteelBlue", 14)) {
	strncpy (input, ".560784 .560784 .737255\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "LimeGreen", 9)) {
	strncpy (input, ".196078 .8 .196078\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Magenta", 7)) {
	strncpy (input, "1 0 1\0", 6);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Maroon", 6)) {
	strncpy (input, ".556863 .137255 .419608\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "MediumAquamarine", 16)) {
	strncpy (input, ".196078 .8 .6\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "MediumBlue", 10)) {
	strncpy (input, ".196078 .196078 .8\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "MediumForestGreen", 17)) {
	strncpy (input, ".419608 .556863 .137255\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "MediumGoldenrod", 15)) {
	strncpy (input, ".917647 .917647 .678431\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "MediumOrchid", 12)) {
	strncpy (input, ".576471 .439216 .858824\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "MediumSeaGreen", 14)) {
	strncpy (input, ".258824 .435294 .258824\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "MediumSlateBlue", 15)) {
	strncpy (input, ".498039 0.  1.0\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "MediumSpringGreen", 17)) {
	strncpy (input, ".498039 1.0 0\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "MediumTurquoise", 15)) {
	strncpy (input, ".439216 .858824 .858824\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "MediumVioletRed", 15)) {
	strncpy (input, ".858824 .439216 .576471\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "MidnightBlue", 12)) {
	strncpy (input, ".184314 .184314 .309804\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "NavyBlue", 8)) {
	strncpy (input, ".137255 .137255 .556863\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Navy", 4)) {
	strncpy (input, ".137255 .137255 .556863\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "OrangeRed", 9)) {
	strncpy (input, "1.0 .498039 0.0\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Orange", 6)) {
	strncpy (input, "1 .5 .0\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Orchid", 6)) {
	strncpy (input, ".858824 .439216 .858824\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "PaleGreen", 9)) {
	strncpy (input, ".560784 .737255 .560784\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Pink", 4)) {
	strncpy (input, ".737255 .560784 .560784\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Plum", 4)) {
	strncpy (input, ".917647 .678431 .917647\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Red", 3)) {
	strncpy (input, "1 0 0\0", 6);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "RichBlue", 8)) {
	strncpy (input, ".35 .35 .67\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Salmon", 6)) {
	strncpy (input, ".435294 .258824 .258824\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "SeaGreen", 8)) {
	strncpy (input, ".137255 .556863 .419608\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Sienna", 6)) {
	strncpy (input, ".556863 .419608 .137255\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "SkyBlue", 7)) {
	strncpy (input, ".196078 .6 .8\0", 14);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "SlateBlue", 9)) {
	strncpy (input, "0 .498039 1\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "SpringGreen", 11)) {
	strncpy (input, "0 1 .498039\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "SteelBlue", 9)) {
	strncpy (input, ".137255 .419608 .556863\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "SummerSky", 9)) {
	strncpy (input, ".22 .69 .87\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Tan", 3)) {
	strncpy (input, ".858824 .576471 .439216\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Thistle", 7)) {
	strncpy (input, ".847059 .74902 .847059\0", 23);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Turquoise", 9)) {
	strncpy (input, ".678431 .917647 .917647\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "VioletRed", 9)) {
	strncpy (input, ".8 .196078 .6\0", 14);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Violet", 6)) {
	strncpy (input, ".309804 .184314 .309804\0", 24);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "VLightGray", 10)) {
	strncpy (input, ".8 .8 .8\0", 9);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "VLightGrey", 10)) {
	strncpy (input, ".8 .8 .8\0", 9);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Wheat", 5)) {
	strncpy (input, ".847059 .847059 .74902\0", 23);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Silver", 6)) {
	strncpy (input, ".90 .91 .98\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "BrightGold", 10)) {
	strncpy (input, ".85 .85 .10\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "OldGold", 7)) {
	strncpy (input, ".81 .71 .23\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Feldspar", 8)) {
	strncpy (input, ".82 .57 .46\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Quartz", 6)) {
	strncpy (input, ".85 .85 .95\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Mica", 4)) {
	strncpy (input, "0 0 0\0", 6);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "NeonPink", 8)) {
	strncpy (input, "1 .43 .78\0", 10);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "DarkPurple", 10)) {
	strncpy (input, ".53 .12 .47\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "NeonBlue", 8)) {
	strncpy (input, ".3 .3 1\0", 8);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "CoolCopper", 10)) {
	strncpy (input, ".85 .53 .10\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "MandarinOrange", 14)) {
	strncpy (input, ".89 .47 .20\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "LightWood", 9)) {
	strncpy (input, ".91 .76 .65\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "MediumWood", 10)) {
	strncpy (input, ".65 .50 .39\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "DarkWood", 8)) {
	strncpy (input, ".52 .37 .26\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "SpicyPink", 9)) {
	strncpy (input, "1 .11 .68\0", 10);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "SemiSweetChoc", 13)) {
	strncpy (input, ".42 .26 .15\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "BakersChoc", 10)) {
	strncpy (input, ".36 .20 .09\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Flesh", 5)) {
	strncpy (input, ".96 .80 .69\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "NewTan", 6)) {
	strncpy (input, ".92 .78 .62\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "NewMidnightBlue", 15)) {
	strncpy (input, "0 0 .61\0", 8);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "VeryDarkBrown", 13)) {
	strncpy (input, ".35 .16 .14\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "DarkBrown", 9)) {
	strncpy (input, ".36 .25 .20\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "DarkTan", 7)) {
	strncpy (input, ".59 .41 .31\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "GreenCopper", 11)) {
	strncpy (input, ".32 .49 .46\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "DkGreenCopper", 13)) {
	strncpy (input, ".29 .46 .43\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "DustyRose", 9)) {
	strncpy (input, ".52 .39 .39\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "HuntersGreen", 12)) {
	strncpy (input, ".13 .37 .31\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Scarlet", 7)) {
	strncpy (input, ".55 .09 .09\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Med_Purple", 10)) {
	strncpy (input, ".73 .16 .96\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "White", 5)) {
	strncpy (input, "1 1 1\0", 6);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Very_Light_Purple", 17)) {
	strncpy (input, ".94 .81 .99\0", 12);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "YellowGreen", 11)) {
	strncpy (input, ".6 .8 .196078\0", 14);
	strcat (input, transp);
	return;
    }
    if (!strncmp (input, "Yellow", 6)) {
	strncpy (input, "1 1 0\0", 6);
	strcat (input, transp);
	return;
    }

/* try to find it in colors.inc, if povray options include +L path or if POV_Include is non-blank */
    colinc = NULL;
    if (strlen (drvui->POV_Include) > 10) {
	colinc = fopen (drvui->POV_Include, "r");
    }
    if (!colinc) {
	strncpy (b, drvui->POV_Options, 29);
	b[29] = '\0';
	a = strstr (b, "+L");
	if (a != NULL) {
	    sscanf (a + 2, "%s", POV_incpath);
	    if (POV_incpath != NULL) {
		strcat (POV_incpath, "/colors.inc");
		colinc = fopen (POV_incpath, "r");
		if (colinc) {
		    strcpy (drvui->POV_Include, POV_incpath);
		    WriteConfig ();
		}
	    }
	}
    }
    if (colinc != NULL) {
	if (transp_start != NULL) 
	    *(transp_start - 1) = '\0';
	while (!feof (colinc)) {
	    if (fgets (line, 80, colinc) == NULL)
		continue;
	    sscanf (line, "%s %s = color red %f green %f blue %f", keyword, thecolor,
		    &thered, &thegreen, &theblue);
	    if (strcmp (keyword, "#declare"))
		continue;
	    if (!strcmp (thecolor, input)) {
		snprintf (input, 30, "%3.2f %3.2f %3.2f%c", thered, thegreen, theblue, 0);
		if (transp_start == NULL) {
		    transp[0] = '\0';
		    transp_start = strstr (line, "filter");
		    if (transp_start != NULL) {
			fval = (float) atof (transp_start + 6);
			sprintf (transp, " transparency %5.2f", fval / 2.);
		    }
		}
		strcat (input, transp);
		fclose (colinc);
		return;
	    }
	}
	fclose (colinc);
    }

    if (Color_Warning <= 5) {
	char string[256];

	sprintf (string, " **** Warning - Non-standard color \"%s\" encountered. ****",
		 input);
	Error_Box (string);
	Color_Warning++;
    }

    if (sscanf (input, " %f %f %f", &thered, &thegreen, &theblue) < 3) {
	strcpy (input, "1. 1. 1.");	// substitute White if no color triplet
    } else {
	snprintf (input, 30, "%3.2f %3.2f %3.2f%c", thered, thegreen, theblue, 0);
	strcat (input, transp);
    }
}

/* ************************************************************** */
/* ************************************************************** */

void
trim_string (char string[], int len)
{
// trim string - remove trailing ^J (Windows) and spaces. 'len' is maximum length
    int i;

    string[len - 1] = 0;
    if (strlen (string) == 0) {
//        strcpy(string,"White");
	return;
    }
    for (i = 0; i < (int) strlen (string); i++) {
	if ((unsigned char) string[i] < ' ') {
	    string[i] = 0;
	    break;
	}
    }
    for (i = strlen (string); i > 0; --i) {
	if (string[i - 1] <= ' ') {
	    string[i - 1] = 0;
	} else
	    break;
    }
}

/* ************************************************************** */
/* ************************************************************** */

int
Unique_Atom (void)
// routine to check if atom at position 'natom' is unique. Returns 1 if it
//   is. If the position is a duplicate, returns 0.
{
    int j, k, test;

    for (k = 0; k < natom - 1; k++) {	/* loop through atoms */
	if (drvui->atoms[k].atom_fn != drvui->frame_no)
	    continue;
	if (drvui->atoms[k].atom_n == drvui->atoms[natom].atom_n) {
	    test = 0;
	    for (j = 0; j < 4; j++) {
		if (drvui->atoms[k].atom_l[j] != drvui->atoms[natom].atom_l[j])
		    test = 1;
	    }
	    if (!test) {	/* name and number the same - check coordinates */
		expand_atom (natom);	/* find all of test atoms in unit cell */
		for (j = 0; j < ncell; j++) {
		    if ((fabs (drvui->atoms[k].atom_xyz[1] - drvui->cell_xyz[j][1]) <
			 0.0001)
			&& (fabs (drvui->atoms[k].atom_xyz[2] - drvui->cell_xyz[j][2]) <
			    0.0001)
			&& (fabs (drvui->atoms[k].atom_xyz[0] - drvui->cell_xyz[j][0]) <
			    0.0001))
			return (0);	/* duplicate */
		}
	    }
	}
    }
    return (1);			/* unique */
}

/* ************************************************************** */
/* ************************************************************** */

int
vec_dif (int n1, float v1[3], int n2, float v2[3], int n3, float v3[3], float v[3])
{
/* routine to determine if vector v == n1*v1+n2*v2+n3*v3, returns 1 if equal, 0 if not */
    int i;

    float d[3];

    for (i = 0; i < 3; i++)
	d[i] = v[i] - n1 * v1[i] - n2 * v2[i] - n3 * v3[i];
    return d[0] * d[0] + d[1] * d[1] + d[2] * d[2] < 1.0e-4 ? 1 : 0;
}
