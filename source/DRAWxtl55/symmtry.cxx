// $Id: symmtry.cxx 900 2009-08-13 20:00:45Z larry $
//
// module symmetry.cxx - part of DRAWxtl V5.5
//
//     Larry W. Finger, Martin Kroeker and Brian Toby
//
// symmetry generator codes
//
// routines contained within this file:
//
// symop - interprets Hermann-Mauguin Space Group Symbol
// getsys - determins crystal system (cubic, tetragonal, etc.) from space group symbol
// monoclinic - interprets monoclinic space group symbol
// orthorhombic - interprets orthorhombic space group symbol
// tetragonal - interprets tetragonal space group symbol
// hexagonal - interprets hexagonal and trigonal space group symbols
// cubic - interprets cubic space group symbol
//
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "DRAWxtlViewUI.h"

void symop (char *);

int getsys (int);

void monoclinic (char[3][4]);

void orthorhombic (char[3][4], int);

void tetragonal (char[3][4], int, int);

void hexagonal (char[3][4], int);

void cubic (char[3][4], int, int);

extern DRAWxtlViewUI *drvui;


/* ************************************************************** */
/* ************************************************************** */

void
symop (char input[100])
/* routine to interpret Hermann-Mauguin Space Group Symbol --
   Program adapted from a routine supplied by Prof. Burzlaff, Univ.
   Erlangen, Germany.  This version obtained from the LAZY-PULVERIX
   source code in FORTRAN and converted to C by L.W. Finger, 19-Apr-1992 */
/* 'sys' is crystal system code: 1-triclinic, 2-monoclinic, 3-orthorhombic,
   4-tetragonal, 5-hexagonal, 6-cubic

   'acentric' is centric indicator: 0-centric, 1-acentric

   'nbr' is Bravais lattice code

   'ng' is number of rotation operators in system

*/
{
    char ibra[8];

    int nbl, i, j, k, ic, nlq, iflg;

    strcpy (ibra, "PABCFIR");
    nbl = -1;
    drvui->acentric = iflg = 1;
    ic = nlq = 0;
/* initialize arrays */
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
    drvui->nbr = 0;
    drvui->nlat = 0;


/* Count number of fields in space group symbol */

    for (i = 1; i <= 36; ++i) {
	if (iflg == 0) {
	    if ((input[i - 1] == '\n') || (input[i - 1] == '\0'))
		for (j = i; j <= 36; ++j)
		    input[j - 1] = ' ';
	    if (input[i - 1] == ' ') {
		nlq = 0;
		ic = 0;
	    } else {
		if (nbl < 0) {
		    for (j = 1; j <= 7; ++j) {
			if (toupper (input[i - 1]) == ibra[j - 1])
			    drvui->nbr = j;
		    }
		    switch (drvui->nbr) {	/* Fill in lattice positions */
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
			drvui->lat_pos[1][0] = 0.5;
			drvui->lat_pos[1][1] = 0.0;
			drvui->lat_pos[1][2] = 0.5;
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
			drvui->lat_pos[1][0] = 0.5;
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
		    nbl = 0;
		} else {
		    if (nlq == 0)
			nbl = nbl + 1;
		    ic = ic + 1;
		    drvui->spg[nbl - 1][ic - 1] = input[i - 1];
		    nlq = 1;
		    if (input[i - 1] == '/')
			drvui->acentric = 0;
		}
	    }
	} else {
	    if (input[i - 1] == ' ')
		iflg = 0;
	}
    }
/* get system code from routine getsys */
    drvui->sys = getsys (nbl);
    switch (drvui->sys) {
    case 1:			/* triclinic */
	if (drvui->spg[0][0] == '-')
	    drvui->acentric = 0;
	drvui->ng = 1;
	break;
    case 2:
	monoclinic (drvui->spg);
	drvui->ng = 2;
	break;
    case 3:
	orthorhombic (drvui->spg, drvui->nbr);
	drvui->ng = 4;
	break;
    case 4:
	tetragonal (drvui->spg, nbl, drvui->nbr);
	break;
    case 5:
	hexagonal (drvui->spg, nbl);
	break;
    case 6:
	cubic (drvui->spg, nbl, drvui->nbr);
	break;
    default:
	drvui->ng = 0;		/* Invalid System  */
	return;
    }
/* get translation operators in range from 0 to 1 */
    for (i = 1; i <= drvui->ng; ++i) {
	for (j = 1; j <= 3; ++j) {
	    while (drvui->ts[i - 1][j - 1] >= 1.0)
		drvui->ts[i - 1][j - 1] = drvui->ts[i - 1][j - 1] - 1.0f;
	    while (drvui->ts[i - 1][j - 1] < 0.0)
		drvui->ts[i - 1][j - 1] = drvui->ts[i - 1][j - 1] + 1.0f;
	}
    }
}

/* ************************************************************** */
/* ************************************************************** */

int
getsys (int nbl)
/* Routine to take partially disassembled space group symbol in 'spg'
   and return the number of the system.  The number of blanks
   in the original symbol is given in 'nbl'.  */
{

    int i, system;		/* loop counter */

    for (i = 1; i <= 4; ++i) {
	if (drvui->spg[1][i - 1] == '3') {	/* '3' in second place => cubic */
	    return (6);
	}
    }
    for (i = 1; i <= 4; ++i) {
	if ((drvui->spg[0][i - 1] == '3') || (drvui->spg[0][i - 1] == '6')) {
/* '3' or '6' in first position => hexagonal */
	    return (5);
	}
	if (drvui->spg[0][i - 1] == '4') {	/* '4' in first position => tetragonal */
	    return (4);
	}
    }
    if (nbl <= 1) {
	if ((drvui->spg[0][0] == '1') || (drvui->spg[0][0] == '-')) {
	    return (1);		/* triclinic */
	}
	system = 2;
	for (i = 1; i <= 4; ++i)
	    drvui->spg[1][i - 1] = drvui->spg[0][i - 1];
	for (i = 0; i <= 3; ++i)
	    drvui->spg[0][i] = ' ';
	drvui->spg[0][0] = '1';
	drvui->spg[2][0] = '1';
    }
    system = 3;
    if ((drvui->spg[0][0] == '1') || (drvui->spg[1][0] == '1'))
	system = 2;
    return (system);
}

/* ************************************************************** */
/* ************************************************************** */

void
monoclinic (char spg[3][4])
/* routine to interpret a monoclinic space group symbol */
{
    int i, j, k, id, ind = 0;

    for (i = 1; i <= 3; ++i) {
	if (spg[i - 1][0] != '1')
	    ind = i;
    }
    id = 1;
    if (spg[ind - 1][0] == '2')
	id = -1;

    if (drvui->nbr == 3)
	ind++;			// support non-standard B settings

    for (i = 1; i <= 3; ++i)
	drvui->ss[1][i - 1][i - 1] = id * drvui->ss[0][i - 1][i - 1];
    drvui->ss[1][ind - 1][ind - 1] = -drvui->ss[1][ind - 1][ind - 1];
    for (i = 1; i <= 3; ++i) {
	if ((spg[i - 1][0] == '2') && (spg[i - 1][1] == '1'))
	    drvui->ts[1][i - 1] = 0.5;
	for (j = 1; j <= 4; ++j) {
	    if (spg[i - 1][j - 1] == 'a')
		drvui->ts[1][0] = 0.5;
	    if (spg[i - 1][j - 1] == 'b')
		drvui->ts[1][1] = 0.5;
	    if (spg[i - 1][j - 1] == 'c')
		drvui->ts[1][2] = 0.5;
	    if (spg[i - 1][j - 1] == 'n') {
		k = i + 1;
		if (k > 3)
		    k = k - 3;
		drvui->ts[1][5 - k - i] = 0.5;
		if (drvui->nbr == 3)
		    k--;
		drvui->ts[1][k - 1] = 0.5;
	    }
	}
    }
}

/* ************************************************************** */
/* ************************************************************** */

void
orthorhombic (char spg[3][4], int nbr)
/* module to generate space group operators for orthorhombic space groups */
{

    float sh[3], tc;

    int i, j, k, l, m, id = 0, jnd, ind, ic, ma, nma[3];

    sh[0] = 0.25;
    sh[1] = 0.25;
    sh[2] = 0.25;
    ind = 1;
    ic = 0;
    if ((spg[0][0] != '2') && (spg[1][0] != '2') && (spg[2][0] != '2')) {
	ind = -1;
	drvui->acentric = 0;
    }
// Convert names of space groups with double glide plane to old symbol
// (mapping Aem2->Abm2, Aea2->Aba2, Cmce->Cmca,Cmme->Cmma,Ccce->Ccca)
    if (spg[0][0] == 'e')
	spg[0][0] = 'b';
    if (spg[2][0] == 'e')
	spg[2][0] = 'a';

    for (i = 1; i <= 3; ++i) {
	id = 1;
	if (spg[i - 1][0] == '2')
	    id = -1;
	for (j = 1; j <= 3; ++j)
	    drvui->ss[i][j - 1][j - 1] = id * ind * drvui->ss[0][j - 1][j - 1];
	drvui->ss[i][i - 1][i - 1] = -drvui->ss[i][i - 1][i - 1];
    }
    for (i = 1; i <= 3; ++i) {
	if ((spg[i - 1][0] == '2') && (spg[i - 1][1] == '1')) {
	    drvui->ts[i][i - 1] = 0.5;
	    k = (i % 3) + 1;	/* Fix BUG in original code for P m n 21  */
	    l = (k % 3) + 1;
	    jnd = 0;
	    if ((spg[k - 1][0] == 'm') && (spg[l - 1][0] == 'n')) {
		drvui->ts[i][k - 1] = 0.5;
		jnd = 1;
	    }
	    if ((spg[k - 1][0] == 'n') && (spg[l - 1][0] == 'm')) {
		drvui->ts[i][l - 1] = 0.5;
		jnd = 1;
	    }
	    if (jnd == 0) {
		if ((spg[k - 1][0] == 'n') || (spg[l - 1][0] == 'n')) {
		    drvui->ts[k][k - 1] = 0.5;	/* Fix BUG in original code for P n a 21 */
		    drvui->ts[l][l - 1] = 0.5;
		}
	    }
	}
	for (j = 1; j <= 4; ++j) {
	    if (spg[i - 1][j - 1] == 'a')
		drvui->ts[i][0] = 0.5;
	    if (spg[i - 1][j - 1] == 'b')
		drvui->ts[i][1] = 0.5;
	    if (spg[i - 1][j - 1] == 'c')
		drvui->ts[i][2] = 0.5;
	    if ((spg[i - 1][j - 1] == 'n') || (spg[i - 1][j - 1] == 'd')) {
		k = i + 1;
		if (k > 3)
		    k = k - 3;
		if (spg[i - 1][j - 1] != 'd') {
		    drvui->ts[i][k - 1] = 0.5;
		    drvui->ts[i][5 - k - i] = 0.5;
		} else {
		    ic = 1;
		    if (drvui->acentric == 1) {
			drvui->ts[i][0] = 0.25;
			drvui->ts[i][1] = 0.25;
			drvui->ts[i][2] = .25;
		    } else {
			drvui->ts[i][k - 1] = 0.25;
			drvui->ts[i][5 - k - i] = 0.25;
		    }
		}
	    }
	}
    }

// Check for space groups with multiple choices for the origin

    if (spg[0][0] == 'd' && spg[1][0] == 'd' && spg[2][0] == 'd') {
	drvui->origin1_flag = 1;	//  F d d d
	for (k = 0; k < 3; k++)
	    drvui->origin_offset[k] = -0.125f;
    }
    if (spg[0][0] == 'n' && spg[1][0] == 'n' && spg[2][0] == 'n') {
	drvui->origin1_flag = 1;	//  P n n n
	for (k = 0; k < 3; k++)
	    drvui->origin_offset[k] = -0.25f;
    }
    if (spg[0][0] == 'b' && spg[1][0] == 'a' && spg[2][0] == 'n') {
	drvui->origin1_flag = 1;	// P b a n
	drvui->origin_offset[0] = -0.25f;
	drvui->origin_offset[1] = -0.25f;
	drvui->origin_offset[2] = 0.0;
    }
    if (spg[0][0] == 'n' && spg[1][0] == 'n' && spg[2][0] == 'b') {
	drvui->origin1_flag = 1;	// P n c b
	drvui->origin_offset[0] = 0.0f;
	drvui->origin_offset[1] = -0.25f;
	drvui->origin_offset[2] = -0.25f;
    }
    if (spg[0][0] == 'c' && spg[1][0] == 'n' && spg[2][0] == 'a') {
	drvui->origin1_flag = 1;	// P c n a
	drvui->origin_offset[0] = -0.25f;
	drvui->origin_offset[1] = 0.0f;
	drvui->origin_offset[2] = -0.25f;
    }
    if (spg[0][0] == 'm' && spg[1][0] == 'm' && spg[2][0] == 'n') {
	drvui->origin1_flag = 1;	// P m m n
	drvui->origin_offset[0] = -0.25f;
	drvui->origin_offset[1] = -0.25f;
	drvui->origin_offset[2] = 0.0f;
    }
    if (spg[0][0] == 'm' && spg[1][0] == 'n' && spg[2][0] == 'm') {
	drvui->origin1_flag = 1;	// P m n m
	drvui->origin_offset[0] = -0.25f;
	drvui->origin_offset[1] = 0.0f;
	drvui->origin_offset[2] = -0.25f;
    }
    if (spg[0][0] == 'n' && spg[1][0] == 'm' && spg[2][0] == 'm') {
	drvui->origin1_flag = 1;	// P n m m
	drvui->origin_offset[0] = 0.0f;
	drvui->origin_offset[1] = -0.25f;
	drvui->origin_offset[2] = -0.25f;
    }
    if (spg[0][0] == 'c' && spg[1][0] == 'c' && spg[2][0] == 'a') {
	drvui->origin1_flag = 1;	// C c c a
	drvui->origin_offset[0] = 0.0f;
	drvui->origin_offset[1] = -0.25f;
	drvui->origin_offset[2] = -0.25f;
    }
    if (spg[0][0] == 'c' && spg[1][0] == 'c' && spg[2][0] == 'b') {
	drvui->origin1_flag = 1;	// C c c b
	drvui->origin_offset[0] = -0.25f;
	drvui->origin_offset[1] = 0.0f;
	drvui->origin_offset[2] = -0.25f;
    }
    if (spg[0][0] == 'b' && spg[1][0] == 'a' && spg[2][0] == 'a') {
	drvui->origin1_flag = 1;	// A b a a
	drvui->origin_offset[0] = -0.25f;
	drvui->origin_offset[1] = 0.0f;
	drvui->origin_offset[2] = -0.25f;
    }
    if (spg[0][0] == 'c' && spg[1][0] == 'a' && spg[2][0] == 'a') {
	drvui->origin1_flag = 1;	// A c a a
	drvui->origin_offset[0] = -0.25f;
	drvui->origin_offset[1] = -0.25f;
	drvui->origin_offset[2] = 0.0f;
    }
    if (spg[0][0] == 'b' && spg[1][0] == 'c' && spg[2][0] == 'b') {
	drvui->origin1_flag = 1;	// B b c b
	drvui->origin_offset[0] = -0.25f;
	drvui->origin_offset[1] = -0.25f;
	drvui->origin_offset[2] = 0.0f;
    }
    if (spg[0][0] == 'b' && spg[1][0] == 'a' && spg[2][0] == 'b') {
	drvui->origin1_flag = 1;	// B b a b
	drvui->origin_offset[0] = 0.0f;
	drvui->origin_offset[1] = -0.25f;
	drvui->origin_offset[2] = -0.25f;
    }
    if (ic == 1)
	return;
    if (drvui->acentric != 1) {
	for (i = 1; i <= 3; ++i) {
	    k = 1 + i;
	    if (k > 3)
		k = k - 3;
	    tc = drvui->ts[k][i - 1] + drvui->ts[6 - i - k][i - 1];
	    if (tc == 1.0)
		tc = 0.0;
	    drvui->ts[i][i - 1] = tc;
	}
	if ((nbr == 1) || (nbr == 5))
	    return;
/* Special treatment of Cmma, Cmca, Imma */
	ma = 0;
	for (i = 1; i <= 3; ++i) {
	    nma[i - 1] = 0;
	    for (j = 1; j <= 4; ++j)
		if (spg[i - 1][j - 1] == 'm')
		    nma[i - 1] = 1;
	    ma = ma + nma[i - 1];
	}
	if ((nbr != 6) || (ma != 2)) {
	    if ((ma == 0) || (ma == 3) || (nbr == 6))
		return;
	    if (nma[nbr - 2] == 1)
		return;
	    sh[nbr - 2] = 0.0;
	}
/* Origin shift */
	for (i = 1; i <= 4; ++i) {
	    for (j = 1; j <= 3; ++j) {
		for (k = 1; k <= 3; ++k) {
		    id = 1;
		    if (j != k)
			id = 0;
		    drvui->ts[i - 1][j - 1] =
			drvui->ts[i - 1][j - 1] + (id -
						   drvui->ss[i - 1][j - 1][k -
									   1]) * sh[k -
										    1];
		}
		if (drvui->ts[i - 1][j - 1] >= 1.0)
		    drvui->ts[i - 1][j - 1] -= 1.0f;
		if (drvui->ts[i - 1][j - 1] < 1.0)
		    drvui->ts[i - 1][j - 1] = drvui->ts[i - 1][j - 1] + 1.0f;
	    }
	}
	return;
    } else {
	ic = 0;
	for (i = 1; i <= 3; ++i)
	    if (drvui->ss[i][0][0] * drvui->ss[i][1][1] * drvui->ss[i][2][2] == -1)
		ic = 1;
	if (ic != 1) {
	    tc = drvui->ts[1][0] + drvui->ts[2][1] + drvui->ts[3][2];
	    if (tc == 0)
		return;
	    for (i = 1; i <= 3; ++i) {	/* do 40 */
		k = i + 1;
		if (k > 3)
		    k = k - 3;
		if (tc <= 0.5) {
		    if (drvui->ts[i][i - 1] != 0.0) {
			m = i - 1;
			if (m == 0)
			    m = 3;
			drvui->ts[m][i - 1] = 0.5;
		    }
		} else {
		    if (tc <= 1.0) {
			if (drvui->ts[i][i - 1] == 0.0) {
			    l = k + 1;
			    if (l > 3)
				l = l - 3;
			    drvui->ts[k][l - 1] = 0.5;
			    drvui->ts[l][k - 1] = 0.5;
			}

		    } else {
			drvui->ts[i][k - 1] = 0.5;
		    }
		}
	    }
	    return;
	} else {
	    for (i = 1; i <= 3; ++i)	/* do 42 */
		if (drvui->ss[i][0][0] * drvui->ss[i][1][1] * drvui->ss[i][2][2] == 1)
		    id = i;
	    for (i = 1; i <= 3; ++i) {	/* do 45 */
		tc = drvui->ts[1][i - 1] + drvui->ts[2][i - 1] + drvui->ts[3][i - 1];
		if ((tc != 0.0) && (tc != 1.0)) {
		    if (((spg[0][0] == 'm') && (spg[1][0] == 'n'))
			|| ((spg[1][0] == 'm') && (spg[2][0] == 'n'))
			|| ((spg[2][0] == 'm') && (spg[0][0] == 'n'))) {
			k = i - 1;
			if (k == 0)
			    k = 3;
			drvui->ts[k][i - 1] = 0.5;
		    } else {
			for (j = 1; j <= 3; ++j)	/* do 43 */
			    if (id != j)
				drvui->ts[j][i - 1] = 0.5;
		    }
		}
	    }
	}
    }
}

/* ************************************************************** */
/* ************************************************************** */

void
tetragonal (char spg[3][4], int nbl, int nbr)
{
    int i, j, k, l;		/* local loop variables */

    int e[3][3], ne;

    float te[3];

    drvui->ng = 4;
    if (nbl == 3)
	drvui->ng = 8;
    drvui->ss[1][0][1] = -1;
    drvui->ss[1][1][0] = 1;
    drvui->ss[1][2][2] = 1;
    if (spg[0][0] == '-') {
	for (i = 0; i <= 2; ++i)
	    for (j = 0; j <= 2; ++j)
		drvui->ss[1][i][j] = -drvui->ss[1][i][j];
    } else {
	if (spg[0][1] == '1')
	    drvui->ts[1][2] = 0.25;
	if (spg[0][1] == '3')
	    drvui->ts[1][2] = 0.75;
	if (spg[0][1] == '2')
	    drvui->ts[1][2] = 0.5;
	if ((spg[0][2] == 'n') || ((spg[0][3] == 'n') && (nbl == 3)))
	    drvui->ts[1][0] = 0.5;
	if (((spg[0][3] == 'n') && (nbl == 1)) || ((spg[0][1] == '1') &&
						   (drvui->acentric == 1) && (nbr == 6)))
	    drvui->ts[1][1] = 0.5;
	if ((spg[0][1] == '1') && (!drvui->acentric) && (nbr == 6)) {
	    if (nbl == 1) {
		drvui->ts[1][0] = 0.75;
		drvui->ts[1][1] = 0.25;
	    } else {
		drvui->ts[1][0] = 0.25;
		drvui->ts[1][1] = 0.75;
	    }
	} else {
	    if ((spg[1][1] == '1') || (spg[1][1] == '3') || ((spg[0][3] != 'n')
							     && (spg[1][0] == 'n')
							     && (spg[2][0] == 'm'))) {
		drvui->ts[1][0] = 0.5;
		drvui->ts[1][1] = 0.5;
	    }
	}
    }
    if (spg[0][2] == 'n') {
	drvui->origin1_flag = 1;	// P4/n**
	drvui->origin_offset[0] = -0.25f;
	drvui->origin_offset[1] = -0.25f;
	if (spg[1][0] == 'n')
	    drvui->origin_offset[2] = -0.25f;	// P4/n n c
	else
	    drvui->origin_offset[2] = 0.0f;	// Other P4/n**
    }
    if (spg[0][1] == '2' && spg[0][3] == 'n') {
	drvui->origin1_flag = 1;	// P42/n * * space groups
	for (k = 0; k < 3; k++)
	    drvui->origin_offset[k] = -0.25f;
    }
    if (spg[0][1] == '1' && spg[0][3] == 'a') {
	drvui->origin1_flag = 1;	// I41/a * * space groups
	drvui->origin_offset[0] = 0.0f;
	drvui->origin_offset[1] = -0.25f;
	drvui->origin_offset[2] = -0.125f;
    }
    drvui->ss[2][0][0] = -1;	/*50 */
    drvui->ss[2][1][1] = -1;
    drvui->ss[2][2][2] = 1;
    drvui->ts[2][0] = drvui->ss[1][0][1] * drvui->ts[1][1] + drvui->ts[1][0];
    drvui->ts[2][1] = drvui->ss[1][1][0] * drvui->ts[1][0] + drvui->ts[1][1];
    drvui->ts[2][2] = drvui->ss[1][2][2] * drvui->ts[1][2] + drvui->ts[1][2];
    if ((nbr == 6) && (drvui->acentric == 1) && (spg[0][1] == '1'))
	drvui->ts[2][2] = 0.5;
    for (i = 0; i <= 2; ++i) {
	drvui->ts[3][i] = drvui->ts[1][i];
	for (j = 0; j <= 2; ++j) {
	    drvui->ts[3][i] = drvui->ts[3][i] + drvui->ss[1][i][j] * drvui->ts[2][j];
	    for (k = 0; k <= 2; ++k)
		drvui->ss[3][i][j] =
		    drvui->ss[3][i][j] + drvui->ss[1][i][k] * drvui->ss[2][k][j];
	}
    }
    if (nbl == 1)
	return;
    for (i = 0; i <= 2; ++i)
	for (j = 0; j <= 2; ++j)
	    e[j][i] = 0;
    for (i = 0; i <= 2; ++i)
	te[i] = 0.0;
    if (!drvui->acentric) {
	e[0][0] = -1;
	e[1][1] = 1;
	e[2][2] = 1;
	if ((spg[1][0] == 'c') || (spg[1][0] == 'n'))
	    te[2] = 0.5;
	if ((spg[1][0] == 'b') || (spg[1][0] == 'n')) {
	    te[1] = 0.5;
	    te[0] = 0.5;
	}
	if ((spg[0][2] == 'n') || (spg[0][3] == 'n'))
	    te[0] = te[0] + 0.5f;
    } else {
	if ((spg[1][0] != '2') && (spg[2][0] != '2')) {	/*55 */
	    e[0][0] = -1;
	    e[1][1] = 1;
	    e[2][2] = 1;
	    if ((spg[1][0] == 'c') || (spg[1][0] == 'n'))
		te[2] = 0.5;
	    if ((spg[1][0] == 'n') || (spg[1][0] == 'b')) {
		te[0] = 0.5;
		te[1] = 0.5;
	    }
	} else {
	    if ((spg[1][0] == '2') && (spg[2][0] == '2')) {	/*54 */
		e[0][1] = 1;
		e[1][0] = 1;
		e[2][2] = -1;
		if ((spg[1][1] == ' ') && (nbr != 6) && (spg[0][1] != ' ')) {
		    if (spg[0][1] == '1')
			te[2] = 0.75;
		    if (spg[0][1] == '2')
			te[2] = 0.5;
		    if (spg[0][1] == '3')
			te[2] = 0.25;
		}
	    } else {
		if (spg[1][0] == '2') {	/*53 */
		    e[0][0] = 1;
		    e[1][1] = -1;
		    e[2][2] = -1;
		    if (spg[2][0] == 'c')
			te[2] = 0.5;
		    if (spg[2][0] == 'd') {
			te[2] = 0.25;
			te[1] = 0.5;
		    }
		    if (spg[1][1] == '1') {
			te[0] = 0.5;
			te[1] = 0.5;
		    }
		} else {
		    if ((spg[1][0] == 'c') || (spg[1][0] == 'n'))
			te[2] = 0.5;
		    e[0][0] = -1;
		    e[1][1] = 1;
		    e[2][2] = 1;
		    if ((spg[1][0] == 'n') || (spg[1][0] == 'b')) {
			te[0] = 0.5;
			te[1] = 0.5;
		    }
		}
	    }
	}
    }
    ne = 4;
    for (i = 0; i < ne; ++i) {
	for (j = 0; j <= 2; ++j) {
	    drvui->ts[ne + i][j] = te[j];
	    for (k = 0; k <= 2; ++k) {
		drvui->ts[ne + i][j] = drvui->ts[ne + i][j] + e[j][k] * drvui->ts[i][k];
		for (l = 0; l <= 2; ++l)
		    drvui->ss[ne + i][j][k] =
			drvui->ss[ne + i][j][k] + e[j][l] * drvui->ss[i][l][k];
	    }
	}
    }
}

/* ************************************************************** */
/* ************************************************************** */

void
hexagonal (char spg[3][4], int nbl)
{

    int e[3][3], i, j, k, l, ne;	/* local variables */

    float te[3];

    drvui->ng = 3;
    if ((spg[0][0] == '-') && (spg[0][1] == '3'))
	drvui->acentric = 0;	/* -3 ... */
    if (spg[0][0] != '6') {
	drvui->ss[1][0][1] = -1;
	drvui->ss[1][1][0] = 1;
	drvui->ss[1][1][1] = -1;
	drvui->ss[1][2][2] = 1;
	if (spg[0][1] == '1')
	    drvui->ts[1][2] = (1.f / 3.f);
	if (spg[0][1] == '2')
	    drvui->ts[1][2] = (2.f / 3.f);
	drvui->ss[2][0][0] = -1;
	drvui->ss[2][1][0] = -1;
	drvui->ss[2][0][1] = 1;
	drvui->ss[2][2][2] = 1;
	drvui->ts[2][2] = 2.0f * drvui->ts[1][2];
	if (drvui->ts[2][2] >= 1.0)
	    drvui->ts[2][2] = drvui->ts[2][2] - 1.0f;
	if ((nbl == 1) && (spg[0][1] != '6'))
	    return;
	if (spg[0][1] == '6') {
	    drvui->ng = 2 * drvui->ng;
	    for (i = 0; i <= 2; ++i)
		for (j = 0; j <= 2; ++j)
		    for (k = 0; k <= 2; ++k) {
			drvui->ss[3 + i][j][k] = drvui->ss[i][j][k];
			drvui->ss[3 + i][2][2] = -1;
		    }
	}
	if (nbl == 1)
	    return;
	if ((spg[1][0] == 'c') || (spg[2][0] == 'c')) {
	    drvui->ts[3][2] = 0.5;
	    drvui->ts[4][2] = 0.5;
	    drvui->ts[5][2] = 0.5;
	}
    } else {			/*61 */
	drvui->ng = 2 * drvui->ng;
	drvui->ss[1][0][0] = 1;
	drvui->ss[1][0][1] = -1;
	drvui->ss[1][1][0] = 1;
	drvui->ss[1][2][2] = 1;
	if (spg[0][1] == '1')
	    drvui->ts[1][2] = (1.0f / 6.0f);
	if (spg[0][1] == '2')
	    drvui->ts[1][2] = (2.0f / 6.0f);
	if (spg[0][1] == '3')
	    drvui->ts[1][2] = (3.0f / 6.0f);
	if (spg[0][1] == '4')
	    drvui->ts[1][2] = (4.0f / 6.0f);
	if (spg[0][1] == '5')
	    drvui->ts[1][2] = (5.0f / 6.0f);
	for (i = 0; i <= 3; ++i) {
	    for (j = 0; j <= 2; ++j) {
		drvui->ts[2 + i][j] = drvui->ts[1][j];
		for (k = 0; k <= 2; ++k) {
		    drvui->ts[2 + i][j] =
			drvui->ts[2 + i][j] + drvui->ss[1][j][k] * drvui->ts[i + 1][k];
		    if (drvui->ts[2 + i][j] >= 1.0)
			drvui->ts[2 + i][j] = drvui->ts[2 + i][j] - 1.0f;
		    for (l = 0; l <= 2; ++l)
			drvui->ss[2 + i][j][k] =
			    drvui->ss[2 + i][j][k] + drvui->ss[1][j][l] * drvui->ss[1 +
										    i][l]
			    [k];
		}
	    }
	}
	if (nbl == 1)
	    return;
    }
    for (i = 0; i <= 2; ++i) {
	te[i] = 0.0;
	for (j = 0; j <= 2; ++j)
	    e[j][i] = 0;
    }
    drvui->ng = 2 * drvui->ng;
    if (spg[1][0] != '1') {
	if (spg[1][0] != '2') {
	    e[0][1] = -1;
	    e[1][0] = -1;
	    e[2][2] = 1;
	    if (spg[1][0] == 'c')
		te[2] = 0.5;
	} else {		/*64 */
	    e[0][1] = 1;
	    e[1][0] = 1;
	    e[2][2] = -1;
	    te[2] = 2.0f * drvui->ts[1][2];
	    if ((spg[0][0] == '3') && ((spg[0][1] == '1') || (spg[0][1] == '2')))
		te[2] = 0.0;
	}
    } else {			/*65 */
	if (spg[2][0] != '2') {
	    e[0][1] = 1;
	    e[1][0] = 1;
	    e[2][2] = 1;
	    if (spg[2][0] == 'c')
		te[2] = 0.5f;
	} else {		/*66 */
	    e[0][1] = -1;
	    e[1][0] = -1;
	    e[2][2] = -1;
	    te[2] = 2.0f * drvui->ts[1][2];
	    if (te[2] > 1.0)
		te[2] = te[2] - 1.0f;
	}
    }
    ne = 6;
    if ((spg[0][0] == '3') || ((spg[0][1] == '3') && (spg[0][0] == '-')))
	ne = 3;
    for (i = 0; i < ne; ++i) {
	for (j = 0; j <= 2; ++j) {
	    drvui->ts[ne + i][j] = te[j];
	    for (k = 0; k <= 2; ++k) {
		drvui->ts[ne + i][j] = drvui->ts[ne + i][j] + e[j][k] * drvui->ts[i][k];
		for (l = 0; l <= 2; ++l)
		    drvui->ss[ne + i][j][k] =
			drvui->ss[ne + i][j][k] + e[j][l] * drvui->ss[i][l][k];
	    }
	}
    }
}

/* ************************************************************** */
/* ************************************************************** */

void
cubic (char spg[3][4], int nbl, int nbr)
{
    int i, j, k, l, m, ne, e[3][3];	/* local scratch */

    float te[3];

    if (nbl == 3)
	drvui->ng = 24;
    else
	drvui->ng = 12;
    memset (te, 0, sizeof (te));
    if ((spg[0][0] != '2') && (spg[0][0] != '4') && (spg[0][0] != '-'))
	drvui->acentric = 0;
    for (i = 0; i <= 2; ++i) {
	for (j = 0; j <= 2; ++j) {
	    drvui->ss[i + 1][j][j] = 1;
	    if (j != i) {
		drvui->ss[i + 1][j][j] = -1;
		if (spg[0][0] == 'n') {	// P n * *
		    drvui->ts[i + 1][j] = 0.5;
		    drvui->origin1_flag = 1;
		    for (k = 0; k < 3; k++)
			drvui->origin_offset[k] = -0.25f;
		}
		if (spg[0][0] == 'd') {
		    drvui->ts[i + 1][j] = 0.25;	// F d 3 *
		    drvui->origin1_flag = 1;	// alternate origin for these space groups
		    for (k = 0; k < 3; k++)
			drvui->origin_offset[k] = -0.125f;
		}
	    }
	}
    }
    if (((spg[0][0] == 'a') || (spg[2][0] == 'd') || (spg[0][1] == '3')
	 || (spg[0][1] == '1')) && (nbr != 5)) {
	for (i = 1; i <= 3; ++i) {
	    drvui->ts[i][i - 1] = 0.5;
	    k = i + 1;
	    if (k == 4)
		k = 1;
	    drvui->ts[i][k - 1] = 0.5;
	}
    }
    for (i = 0; i <= 3; ++i) {
	for (j = 0; j <= 2; ++j) {
	    for (k = 0; k <= 2; ++k) {
		l = j + 1;
		if (l == 3)
		    l = 0;
		m = j - 1;
		if (m < 0)
		    m = 2;
		drvui->ss[i + 4][j][k] = drvui->ss[i][l][k];
		drvui->ss[i + 8][j][k] = drvui->ss[i][m][k];
		drvui->ts[i + 4][j] = drvui->ts[i][l];
		drvui->ts[i + 8][j] = drvui->ts[i][m];
	    }
	}
    }
    if (drvui->ng == 12)
	return;
    ne = 12;
    memset (te, 0, sizeof (te));
    memset (e, 0, sizeof (e));
    e[0][1] = 1;
    e[1][0] = 1;
    e[2][2] = 1;
    if (spg[2][0] == '2')
	e[2][2] = -1;
    if (spg[2][0] == 'c')
	te[2] = 0.5;
    for (i = 0; i <= 2; ++i) {	/* do 73 */
	if ((spg[2][0] == 'n') || (spg[0][1] == '2'))
	    te[i] = 0.5;
	if ((spg[0][1] == '3') || (spg[0][1] == '1'))
	    te[i] = 0.25;
	if (spg[2][0] == 'd')
	    te[i] = 0.75;
    }
    if ((spg[0][1] == '1') && (nbr == 1))
	te[0] = 0.75;
/*     if ((spg[0][1] == '3') && (nbr == 1))
	 te[0] = 0.25;  ??????  */
    if (((spg[0][1] == '1') && (nbr == 6))
	|| ((spg[0][1] == '3') && (nbr == 1))) {
	te[1] = 0.75;
	te[2] = 0.75;
    }
    for (i = 0; i < ne; ++i) {
	for (j = 0; j <= 2; ++j) {
	    drvui->ts[ne + i][j] = te[j];
	    for (k = 0; k <= 2; ++k) {
		drvui->ts[ne + i][j] = drvui->ts[ne + i][j] + e[j][k] * drvui->ts[i][k];
		for (l = 0; l <= 2; ++l)
		    drvui->ss[ne + i][j][k] =
			drvui->ss[ne + i][j][k] + e[j][l] * drvui->ss[i][l][k];
	    }
	}
    }
}
