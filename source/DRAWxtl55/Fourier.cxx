// $Id: Fourier.cxx 1102 2011-01-03 15:17:35Z martin $
//
// Fourier.cxx - Source module for DRAWxtl V5.5
// Coded using the FLTK 1.1.6 widget set
//
//     Larry W. Finger, Martin Kroeker and Brian Toby
//
//
// This module contains the following routines:
//
//  fillcube - read a 512 byte record of a dn6 Fourier map
//  read_dn6 - read O Format (dn6) Fourier file
//  read_fcf - read CIF-standard Fo/Fc file and compute density
//  read_grd - read GSAS (grd) Fourier file
//  read_stf - read JANA (stf) Fourier file
//  read_m80 - read JANA (m80) Fourier file
//  read_m81 - read JANA (m81) density file
//  read_flp - read FullProf (flp) Fourier file
//  read_w2k - read WIEN2k calculated charge density file
//  read_exc - read EXCITING calculated charge density or ELF file
//  read_vasp - read VASP calculated charge density file
//  read_aim - read WIEN2k calculated Bader AIM surface file
//  read_xsf - read charge density file in XCrysDen format

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

/* ************************************************************** */
/* ************************************************************** */

int
fillcube (uchar * cube, FILE * mapin)
{
/* routine to read the next record of a dn6 map.
 *
 * The first record is a header. After that, each record contains
 * an 8 x 8 x 8 piece of the map scaled so that each map point
 * occupies a single byte.
 */

    int i, val;

    for (i = 0; i < 512; i++) {
	if ((val = fgetc (mapin)) == EOF)
	    return 0;
	cube[i] = (uchar) val;
    }
    return 1;
}

/* ************************************************************** */
/* ************************************************************** */

void
read_dn6 (char *infile, int Quick)
{
/* routine to read a dn6 or 'O' format Fourier map */
    FILE *mapin;

    float rhomin = 1.0e15f, rhomax = -1.0e15f;

    char string[200];

    uchar cube[512];

    float rhoscal;

    int rhoadd;

    int iscal1, iscal2;

    int i, j, k;

    int ori[3];

    int nga, ngb, ngc;

    int kk, jj, ii;

    int ncube;

    int mappos;

    int nz, ny, nx;

    if ((mapin = fopen (infile, "rb")) == NULL) {
	sprintf (string, "Cannot open Fourier map (.dn6) file %s", infile);
	Error_Box (string);
	return;
    }
    if (!fillcube (cube, mapin)) {
	Error_Box ("Error reading DN6 Fourier File.");
	fclose (mapin);
	return;
    }
/* TODO confirm the description of ori */
    ori[0] = (int) cube[0] * 256 + (int) cube[1];	/* get starting point of map */
    ori[1] = (int) cube[2] * 256 + (int) cube[3];
    ori[2] = (int) cube[4] * 256 + (int) cube[5];
    nga = (int) cube[6] * 256 + (int) cube[7];
    ngb = (int) cube[8] * 256 + (int) cube[9];
    ngc = (int) cube[10] * 256 + (int) cube[11];
    mapstep_a = (int) cube[12] * 256 + (int) cube[13];
    mapstep_b = (int) cube[14] * 256 + (int) cube[15];
    mapstep_c = (int) cube[16] * 256 + (int) cube[17];
    iscal1 = (int) cube[34] * 256 + (int) cube[35];
    map_a = (float) ((int) cube[18] * 256 + (int) cube[19]) / (float) iscal1;
    map_b = (float) ((int) cube[20] * 256 + (int) cube[21]) / (float) iscal1;
    map_c = (float) ((int) cube[22] * 256 + (int) cube[23]) / (float) iscal1;
    map_alpha = (float) ((int) cube[24] * 256 + (int) cube[25]) / (float) iscal1;
    map_beta = (float) ((int) cube[26] * 256 + (int) cube[27]) / (float) iscal1;
    map_gamma = (float) ((int) cube[28] * 256 + (int) cube[29]) / (float) iscal1;
    iscal2 = (int) cube[36] * 256 + (int) cube[37];
    rhoscal = (float) ((int) cube[30] * 256 + (int) cube[31]) / (float) iscal2;
    rhoadd = (int) cube[32] * 256 + (int) cube[33];
    Map_Info.lat_con[0] = map_a;
    Map_Info.lat_con[1] = map_b;
    Map_Info.lat_con[2] = map_c;
    Map_Info.lat_con[3] = map_alpha;
    Map_Info.lat_con[4] = map_beta;
    Map_Info.lat_con[5] = map_gamma;
    Map_Info.map_int[0] = mapstep_a;
    Map_Info.map_int[1] = mapstep_b;
    Map_Info.map_int[2] = mapstep_c;
    Map_Info.xlim[0] = (float) ori[0] / (float) mapstep_a;
    Map_Info.ylim[0] = (float) ori[1] / (float) mapstep_b;
    Map_Info.zlim[0] = (float) ori[2] / (float) mapstep_c;
    Map_Info.xlim[1] = (float) (nga + ori[0] - 1) / (float) mapstep_a++;
    Map_Info.ylim[1] = (float) (ngb + ori[1] - 1) / (float) mapstep_b++;
    Map_Info.zlim[1] = (float) (ngc + ori[2] - 1) / (float) mapstep_c++;

    Map_Info.info_valid = 1;

    if (FourierPt)
	free (FourierPt);

    FourierPt = (float *) zalloc (mapstep_a * mapstep_b * mapstep_c * sizeof (float));
    if (FourierPt == NULL) {
	Error_Box ("ERROR -- Unable to allocate space for map.\n");
	fclose (mapin);
	return;
    }

/* read the Rho values - map file has x fastest, and z slowest,
 * In FourierPt, we store with z fastest */

    for (k = 0; k < ngc; k = k + 8) {
	for (j = 0; j < ngb; j = j + 8) {
	    for (i = 0; i < nga; i = i + 8) {
		if (!fillcube (cube, mapin)) {	/* read the next 8 x 8 x 8 block */
		    Error_Box ("Error reading DN6 Fourier File.");
		    fclose (mapin);
		    return;
		}
		for (kk = 0; kk < min (8, ngc - k); kk++) {	/* unpack the cube */
		    nz = k + kk + ori[2];
		    for (jj = 0; jj < min (8, ngb - j); jj++) {
			ny = (j + jj + ori[1]) * mapstep_c;
			for (ii = 0; ii < min (8, nga - i); ii = ii + 2) {
			    nx = (i + ii + ori[0]) * mapstep_c * mapstep_b;
			    mappos = nx + ny + nz;
			    ncube = 64 * kk + 8 * jj + ii;
			    FourierPt[mappos] =
				(float) (cube[ncube + 1] - rhoadd) / rhoscal;
			    if (FourierPt[mappos] > rhomax)
				rhomax = FourierPt[mappos];
			    if (FourierPt[mappos] < rhomin)
				rhomin = FourierPt[mappos];
			    if (ii + 1 < min (8, nga - i)) {
				FourierPt[mappos + 1] =
				    (float) (cube[ncube] - rhoadd) / rhoscal;
				if (FourierPt[mappos + 1] > rhomax)
				    rhomax = FourierPt[mappos + 1];
				if (FourierPt[mappos + 1] < rhomin)
				    rhomin = FourierPt[mappos + 1];
			    }
			}
		    }
		}
	    }
	}
    }
    if (!Quick) {
	fprintf (drvui->flout, "Reading Fourier map file %s\n"
		 "cell = %f %f %f %f %f %f\n  map steps: a=%d b=%d c=%d\n  "
		 "Range of rho values: %f to %f\n", infile,
		 map_a, map_b, map_c, map_alpha, map_beta, map_gamma,
		 mapstep_a, mapstep_b, mapstep_c, rhomin, rhomax);
    }
    Map_Info.rhomn = rhomin;
    Map_Info.rhomx = rhomax;
    ReadFourMap = 1;
    fclose (mapin);
}

/* ************************************************************** */
/* ************************************************************** */

void
read_fcf (char *infile, int Quick)
{
// read CIF-standard Fo/Fc file and compute density

    FILE *mapin;

    float *fo, *fc, *f2, *ftmp1, *ftmp2, *ftmp3;

    short *ih, *ik, *il, *itmp1, *itmp2, *itmp3;

    float *Cklx;

    float *Dklx;

    float *Elxy;

    float *Flxy;

    double x, y, z;

    char string[200];

    char maptitle[61];

    int squared, phases, ab;

    int nval = 200;

    int nr;

    int i = 0, j, k, l;

    int ijk = 0;

    float rhomin = 1.0e15f, rhomax = -1.0e15f;

    float factor;

    short ih0, ik0, il0;

    short ih1, ik1, il1;

    short im, in, io;

    float fo0, fc0, fc1, f20, f21;

    int kmin, kmax, nk;

    int lmin, lmax, nl;

    float phase;

    int skip;

    int progress;

    int dimension = 3;

    if (FourierPt != NULL) {
	return;
    }

    if ((mapin = fopen (infile, "r")) == NULL) {
	sprintf (string, "Cannot open Fo/Fc (CIF fcf) file %s", infile);
	Error_Box (string);
	return;
    }

    fo = (float *) zalloc (200 * sizeof (float));
    if (fo == NULL) {
	sprintf (string, "Unable to obtain storage for Fo data");
	Error_Box (string);
	fclose (mapin);
	return;
    }
    fc = (float *) zalloc (200 * sizeof (float));
    if (fc == NULL) {
	sprintf (string, "Unable to obtain storage for Fc data");
	Error_Box (string);
	free (fo);
	fclose (mapin);
	return;
    }
    f2 = (float *) zalloc (200 * sizeof (float));
    if (f2 == NULL) {
	sprintf (string, "Unable to obtain storage for B or phase data");
	Error_Box (string);
	free (fo);
	free (fc);
	fclose (mapin);
	return;
    }
    ih = (short *) zalloc (200 * sizeof (short));
    if (ih == NULL) {
	sprintf (string, "Unable to obtain storage for h data");
	Error_Box (string);
	free (fo);
	free (fc);
	free (f2);
	fclose (mapin);
	return;
    }
    ik = (short *) zalloc (200 * sizeof (short));
    if (ik == NULL) {
	sprintf (string, "Unable to obtain storage for k data");
	Error_Box (string);
	free (fo);
	free (fc);
	free (f2);
	free (ih);
	fclose (mapin);
	return;
    }
    il = (short *) zalloc (200 * sizeof (short));
    if (il == NULL) {
	sprintf (string, "Unable to obtain storage for l data");
	Error_Box (string);
	free (fo);
	free (fc);
	free (f2);
	free (ih);
	free (ik);
	fclose (mapin);
	return;
    }

    do {
	fgets (string, 80, mapin);
	if (strstr (string, "title")) {
	    memset (maptitle, 0, 61);
	    sscanf (string, "%*s %60c", maptitle);
	    trim_string (maptitle, 60);
	    strcpy (Map_Info.title, maptitle);
	    Map_Info.info_valid = 1;
	}
    } while (strncmp (string, " _refln", 6) && !feof(mapin) );	// search for start of _refln group

    if (strncmp (string, " _refln_index_h", 14)) {
	sprintf (string, "Unsupported field sequence in fcf file - need h k l first");
	Error_Box (string);
	free (fo);
	free (fc);
	free (f2);
	free (ih);
	free (ik);
	free (il);
	fclose (mapin);
	return;
    }

    squared = 0;
    phases = 0;
    ab = 0;
    while (!strncmp (string, " _refln", 6)) {
	fgets (string, 80, mapin);
	if (strstr (string, "squared"))
	    squared = 1;
	if (strstr (string, "_phase_"))
	    phases = 1;
	if (strstr (string, "A_calc"))
	    ab = 1;
	if (strstr (string, "index_m"))
	    dimension = 4;
	if (strstr (string, "index_n"))
	    dimension = 5;
	if (strstr (string, "index_o"))
	    dimension = 6;
    }

    if (dimension > 3) {
	sprintf (string,
		 "%d-dimensional data encountered - only main reflections will be used",
		 dimension);
	Error_Box (string);
    }

    if (phases == 0 && ab == 0) {
	sprintf (string, "fcf file does not contain phase angle or A,B data");
	Error_Box (string);
	free (fo);
	free (fc);
	free (f2);
	free (ih);
	free (ik);
	free (il);
	fclose (mapin);
	return;
    }

    nr = 0;
    do {
	skip = 0;
	if (dimension > 3) {
	    switch (dimension) {
	    case 4:
		im = 0;
		i = sscanf (string, "%hd %hd %hd %hd %f %f %f", &ih0, &ik0, &il0, &im,
			    &fo0, &fc0, &f20);
		if (im != 0)
		    skip = 1;
		break;
	    case 5:
		im = in = 0;
		i = sscanf (string, "%hd %hd %hd %hd %hd %f %f %f", &ih0, &ik0, &il0, &im,
			    &in, &fo0, &fc0, &f20);
		if (im != 0 || in != 0)
		    skip = 1;
		break;
	    case 6:
		im = in = io = 0;
		i = sscanf (string, "%hd %hd %hd %hd %hd %hd %f %f %f",
			    &ih0, &ik0, &il0, &im, &in, &io, &fo0, &fc0, &f20);
		if (im != 0 || in != 0 || io != 0)
		    skip = 1;
		break;
	    }
	} else {
	    i = sscanf (string, "%hd %hd %hd %f %*f %f %f", &ih0, &ik0, &il0, &fo0, &fc0,
			&f20);
	}
	if (i < dimension + 3)
	    break;
	if (skip == 1) {
	    fgets (string, 80, mapin);
	    continue;
	}

	if (squared == 1)
	    fo0 = (float) sqrt (fo0);
	if (ab) {		// change all data to 'phase' type
	    fc1 = (float) sqrt (fc0 * fc0 + f20 * f20);
	    f20 = (float) atan2 (f20 / fc1, fc0 / fc1);
	    fc0 = fc1;
	} else {
	    f20 *= (float) (PI / 180.0);	// phase in radians
	}

// apply spacegroup symmetry (including xyz - catch duplicates in input list)
	for (k = 0; k < drvui->ng; ++k) {
	    skip = 0;

	    ih1 =
		drvui->ss[k][0][0] * ih0 + drvui->ss[k][1][0] * ik0 +
		drvui->ss[k][2][0] * il0;
	    ik1 =
		drvui->ss[k][0][1] * ih0 + drvui->ss[k][1][1] * ik0 +
		drvui->ss[k][2][1] * il0;
	    il1 =
		drvui->ss[k][0][2] * ih0 + drvui->ss[k][1][2] * ik0 +
		drvui->ss[k][2][2] * il0;
	    phase =
		2.0f * (float) PI *(drvui->ts[k][0] * ih0 + drvui->ts[k][1] * ik0 +
				    drvui->ts[k][2] * il0);
	    f21 = f20 + phase;
	    f21 = (float) atan2 (sin (f21), cos (f21));	// get f21 in range -PI to PI

	    for (j = 0; j < nr; j++) {
		if (ih[j] == ih1 && ik[j] == ik1 && il[j] == il1) {
		    skip = 1;
		    break;
		}
	    }
	    if (skip == 0) {
		ih[nr] = ih1;
		ik[nr] = ik1;
		il[nr] = il1;
		fo[nr] = fo0;
		fc[nr] = fc0;
		f2[nr] = f21;
		nr++;

	    }
	    ih1 *= -1;		// add Friedel opposite
	    ik1 *= -1;
	    il1 *= -1;
	    f21 *= -1.0;
	    skip = 0;

	    for (j = 0; j < nr; j++) {
		if (ih[j] == ih1 && ik[j] == ik1 && il[j] == il1) {
		    skip = 1;
		    break;
		}
	    }
	    if (skip == 0) {
		ih[nr] = ih1;
		ik[nr] = ik1;
		il[nr] = il1;
		fo[nr] = fo0;
		fc[nr] = fc0;
		f2[nr] = f21;
		nr++;
	    }

	    if (nr > nval - 10) {
		nval += 100;
		ftmp1 = (float *) realloc (fo, nval * sizeof (float));
		ftmp2 = (float *) realloc (fc, nval * sizeof (float));
		ftmp3 = (float *) realloc (f2, nval * sizeof (float));
		itmp1 = (short *) realloc (ih, nval * sizeof (short));
		itmp2 = (short *) realloc (ik, nval * sizeof (short));
		itmp3 = (short *) realloc (il, nval * sizeof (short));
		if (!ftmp1 || !ftmp2 || !ftmp3 || !itmp1 || !itmp2 || !itmp3) {
		    sprintf (string,
			     "Unable to expand storage for h, k, l, fo, fc or phase data");
		    free (fo);
		    free (fc);
		    free (f2);
		    free (ih);
		    free (ik);
		    free (il);
		    if (ftmp1)
			free (ftmp1);
		    if (ftmp2)
			free (ftmp2);
		    if (ftmp3)
			free (ftmp3);
		    if (itmp1)
			free (itmp1);
		    if (itmp2)
			free (itmp2);
		    if (itmp3)
			free (itmp3);
		    Error_Box (string);
		    fclose (mapin);
		    return;
		} else {
		    fo = ftmp1;
		    fc = ftmp2;
		    f2 = ftmp3;
		    ih = itmp1;
		    ik = itmp2;
		    il = itmp3;
		}
	    }			// if we need to realloc
	}			// loop over all symmetry operators
	fgets (string, 80, mapin);
    } while (i > 0 && !feof (mapin));

    fclose (mapin);


    factor = 1.0f / (drvui->lat_con[0] * drvui->lat_con[1] * drvui->lat_con[2]
		     * (float) sqrt (1.0 - cos (drvui->lat_con[3] * PI / 180.0) *
				     cos (drvui->lat_con[3] * PI / 180.0) *
				     cos (drvui->lat_con[4] * PI / 180.0) *
				     cos (drvui->lat_con[4] * PI / 180.0) *
				     cos (drvui->lat_con[5] * PI / 180.0) *
				     cos (drvui->lat_con[5] * PI / 180.0)));
    kmin = kmax = ik[0];
    lmin = lmax = il[0];

    for (l = 0; l < nr; l++) {
	if (ik[l] < kmin)
	    kmin = ik[l];	// find range of k and l
	if (ik[l] > kmax)
	    kmax = ik[l];
	if (il[l] < lmin)
	    lmin = il[l];
	if (il[l] > lmax)
	    lmax = il[l];
	switch (Map_Info.map_type) {	// get map coefficients/V
	default:
	case 0:		// Fo map
	    break;
	case 1:		// Fc map
	    fo[l] = fc[l];
	    break;
	case 2:		// Fo - Fc map
	    fo[l] = fo[l] - fc[l];
	    break;
	case 3:		// 2Fo - Fc map
	    fo[l] = 2.0f * fo[l] - fc[l];
	    break;
	case 4:		// Fo2 (Patterson)
	    fo[l] = fo[l] * fo[l];
	    f2[l] = 0.0f;
	    break;
	}
	fo[l] *= factor;	// update coefficients times 1/V
    }

    nk = kmax - kmin + 1;
    nl = lmax - lmin + 1;
    Cklx = (float *) zalloc (nk * nl * sizeof (float));
    Dklx = (float *) zalloc (nk * nl * sizeof (float));
    Elxy = (float *) zalloc (nl * sizeof (float));
    Flxy = (float *) zalloc (nl * sizeof (float));

    if (!Cklx || !Dklx || !Elxy || !Flxy) {
	Error_Box ("Unable to allocate Beevers-Lipson arrays.");
	if (Cklx) free (Cklx);
	if (Dklx) free (Dklx);
	if (Elxy) free (Elxy);
	if (Flxy) free (Flxy);
	return;
    }
    mapstep_a = (int) ((float)Map_Info.res * drvui->lat_con[0] + 0.5f);
    mapstep_b = (int) ((float)Map_Info.res * drvui->lat_con[1] + 0.5f);
    mapstep_c = (int) ((float)Map_Info.res * drvui->lat_con[2] + 0.5f);

    FourierPt = (float *) zalloc (mapstep_a * mapstep_b * mapstep_c * sizeof (float));
    if (FourierPt == NULL) {
	Error_Box ("ERROR -- Unable to allocate space for map.");
	return;
    }
// start the Beevers-Lipson expansion

    Progress_Window (-1, "Computing Density", (float) mapstep_a * mapstep_b);

    progress = 0;
    for (i = 0; i < mapstep_a; i++) {
	x = (float) i / (float) mapstep_a;
	memset (Cklx, 0, nk * nl * sizeof (float));
	memset (Dklx, 0, nk * nl * sizeof (float));
	for (l = 0; l < nr; l++) {
	    ijk = (il[l] - lmin) * nk + ik[l] - kmin;
	    Cklx[ijk] += fo[l] * (float) cos (2.0 * PI * ih[l] * x - f2[l]);
	    Dklx[ijk] -= fo[l] * (float) sin (2.0 * PI * ih[l] * x - f2[l]);
	}
	for (j = 0; j < mapstep_b; j++) {
	    Progress_Window (0, NULL, (float) ++progress);
	    y = (float) j / (float) mapstep_b;
	    memset (Elxy, 0, nl * sizeof (float));
	    memset (Flxy, 0, nl * sizeof (float));
	    for (l = 0; l < nl; l++) {
		for (k = 0; k < nk; k++) {
		    ijk = l * nk + k;
		    Elxy[l] += Cklx[ijk] * (float) cos (2.0 * PI * (k + kmin) * y)
			+ Dklx[ijk] * (float) sin (2.0 * PI * (k + kmin) * y);
		    Flxy[l] += Dklx[ijk] * (float) cos (2.0 * PI * (k + kmin) * y)
			- Cklx[ijk] * (float) sin (2.0 * PI * (k + kmin) * y);
		}
	    }
	    for (k = 0; k < mapstep_c; k++) {
		z = (float) k / (float) mapstep_c;
		ijk = i * mapstep_b * mapstep_c + j * mapstep_c + k;

		for (l = 0; l < nl; l++) {
		    FourierPt[ijk] += Elxy[l] * (float) cos (2.0 * PI * (l + lmin) * z)
			+ Flxy[l] * (float) sin (2.0 * PI * (l + lmin) * z);
		}
		if (FourierPt[ijk] < rhomin)
		    rhomin = FourierPt[ijk];
		if (FourierPt[ijk] > rhomax)
		    rhomax = FourierPt[ijk];
	    }
	}
    }
    if (!Quick)
	fprintf (drvui->fcns, "Reading Fourier map file %s\n  "
		 "map steps: a=%d b=%d c=%d\n  "
		 "Range of rho values: %f to %f, gridsize %d\n", infile,
		 mapstep_a, mapstep_b, mapstep_c, rhomin, rhomax, ijk);


    Progress_Window (-2, NULL, 0.0f);
    ReadFourMap = 1;
    Map_Info.rhomn = rhomin;
    Map_Info.rhomx = rhomax;
    Map_Info.map_int[0] = mapstep_a;
    Map_Info.map_int[1] = mapstep_b;
    Map_Info.map_int[2] = mapstep_c;
    Map_Info.lat_con[0] = drvui->lat_con[0];
    Map_Info.lat_con[1] = drvui->lat_con[1];
    Map_Info.lat_con[2] = drvui->lat_con[2];
    Map_Info.lat_con[3] = drvui->lat_con[3];
    Map_Info.lat_con[4] = drvui->lat_con[4];
    Map_Info.lat_con[5] = drvui->lat_con[5];
    Map_Info.xlim[0] = 0.0f;
    Map_Info.xlim[1] = 1.0f;
    Map_Info.ylim[0] = 0.0f;
    Map_Info.ylim[1] = 1.0f;
    Map_Info.zlim[0] = 0.0f;
    Map_Info.zlim[1] = 1.0f;

    free (Cklx);
    free (Dklx);
    free (Elxy);
    free (Flxy);
    free (fo);
    free (fc);
    free (f2);
    free (ih);
    free (ik);
    free (il);
}

/* ************************************************************** */
/* ************************************************************** */

void
read_m81 (char *infile, int Quick)
{
// read JANA binary .m81 Fo or jpdf file 
    FILE *mapin;

    char string[200];

    int iheaders[6];

    float fheaders[6];

    int i, j, k, l, ijk;

    int nskip, nskip5, nskip6;

    int nrec, nmap, maptype;

    int astart, aend, bstart, bend, cstart, cend;

    int step1, step2, step3, step4, step5, step6;

    int new_mapstep_a, new_mapstep_b, new_mapstep_c;

    int ireslt = 0;

    int fullcell;

    float min1, max1, min2, max2, min3, max3;

    float min4, max4, min5, max5, min6, max6;

    float xmin, xmax, ymin, ymax, zmin, zmax;

    float x4min = 0.f, x4max = 0.f, x5min = 0.f, x5max = 0.f, x6min = 0.f, x6max = 0.f;

    float xstep4, xstep5, xstep6;

    float rhomin, rhomax;

    float *ftmp;

    int axis[6], axes, modaxes;
    const char *theaxes[] = { " ", "x", "y", "z", "x4", "x5", "x6" };
    const char *maptypes[] = { "Patterson", "checking Patterson", "difference Patterson",
	"Fourier", "checking Fourier", "difference Fourier", "shape function"
    };

    if (FourierPt)
	free (FourierPt);

    if ((mapin = fopen (infile, "rb")) == NULL) {
	sprintf (string, "Cannot open JANA m81 file %s", infile);
	Error_Box (string);
	return;
    }
    i = fread (iheaders, 4, 6, mapin);	// nx ny nz nx4 nx5 nx6
    step1 = iheaders[0];
    step2 = iheaders[1];
    step3 = iheaders[2];
    step4 = iheaders[3];
    step5 = iheaders[4];
    step6 = iheaders[5];
    i = fread (iheaders, 4, 2, mapin);	// nx*ny  nmaps
    nrec = iheaders[0];
    nmap = iheaders[1];
    i = fread (fheaders, 4, 6, mapin);	//xmin..xmax...zmin..zmax
    min1 = fheaders[0];
    max1 = fheaders[1];
    min2 = fheaders[2];
    max2 = fheaders[3];
    min3 = fheaders[4];
    max3 = fheaders[5];
    i = fread (fheaders, 4, 6, mapin);	// qmin..qmax...smin..smax (ignored)
    min4 = fheaders[0];
    max4 = fheaders[1];
    min5 = fheaders[2];
    max5 = fheaders[3];
    min6 = fheaders[4];
    max6 = fheaders[5];
    i = fread (fheaders, 4, 6, mapin);	// step sizes (ignored)
    xstep4 = fheaders[3];
    xstep5 = fheaders[4];
    xstep6 = fheaders[5];
    if (max4 == 0.0f && xstep4 == 1.0f)
	xstep4 = 0.0f;
    if (max5 == 0.0f && xstep5 == 1.0f)
	xstep5 = 0.0f;
    if (max6 == 0.0f && xstep6 == 1.0f)
	xstep6 = 0.0f;
    i = fread (iheaders, 4, 6, mapin);	// axis sequence, e.g. 2 1 3 
    axis[0] = iheaders[0];
    axis[1] = iheaders[1];
    axis[2] = iheaders[2];
    axis[3] = iheaders[3];
    axis[4] = iheaders[4];
    axis[5] = iheaders[5];
    axes = 100 * axis[0] + 10 * axis[1] + axis[2];
    modaxes = 100 * axis[3] + 10 * axis[4] + axis[5];
    i = fread (iheaders, 4, 2, mapin);	// composite number, maptype
    maptype = iheaders[0] - 1;

    ftmp = (float *) zalloc (nrec * sizeof (float));
    if (ftmp == NULL) {
	Error_Box ("ERROR -- Unable to allocate space for map.\n");
	fclose (mapin);
	return;
    }

    switch (axes) {
    case 123:
	mapstep_a = step1;
	mapstep_b = step2;
	mapstep_c = step3;
	xmin = min1;
	xmax = max1;
	ymin = min2;
	ymax = max2;
	zmin = min3;
	zmax = max3;
	break;
    case 132:
	mapstep_a = step1;
	mapstep_b = step3;
	mapstep_c = step2;
	xmin = min1;
	xmax = max1;
	ymin = min3;
	ymax = max3;
	zmin = min2;
	zmax = max2;
	break;
    case 312:
	mapstep_a = step2;
	mapstep_b = step3;
	mapstep_c = step1;
	xmin = min2;
	xmax = max2;
	ymin = min3;
	ymax = max3;
	zmin = min1;
	zmax = max1;
	break;
    case 321:
	mapstep_a = step3;
	mapstep_b = step2;
	mapstep_c = step1;
	xmin = min3;
	xmax = max3;
	ymin = min2;
	ymax = max2;
	zmin = min1;
	zmax = max1;
	break;
    case 213:
	mapstep_a = step2;
	mapstep_b = step1;
	mapstep_c = step3;
	xmin = min2;
	xmax = max2;
	ymin = min1;
	ymax = max1;
	zmin = min3;
	zmax = max3;
	break;
    case 231:
	mapstep_a = step3;
	mapstep_b = step1;
	mapstep_c = step2;
	xmin = min3;
	xmax = max3;
	ymin = min1;
	ymax = max1;
	zmin = min2;
	zmax = max2;
	break;
    default:
	sprintf (string, "Primary axes must be some permutation of x y z,\n"
		 "while this file has %s %s %s %s %s %s.\n", theaxes[axis[0]],
		 theaxes[axis[1]], theaxes[axis[2]], theaxes[axis[3]], theaxes[axis[4]],
		 theaxes[axis[5]]);
	Error_Box (string);
	free (ftmp);
	fclose (mapin);
	return;
	break;
    }
    switch (modaxes) {
    case 0:
	break;
    case 400:
	x4step = xstep4;
	x4min = min4;
	x4max = max4;
	x5min = 0.;
	x5max = 0.;
	x6min = 0.;
	x6max = 0.;
	break;
    case 450:
	x4step = xstep4;
	x4min = min4;
	x4max = max4;
	x5step = xstep5;
	x5min = min5;
	x5max = max5;
	x6min = 0.;
	x6max = 0.;
	break;
    case 540:
	x4step = xstep5;
	x4min = min5;
	x4max = max5;
	x5step = xstep4;
	x5min = min4;
	x5max = max4;
	x6min = 0.;
	x6max = 0.;
	break;
    case 456:
	x4step = xstep4;
	x4min = min4;
	x4max = max4;
	x5step = xstep5;
	x5min = min5;
	x5max = max5;
	x6step = xstep6;
	x6min = min6;
	x6max = max6;
	break;
    case 465:
	x4step = xstep4;
	x4min = min4;
	x4max = max4;
	x5step = xstep6;
	x5min = min6;
	x5max = max6;
	x6step = xstep5;
	x6min = min5;
	x6max = max5;
	break;
    case 546:
	x4step = xstep5;
	x4min = min5;
	x4max = max5;
	x5step = xstep4;
	x5min = min4;
	x5max = max4;
	x6step = xstep6;
	x6min = min6;
	x6max = max6;
	break;
    case 564:
	x4step = xstep6;
	x4min = min6;
	x4max = max6;
	x5step = xstep4;
	x5min = min4;
	x5max = max4;
	x6step = xstep5;
	x6min = min5;
	x6max = max5;
	break;
    case 654:
	x4step = xstep6;
	x4min = min6;
	x4max = max6;
	x5step = xstep5;
	x5min = min5;
	x5max = max5;
	x6step = xstep4;
	x6min = min4;
	x6max = max4;
	break;
    case 645:
	x4step = xstep5;
	x4min = min5;
	x4max = max5;
	x5step = xstep6;
	x5min = min6;
	x5max = max6;
	x6step = xstep4;
	x6min = min4;
	x6max = max4;
	break;
    default:
	Error_Box ("Unhandled sequence of modulated axes.\n");
	free (ftmp);
	fclose (mapin);
	return;
	break;
    }

    // skip to end of first record - record length is max(nrec,34)*4 

    nskip = 34;
    if (nrec > nskip)
	nskip = nrec;
    nskip -= 34;		// header values already read
    i = fread (ftmp, 4, nskip, mapin);

    Map_Info.xlim[0] = xmin;
    Map_Info.xlim[1] = xmax;
    Map_Info.ylim[0] = ymin;
    Map_Info.ylim[1] = ymax;
    Map_Info.zlim[0] = zmin;
    Map_Info.zlim[1] = zmax;

    fullcell = 1;
    if ((fabs (Map_Info.xlim[0]) > 0.001) || (fabs (Map_Info.xlim[1] - 1.0f) > 0.001)) {
	// not full cell along x
	new_mapstep_a =
	    (int) ((1.0f / (Map_Info.xlim[1] - Map_Info.xlim[0])) * (float) mapstep_a +
		   0.5f);
	astart = (int) (Map_Info.xlim[0] * new_mapstep_a + 0.5f);
	aend = (int) (Map_Info.xlim[1] * new_mapstep_a + 0.5f);
	fullcell = 0;
	if (astart < 0 && aend - astart != mapstep_a)
	    aend++;
	if (astart > 0 && aend - astart != mapstep_a)
	    aend--;
	if (aend - astart != mapstep_a)
	    Error_Box ("Error recalculating map a limits for full cell");
    } else {
	new_mapstep_a = mapstep_a;
	astart = 0;
	aend = mapstep_a;
    }
    if ((fabs (Map_Info.ylim[0]) > 0.001) || (fabs (Map_Info.ylim[1] - 1.0f) > 0.001)) {
	// not full cell along y
	new_mapstep_b =
	    (int) ((1.0f / (Map_Info.ylim[1] - Map_Info.ylim[0])) * (float) mapstep_b +
		   0.5f);
	bstart = (int) (Map_Info.ylim[0] * new_mapstep_b + 0.5f);
	bend = (int) (Map_Info.ylim[1] * new_mapstep_b + 0.5f);
	fullcell = 0;
	if (bstart < 0 && bend - bstart != mapstep_b)
	    bend++;
	if (bstart > 0 && bend - bstart != mapstep_b)
	    bend--;
	if (bend - bstart != mapstep_b)
	    Error_Box ("Error recalculating map b limits for full cell");
    } else {
	new_mapstep_b = mapstep_b;
	bstart = 0;
	bend = mapstep_b;
    }
    if ((fabs (Map_Info.zlim[0]) > 0.001) || (fabs (Map_Info.zlim[1] - 1.0f) > 0.001)) {
	// not full cell along z
	new_mapstep_c =
	    (int) ((1.0f / (Map_Info.zlim[1] - Map_Info.zlim[0])) * (float) mapstep_c +
		   0.5f);
	cstart = (int) (Map_Info.zlim[0] * new_mapstep_c + 0.5f);
	cend = (int) (Map_Info.zlim[1] * new_mapstep_c + 0.5f);
	fullcell = 0;
	if (cstart < 0 && cend - cstart != mapstep_c)
	    cend++;
	if (cstart > 0 && cend - cstart != mapstep_c)
	    cend--;
	if (cend - cstart != mapstep_c)
	    Error_Box ("Error recalculating map c limits for full cell");
    } else {
	new_mapstep_c = mapstep_c;
	cstart = 0;
	cend = mapstep_c;
    }
//    if (FourierPt) free(FourierPt);
    FourierPt =
	(float *) zalloc (new_mapstep_a * new_mapstep_b * new_mapstep_c * sizeof (float));
    if (FourierPt == NULL) {
	ireslt = 0;
	Error_Box ("ERROR -- Unable to allocate space for map.");
	free (ftmp);
	fclose (mapin);
	return;
    }

    switch (modaxes) {
    case 0:
	break;
    case 400:
	if (x4Val > 0.0f && x4Val <= x4max) {
	    nskip = (int) ((x4Val - x4min) / x4step + .5);
	    for (i = 0; i < nskip; i++) {
		for (j = 0; j < step3; j++) {
		    ireslt = fread (ftmp, 4, nrec, mapin);
		    if (ireslt <= 0) {
//			fprintf (stderr, "error parsing superspace map\n");
			Error_Box ("Error seeking section in superspace map.");
			free (ftmp);
			fclose (mapin);
			return;
		    }
		}
	    }
	}
	break;
    case 450:
	if (x5Val > 0.0f && x5Val <= x5max) {
	    nskip5 = (int) ((x5Val - x5min) / x5step + .5);
	    for (i = 0; i < nskip5; i++) {
		for (j = 0; j < step4; j++) {
		    for (k = 0; k < step3; k++) {
			ireslt = fread (ftmp, 4, nrec, mapin);
			if (ireslt <= 0) {
//			    fprintf (stderr, "error parsing superspace map\n");
			    Error_Box ("Error seeking section in superspace map.");
			    free (ftmp);
			    fclose (mapin);
			    return;
			}
		    }
		}
	    }
	}
	if (x4Val > 0.0f && x4Val <= x4max) {
	    nskip = (int) ((x4Val - x4min) / x4step + .5);
	    for (i = 0; i < nskip; i++) {
		for (j = 0; j < step3; j++) {
		    ireslt = fread (ftmp, 4, nrec, mapin);
		    if (ireslt <= 0) {
//			fprintf (stderr, "error parsing superspace map\n");
			Error_Box ("Error seeking section in superspace map.");
			free (ftmp);
			fclose (mapin);
			return;
		    }
		}
	    }
	}
	break;
    case 540:
	if (x4Val > 0.0f && x4Val <= x4max) {
	    nskip = (int) ((x4Val - x4min) / x4step + .5);
	    for (i = 0; i < nskip; i++) {
		for (j = 0; j < step5; j++) {
		    for (k = 0; k < step3; k++) {
			ireslt = fread (ftmp, 4, nrec, mapin);
			if (ireslt <= 0) {
//			    fprintf (stderr, "error parsing superspace map\n");
			    Error_Box ("Error seeking section in superspace map.");
			    free (ftmp);
			    fclose (mapin);
			    return;
			}
		    }
		}
	    }
	}
	if (x5Val > 0.0f && x5Val <= x5max) {
	    nskip5 = (int) ((x5Val - x5min) / x5step + .5);
	    for (i = 0; i < nskip5; i++) {
		for (j = 0; j < step3; j++) {
		    ireslt = fread (ftmp, 4, nrec, mapin);
		    if (ireslt <= 0) {
//			fprintf (stderr, "error parsing superspace map\n");
			Error_Box ("Error seeking section in superspace map.");
			free (ftmp);
			fclose (mapin);
			return;
		    }
		}
	    }
	}
	break;
    case 456:
	if (x6Val > 0.0f && x6Val <= x6max) {
	    nskip6 = (int) ((x6Val - x6min) / x6step + .5);
	    for (i = 0; i < nskip6; i++) {
		for (j = 0; j < step5; j++) {
		    for (k = 0; k < step4; k++) {
			for (l = 0; l < step3; l++) {
			    ireslt = fread (ftmp, 4, nrec, mapin);
			    if (ireslt <= 0) {
//				fprintf (stderr, "error parsing superspace map\n");
				Error_Box ("Error seeking section in superspace map.");
				free (ftmp);
				fclose (mapin);
				return;
			    }
			}
		    }
		}
	    }
	}
	if (x5Val > 0.0f && x5Val <= x5max) {
	    nskip5 = (int) ((x5Val - x5min) / x5step + .5);
	    for (i = 0; i < nskip5; i++) {
		for (j = 0; j < step4; j++) {
		    for (k = 0; k < step3; k++) {
			ireslt = fread (ftmp, 4, nrec, mapin);
			if (ireslt <= 0) {
//			    fprintf (stderr, "error parsing superspace map\n");
			    Error_Box ("Error seeking section in superspace map.");
			    free (ftmp);
			    fclose (mapin);
			    return;
			}
		    }
		}
	    }
	}
	if (x4Val > 0.0f && x4Val <= x4max) {
	    nskip = (int) ((x4Val - x4min) / x4step + .5);
	    for (i = 0; i < nskip; i++) {
		for (j = 0; j < step3; j++) {
		    ireslt = fread (ftmp, 4, nrec, mapin);
		    if (ireslt <= 0) {
//			fprintf (stderr, "error parsing superspace map\n");
			Error_Box ("Error seeking section in superspace map.");
			free (ftmp);
			fclose (mapin);
			return;
		    }
		}
	    }
	}
	break;
    case 465:
	if (x5Val > 0.0f && x5Val <= x5max) {
	    nskip5 = (int) ((x5Val - x5min) / x5step + .5);
	    for (i = 0; i < nskip5; i++) {
		for (j = 0; j < step6; j++) {
		    for (k = 0; k < step4; k++) {
			for (l = 0; l < step3; l++) {
			    ireslt = fread (ftmp, 4, nrec, mapin);
			    if (ireslt <= 0) {
//				fprintf (stderr, "error parsing superspace map\n");
				Error_Box ("Error seeking section in superspace map.");
				free (ftmp);
				fclose (mapin);
				return;
			    }
			}
		    }
		}
	    }
	}
	if (x6Val > 0.0f && x6Val <= x6max) {
	    nskip6 = (int) ((x6Val - x6min) / x6step + .5);
	    for (i = 0; i < nskip6; i++) {
		for (j = 0; j < step4; j++) {
		    for (k = 0; k < step3; k++) {
			ireslt = fread (ftmp, 4, nrec, mapin);
			if (ireslt <= 0) {
//			    fprintf (stderr, "error parsing superspace map\n");
			    Error_Box ("Error seeking section in superspace map.");
			    free (ftmp);
			    fclose (mapin);
			    return;
			}
		    }
		}
	    }
	}
	if (x4Val > 0.0f && x4Val <= x4max) {
	    nskip = (int) ((x4Val - x4min) / x4step + .5);
	    for (i = 0; i < nskip; i++) {
		for (j = 0; j < step3; j++) {
		    ireslt = fread (ftmp, 4, nrec, mapin);
		    if (ireslt <= 0) {
//			fprintf (stderr, "error parsing superspace map\n");
			Error_Box ("Error seeking section in superspace map.");
			free (ftmp);
			fclose (mapin);
			return;
		    }
		}
	    }
	}
	break;
    case 546:
	if (x6Val > 0.0f && x6Val <= x6max) {
	    nskip6 = (int) ((x6Val - x6min) / x6step + .5);
	    for (i = 0; i < nskip6; i++) {
		for (j = 0; j < step4; j++) {
		    for (k = 0; k < step5; k++) {
			for (l = 0; l < step3; l++) {
			    ireslt = fread (ftmp, 4, nrec, mapin);
			    if (ireslt <= 0) {
//				fprintf (stderr, "error parsing superspace map\n");
				Error_Box ("Error seeking section in superspace map.");
				free (ftmp);
				fclose (mapin);
				return;
			    }
			}
		    }
		}
	    }
	}
	if (x4Val > 0.0f && x4Val <= x4max) {
	    nskip = (int) ((x4Val - x4min) / x4step + .5);
	    for (i = 0; i < nskip; i++) {
		for (j = 0; j < step5; j++) {
		    for (k = 0; k < step3; k++) {
			ireslt = fread (ftmp, 4, nrec, mapin);
			if (ireslt <= 0) {
//			    fprintf (stderr, "error parsing superspace map\n");
			    Error_Box ("Error seeking section in superspace map.");
			    free (ftmp);
			    fclose (mapin);
			    return;
			}
		    }
		}
	    }
	}
	if (x5Val > 0.0f && x5Val <= x5max) {
	    nskip5 = (int) ((x5Val - x5min) / x5step + .5);
	    for (i = 0; i < nskip5; i++) {
		for (j = 0; j < step3; j++) {
		    ireslt = fread (ftmp, 4, nrec, mapin);
		    if (ireslt <= 0) {
//			fprintf (stderr, "error parsing superspace map\n");
			Error_Box ("Error seeking section in superspace map.");
			free (ftmp);
			fclose (mapin);
			return;
		    }
		}
	    }
	}
	break;
    case 564:
	if (x4Val > 0.0f && x4Val <= x4max) {
	    nskip = (int) ((x4Val - x4min) / x4step + .5);
	    for (i = 0; i < nskip; i++) {
		for (j = 0; j < step6; j++) {
		    for (k = 0; k < step5; k++) {
			for (l = 0; l < step3; l++) {
			    ireslt = fread (ftmp, 4, nrec, mapin);
			    if (ireslt <= 0) {
//				fprintf (stderr, "error parsing superspace map\n");
				Error_Box ("Error seeking section in superspace map.");
				free (ftmp);
				fclose (mapin);
				return;
			    }
			}
		    }
		}
	    }
	}
	if (x6Val > 0.0f && x6Val <= x6max) {
	    nskip6 = (int) ((x6Val - x6min) / x6step + .5);
	    for (i = 0; i < nskip6; i++) {
		for (j = 0; j < step5; j++) {
		    for (k = 0; k < step3; k++) {
			ireslt = fread (ftmp, 4, nrec, mapin);
			if (ireslt <= 0) {
//			    fprintf (stderr, "error parsing superspace map\n");
			    Error_Box ("Error seeking section in superspace map.");
			    free (ftmp);
			    fclose (mapin);
			    return;
			}
		    }
		}
	    }
	}
	if (x5Val > 0.0f && x5Val <= x5max) {
	    nskip5 = (int) ((x5Val - x5min) / x5step + .5);
	    for (i = 0; i < nskip5; i++) {
		for (j = 0; j < step3; j++) {
		    ireslt = fread (ftmp, 4, nrec, mapin);
		    if (ireslt <= 0) {
//			fprintf (stderr, "error parsing superspace map\n");
			Error_Box ("Error seeking section in superspace map.");
			free (ftmp);
			fclose (mapin);
			return;
		    }
		}
	    }
	}
	break;
    case 654:
	if (x4Val > 0.0f && x4Val <= x4max) {
	    nskip = (int) ((x4Val - x4min) / x4step + .5);
	    for (i = 0; i < nskip; i++) {
		for (j = 0; j < step5; j++) {
		    for (k = 0; k < step6; k++) {
			for (l = 0; l < step3; l++) {
			    ireslt = fread (ftmp, 4, nrec, mapin);
			    if (ireslt <= 0) {
//				fprintf (stderr, "error parsing superspace map\n");
				Error_Box ("Error seeking section in superspace map.");
				free (ftmp);
				fclose (mapin);
				return;
			    }
			}
		    }
		}
	    }
	}
	if (x5Val > 0.0f && x5Val <= x5max) {
	    nskip5 = (int) ((x5Val - x5min) / x5step + .5);
	    for (i = 0; i < nskip5; i++) {
		for (j = 0; j < step6; j++) {
		    for (k = 0; k < step3; k++) {
			ireslt = fread (ftmp, 4, nrec, mapin);
			if (ireslt <= 0) {
//			    fprintf (stderr, "error parsing superspace map\n");
			    Error_Box ("Error seeking section in superspace map.");
			    free (ftmp);
			    fclose (mapin);
			    return;
			}
		    }
		}
	    }
	}
	if (x6Val > 0.0f && x6Val <= x6max) {
	    nskip6 = (int) ((x6Val - x6min) / x6step + .5);
	    for (i = 0; i < nskip6; i++) {
		for (j = 0; j < step3; j++) {
		    ireslt = fread (ftmp, 4, nrec, mapin);
		    if (ireslt <= 0) {
//			fprintf (stderr, "error parsing superspace map\n");
			Error_Box ("Error seeking section in superspace map.");
			free (ftmp);
			fclose (mapin);
			return;
		    }
		}
	    }
	}
	break;
    case 645:
	if (x5Val > 0.0f && x5Val <= x5max) {
	    nskip5 = (int) ((x5Val - x5min) / x5step + .5);
	    for (i = 0; i < nskip5; i++) {
		for (j = 0; j < step4; j++) {
		    for (k = 0; k < step6; k++) {
			for (l = 0; l < step3; l++) {
			    ireslt = fread (ftmp, 4, nrec, mapin);
			    if (ireslt <= 0) {
//				fprintf (stderr, "error parsing superspace map\n");
				Error_Box ("Error seeking section in superspace map.");
				free (ftmp);
				fclose (mapin);
				return;
			    }
			}
		    }
		}
	    }
	}
	if (x4Val > 0.0f && x4Val <= x4max) {
	    nskip = (int) ((x4Val - x4min) / x4step + .5);
	    for (i = 0; i < nskip; i++) {
		for (j = 0; j < step6; j++) {
		    for (k = 0; k < step3; k++) {
			ireslt = fread (ftmp, 4, nrec, mapin);
			if (ireslt <= 0) {
//			    fprintf (stderr, "error parsing superspace map\n");
			    Error_Box ("Error seeking section in superspace map.");
			    free (ftmp);
			    fclose (mapin);
			    return;
			}
		    }
		}
	    }
	}
	if (x6Val > 0.0f && x6Val <= x6max) {
	    nskip6 = (int) ((x6Val - x6min) / x6step + .5);
	    for (i = 0; i < nskip6; i++) {
		for (j = 0; j < step3; j++) {
		    ireslt = fread (ftmp, 4, nrec, mapin);
		    if (ireslt <= 0) {
//			fprintf (stderr, "error parsing superspace map\n");
			Error_Box ("Error seeking section in superspace map.");
			free (ftmp);
			fclose (mapin);
			return;
		    }
		}
	    }
	}
	break;

    default:
	Error_Box
	    ("ERROR -- Unable to select desired map section (unsupported axis numbering");
	free (ftmp);
	fclose (mapin);
	return;
	break;
    }

// read the Rho values
    rhomin = 100000.0f;
    rhomax = 0.000001f;
    if (fullcell == 1) {	//need to skip redundant information at end
	new_mapstep_a--;
	new_mapstep_b--;
	new_mapstep_c--;
	switch (axes) {
	case 123:		// X Y Z order  
	    for (i = cstart; i < cend; i++) {	//c  
		for (j = bstart; j < bend - 1; j++) {	//b
		    for (k = astart; k < aend - 1; k++) {
			ireslt = fread (ftmp, 4, 1, mapin);
			ijk = k * new_mapstep_c * new_mapstep_b + j * new_mapstep_c + i;
			FourierPt[ijk] = ftmp[0];
			if (FourierPt[ijk] < rhomin)
			    rhomin = FourierPt[ijk];
			if (FourierPt[ijk] > rhomax)
			    rhomax = FourierPt[ijk];
		    }
		    fread (ftmp, 4, 1, mapin);
		}
		for (k = astart; k < aend; k++)
		    fread (ftmp, 4, 1, mapin);
	    }
	    break;
	case 312:		// ZXY      
	    for (i = bstart; i < bend - 1; i++) {
		for (j = astart; j < aend - 1; j++) {
		    for (k = cstart; k < cend - 1; k++) {
			fread (ftmp, 4, 1, mapin);
			ijk = j * new_mapstep_c * new_mapstep_b + i * new_mapstep_c + k;
//                if (wrong_end) dens[k] = end_flip_real(dens[k]);
			FourierPt[ijk] = ftmp[0];
			if (FourierPt[ijk] < rhomin)
			    rhomin = FourierPt[ijk];
			if (FourierPt[ijk] > rhomax)
			    rhomax = FourierPt[ijk];
		    }
		    fread (ftmp, 4, 1, mapin);
		}
		for (k = cstart; k < cend; k++)
		    fread (ftmp, 4, 1, mapin);
	    }
	    break;
	case 321:		//ZYX  
	    for (i = astart; i < aend - 1; i++) {
		for (j = cstart; j < cend - 1; j++) {
		    for (k = bstart; k < bend - 1; k++) {
			fread (ftmp, 4, 1, mapin);
			ijk = i * new_mapstep_c * new_mapstep_b + j * new_mapstep_b + k;
//                if (wrong_end) dens[k] = end_flip_real(dens[k]);
			FourierPt[ijk] = ftmp[0];
			if (FourierPt[ijk] < rhomin)
			    rhomin = FourierPt[ijk];
			if (FourierPt[ijk] > rhomax)
			    rhomax = FourierPt[ijk];
		    }
		    fread (ftmp, 4, 1, mapin);
		}
		for (k = bstart; k < bend; k++)
		    fread (ftmp, 4, 1, mapin);
	    }
	    break;
	case 213:		//YXZ 
	    for (i = cstart; i < cend - 1; i++) {
		for (j = astart; j < aend - 1; j++) {
		    for (k = bstart; k < bend - 1; k++) {
			fread (ftmp, 4, 1, mapin);
			ijk = j * new_mapstep_c * new_mapstep_b + k * new_mapstep_c + i;
//                if (wrong_end) dens[k] = end_flip_real(dens[k]);
			FourierPt[ijk] = ftmp[0];
			if (FourierPt[ijk] < rhomin)
			    rhomin = FourierPt[ijk];
			if (FourierPt[ijk] > rhomax)
			    rhomax = FourierPt[ijk];
		    }
		    fread (ftmp, 4, 1, mapin);
		}
		for (k = bstart; k < bend; k++)
		    fread (ftmp, 4, 1, mapin);
	    }
	    break;
	case 231:		//YZX     
	    for (i = astart; i < aend - 1; i++) {
		for (j = cstart; j < cend - 1; j++) {
		    for (k = bstart; k < bend - 1; k++) {
			fread (ftmp, 4, 1, mapin);
			ijk = i * new_mapstep_c * new_mapstep_b + k * new_mapstep_c + j;
//                if (wrong_end) dens[k] = end_flip_real(dens[k]);
			FourierPt[ijk] = ftmp[0];
			if (FourierPt[ijk] < rhomin)
			    rhomin = FourierPt[ijk];
			if (FourierPt[ijk] > rhomax)
			    rhomax = FourierPt[ijk];
		    }
		    fread (ftmp, 4, 1, mapin);
		}
		for (k = bstart; k < bend; k++)
		    fread (ftmp, 4, 1, mapin);
	    }
	    break;
	case 132:		//XZY      
	    for (i = bstart; i < bend - 1; i++) {
		for (j = cstart; j < cend - 1; j++) {
		    for (k = astart; k < aend - 1; k++) {
			fread (ftmp, 4, 1, mapin);
			ijk = k * new_mapstep_c * new_mapstep_b + i * new_mapstep_c + j;
//                if (wrong_end) dens[k] = end_flip_real(dens[k]);
			FourierPt[ijk] = ftmp[0];
			if (FourierPt[ijk] < rhomin)
			    rhomin = FourierPt[ijk];
			if (FourierPt[ijk] > rhomax)
			    rhomax = FourierPt[ijk];
		    }
		    fread (ftmp, 4, 1, mapin);
		}
		for (k = astart; k < aend; k++)
		    fread (ftmp, 4, 1, mapin);
	    }
	    break;
	default:
	    Error_Box ("Primary axes must be some permutation of X,Y,Z .\n");
	    free (ftmp);
	    fclose (mapin);
	    return;
	    break;
	}
//      new_mapstep_a-=1;
//      new_mapstep_b-=1;
//      new_mapstep_c-=1;
    } else {
	switch (axes) {
	case 123:		// X Y Z order  
	    for (i = cstart; i < cend; i++) {	//c  
		for (j = bstart; j < bend; j++) {	//b
		    for (k = astart; k < aend; k++) {
			ireslt = fread (ftmp, 4, 1, mapin);
			ijk = k * new_mapstep_c * new_mapstep_b + j * new_mapstep_c + i;
//                if (wrong_end) dens[k] = end_flip_real(dens[k]);
			FourierPt[ijk] = ftmp[0];
			if (FourierPt[ijk] < rhomin)
			    rhomin = FourierPt[ijk];
			if (FourierPt[ijk] > rhomax)
			    rhomax = FourierPt[ijk];
		    }
		}
	    }
	    break;
	case 312:		// ZXY      
	    for (i = bstart; i < bend; i++) {
		for (j = astart; j < aend; j++) {
		    for (k = cstart; k < cend; k++) {
			fread (ftmp, 4, 1, mapin);
			ijk = j * new_mapstep_c * new_mapstep_b + i * new_mapstep_c + k;
//                if (wrong_end) dens[k] = end_flip_real(dens[k]);
			FourierPt[ijk] = ftmp[0];
			if (FourierPt[ijk] < rhomin)
			    rhomin = FourierPt[ijk];
			if (FourierPt[ijk] > rhomax)
			    rhomax = FourierPt[ijk];
		    }
		}
	    }
	    break;
	case 321:		//ZYX  
	    for (i = astart; i < aend; i++) {
		for (j = cstart; j < cend; j++) {
		    for (k = bstart; k < bend; k++) {
			fread (ftmp, 4, 1, mapin);
			ijk = i * new_mapstep_c * new_mapstep_b + j * new_mapstep_b + k;
//                if (wrong_end) dens[k] = end_flip_real(dens[k]);
			FourierPt[ijk] = ftmp[0];
			if (FourierPt[ijk] < rhomin)
			    rhomin = FourierPt[ijk];
			if (FourierPt[ijk] > rhomax)
			    rhomax = FourierPt[ijk];
		    }
		}
	    }
	    break;
	case 213:		//YXZ 
	    for (i = cstart; i < cend; i++) {
		for (j = astart; j < aend; j++) {
		    for (k = bstart; k < bend; k++) {
			fread (ftmp, 4, 1, mapin);
			ijk = j * new_mapstep_c * new_mapstep_b + k * new_mapstep_c + i;
//                if (wrong_end) dens[k] = end_flip_real(dens[k]);
			FourierPt[ijk] = ftmp[0];
			if (FourierPt[ijk] < rhomin)
			    rhomin = FourierPt[ijk];
			if (FourierPt[ijk] > rhomax)
			    rhomax = FourierPt[ijk];
		    }
		}
	    }
	    break;
	case 231:		//YZX     
	    for (i = astart; i < aend; i++) {
		for (j = cstart; j < cend; j++) {
		    for (k = bstart; k < bend; k++) {
			fread (ftmp, 4, 1, mapin);
			ijk = i * new_mapstep_c * new_mapstep_b + k * new_mapstep_c + j;
//                if (wrong_end) dens[k] = end_flip_real(dens[k]);
			FourierPt[ijk] = ftmp[0];
			if (FourierPt[ijk] < rhomin)
			    rhomin = FourierPt[ijk];
			if (FourierPt[ijk] > rhomax)
			    rhomax = FourierPt[ijk];
		    }
		}
	    }
	    break;
	case 132:		//XZY      
	    for (i = bstart; i < bend; i++) {
		for (j = cstart; j < cend; j++) {
		    for (k = astart; k < aend; k++) {
			fread (ftmp, 4, 1, mapin);
			ijk = k * new_mapstep_c * new_mapstep_b + i * new_mapstep_c + j;
//                if (wrong_end) dens[k] = end_flip_real(dens[k]);
			FourierPt[ijk] = ftmp[0];
			if (FourierPt[ijk] < rhomin)
			    rhomin = FourierPt[ijk];
			if (FourierPt[ijk] > rhomax)
			    rhomax = FourierPt[ijk];
		    }
		}
	    }
	    break;
	default:
	    Error_Box ("Primary axes must be some permutation of X,Y,Z .\n");
	    free (ftmp);
	    fclose (mapin);
	    return;
	    break;
	}
    }

    mapstep_a = new_mapstep_a;
    mapstep_b = new_mapstep_b;
    mapstep_c = new_mapstep_c;
    fclose (mapin);
    free (ftmp);

    ReadFourMap = 1;
    sprintf (Map_Info.title, "%s map with axes %s %s %s %s %s %s\n", maptypes[maptype],
	     theaxes[axis[0]], theaxes[axis[1]], theaxes[axis[2]],
	     theaxes[axis[3]], theaxes[axis[4]], theaxes[axis[5]]);
    Map_Info.rhomn = rhomin;
    Map_Info.rhomx = rhomax;
    Map_Info.map_int[0] = mapstep_a - 1;
    Map_Info.map_int[1] = mapstep_b - 1;
    Map_Info.map_int[2] = mapstep_c - 1;
    Map_Info.lat_con[0] = drvui->lat_con[0];
    Map_Info.lat_con[1] = drvui->lat_con[1];
    Map_Info.lat_con[2] = drvui->lat_con[2];
    Map_Info.lat_con[3] = drvui->lat_con[3];
    Map_Info.lat_con[4] = drvui->lat_con[4];
    Map_Info.lat_con[5] = drvui->lat_con[5];
    Map_Info.info_valid = 1;

    Map_Info.x4lim[0] = x4min;
    Map_Info.x4lim[1] = x4max;
    Map_Info.x5lim[0] = x5min;
    Map_Info.x5lim[1] = x5max;
    Map_Info.x6lim[0] = x6min;
    Map_Info.x6lim[1] = x6max;
#if 1
    xMin = Map_Info.xlim[0] = xmin;
    xMax = Map_Info.xlim[1] = xmax;
    yMin = Map_Info.ylim[0] = ymin;
    yMax = Map_Info.ylim[1] = ymax;
    zMin = Map_Info.zlim[0] = zmin;
    zMax = Map_Info.zlim[1] = zmax;
#endif

    if (!Quick) {
	fprintf (drvui->flout, "Reading %s map file %s\n"
		 "axis sequence = %s %s %s %s %s %s\n  map steps: a=%d b=%d c=%d\n  "
		 "Range of rho values: %f to %f\n", maptypes[maptype], infile,
		 theaxes[axis[0]], theaxes[axis[1]], theaxes[axis[2]],
		 theaxes[axis[3]], theaxes[axis[4]], theaxes[axis[5]],
		 mapstep_a, mapstep_b, mapstep_c, rhomin, rhomax);
	if (axis[3] != 0)
	    fprintf (drvui->flout, "X4 section %2f to %2f stepsize %2f\n",
		     x4min, x4max, x4step);
	if (axis[4] != 0)
	    fprintf (drvui->flout, "X5 section %2f to %2f stepsize %2f\n",
		     x5min, x5max, x5step);
	if (axis[5] != 0)
	    fprintf (drvui->flout, "X6 section %2f to %2f stepsize %2f\n",
		     x6min, x6max, x6step);
    }

}

void
read_m80 (char *infile, int Quick)
{
// read JANA .m80 Fo/Fc file and compute density

    FILE *mapin;

    float *fo, *fc, *f2, *ftmp1, *ftmp2, *ftmp3;

    short *ih, *ik, *il, *itmp1, *itmp2, *itmp3;

    float *Cklx;

    float *Dklx;

    float *Elxy;

    float *Flxy;

    double x, y, z;

    char string[200];

    int nval = 200;

    int nr;

    int i, j, k, l;

    int ijk = 0;

    float rhomin, rhomax;

    float factor;

    short ih0, ik0, il0;

    short ih1, ik1, il1;

    short im, in, io, is;

    float fo0, fc0, fc1, f20, f21;

    int kmin, kmax, nk;

    int lmin, lmax, nl;

    float phase;

    int skip;

    int progress;

    float foj;

    int warnings;

    int dimension;

    if (FourierPt != NULL) {
	return;
    }

    if ((mapin = fopen (infile, "r")) == NULL) {
	sprintf (string, "Cannot open Fo/Fc (JANA m80) file %s", infile);
	Error_Box (string);
	return;
    }

    fo = (float *) zalloc (200 * sizeof (float));
    if (fo == NULL) {
	sprintf (string, "Unable to obtain storage for Fo data");
	Error_Box (string);
	fclose (mapin);
	return;
    }
    fc = (float *) zalloc (200 * sizeof (float));
    if (fc == NULL) {
	sprintf (string, "Unable to obtain storage for Fc data");
	Error_Box (string);
	free (fo);
	fclose (mapin);
	return;
    }
    f2 = (float *) zalloc (200 * sizeof (float));
    if (f2 == NULL) {
	sprintf (string, "Unable to obtain storage for B or phase data");
	Error_Box (string);
	free (fo);
	free (fc);
	fclose (mapin);
	return;
    }
    ih = (short *) zalloc (200 * sizeof (short));
    if (ih == NULL) {
	sprintf (string, "Unable to obtain storage for h data");
	Error_Box (string);
	free (fo);
	free (fc);
	free (f2);
	fclose (mapin);
	return;
    }
    ik = (short *) zalloc (200 * sizeof (short));
    if (ik == NULL) {
	sprintf (string, "Unable to obtain storage for k data");
	Error_Box (string);
	free (fo);
	free (fc);
	free (f2);
	free (ih);
	fclose (mapin);
	return;
    }
    il = (short *) zalloc (200 * sizeof (short));
    if (il == NULL) {
	sprintf (string, "Unable to obtain storage for l data");
	Error_Box (string);
	free (fo);
	free (fc);
	free (f2);
	free (ih);
	free (ik);
	fclose (mapin);
	return;
    }

    if (fgets (string, 199, mapin)) {
	i = sscanf (string, "%hd %hd %hd %hd %hd %hd %hd %f", &ih0, &ik0, &il0, &im, &in,
		    &io, &is, &fo0);
	if (i > 4)
	    Map_Info.info_valid = 1;
	dimension = i - 2;
    } else {
	sprintf (string, "Error reading M80 file");
	Error_Box (string);
	free (fo);
	free (fc);
	free (f2);
	free (ih);
	free (ik);
	free (il);
	fclose (mapin);
	return;
    }

    if (dimension < 3) {
	sprintf (string, "Error reading M80 file");
	Error_Box (string);
	free (fo);
	free (fc);
	free (f2);
	free (ih);
	free (ik);
	free (il);
	fclose (mapin);
	return;
    }

    nr = 0;
    warnings = 0;

    do {
	skip = 0;
	if (dimension > 3) {
	    switch (dimension) {
	    case 4:
		im = 0;
		i = sscanf (string, "%hd %hd %hd %hd %*d %f %*f %*f %f %f", &ih0, &ik0,
			    &il0, &im, &fo0, &fc0, &f20);
		if (im != 0)
		    skip = 1;
		break;
	    case 5:
		im = in = 0;
		i = sscanf (string, "%hd %hd %hd %hd %hd %*d %f %*f %*f %f %f", &ih0,
			    &ik0, &il0, &im, &in, &fo0, &fc0, &f20);
		if (im != 0 || in != 0)
		    skip = 1;
		break;
	    case 6:
		im = in = io = 0;
		i = sscanf (string, "%hd %hd %hd %hd %hd %hd %*d %f %*f %*f %f %f",
			    &ih0, &ik0, &il0, &im, &in, &io, &fo0, &fc0, &f20);
		if (im != 0 || in != 0 || io != 0)
		    skip = 1;
		break;
	    }
	} else {
	    i = sscanf (string, "%hd %hd %hd %*d %f %*f %*f %f %f", &ih0, &ik0, &il0,
			&fo0, &fc0, &f20);
	}
	if (i < dimension + 3)
	    break;

	if (skip == 1) {	// satellite reflection of modulated structure
	    fgets (string, 199, mapin);
	    continue;
	}
	// change all data to 'phase' type
	fc1 = (float) sqrt (fc0 * fc0 + f20 * f20);
	if (fc1 == 0.)		// happens with data generated by superflip
	    f20 = 0.;
	else
	    f20 = (float) atan2 (f20 / fc1, fc0 / fc1);
	fc0 = fc1;

// apply spacegroup symmetry (including xyz - catch duplicates in input list)
	for (k = 0; k < drvui->ng; ++k) {
	    skip = 0;

	    ih1 =
		drvui->ss[k][0][0] * ih0 + drvui->ss[k][1][0] * ik0 +
		drvui->ss[k][2][0] * il0;
	    ik1 =
		drvui->ss[k][0][1] * ih0 + drvui->ss[k][1][1] * ik0 +
		drvui->ss[k][2][1] * il0;
	    il1 =
		drvui->ss[k][0][2] * ih0 + drvui->ss[k][1][2] * ik0 +
		drvui->ss[k][2][2] * il0;
	    phase =
		2.0f * (float) PI *(drvui->ts[k][0] * ih0 + drvui->ts[k][1] * ik0 +
				    drvui->ts[k][2] * il0);
	    f21 = f20 + phase;
	    f21 = (float) atan2 (sin (f21), cos (f21));	// get f21 in range -PI to PI

	    for (j = 0; j < nr; j++) {
		if (ih[j] == ih1 && ik[j] == ik1 && il[j] == il1) {
		    skip = 1;
		    foj = fo[j];
		    break;
		}
	    }
	    if (skip == 0) {
		ih[nr] = ih1;
		ik[nr] = ik1;
		il[nr] = il1;
		fo[nr] = fo0;
		fc[nr] = fc0;
		f2[nr] = f21;
		nr++;
	    } else {
		if (fabs (foj - fo0) > 1.e-3) {
		    if (!Quick && warnings < 100)
			fprintf (drvui->flout,
				 "Warning - different Fo for %d %d %d and %d %d %d\n",
				 ih0, ik0, il0, ih1, ik1, il1);
		    warnings++;
		}
	    }

	    ih1 *= -1;		// add Friedel opposite
	    ik1 *= -1;
	    il1 *= -1;
	    f21 *= -1.0;
	    skip = 0;

	    for (j = 0; j < nr; j++) {
		if (ih[j] == ih1 && ik[j] == ik1 && il[j] == il1) {
		    skip = 1;
		    foj = fo[j];
		    break;
		}
	    }
	    if (skip == 0) {
		ih[nr] = ih1;
		ik[nr] = ik1;
		il[nr] = il1;
		fo[nr] = fo0;
		fc[nr] = fc0;
		f2[nr] = f21;
		nr++;
	    } else {
		if (fabs (foj - fo0) > 1.e-3) {
		    if (!Quick && warnings < 100)
			fprintf (drvui->flout,
				 "Warning - different Fo for %d %d %d and %d %d %d\n",
				 ih0, ik0, il0, ih1, ik1, il1);
		    warnings++;
		}
	    }

	    if (nr > nval - 10) {
		nval += 100;
		ftmp1 = (float *) realloc (fo, nval * sizeof (float));
		ftmp2 = (float *) realloc (fc, nval * sizeof (float));
		ftmp3 = (float *) realloc (f2, nval * sizeof (float));
		itmp1 = (short *) realloc (ih, nval * sizeof (short));
		itmp2 = (short *) realloc (ik, nval * sizeof (short));
		itmp3 = (short *) realloc (il, nval * sizeof (short));
		if (!ftmp1 || !ftmp2 || !ftmp3 || !itmp1 || !itmp2 || !itmp3) {
		    sprintf (string,
			     "Unable to expand storage for h, k, l, fo, fc or phase data");
		    free (fo);
		    free (fc);
		    free (f2);
		    free (ih);
		    free (ik);
		    free (il);
		    if (ftmp1)
			free (ftmp1);
		    if (ftmp2)
			free (ftmp2);
		    if (ftmp3)
			free (ftmp3);
		    if (itmp1)
			free (itmp1);
		    if (itmp2)
			free (itmp2);
		    if (itmp3)
			free (itmp3);
		    Error_Box (string);
		    fclose (mapin);
		    return;
		} else {
		    fo = ftmp1;
		    fc = ftmp2;
		    f2 = ftmp3;
		    ih = itmp1;
		    ik = itmp2;
		    il = itmp3;
		}
	    }			// if we need to realloc
	}			// loop over all symmetry operators
	fgets (string, 199, mapin);
    } while (i > 0 && !feof (mapin));

    fclose (mapin);

    if (dimension > 3) {
	sprintf (string,
		 "%d-dimensional data encountered - only main reflections will be used",
		 dimension);
	Error_Box (string);
    }
    if (warnings > 0) {
	sprintf (string,
		 "Symmetry problem - %d mismatches in symmetry-equivalent Fo values\n",
		 warnings);
	Error_Box (string);
    }
    factor = 1.0f / (drvui->lat_con[0] * drvui->lat_con[1] * drvui->lat_con[2]
		     * (float) sqrt (1.0 - cos (drvui->lat_con[3] * PI / 180.0) *
				     cos (drvui->lat_con[3] * PI / 180.0) *
				     cos (drvui->lat_con[4] * PI / 180.0) *
				     cos (drvui->lat_con[4] * PI / 180.0) *
				     cos (drvui->lat_con[5] * PI / 180.0) *
				     cos (drvui->lat_con[5] * PI / 180.0)));
    kmin = kmax = ik[0];
    lmin = lmax = il[0];

    for (l = 0; l < nr; l++) {
	if (ik[l] < kmin)
	    kmin = ik[l];	// find range of k and l
	if (ik[l] > kmax)
	    kmax = ik[l];
	if (il[l] < lmin)
	    lmin = il[l];
	if (il[l] > lmax)
	    lmax = il[l];
	switch (Map_Info.map_type) {	// get map coefficients/V
	default:
	case 0:		// Fo map
	    break;
	case 1:		// Fc map
	    fo[l] = fc[l];
	    break;
	case 2:		// Fo - Fc map
	    fo[l] = fo[l] - fc[l];
	    break;
	case 3:		// 2Fo - Fc map
	    fo[l] = 2.0f * fo[l] - fc[l];
	    break;
	case 4:		// Fo2 (Patterson)
	    fo[l] = fo[l] * fo[l];
	    f2[l] = 0.0f;
	    break;
	}
	fo[l] *= factor;	// update coefficients times 1/V
    }

    nk = kmax - kmin + 1;
    nl = lmax - lmin + 1;
    Cklx = (float *) zalloc (nk * nl * sizeof (float));
    Dklx = (float *) zalloc (nk * nl * sizeof (float));
    Elxy = (float *) zalloc (nl * sizeof (float));
    Flxy = (float *) zalloc (nl * sizeof (float));

    if (!Cklx || !Dklx || !Elxy || !Flxy) {
	Error_Box ("Unable to allocate Beevers-Lipson arrays.");
	if (Cklx) free (Cklx);
	if (Dklx) free (Dklx);
	if (Elxy) free (Elxy);
	if (Flxy) free (Flxy);
	free (fo);
	free (fc);
	free (f2);
	free (ih);
	free (ik);
	free (il);
	return;
    }
    mapstep_a = (int) (4.0f * drvui->lat_con[0] + 0.5f);
    mapstep_b = (int) (4.0f * drvui->lat_con[1] + 0.5f);
    mapstep_c = (int) (4.0f * drvui->lat_con[2] + 0.5f);

    FourierPt = (float *) zalloc (mapstep_a * mapstep_b * mapstep_c * sizeof (float));
    if (FourierPt == NULL) {
	Error_Box ("ERROR -- Unable to allocate space for map.");
	return;
    }

    rhomin = 1.0e15f, rhomax = -1.0e15f;

// start the Beevers-Lipson expansion

    Progress_Window (-1, "Computing Density", (float) mapstep_a * mapstep_b);

    progress = 0;
    for (i = 0; i < mapstep_a; i++) {
	x = (float) i / (float) mapstep_a;
	memset (Cklx, 0, nk * nl * sizeof (float));
	memset (Dklx, 0, nk * nl * sizeof (float));
	for (l = 0; l < nr; l++) {
	    ijk = (il[l] - lmin) * nk + ik[l] - kmin;
	    Cklx[ijk] += fo[l] * (float) cos (2.0 * PI * ih[l] * x - f2[l]);
	    Dklx[ijk] -= fo[l] * (float) sin (2.0 * PI * ih[l] * x - f2[l]);
	}
	for (j = 0; j < mapstep_b; j++) {
	    Progress_Window (0, NULL, (float) ++progress);
	    y = (float) j / (float) mapstep_b;
	    memset (Elxy, 0, nl * sizeof (float));
	    memset (Flxy, 0, nl * sizeof (float));
	    for (l = 0; l < nl; l++) {
		for (k = 0; k < nk; k++) {
		    ijk = l * nk + k;
		    Elxy[l] += Cklx[ijk] * (float) cos (2.0 * PI * (k + kmin) * y)
			+ Dklx[ijk] * (float) sin (2.0 * PI * (k + kmin) * y);
		    Flxy[l] += Dklx[ijk] * (float) cos (2.0 * PI * (k + kmin) * y)
			- Cklx[ijk] * (float) sin (2.0 * PI * (k + kmin) * y);
		}
	    }
	    for (k = 0; k < mapstep_c; k++) {
		z = (float) k / (float) mapstep_c;
		ijk = i * mapstep_b * mapstep_c + j * mapstep_c + k;

		for (l = 0; l < nl; l++) {
		    FourierPt[ijk] += Elxy[l] * (float) cos (2.0 * PI * (l + lmin) * z)
			+ Flxy[l] * (float) sin (2.0 * PI * (l + lmin) * z);
		}
		if (FourierPt[ijk] < rhomin)
		    rhomin = FourierPt[ijk];
		if (FourierPt[ijk] > rhomax)
		    rhomax = FourierPt[ijk];
	    }
	}
    }
    if (!Quick)
	fprintf (drvui->fcns, "Reading Fourier map file %s\n  "
		 "map steps: a=%d b=%d c=%d\n  "
		 "Range of rho values: %f to %f, gridsize %d\n", infile,
		 mapstep_a, mapstep_b, mapstep_c, rhomin, rhomax, ijk);


    Progress_Window (-2, NULL, 0.0f);
    ReadFourMap = 1;
    Map_Info.rhomn = rhomin;
    Map_Info.rhomx = rhomax;
    Map_Info.map_int[0] = mapstep_a;
    Map_Info.map_int[1] = mapstep_b;
    Map_Info.map_int[2] = mapstep_c;
    Map_Info.lat_con[0] = drvui->lat_con[0];
    Map_Info.lat_con[1] = drvui->lat_con[1];
    Map_Info.lat_con[2] = drvui->lat_con[2];
    Map_Info.lat_con[3] = drvui->lat_con[3];
    Map_Info.lat_con[4] = drvui->lat_con[4];
    Map_Info.lat_con[5] = drvui->lat_con[5];
    Map_Info.xlim[0] = 0.0f;
    Map_Info.xlim[1] = 1.0f;
    Map_Info.ylim[0] = 0.0f;
    Map_Info.ylim[1] = 1.0f;
    Map_Info.zlim[0] = 0.0f;
    Map_Info.zlim[1] = 1.0f;

    free (Cklx);
    free (Dklx);
    free (Elxy);
    free (Flxy);
    free (fo);
    free (fc);
    free (f2);
    free (ih);
    free (ik);
    free (il);
}

/* ************************************************************** */
/* ************************************************************** */

void
read_flp (char *infile, int Quick)
{
// read FullProf binary map file

#include "read_flp.h"		// define all the struct's

    char titulo[81];

    char version[6];

    struct crystal_cell_type celda;

    FILE *mapin;

    float rhomin = 1.0e+15f, rhomax = -1.0e+15f;

    int ireslt = 0;

    char string[200];

    int icount;

    struct space_group spg;

    char temp[16384];

    int i, j, k, ijk = 0;

    int natoms;

    struct FFT_param map_param;

    struct Atom_Type atm[200];

    float dens[200];

    int new_mapstep_a, astart, aend;

    int new_mapstep_b, bstart, bend;

    int new_mapstep_c, cstart, cend;

    int wrong_end = 0;

    memset (version, 0, 6);
    memset (titulo, 0, 81);
    if ((mapin = fopen (infile, "r+b")) == NULL) {
	sprintf (string, "Cannot open Fourier map (FulProf) file %s", infile);
	Error_Box (string);
	return;
    }
    ireslt = fread (&icount, sizeof (int), 1, mapin);	// read record length
    if (!ireslt) {
	char string[200];

	sprintf (string, "Error reading Fourier map file %s", infile);
	Error_Box (string);
	fclose (mapin);
	return;
    }
    if (icount != 80) {
	icount = end_flip (icount);
	wrong_end = 1;
    }
    if (icount != 80) {
	Error_Box ("Map does not have the correct format\nJob aborted.");
	fclose (mapin);
	return;
    }
    fread (titulo, sizeof (char), icount, mapin);	// read title
    fread (&icount, sizeof (int), 1, mapin);	// skip trailer
    trim_string (titulo, 80);
    strcpy (Map_Info.title, titulo);

    fread (&icount, sizeof (int), 1, mapin);	// read header for next record
    if (wrong_end)
	icount = end_flip (icount);
    fread (version, sizeof (char), icount, mapin);
    fread (&icount, sizeof (int), 1, mapin);	// skip trailer

    fread (&icount, sizeof (int), 1, mapin);	// read header for next record
    if (wrong_end)
	icount = end_flip (icount);
    fread (&celda, sizeof (char), icount, mapin);
    if (wrong_end) {
	for (i = 0; i < 3; i++) {
	    celda.ang[i] = end_flip_real (celda.ang[i]);
	    celda.cell[i] = end_flip_real (celda.cell[i]);
	}
    }
    fread (&icount, sizeof (int), 1, mapin);	// skip trailer
    for (i = 0; i < 3; i++) {
	Map_Info.lat_con[i] = celda.cell[i];
	Map_Info.lat_con[i + 3] = celda.ang[i];
    }

    fread (&icount, sizeof (int), 1, mapin);	// read header for next record
    if (wrong_end)
	icount = end_flip (icount);
    fread (&spg, sizeof (char), icount, mapin);
    fread (&icount, sizeof (int), 1, mapin);	// skip trailer

    fread (&icount, sizeof (int), 1, mapin);	// read header for next record
    if (wrong_end)
	icount = end_flip (icount);
    fread (&natoms, sizeof (char), icount, mapin);
    fread (&icount, sizeof (int), 1, mapin);	// skip trailer

    fread (&icount, sizeof (int), 1, mapin);	// read header for next record
    if (wrong_end)
	icount = end_flip (icount);
    fread (&atm, sizeof (char), icount, mapin);
    fread (&icount, sizeof (int), 1, mapin);	// skip trailer

    for (i = 0; i < 2; i++) {	// skip next 2 records
	fread (&icount, sizeof (int), 1, mapin);	// read header for next record
	if (wrong_end)
	    icount = end_flip (icount);
	fread (&temp, sizeof (char), icount, mapin);
	fread (&icount, sizeof (int), 1, mapin);	// skip trailer
    }

    fread (&icount, sizeof (int), 1, mapin);	// read header for next record
    if (wrong_end)
	icount = end_flip (icount);
    fread (&map_param, sizeof (char), icount, mapin);
    fread (&icount, sizeof (int), 1, mapin);	// skip trailer
    if (wrong_end) {
	for (i = 0; i < 3; i++) {
	    map_param.ngrid[i] = end_flip (map_param.ngrid[i]);
	}
    }
    Map_Info.map_int[0] = mapstep_a = map_param.ngrid[0];
    Map_Info.map_int[1] = mapstep_b = map_param.ngrid[1];
    Map_Info.map_int[2] = mapstep_c = map_param.ngrid[2];
    for (i = 0; i < 2; i++) {
	if (wrong_end) {
	    map_param.xlim[i] = end_flip_real (map_param.xlim[i]);
	    map_param.ylim[i] = end_flip_real (map_param.ylim[i]);
	    map_param.zlim[i] = end_flip_real (map_param.zlim[i]);
	}
	Map_Info.xlim[i] = map_param.xlim[i];
	Map_Info.ylim[i] = map_param.ylim[i];
	Map_Info.zlim[i] = map_param.zlim[i];
    }
    Map_Info.info_valid = 1;

    fread (&icount, sizeof (int), 1, mapin);	// read header for next record
    if (wrong_end)
	icount = end_flip (icount);
    fread (&temp, sizeof (char), icount, mapin);
    fread (&icount, sizeof (int), 1, mapin);	// skip trailer

    if ((fabs (Map_Info.xlim[0]) > 0.001) || (fabs (Map_Info.xlim[1] - 1.0f) > 0.001)) {
// not full cell along x
	new_mapstep_a =
	    (int) (1.0f / (Map_Info.xlim[1] - Map_Info.xlim[0]) + 0.5f) * mapstep_a;
	astart = (int) (Map_Info.xlim[0] * new_mapstep_a);
	aend = (int) (Map_Info.xlim[1] * new_mapstep_a);
    } else {
	new_mapstep_a = mapstep_a;
	astart = 0;
	aend = mapstep_a;
    }
    if ((fabs (Map_Info.ylim[0]) > 0.001) || (fabs (Map_Info.ylim[1] - 1.0f) > 0.001)) {
// not full cell along y
	new_mapstep_b =
	    (int) (1.0f / (Map_Info.ylim[1] - Map_Info.ylim[0]) + 0.5f) * mapstep_b;
	bstart = (int) (Map_Info.ylim[0] * new_mapstep_b);
	bend = (int) (Map_Info.ylim[1] * new_mapstep_b);
    } else {
	new_mapstep_b = mapstep_b;
	bstart = 0;
	bend = mapstep_b;
    }
    if ((fabs (Map_Info.zlim[0]) > 0.001) || (fabs (Map_Info.zlim[1] - 1.0f) > 0.001)) {
// not full cell along z
	new_mapstep_c =
	    (int) (1.0f / (Map_Info.zlim[1] - Map_Info.zlim[0]) + 0.5f) * mapstep_c;
	cstart = (int) (Map_Info.zlim[0] * new_mapstep_c);
	cend = (int) (Map_Info.zlim[1] * new_mapstep_c);
    } else {
	new_mapstep_c = mapstep_c;
	cstart = 0;
	cend = mapstep_c;
    }
    if (FourierPt)
	free (FourierPt);
    FourierPt =
	(float *) zalloc (new_mapstep_a * new_mapstep_b * new_mapstep_c * sizeof (float));
    if (FourierPt == NULL) {
	ireslt = 0;
	Error_Box ("ERROR -- Unable to allocate space for map.");
	fclose (mapin);
	return;
    }
// read the Rho values

    for (i = astart; i < aend; i++) {
	for (j = bstart; j < bend; j++) {
	    fread (&icount, sizeof (int), 1, mapin);	// read header for next record
	    if (wrong_end)
		icount = end_flip (icount);
	    fread (&dens, sizeof (char), icount, mapin);
	    fread (&icount, sizeof (int), 1, mapin);	// skip trailer
	    for (k = cstart; k < cend; k++) {
		ijk = i * new_mapstep_c * new_mapstep_b + j * new_mapstep_c + k;
		if (wrong_end)
		    dens[k] = end_flip_real (dens[k]);
		FourierPt[ijk] = dens[k];
		if (FourierPt[ijk] < rhomin)
		    rhomin = FourierPt[ijk];
		if (FourierPt[ijk] > rhomax)
		    rhomax = FourierPt[ijk];
	    }
	}
    }
    fclose (mapin);
    if (!Quick)
	fprintf (drvui->fcns, "Reading Fourier map file %s\n  "
		 "cell = %f %f %f %f %f %f\n  map steps: a=%d b=%d c=%d\n  "
		 "Range of rho values: %f to %f, gridsize %d\n", infile,
		 celda.cell[0], celda.cell[1], celda.cell[2],
		 celda.ang[0], celda.ang[1], celda.ang[2],
		 mapstep_a, mapstep_b, mapstep_c, rhomin, rhomax, ijk);
    ReadFourMap = 1;
    Map_Info.rhomn = rhomin;
    Map_Info.rhomx = rhomax;
    mapstep_a = new_mapstep_a;
    mapstep_b = new_mapstep_b;
    mapstep_c = new_mapstep_c;
    Map_Info.map_int[0] = mapstep_a;
    Map_Info.map_int[1] = mapstep_b;
    Map_Info.map_int[2] = mapstep_c;
}

/* ************************************************************** */
/* ************************************************************** */

void
read_grd (char *infile, int Quick)
{
    FILE *mapin;

    char line[81];

    char *reslt;

    float rhomin = 1.0e15f, rhomax = -1.0e15f;

    int ireslt = 0;

    char string[200];

    int new_mapstep_a = 0, astart, aend;

    int new_mapstep_b = 0, bstart, bend;

    int new_mapstep_c = 0, cstart, cend;

    if ((mapin = fopen (infile, "r")) == NULL) {
	sprintf (string, "Cannot open Fourier map (.grd) file %s", infile);
	Error_Box (string);
	return;
    }
    reslt = fgets (line, sizeof (line), mapin);
    trim_string (line, 81);
    strcpy (Map_Info.title, line);
    Map_Info.info_valid = 1;
    if (reslt != NULL) {
	ireslt = fscanf (mapin, "%f %f %f %f %f %f",
			 &map_a, &map_b, &map_c, &map_alpha, &map_beta, &map_gamma);
	Map_Info.lat_con[0] = map_a;
	Map_Info.lat_con[1] = map_b;
	Map_Info.lat_con[2] = map_c;
	Map_Info.lat_con[3] = map_alpha;
	Map_Info.lat_con[4] = map_beta;
	Map_Info.lat_con[5] = map_gamma;
    }
    if (ireslt > 0)
	ireslt = fscanf (mapin, "%d %d %d", &mapstep_a, &mapstep_b, &mapstep_c);
    if (ireslt > 0) {
	int i, j, k, ijk;

	if ((fabs (Map_Info.xlim[0]) > 0.001) || (fabs (Map_Info.xlim[1] - 1.0f) > 0.001)) {
// not full cell along x
	    new_mapstep_a =
		(int) (1.0f / (Map_Info.xlim[1] - Map_Info.xlim[0]) + 0.5f) * mapstep_a;
	    astart = (int) (Map_Info.xlim[0] * new_mapstep_a);
	    aend = (int) (Map_Info.xlim[1] * new_mapstep_a);
	} else {
	    new_mapstep_a = mapstep_a;
	    astart = 0;
	    aend = mapstep_a;
	}
	if ((fabs (Map_Info.ylim[0]) > 0.001) || (fabs (Map_Info.ylim[1] - 1.0f) > 0.001)) {
// not full cell along y
	    new_mapstep_b =
		(int) (1.0f / (Map_Info.ylim[1] - Map_Info.ylim[0]) + 0.5f) * mapstep_b;
	    bstart = (int) (Map_Info.ylim[0] * new_mapstep_b);
	    bend = (int) (Map_Info.ylim[1] * new_mapstep_b);
	} else {
	    new_mapstep_b = mapstep_b;
	    bstart = 0;
	    bend = mapstep_b;
	}
	if ((fabs (Map_Info.zlim[0]) > 0.001) || (fabs (Map_Info.zlim[1] - 1.0f) > 0.001)) {
// not full cell along z
	    new_mapstep_c =
		(int) (1.0f / (Map_Info.zlim[1] - Map_Info.zlim[0]) + 0.5f) * mapstep_c;
	    cstart = (int) (Map_Info.zlim[0] * new_mapstep_c);
	    cend = (int) (Map_Info.zlim[1] * new_mapstep_c);
	} else {
	    new_mapstep_c = mapstep_c;
	    cstart = 0;
	    cend = mapstep_c;
	}

	if (FourierPt)
	    free (FourierPt);

	FourierPt =
	    (float *) zalloc (new_mapstep_a * new_mapstep_b * new_mapstep_c *
			      sizeof (float));
	if (FourierPt == NULL) {
	    ireslt = 0;
	    Error_Box ("ERROR -- Unable to allocate space for map.\n");
	    fclose (mapin);
	    return;
	}

/* read the Rho values */
	for (i = astart; i < aend; i++) {
	    for (j = bstart; j < bend; j++) {
		for (k = cstart; k < cend; k++) {
		    ijk = i * new_mapstep_b * new_mapstep_c + j * new_mapstep_c + k;
		    ireslt = fscanf (mapin, "%f", &FourierPt[ijk]);
		    if (FourierPt[ijk] < rhomin)
			rhomin = FourierPt[ijk];
		    if (FourierPt[ijk] > rhomax)
			rhomax = FourierPt[ijk];
		}
	    }
	}
    }
    Map_Info.rhomn = rhomin;
    Map_Info.rhomx = rhomax;
    Map_Info.map_int[0] = mapstep_a = new_mapstep_a;
    Map_Info.map_int[1] = mapstep_b = new_mapstep_b;
    Map_Info.map_int[2] = mapstep_c = new_mapstep_c;

    if (ireslt > 0) {
	if (!Quick) {
	    fprintf (drvui->flout, "Reading Fourier map file %s\n\tTitle= %s\n  "
		     "cell = %f %f %f %f %f %f\n  map steps: a=%d b=%d c=%d\n  "
		     "Range of rho values: %f to %f\n", infile, line,
		     map_a, map_b, map_c, map_alpha, map_beta, map_gamma,
		     mapstep_a, mapstep_b, mapstep_c, rhomin, rhomax);
	}
	ReadFourMap = 1;
    } else {
	char string[200];

	sprintf (string, "Error reading Fourier map file %s", infile);
	Error_Box (string);
    }
    fclose (mapin);
}

/* ************************************************************** */
/* ************************************************************** */

void
read_stf (char *infile, int Quick)
{
    FILE *mapin;

    char line[81];

    char *reslt;

    float rhomin = 1.0e15f, rhomax = -1.0e15f;

    int ireslt = 0;

    int ijk = 0;

    int i, j, k;

    int ii, jj, kk;

    int fullcell = 1;

    float bx1, bx2, by1, by2, bz1, bz2;

    char string[200];

    int new_mapstep_a, astart, aend;

    int new_mapstep_b, bstart, bend;

    int new_mapstep_c, cstart, cend;

    if ((mapin = fopen (infile, "r")) == NULL) {
	sprintf (string, "Cannot open Fourier map (.stf) file %s\n", infile);
	Error_Box (string);
	return;
    }
    reslt = fgets (line, sizeof (line), mapin);
    reslt = fgets (line, sizeof (line), mapin);
    reslt = fgets (line, sizeof (line), mapin);
    if (reslt) {
	ireslt = sscanf (line, "%*s %d %d %d", &mapstep_a, &mapstep_b, &mapstep_c);
    }
    if (ireslt > 0) {
	reslt = fgets (line, sizeof (line), mapin);
	sscanf (line, "%*s %f %f %f %f %f %f", &bx1, &bx2, &by1, &by2, &bz1, &bz2);
	Map_Info.xlim[0] = bx1;
	Map_Info.xlim[1] = bx2;
	Map_Info.ylim[0] = by1;
	Map_Info.ylim[1] = by2;
	Map_Info.zlim[0] = bz1;
	Map_Info.zlim[1] = bz2;
	Map_Info.lat_con[0] = map_a = drvui->lat_con[0];
	Map_Info.lat_con[1] = map_b = drvui->lat_con[1];
	Map_Info.lat_con[2] = map_c = drvui->lat_con[2];
	Map_Info.lat_con[3] = map_alpha = drvui->lat_con[3];
	Map_Info.lat_con[4] = map_beta = drvui->lat_con[4];
	Map_Info.lat_con[5] = map_gamma = drvui->lat_con[5];

	if ((fabs (Map_Info.xlim[0]) > 0.001) || (fabs (Map_Info.xlim[1] - 1.0f) > 0.001)) {
// not full cell along x
	    new_mapstep_a =
		(int) ((1.0f / (Map_Info.xlim[1] - Map_Info.xlim[0])) *
		       (float) mapstep_a + 0.5f);
	    astart = (int) (Map_Info.xlim[0] * new_mapstep_a + 0.5f);
	    aend = (int) (Map_Info.xlim[1] * new_mapstep_a + 0.5f);
	    fullcell = 0;
	    if (astart < 0 && aend - astart != mapstep_a)
		aend++;
	    if (astart > 0 && aend - astart != mapstep_a)
		aend--;
	    if (aend - astart != mapstep_a)
		Error_Box ("Error recalculating map a limits for full cell");
	} else {
	    new_mapstep_a = mapstep_a;
	    astart = 0;
	    aend = mapstep_a;
	}
	if ((fabs (Map_Info.ylim[0]) > 0.001) || (fabs (Map_Info.ylim[1] - 1.0f) > 0.001)) {
// not full cell along y
	    new_mapstep_b =
		(int) ((1.0f / (Map_Info.ylim[1] - Map_Info.ylim[0])) *
		       (float) mapstep_b + 0.5f);
	    bstart = (int) (Map_Info.ylim[0] * new_mapstep_b + 0.5f);
	    bend = (int) (Map_Info.ylim[1] * new_mapstep_b + 0.5f);
	    fullcell = 0;
	    if (bstart < 0 && bend - bstart != mapstep_b)
		bend++;
	    if (bstart > 0 && bend - bstart != mapstep_b)
		bend--;
	    if (bend - bstart != mapstep_b)
		Error_Box ("Error recalculating map b limits for full cell");
	} else {
	    new_mapstep_b = mapstep_b;
	    bstart = 0;
	    bend = mapstep_b;
	}
	if ((fabs (Map_Info.zlim[0]) > 0.001) || (fabs (Map_Info.zlim[1] - 1.0f) > 0.001)) {
// not full cell along z
	    new_mapstep_c =
		(int) ((1.0f / (Map_Info.zlim[1] - Map_Info.zlim[0])) *
		       (float) mapstep_c + 0.5f);
	    cstart = (int) (Map_Info.zlim[0] * new_mapstep_c + 0.5f);
	    cend = (int) (Map_Info.zlim[1] * new_mapstep_c + 0.5f);
	    fullcell = 0;
	    if (cstart < 0 && cend - cstart != mapstep_c)
		cend++;
	    if (cstart > 0 && cend - cstart != mapstep_c)
		cend--;
	    if (cend - cstart != mapstep_c)
		Error_Box ("Error recalculating map c limits for full cell");
	} else {
	    new_mapstep_c = mapstep_c;
	    cstart = 0;
	    cend = mapstep_c;
	}

	if (FourierPt)
	    free (FourierPt);

	FourierPt =
	    (float *) zalloc (new_mapstep_a * new_mapstep_b * new_mapstep_c *
			      sizeof (float));
	if (FourierPt == NULL) {
	    ireslt = 0;
	    Error_Box ("ERROR -- Unable to allocate space for map.\n");
	    fclose (mapin);
	    return;
	}
	reslt = fgets (line, sizeof (line), mapin);
	reslt = fgets (line, sizeof (line), mapin);
	reslt = fgets (line, sizeof (line), mapin);

	for (i = astart; i < aend; i++) {
	    for (j = bstart; j < bend - 1; j++) {
		for (k = cstart; k < cend - 1; k++) {
		    ii = i + (i < 0 ? new_mapstep_a : 0);
		    jj = j + (j < 0 ? new_mapstep_b : 0);
		    kk = k + (k < 0 ? new_mapstep_c : 0);
		    ijk =
			kk * (new_mapstep_a - 1) * (new_mapstep_b - 1) +
			jj * (new_mapstep_a - 1) + ii;
		    ireslt = fscanf (mapin, "%E", &FourierPt[ijk]);
		    if (FourierPt[ijk] < rhomin)
			rhomin = FourierPt[ijk];
			if (FourierPt[ijk] > rhomax)
			    rhomax = FourierPt[ijk];
		}
		fscanf (mapin, "%*f");
	    }
	    for (k = 0; k < mapstep_c; k++)
		fscanf (mapin, "%*f");
	}
	new_mapstep_a--;
	new_mapstep_b--;
	new_mapstep_c--;

	Map_Info.map_int[0] = mapstep_a = new_mapstep_a;
	Map_Info.map_int[1] = mapstep_b = new_mapstep_b;
	Map_Info.map_int[2] = mapstep_c = new_mapstep_c;
	Map_Info.rhomn = rhomin;
	Map_Info.rhomx = rhomax;
	Map_Info.info_valid = 1;

    }
    /* if mapstep line present */
    if (ijk > 0) {
	if (!Quick)
	    fprintf (drvui->flout, "Reading Fourier map file %s\n  "
		     "cell = %f %f %f %f %f %f\n  map steps: a=%d b=%d c=%d\n  "
		     "Range of rho values: %f to %f, gridsize %d\n", infile,
		     map_a, map_b, map_c, map_alpha, map_beta, map_gamma,
		     mapstep_a, mapstep_b, mapstep_c, rhomin, rhomax, ijk);
	ReadFourMap = 1;
    } else {
	sprintf (string, "Error reading Fourier map file %s", infile);
	Error_Box (string);
    }
    fclose (mapin);
}

/* ************************************************************** */
/* ************************************************************** */

void
read_vasp (char *infile, int Quick)
{
    FILE *mapin;

    char line[200];

    char *reslt;

    float rhomin = 1.0e15f, rhomax = -1.0e15f;

    int ireslt = 0;

    int ijk = 0;

    int i, j, k, l, m;

    float tmp_rho[10];

    char string[200];

    float scale;

    float cmat[3][3];

    int new_mapstep_a = 0, astart, aend;

    int new_mapstep_b = 0, bstart, bend;

    int new_mapstep_c = 0, cstart, cend;

    map_a = drvui->lat_con[0];
    map_b = drvui->lat_con[1];
    map_c = drvui->lat_con[2];
    map_alpha = drvui->lat_con[3];
    map_beta = drvui->lat_con[4];
    map_gamma = drvui->lat_con[5];
    if ((mapin = fopen (infile, "r")) == NULL) {
	sprintf (string, "Cannot open VASP file %s\n", infile);
	Error_Box (string);
	return;
    }
    fgets (line, sizeof (line), mapin);
    trim_string (line, 80);
    strcpy (Map_Info.title, line);
    fgets (line, sizeof (line), mapin);
    sscanf (line, "%f", &scale);
    for (i = 0; i < 3; i++) {	// read the cell matrix
	if (!(reslt = fgets (line, sizeof (line), mapin))) {
	    sprintf (string, "Error reading VASP file.");
	    Error_Box (string);
	    fclose (mapin);
	    return;
	}
	sscanf (line, "%f %f %f", &cmat[i][0], &cmat[i][1], &cmat[i][2]);
    }
    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    cmat[i][j] *= scale;
    fgets (line, sizeof (line), mapin);
    Map_Info.lat_con[0] = (float) sqrt (cmat[0][0]);	// this section needs fixing
    Map_Info.lat_con[1] = (float) sqrt (cmat[1][1]);
    Map_Info.lat_con[2] = (float) sqrt (cmat[2][2]);
    Map_Info.lat_con[3] = 90.0f;
    Map_Info.lat_con[4] = 90.0f;
    Map_Info.lat_con[5] = 90.0f;
    j = 0;
    i = strlen (line) / 4;
    for (k = 0; k < i; k++) {
	sscanf (line + 4 * k, "%4d", &l);
	j += l;
    }
    for (i = 0; i < j + 2; i++) {	// set to skip k-point stuff
	if (!(reslt = fgets (line, sizeof (line), mapin))) {
	    sprintf (string, "Error reading VASP file.");
	    Error_Box (string);
	    fclose (mapin);
	    return;
	}
    }
    fgets (line, sizeof (line), mapin);
    ireslt = sscanf (line, "%d %d %d", &mapstep_a, &mapstep_b, &mapstep_c);
    if (ireslt > 0) {
	if ((fabs (Map_Info.xlim[0]) > 0.001) || (fabs (Map_Info.xlim[1] - 1.0f) > 0.001)) {
// not full cell along x
	    new_mapstep_a =
		(int) (1.0f / (Map_Info.xlim[1] - Map_Info.xlim[0]) + 0.5f) * mapstep_a;
	    astart = (int) (Map_Info.xlim[0] * new_mapstep_a);
	    aend = (int) (Map_Info.xlim[1] * new_mapstep_a);
	} else {
	    new_mapstep_a = mapstep_a;
	    astart = 0;
	    aend = mapstep_a;
	}
	if ((fabs (Map_Info.ylim[0]) > 0.001) || (fabs (Map_Info.ylim[1] - 1.0f) > 0.001)) {
// not full cell along y
	    new_mapstep_b =
		(int) (1.0f / (Map_Info.ylim[1] - Map_Info.ylim[0]) + 0.5f) * mapstep_b;
	    bstart = (int) (Map_Info.ylim[0] * new_mapstep_b);
	    bend = (int) (Map_Info.ylim[1] * new_mapstep_b);
	} else {
	    new_mapstep_b = mapstep_b;
	    bstart = 0;
	    bend = mapstep_b;
	}
	if ((fabs (Map_Info.zlim[0]) > 0.001) || (fabs (Map_Info.zlim[1] - 1.0f) > 0.001)) {
// not full cell along z
	    new_mapstep_c =
		(int) (1.0f / (Map_Info.zlim[1] - Map_Info.zlim[0]) + 0.5f) * mapstep_c;
	    cstart = (int) (Map_Info.zlim[0] * new_mapstep_c);
	    cend = (int) (Map_Info.zlim[1] * new_mapstep_c);
	} else {
	    new_mapstep_c = mapstep_c;
	    cstart = 0;
	    cend = mapstep_c;
	}

	if (FourierPt)
	    free (FourierPt);

	FourierPt =
	    (float *) zalloc (new_mapstep_a * new_mapstep_b * new_mapstep_c *
			      sizeof (float));
	if (FourierPt == NULL) {
	    ireslt = 0;
	    Error_Box ("ERROR -- Unable to allocate space for map.\n");
	    fclose (mapin);
	    return;
	}
// read the Rho values

	m = 11;
	for (i = cstart; i < cend; i++) {
	    for (j = bstart; j < bend; j++) {
		for (k = astart; k < aend; k++) {
		    if (m > 9) {
			reslt = fgets (line, sizeof (line), mapin);
			sscanf (line, "%f %f %f %f %f %f %f %f %f %f", &tmp_rho[0],
				&tmp_rho[1], &tmp_rho[2], &tmp_rho[3], &tmp_rho[4],
				&tmp_rho[5], &tmp_rho[6], &tmp_rho[7], &tmp_rho[8],
				&tmp_rho[9]);
			m = 0;
		    }
		    ijk = k * new_mapstep_b * new_mapstep_c + j * new_mapstep_c + i;
		    FourierPt[ijk] = tmp_rho[m++];
		    if (FourierPt[ijk] < rhomin)
			rhomin = FourierPt[ijk];
		    if (FourierPt[ijk] > rhomax)
			rhomax = FourierPt[ijk];
		}
	    }
	}
    }

    if (reslt) {
	ReadFourMap = 1;
    } else {
	char string[200];

	sprintf (string, "Error reading Fourier map file %s", infile);
	Error_Box (string);
    }
    Map_Info.rhomn = rhomin;
    Map_Info.rhomx = rhomax;
    Map_Info.map_int[0] = mapstep_a = new_mapstep_a;
    Map_Info.map_int[1] = mapstep_b = new_mapstep_b;
    Map_Info.map_int[2] = mapstep_c = new_mapstep_c;
    if (!Quick) {
	fprintf (drvui->flout, "Reading Fourier map file %s\n"
		 "cell = %f %f %f %f %f %f\n  map steps: a=%d b=%d c=%d\n  "
		 "Range of rho values: %f to %f\n", infile,
		 map_a, map_b, map_c, map_alpha, map_beta, map_gamma,
		 mapstep_a, mapstep_b, mapstep_c, rhomin, rhomax);
    }
    Map_Info.info_valid = 1;
    fclose (mapin);

}

/* ************************************************************** */
/* ************************************************************** */

void
read_w2k (char *infile, int Quick)
{
    FILE *mapin;

    char line[81];

    char *reslt;

    float rhomin = 1.0e15f, rhomax = -1.0e15f;

    int ireslt = 0;

    int dimen = 0;

    int ijk = 0;

    int i, j, k;

    int ii, jj, kk;

    char string[200];

    int new_mapstep_a = 0, astart, aend;

    int new_mapstep_b = 0, bstart, bend;

    int new_mapstep_c = 0, cstart, cend;

    Map_Info.lat_con[0] = map_a = drvui->lat_con[0];
    Map_Info.lat_con[1] = map_b = drvui->lat_con[1];
    Map_Info.lat_con[2] = map_c = drvui->lat_con[2];
    Map_Info.lat_con[3] = map_alpha = drvui->lat_con[3];
    Map_Info.lat_con[4] = map_beta = drvui->lat_con[4];
    Map_Info.lat_con[5] = map_gamma = drvui->lat_con[5];

    if ((mapin = fopen (infile, "r")) == NULL) {
	sprintf (string, "Cannot open Fourier map (.w2k) file %s\n", infile);
	Error_Box (string);
	return;
    }

    do {
	reslt = fgets (line, sizeof (line), mapin);
	if (!strncmp (line, "(@3 [3]", 7))
	    dimen = 1;
    } while (!feof (mapin) && dimen == 0);

    if (dimen == 0) {		// not in lapw5_3d format, probably from wien2venus.py
	rewind (mapin);
	reslt = fgets (line, sizeof (line), mapin);
	if (!strncmp (line, "cell", 4))
	    dimen = 2;
    }

    if (dimen == 0) {
	sprintf (string, "Cannot read Fourier map (.w2k) file %s\n", infile);
	Error_Box (string);
	fclose (mapin);
	return;
    }

    if (dimen == 2) {		// file generated by wien2venus.py
	reslt = fgets (line, sizeof (line), mapin);
	reslt = fgets (line, sizeof (line), mapin);
	reslt = fgets (line, sizeof (line), mapin);
	ireslt = sscanf (line, "%d %d %d", &mapstep_a, &mapstep_b, &mapstep_c);
    } else {
	ireslt = sscanf (line, "%*s %*s %d %d %d", &mapstep_a, &mapstep_b, &mapstep_c);
    }

    if (ireslt > 0) {
	if ((fabs (Map_Info.xlim[0]) > 0.001) || (fabs (Map_Info.xlim[1] - 1.0f) > 0.001)) {
// not full cell along x
	    new_mapstep_a =
		(int) (1.0f / (Map_Info.xlim[1] - Map_Info.xlim[0]) + 0.5f) * mapstep_a;
	    astart = (int) (Map_Info.xlim[0] * new_mapstep_a);
	    aend = (int) (Map_Info.xlim[1] * new_mapstep_a);
	} else {
	    new_mapstep_a = mapstep_a;
	    astart = 0;
	    aend = mapstep_a;
	}
	if ((fabs (Map_Info.ylim[0]) > 0.001) || (fabs (Map_Info.ylim[1] - 1.0f) > 0.001)) {
// not full cell along y
	    new_mapstep_b =
		(int) (1.0f / (Map_Info.ylim[1] - Map_Info.ylim[0]) + 0.5f) * mapstep_b;
	    bstart = (int) (Map_Info.ylim[0] * new_mapstep_b);
	    bend = (int) (Map_Info.ylim[1] * new_mapstep_b);
	} else {
	    new_mapstep_b = mapstep_b;
	    bstart = 0;
	    bend = mapstep_b;
	}
	if ((fabs (Map_Info.zlim[0]) > 0.001) || (fabs (Map_Info.zlim[1] - 1.0f) > 0.001)) {
// not full cell along z
	    new_mapstep_c =
		(int) (1.0f / (Map_Info.zlim[1] - Map_Info.zlim[0]) + 0.5f) * mapstep_c;
	    cstart = (int) (Map_Info.zlim[0] * new_mapstep_c);
	    cend = (int) (Map_Info.zlim[1] * new_mapstep_c);
	} else {
	    new_mapstep_c = mapstep_c;
	    cstart = 0;
	    cend = mapstep_c;
	}

	if (FourierPt)
	    free (FourierPt);

	FourierPt =
	    (float *) zalloc (new_mapstep_a * new_mapstep_b * new_mapstep_c *
			      sizeof (float));
	if (FourierPt == NULL) {
	    ireslt = 0;
	    Error_Box ("ERROR -- Unable to allocate space for map.\n");
	    fclose (mapin);
	    return;
	}

	if (dimen == 1) {
	    reslt = fgets (line, sizeof (line), mapin);

	    for (k = cstart; k < cend; k++) {
		for (j = bstart; j < bend - 1; j++) {
		    for (i = astart; i < aend - 1; i++) {
			ijk =
			    j * (new_mapstep_a - 1) * (new_mapstep_b - 1) +
			    k * (new_mapstep_a - 1) + i;
			ireslt = fscanf (mapin, "%f", &FourierPt[ijk]);
			if (FourierPt[ijk] < rhomin)
			    rhomin = FourierPt[ijk];
			if (FourierPt[ijk] > rhomax)
			    rhomax = FourierPt[ijk];
		    }
		    fscanf (mapin, "%*f");
		}
		for (i = 0; i < mapstep_a; i++)
		    fscanf (mapin, "%*f");
	    }
	} else {
	    ijk = -1;
	    for (k = astart; k < aend; k++) {
		for (j = bstart; j < bend - 1; j++) {
		    for (i = cstart; i < cend - 1; i++) {
			ijk++;

			if (drvui->sys == 5) {	/* convert from rhombohedral setting */
			    float xi = (float) k / (float) aend;

			    float yi = (float) j / (float) bend;

			    float zi = (float) i / (float) cend;

			    float xn = xi * 2.f / 3.f - yi * 1.f / 3.f - zi * 1.f / 3.f;

			    if (xn < 0.)
				xn += 1.;
			    if (xn > 1.)
				xn -= 1.;
			    float yn = xi * 1.f / 3.f + yi * 1.f / 3.f - zi * 2.f / 3.f;

			    if (yn < 0.)
				yn += 1.;
			    if (yn > 1.)
				yn -= 1.;
			    float zn = xi * 1.f / 3.f + yi * 1.f / 3.f + zi * 1.f / 3.f;

			    if (zn < 0.)
				zn += 1.;
			    if (zn > 1.)
				zn -= 1.;
			    kk = (int) (xn * aend + .5);
			    jj = (int) (yn * bend + .5);
			    ii = (int) (zn * cend + .5);

			    ijk =
				kk * (new_mapstep_c - 1) * (new_mapstep_b - 1) +
				jj * (new_mapstep_c - 1) + ii;
			}

			ireslt = fscanf (mapin, "%f", &FourierPt[ijk]);
			if (FourierPt[ijk] < rhomin)
			    rhomin = FourierPt[ijk];
			if (FourierPt[ijk] > rhomax)
			    rhomax = FourierPt[ijk];
		    }
		    fscanf (mapin, "%*f");
		}
		for (i = 0; i < mapstep_a; i++)
		    fscanf (mapin, "%*f");
	    }
	}

    }

    if (ijk > 0) {
	if (!Quick)
	    fprintf (drvui->flout, "Reading Fourier map file %s\n  "
		     "cell = %f %f %f %f %f %f\n  map steps: a=%d b=%d c=%d\n  "
		     "Range of rho values: %f to %f, gridsize %d\n", infile,
		     map_a, map_b, map_c, map_alpha, map_beta, map_gamma,
		     mapstep_a, mapstep_b, mapstep_c, rhomin, rhomax, ijk);
	ReadFourMap = 1;
    } else {
	sprintf (string, "Error reading Fourier map file %s", infile);
	Error_Box (string);
    }
    fclose (mapin);
    Map_Info.rhomn = rhomin;
    Map_Info.rhomx = rhomax;
    Map_Info.map_int[0] = mapstep_a = --new_mapstep_a;
    Map_Info.map_int[1] = mapstep_b = --new_mapstep_b;
    Map_Info.map_int[2] = mapstep_c = --new_mapstep_c;
    Map_Info.info_valid = 1;
}

/* ************************************************************** */
/* ************************************************************** */

void
read_exc (char *infile, int Quick)
{
    FILE *mapin;

    char line[81];

    char *reslt;

    float rhomin, rhomax;

    int ireslt = 0;

    int ijk = 0;

    int i, j, k;

    char string[200];

    int new_mapstep_a = 0, astart, aend;

    int new_mapstep_b = 0, bstart, bend;

    int new_mapstep_c = 0, cstart, cend;

    Map_Info.lat_con[0] = map_a = drvui->lat_con[0];
    Map_Info.lat_con[1] = map_b = drvui->lat_con[1];
    Map_Info.lat_con[2] = map_c = drvui->lat_con[2];
    Map_Info.lat_con[3] = map_alpha = drvui->lat_con[3];
    Map_Info.lat_con[4] = map_beta = drvui->lat_con[4];
    Map_Info.lat_con[5] = map_gamma = drvui->lat_con[5];

    rhomin = 1.0e15f;
    rhomax = -1.0e15f;

    if ((mapin = fopen (infile, "r")) == NULL) {
	sprintf (string, "Cannot open Fourier map (.OUT) file %s\n", infile);
	Error_Box (string);
	return;
    }

    ireslt = fscanf (mapin, "%d %d %d", &mapstep_a, &mapstep_b, &mapstep_c);

    if (ireslt > 0) {
	if ((fabs (Map_Info.xlim[0]) > 0.001) || (fabs (Map_Info.xlim[1] - 1.0f) > 0.001)) {
// not full cell along x
	    new_mapstep_a =
		(int) (1.0f / (Map_Info.xlim[1] - Map_Info.xlim[0]) + 0.5f) * mapstep_a;
	    astart = (int) (Map_Info.xlim[0] * new_mapstep_a);
	    aend = (int) (Map_Info.xlim[1] * new_mapstep_a);
	} else {
	    new_mapstep_a = mapstep_a;
	    astart = 0;
	    aend = mapstep_a;
	}
	if ((fabs (Map_Info.ylim[0]) > 0.001) || (fabs (Map_Info.ylim[1] - 1.0f) > 0.001)) {
// not full cell along y
	    new_mapstep_b =
		(int) (1.0f / (Map_Info.ylim[1] - Map_Info.ylim[0]) + 0.5f) * mapstep_b;
	    bstart = (int) (Map_Info.ylim[0] * new_mapstep_b);
	    bend = (int) (Map_Info.ylim[1] * new_mapstep_b);
	} else {
	    new_mapstep_b = mapstep_b;
	    bstart = 0;
	    bend = mapstep_b;
	}
	if ((fabs (Map_Info.zlim[0]) > 0.001) || (fabs (Map_Info.zlim[1] - 1.0f) > 0.001)) {
// not full cell along z
	    new_mapstep_c =
		(int) (1.0f / (Map_Info.zlim[1] - Map_Info.zlim[0]) + 0.5f) * mapstep_c;
	    cstart = (int) (Map_Info.zlim[0] * new_mapstep_c);
	    cend = (int) (Map_Info.zlim[1] * new_mapstep_c);
	} else {
	    new_mapstep_c = mapstep_c;
	    cstart = 0;
	    cend = mapstep_c;
	}

	if (FourierPt)
	    free (FourierPt);

	FourierPt =
	    (float *) zalloc (new_mapstep_a * new_mapstep_b * new_mapstep_c *
			      sizeof (float));
	if (FourierPt == NULL) {
	    ireslt = 0;
	    Error_Box ("ERROR -- Unable to allocate space for map.\n");
	    fclose (mapin);
	    return;
	}

	reslt = fgets (line, sizeof (line), mapin);	// trailing comments from mapstep line
	fgetc (mapin);

	for (i = cstart; i < cend; i++) {
	    for (j = bstart; j < bend; j++) {
		for (k = astart; k < aend; k++) {
		    ijk = k * new_mapstep_b * new_mapstep_c + j * new_mapstep_b + i;
		    reslt = fgets (line, sizeof (line), mapin);
		    ireslt = sscanf (line, "%*f %*f %*f %f", &FourierPt[ijk]);
		    if (FourierPt[ijk] < rhomin)
			rhomin = FourierPt[ijk];
		    if (FourierPt[ijk] > rhomax)
			rhomax = FourierPt[ijk];
		}
	    }
	}
    }

    if (ijk > 0) {
	if (!Quick)
	    fprintf (drvui->flout, "Reading Fourier map file %s\n  "
		     "cell = %f %f %f %f %f %f\n  map steps: a=%d b=%d c=%d\n  "
		     "Range of rho values: %f to %f, gridsize %d\n", infile,
		     map_a, map_b, map_c, map_alpha, map_beta, map_gamma,
		     mapstep_a, mapstep_b, mapstep_c, rhomin, rhomax, ijk);
	ReadFourMap = 1;
    } else {
	sprintf (string, "Error reading Fourier map file %s", infile);
	Error_Box (string);
    }
    fclose (mapin);
    Map_Info.rhomn = rhomin;
    Map_Info.rhomx = rhomax;
    Map_Info.map_int[0] = mapstep_a = new_mapstep_a;
    Map_Info.map_int[1] = mapstep_b = new_mapstep_b;
    Map_Info.map_int[2] = mapstep_c = new_mapstep_c;
    Map_Info.info_valid = 1;
}

void
read_aim (char *infile, int Quick)
{
    FILE *surf;

    char line[81];

    char *reslt;

    int ireslt = 0;

    int i, j, k, m;

    int nphi, ntet;

    int nvrt;

    float phi, theta, r, e;

    double sinphi, sinthe, cosphi, costhe;

    char string[200];

    if ((surf = fopen (infile, "r")) == NULL) {
	sprintf (string, "Cannot open AIM surface (.surf) file %s\n", infile);
	Error_Box (string);
	return;
    }

    m = drvui->nsurf;

    reslt = fgets (line, 80, surf);
    reslt = fgets (line, 80, surf);
    ireslt = sscanf (line, "%d", &ntet);
    drvui->ntet[m] = ntet;
    reslt = fgets (line, 80, surf);
    ireslt = sscanf (line, "%d", &nphi);
    drvui->nphi[m] = nphi;
    if (ireslt <= 0) {
	sprintf (string, "Error reading AIM surface (.surf) file %s\n", infile);
	Error_Box (string);
	fclose (surf);
	return;
    }

    nvrt = ntet * nphi;
    drvui->surfx[m] = (float *) malloc (nvrt * sizeof (float));
    drvui->surfy[m] = (float *) malloc (nvrt * sizeof (float));
    drvui->surfz[m] = (float *) malloc (nvrt * sizeof (float));

    k = 0;
    for (i = 0; i < ntet; ++i) {
	for (j = 0; j < nphi; ++j) {
	    reslt = fgets (line, 80, surf);
	    ireslt = sscanf (line, "%f %f %gD%f %*f", &phi, &theta, &r, &e);
//if (ireslt<3) break;
	    r *= powf (10., e);
	    r *= 0.529f;
	    sinphi = sin (phi);
	    cosphi = cos (phi);
	    sinthe = sin (theta);
	    costhe = cos (theta);

	    drvui->surfx[m][k] = (float) (sinphi * costhe * r);
	    drvui->surfy[m][k] = (float) (sinphi * sinthe * r);
	    drvui->surfz[m][k] = (float) (cosphi * r);
//fprintf(stderr,"%f %f %f\n",drvui->surfx[m][k],drvui->surfy[m][k],drvui->surfz[m][k]);
	    k++;
	}
    }

    if (!Quick)
	fprintf (drvui->flout, "Read %d vertices from AIM surface file %s\n",
		 nvrt, infile);

    fclose (surf);
    return;

}

/* ************************************************************** */
/* ************************************************************** */

void
read_xsf (char *infile, int Quick)
{
    FILE *mapin;

    char line[81];

    char *reslt = NULL;

    float rhomin, rhomax;

    int ireslt = 0;

    int ijk = 0;

    int i, j, k;

    char string[200];

    int new_mapstep_a = 0, astart, aend;

    int new_mapstep_b = 0, bstart, bend;

    int new_mapstep_c = 0, cstart, cend;

    Map_Info.lat_con[0] = map_a = drvui->lat_con[0];
    Map_Info.lat_con[1] = map_b = drvui->lat_con[1];
    Map_Info.lat_con[2] = map_c = drvui->lat_con[2];
    Map_Info.lat_con[3] = map_alpha = drvui->lat_con[3];
    Map_Info.lat_con[4] = map_beta = drvui->lat_con[4];
    Map_Info.lat_con[5] = map_gamma = drvui->lat_con[5];

    rhomin = 1.0e15f;
    rhomax = -1.0e15f;

    if ((mapin = fopen (infile, "r")) == NULL) {
	sprintf (string, "Cannot open Fourier map (.XSF) file %s\n", infile);
	Error_Box (string);
	return;
    }

    while ( !strstr(line,"BEGIN_BLOCK_DATAGRID_3D")) {
	reslt = fgets (line, sizeof (line), mapin);
	if (feof(mapin)) {
	    sprintf (string, "Error reading Fourier map (.XSF) file %s\n", infile);
	    Error_Box (string);
	    fclose (mapin);
	    return;
	}
    }
    reslt = fgets (line, sizeof (line), mapin);	// skip comment
    reslt = fgets (line, sizeof (line), mapin);

    if ( !strstr(line,"BEGIN_DATAGRID_3D_") ) {
	sprintf (string, "Error reading Fourier map (.XSF) file %s\n", infile);
	Error_Box (string);
	fclose (mapin);
	return;
    }

    reslt = fgets (line, sizeof (line), mapin);	// 
    ireslt=sscanf (line, "%d %d %d", &mapstep_a, &mapstep_b, &mapstep_c);

    if (ireslt > 0) {
	if ((fabs (Map_Info.xlim[0]) > 0.001) || (fabs (Map_Info.xlim[1] - 1.0f) > 0.001)) {
// not full cell along x
	    new_mapstep_a =
		(int) (1.0f / (Map_Info.xlim[1] - Map_Info.xlim[0]) + 0.5f) * mapstep_a;
	    astart = (int) (Map_Info.xlim[0] * new_mapstep_a);
	    aend = (int) (Map_Info.xlim[1] * new_mapstep_a);
	} else {
	    new_mapstep_a = mapstep_a;
	    astart = 0;
	    aend = mapstep_a;
	}
	if ((fabs (Map_Info.ylim[0]) > 0.001) || (fabs (Map_Info.ylim[1] - 1.0f) > 0.001)) {
// not full cell along y
	    new_mapstep_b =
		(int) (1.0f / (Map_Info.ylim[1] - Map_Info.ylim[0]) + 0.5f) * mapstep_b;
	    bstart = (int) (Map_Info.ylim[0] * new_mapstep_b);
	    bend = (int) (Map_Info.ylim[1] * new_mapstep_b);
	} else {
	    new_mapstep_b = mapstep_b;
	    bstart = 0;
	    bend = mapstep_b;
	}
	if ((fabs (Map_Info.zlim[0]) > 0.001) || (fabs (Map_Info.zlim[1] - 1.0f) > 0.001)) {
// not full cell along z
	    new_mapstep_c =
		(int) (1.0f / (Map_Info.zlim[1] - Map_Info.zlim[0]) + 0.5f) * mapstep_c;
	    cstart = (int) (Map_Info.zlim[0] * new_mapstep_c);
	    cend = (int) (Map_Info.zlim[1] * new_mapstep_c);
	} else {
	    new_mapstep_c = mapstep_c;
	    cstart = 0;
	    cend = mapstep_c;
	}

	if (FourierPt)
	    free (FourierPt);

	FourierPt =
	    (float *) zalloc (new_mapstep_a * new_mapstep_b * new_mapstep_c *
			      sizeof (float));
	if (FourierPt == NULL) {
	    ireslt = 0;
	    Error_Box ("ERROR -- Unable to allocate space for map.\n");
	    fclose (mapin);
	    return;
	}

	reslt = fgets (line, sizeof (line), mapin);	// ignore cell vectors
	reslt = fgets (line, sizeof (line), mapin);	// for now
	reslt = fgets (line, sizeof (line), mapin);	// 
	reslt = fgets (line, sizeof (line), mapin);	// 

	for (i = cstart; i < cend; i++) {
	    for (j = bstart; j < bend; j++) {
		for (k = astart; k < aend; k++) {
		    ijk = k * new_mapstep_b * new_mapstep_c + j * new_mapstep_c + i;
		    ireslt = fscanf (mapin, "%E", &FourierPt[ijk]);
		    if (FourierPt[ijk] < rhomin)
			rhomin = FourierPt[ijk];
		    if (FourierPt[ijk] > rhomax)
			rhomax = FourierPt[ijk];
		}
	    }
	}
    }

    if (ijk > 0) {
	if (!Quick)
	    fprintf (drvui->flout, "Reading Fourier map file %s\n  "
		     "cell = %f %f %f %f %f %f\n  map steps: a=%d b=%d c=%d\n  "
		     "Range of rho values: %f to %f, gridsize %d\n", infile,
		     map_a, map_b, map_c, map_alpha, map_beta, map_gamma,
		     mapstep_a, mapstep_b, mapstep_c, rhomin, rhomax, ijk);
	ReadFourMap = 1;
    } else {
	sprintf (string, "Error reading Fourier map file %s", infile);
	Error_Box (string);
    }
    fclose (mapin);
    Map_Info.rhomn = rhomin;
    Map_Info.rhomx = rhomax;
    Map_Info.map_int[0] = mapstep_a = new_mapstep_a;
    Map_Info.map_int[1] = mapstep_b = new_mapstep_b;
    Map_Info.map_int[2] = mapstep_c = new_mapstep_c;
    Map_Info.info_valid = 1;
}
