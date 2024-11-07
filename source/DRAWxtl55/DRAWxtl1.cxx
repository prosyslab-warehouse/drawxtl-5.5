// $Id: DRAWxtl1.cxx 1114 2011-02-23 20:29:18Z martin $
//
// module drawxtl1.cxx - part of DRAWxtl V5.5 - the GUI version
// Coded using the FLTK 1.1.6 widget set
//
// This program is copyrighted by Larry W. Finger, Martin Kroeker and Brian Toby
// and is distributed under the GNU General Public License - see the
// accompanying COPYING file for more details.
//
// This module contains the following routines:
//
//  add_to_list - add atom position to unit-cell list
//  add_vert - add vertices (atom positions) to display lists
//  add_vert_nc - add vertices without checking
//  analyze_bonds - print the bond distance tables
//  axeqb - solve matrix equation AX = B for B (A is 3x3 only)
//  build_box_contents - find all atoms in the display parallelopiped
//  Calc_Rot - calculate rotation matrices
//  convert_ellipsoid - handles conversions from Bij or Uij to betaij
//  Conv_Sym_Mat - converts rotational part of symmetry matrix to a Cartesian rotation
//  dist - calculates distance between two vertices
//  dot0_3d - calculated dot product of two 3D vectors
//  eigen - calculates eigen values and vectors for ellipsoids
//  find_all_in_box - finds all atoms of a given type in the display box
//  generate_arrows - generates magnetic vector arrows
//  generate_bonds - Generate bond descriptions for the display lists
//  generate_cones - generate lone-pair cone descriptions for display lists
//  get_atom_id - reworks atom labels
//  get_input - calls actual routine that reads input, generates lattice metric and does preliminary ellipsoid processing
//  Locate_Triple - routine to place unit-cell vector triple
//  make_bmat - calculates lattice metric routines for conversion to cartesian coordinates
//  modulate_parameters - do the adjustments to positions and occupancies for aperiodic xtals
//  modulate_uij - do the adjustments to Uij's for aperiodic crystals
//  polygon_normal_3d - calculates normal vector of a polygon
//  polygon_solid_angle_3d - calculates projected solid angle of a 3d plane polygon
//  not_in_slab - check if a point is inside a given parallelepiped
//  Output_Spheres - adds sphere descriptions to output display lists
//  plot_vrml_poly - generates polyhedral descriptions to VRML output
//  print_sym - print symmetry operators

#include "drawxtl.h"
#include "DRAWxtlViewUI.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <FL/glut.H>
#if defined(__APPLE__)
#  include <openGL/glu.h>
#else
#  include <GL/glu.h>
#endif
#ifdef WIN32
#define cbrt(val) pow((val), 1.0 / 3.0)
#endif
// forward references for routines

#include "DRAWxtl_proto.h"

// global variables

#include "draw_ext.h"

/* ************************************************************** */
/* ************************************************************** */

void
add_to_list (float xp[3], int no, int nos)

/* add atom position to unit-cell list if unique */
/* no - symmetry operator number */
/* nos - signed symmetry operator number to track inversion */
{
    int i, ind;

    for (i = 0; i <= 2; ++i) {	/* get coordinates between 0 and 1 */
	if (xp[i] >= 1.0)
	    xp[i] = xp[i] - 1.0f;
	if (xp[i] < 0.0)
	    xp[i] = xp[i] + 1.0f;
    }
    ind = 0;
    for (i = 0; i < ncell; ++i) {
	if ((fabs (drvui->cell_xyz[i][0] - xp[0]) <= 0.00002) &&
	    (fabs (drvui->cell_xyz[i][1] - xp[1]) <= 0.00002) &&
	    (fabs (drvui->cell_xyz[i][2] - xp[2]) <= 0.00002))
	    ind = 1;		/*set if same */
    }
    if (ind == 0) {		/* true for new position */
	for (i = 0; i <= 2; ++i)	/*copy to list */
	    drvui->cell_xyz[ncell][i] = xp[i];
	drvui->sym_op_signed[ncell] = nos;
	drvui->sym_op_no[ncell++] = no;
    }
}

/* ************************************************************** */
/* ************************************************************** */

void
add_vert (float vert[3], int no, int check, int modulate, int sign)

/* procedure to add the vertex in 'vert' to the current list if the
   vertex is contained within the limits of the box, within the limits
   of a polyhedral center, or needed to complete a molecule. If modulate
   is true (!= 0), the coordinates in vert will be adjusted for modulation
   before testing */
{
    float c_vert[3];		/* vertex in cartesian coords */

    float r_vert[3];		/* vertex in rotated coordinates */

    float vert_saved[3];	/* save the unmodulated vetex */

    int i, j;			/* loop variables */

    float temp;			/* temporary storage */

    int outside;		/* temporary for molecular completion */

    int number;			/* decoded atom number */

    int poly_center;		/* true if this atom is center of a polyhedron */

    int consider_extra;

    float d = 0.0f;

    double occupancy;

    int sym_no;

    if (!check_vert_alloc (NvertM, 0))
	return;			//abort immediately if out of array space
    number = no / 1000;
    sym_no = no - 1000 * number;
    for (i = 0; i < 3; i++)
	vert_saved[i] = vert[i];
    if (modulate & (drvui->modulated > 0)) {
	modulate_parameters (vert, &occupancy, sym_no, number);
	if (occupancy < 0.01)
	    return;
    }
    for (i = 0; i < 3; ++i) {	/* convert vertex coordinates to Cartesian */
	c_vert[i] = 0.0f;
	for (j = 0; j < 3; ++j)
	    c_vert[i] += (float) drvui->b_mat[i][j] * (vert[j] - origin[j]);
    }
    if (number >= 0 && number < drvui->atom_alloc)
	poly_center = (((drvui->atoms[number].atom_n >> 8) & 255) > 0);
    else
	poly_center = 0;
    consider_extra = (domolcomp != 0 || (drvui->npoly > 1 && !poly_center));
    outside = 0;
    if (check == 0) {
	for (i = 0; i <= 2; i++) {	/* check packing range */
	    if ((vert[i] < drvui->frames[drvui->frame_no].cryst_lim[i]) ||
		(vert[i] > drvui->frames[drvui->frame_no].cryst_lim[i + 3]))
		outside = 1;
	}
	if (outside == 0 && slabmode == 1)
	    outside = not_in_slab (c_vert[0], c_vert[1], c_vert[2]);

	if (((outside == 0) && (fabs (c_vert[0]) > boxlim[0])) ||
	    (fabs (c_vert[1]) > boxlim[1]) || (fabs (c_vert[2]) > boxlim[2]))
	    outside = 1;	/* check box limits */
	if (outside != 0 && !consider_extra)
	    return;
    }
    d = 0.0f;
    if (outside == 1) {		/* molecule or polyhedral completion is in progress */
	float dmin, dmax;

	int polyno;		/* number of the polyhedron */

	for (i = 1; i < NvertM; i++) {
	    if (drvui->atoms[drvui->atom_no[i] / 1000].atom_fn != drvui->frame_no)
		continue;

	    dmax = drvui->mol_d * drvui->mol_d;
	    dmin = 0.;
	    d = 0.0f;
	    for (j = 0; j < 3; j++) {
		temp = (s_vert[3 * i + j] - c_vert[j]);
		d += temp * temp;
	    }
	    j = drvui->atom_no[i] / 1000;
	    poly_center = (((drvui->atoms[j].atom_n >> 8) & 255) > 0);
	    if (poly_center) {
		polyno = (drvui->atoms[j].atom_n >> 8) & 255;
		dmax =
		    drvui->polyhedra[polyno].poly_size *
		    drvui->polyhedra[polyno].poly_size;
		dmin =
		    drvui->polyhedra[polyno].poly_min * drvui->polyhedra[polyno].poly_min;
	    }
	    if ((domolcomp != 0 || poly_center) && (d <= dmax && d > 0.02)) {
		outside = 2;
		if (poly_center && d < dmin)
		    outside = 1;
		if ((fabs (c_vert[0]) > 10. * boxlim[0])
		    || (fabs (c_vert[1]) > 10. * boxlim[1])
		    || (fabs (c_vert[2]) > 10. * boxlim[2])) {
		    fprintf (drvui->fcns,
			     "mol'completing beyond 10*lattice constant - aborting run\n");
		    outside = 1;	/* far out of check box limits - combinatorial explosion */
		    break;
		}
	    }
	}
    }
    if (outside == 1)
	return;
    i = (int) (no / 1000);
    for (j = 1; j < NvertM; ++j) {
	if (fabs (s_vert[3 * j] - c_vert[0]) < 0.01
	    && fabs (s_vert[3 * j + 1] - c_vert[1]) < 0.01
	    && fabs (s_vert[3 * j + 2] - c_vert[2]) < 0.01
	    && check_atom_name (drvui->atoms[i].atom_l,
				drvui->atoms[drvui->atom_no[j] / 1000].atom_l)
	    && (number == i))
	    return;		/* eliminate duplicates */
    }
    for (i = 0; i <= 2; ++i) {	/* add atom to Master Lists */
	xypos[3 * NvertM + i] = vert[i];
	xypos_nm[3 * NvertM + i] = vert_saved[i];
	s_vert[3 * NvertM + i] = c_vert[i];
//    if (offset[i] > c_vert[i])
//      offset[i] = c_vert[i];
    }
    i = no - 1000 * (no / 1000);

//if (occupancy != 1.) fprintf (stderr,"occ=%f\n",occupancy);
    drvui->vert_occ[NvertM] = (float) occupancy;
    drvui->atom_so[NvertM] = sign;
    drvui->atom_no[NvertM++] = no;
    if (!check_vert_alloc (NvertM, 0)) {
	if (domolcomp != 0) {
	    Error_Box ("Too many vertices for dimensions:\n"
		       " Did you use a 'molcomp' command incorrectly?\n"
		       " If not, please increase parameter MAX_VERTS.");
	} else {
	    Error_Box
		("Too many vertices for dimensions - please increase parameter MAX_VERTS.");
	}
	return;
    }
    for (i = 0; i <= 2; ++i) {	/* calculate position of point after POV rotation */
	r_vert[i] = 0.0f;
	for (j = 0; j <= 2; ++j)
	    r_vert[i] += (float) G_Rot[j][i] * c_vert[j];
    }
    if (drvui->automation)
	return;			/* if automation in progress, keep old scaling */
    for (i = 0; i <= 2; ++i) {	/* Update Min and Max of output */
	if (POV_Max[i] < r_vert[i])
	    POV_Max[i] = r_vert[i];
	if (POV_Min[i] > r_vert[i])
	    POV_Min[i] = r_vert[i];
    }
}

/* ************************************************************** */
/* ************************************************************** */

void
add_vert_nc (float vert[3])

/* procedure to add the vertices in 'vert' to the current working list
   (NOT Master List ) without checking - NO modulation */
{

    int i, j;			/* loop variables */

    if (!s_vert) {
	if ((s_vert = (float *) zalloc ((unsigned) (3 * (long) (drvui->verts_alloc
								* sizeof (float))))) ==
	    NULL) {
	    Error_Box ("Unable to get initial s_vert allocation");
	    return;
	}
    }
    for (i = 0; i <= 2; ++i) {	/* convert vertex coordinates to Cartesian */
	s_vert[3 * nvert + i] = 0.0f;
	for (j = 0; j <= 2; ++j)
	    s_vert[3 * nvert + i] += (float) drvui->b_mat[i][j] * (vert[j] - origin[j]);
    }
    nvert++;
}

/* ************************************************************** */
/* ************************************************************** */

void
analyze_bonds (void)

/* routine to print bond distances up to a maximum of 'printdist' input units */
{
    int i, j, k, l, nvert1, *itype;

    float *d, dmin;

    double occ;

    if (!(itype = (int *) zalloc ((unsigned) (drvui->verts_alloc * sizeof (int))))) {
	Error_Box ("Unable to assign space for bond analysis.");
	return;
    }
    fprintf (drvui->flout, "\n\n Bond Distance Analysis\n");
/* Find all atoms of this type in master list */
    l = nvert = nvert1 = 1;
    for (j = 0; j < natom; ++j) {
	if (drvui->atoms[j].atom_fn != drvui->frame_no)
	    continue;
	l = nvert;
	find_all_in_box (j);
	for (k = l; k < nvert; ++k)
	    itype[k] = j;	// itype[k] has atom number j
    }
    for (i = 0; i < natom; ++i) {	/* loop through atom types */
	if (drvui->atoms[i].atom_fn != drvui->frame_no)
	    continue;
/* add the base entry to the beginning of the atom list  - picked from vertex list
   because position originating from x,y,z symmetry operator may be absent due to
   crenel exclusion */
	k = -1;
	for (l = 1; l < nvert; l++) {	/* find an entry for atom i */
	    if (itype[l] == i) {
		for (j = 0; j < 3; ++j) {
		    o_vert[j] = o_vert[3 * l + j];
		    o_vert_nm[j] = o_vert_nm[3 * l + j];
		    s_vert[j] = s_vert[3 * l + j];
		}
		k = l;
		break;
	    }
	}
	if (k < 0) {
/* There is no entry for atom i in the list due to limits, etc. Add the base one. */
	    for (j = 0; j < 3; j++) {
		s_vert[j] = 0;
		o_vert[j] = drvui->atoms[i].atom_xyz[j];
		o_vert_nm[j] = o_vert[j];
		for (l = 0; l < 3; l++)
		    s_vert[j] +=
			(float) drvui->b_mat[l][j] * (drvui->atoms[i].atom_xyz[l] -
						      origin[l]);
	    }
	}
	if ((d = (float *) zalloc ((unsigned) (nvert * sizeof (float)))) == NULL) {
	    free (itype);
	    Error_Box ("Unable to get bond distance allocation");
	    return;
	}
	for (j = 1; j < nvert; ++j)	/* calculate distances */
	    d[j] = dist (0, j);
	fprintf (drvui->flout,
		 "\nDistances from %c%c%c%c%3d at %8.5f%8.5f%8.5f    to\n\n",
		 drvui->atoms[i].atom_l[0], drvui->atoms[i].atom_l[1],
		 drvui->atoms[i].atom_l[2], drvui->atoms[i].atom_l[3],
		 drvui->atoms[i].sv_atom_n, o_vert[0], o_vert[1], o_vert[2]);
	for (k = nvert1; k < nvert; ++k) {
	    dmin = 999999.9f;
	    for (j = nvert1; j < nvert; ++j) {	/* loop through distances */
		if (d[j] < dmin) {
		    l = j;
		    dmin = d[j];
		}
	    }
	    if (dmin == 0.0 && itype[l] == i) {
		d[l] = 999999.9f;
	    } else {
		if (dmin <= printdist) {
		    fprintf (drvui->flout, " %c%c%c%c%3d at %9.5f%9.5f%9.5f  %8.3f",
			     drvui->atoms[itype[l]].atom_l[0],
			     drvui->atoms[itype[l]].atom_l[1],
			     drvui->atoms[itype[l]].atom_l[2],
			     drvui->atoms[itype[l]].atom_l[3],
			     drvui->atoms[itype[l]].sv_atom_n, o_vert[3 * l],
			     o_vert[3 * l + 1], o_vert[3 * l + 2], dmin);
		    d[l] = 999999.9f;
		    if (drvui->atoms[i].atom_ismod != 0 || drvui->atoms[itype[l]].atom_ismod != 0) {	/* one or both ends are modulated */
			float temp_pa[3];

			float vert1[3], vert2[3], vert3[3];

			float dsum = 0.0f, dmax = 0.0f, thisd;

			int m, mx[3], m1, m2, number = 0;

			dmin = 999999.0f;
			for (m = 0; m < 3; m++) {
			    temp_pa[m] = drvui->phaseshift[m];
			    mx[m] = 50;
			    if (drvui->modulated <= m)
				mx[m] = 1;
			}
			for (m = 0; m < mx[0]; m++) {	/* sum over values of initial phaseshift for first modulation */
			    drvui->phaseshift[0] = m * 0.02f;
			    for (m1 = 0; m1 < mx[1]; m1++) {
				drvui->phaseshift[1] = m1 * 0.02f;	/* second modulation */
				for (m2 = 0; m2 < mx[2]; m2++) {
				    drvui->phaseshift[2] = m2 * 0.02f;	/* third modulation */
				    for (j = 0; j < 3; j++) {
					vert1[j] = o_vert_nm[j];
					vert2[j] = o_vert_nm[3 * l + j];
				    }
				    modulate_parameters (vert1, &occ, 0, i);
				    modulate_parameters (vert2, &occ, vert_sym_no[l],
							 itype[l]);
				    for (j = 0; j < 3; ++j) {	/* convert vertex coordinates to Cartesian */
					vert3[j] = 0.0f;
					for (k = 0; k <= 2; ++k)
					    vert3[j] +=
						(float) drvui->b_mat[j][k] * (vert2[k] -
									      vert1[k]);
				    }
				    thisd =
					(float) sqrt (vert3[0] * vert3[0] +
						      vert3[1] * vert3[1] +
						      vert3[2] * vert3[2]);
				    dsum += thisd;
				    dmin = min (thisd, dmin);
				    dmax = max (thisd, dmax);
				    number++;
				}
			    }
			}
			dsum /= (float) number;	/* get mean distance */
			for (m = 0; m < 3; m++)
			    drvui->phaseshift[m] = temp_pa[m];
			fprintf (drvui->flout, "%8.3f %8.3f %8.3f", dsum, dmin, dmax);
		    }
		    fprintf (drvui->flout, "\n");
		} else
		    break;
	    }
	}
	free (d);
    }
    free (itype);
    fprintf (drvui->flout, "\n\n");
}

/* ************************************************************** */
/* ************************************************************** */

void
axeqb (double a1[3][3], double x[3], double b[3])

/* solve matrix equation ax = b for x, where a is 3x3 and x and b are
   vectors */
{
    int i, j;

    double a[3][3], det;

    det = determinant (a1);
    if (fabs (det) < 1.0e-6)
	det = 1.0e-6;
    a[0][0] = (a1[1][1] * a1[2][2] - a1[1][2] * a1[2][1]) / det;
    a[1][0] = -(a1[1][0] * a1[2][2] - a1[1][2] * a1[2][0]) / det;
    a[2][0] = (a1[1][0] * a1[2][1] - a1[1][1] * a1[2][0]) / det;
    a[0][1] = -(a1[0][1] * a1[2][2] - a1[0][2] * a1[2][1]) / det;
    a[1][1] = (a1[0][0] * a1[2][2] - a1[0][2] * a1[2][0]) / det;
    a[2][1] = -(a1[0][0] * a1[2][1] - a1[0][1] * a1[2][0]) / det;
    a[0][2] = (a1[0][1] * a1[1][2] - a1[0][2] * a1[1][1]) / det;
    a[1][2] = -(a1[0][0] * a1[1][2] - a1[0][2] * a1[1][0]) / det;
    a[2][2] = (a1[0][0] * a1[1][1] - a1[0][1] * a1[1][0]) / det;
    for (i = 0; i < 3; i++) {
	x[i] = 0.0;
	for (j = 0; j < 3; j++)
	    x[i] += a[i][j] * b[j];
    }
}

/* ************************************************************** */
/* ************************************************************** */

void
build_box_contents (void)

/* routine to locate all atoms in the display parallelopiped plus
   any atoms needed to complete molecules or polyhedra */
{
    float vert[3];		/* storage for center */

    int ix[3];			/* array for lattice offset */

    int j, k, saved_nvert;

    int atoms, loop;

    loop = 0;
    while (1) {			// infinite loop here
	saved_nvert = NvertM;	// save number currently in list
	for (atoms = 0; atoms < natom; atoms++) {	// loop through atoms
	    if (drvui->atoms[atoms].atom_fn != drvui->frame_no)
		continue;
	    for (j = 0; j < 3; j++) {	// put the base atom in the range 0-1
		drvui->atoms[atoms].saved_xyz[j] = drvui->atoms[atoms].atom_xyz[j];
		if (drvui->atoms[atoms].atom_xyz[j] < 0.0f)
		    drvui->atoms[atoms].atom_xyz[j] += 1.0f;
		if (drvui->atoms[atoms].atom_xyz[j] > 1.0f)
		    drvui->atoms[atoms].atom_xyz[j] -= 1.0f;
	    }
	    expand_atom (atoms);	// generate all atoms of this flavor in unit cell
	    for (j = 0; j < ncell; ++j) {	// atoms in cell
		for (k = 0; k <= 2; ++k)
		    vert[k] = drvui->cell_xyz[j][k];
		add_vert (vert, 1000 * atoms + drvui->sym_op_no[j], 0, 1, drvui->sym_op_signed[j]);	// add point to cell if inside box
	    }
	}
	if (saved_nvert == NvertM)
	    break;
    }
    while (1) {			/* infinite loop here */
	saved_nvert = NvertM;	/* save number currently in list */
	for (atoms = 0; atoms < natom; atoms++) {	/* loop through atoms */
	    if (drvui->atoms[atoms].atom_fn != drvui->frame_no)
		continue;
	    expand_atom (atoms);	/* generate all atoms of this flavor in unit cell */
	    for (j = 0; j < ncell; ++j) {	/*atoms in cell */
		for (k = 0; k <= 2; ++k)
		    vert[k] = drvui->cell_xyz[j][k];
		add_vert (vert, 1000 * atoms + drvui->sym_op_no[j], 0, 1, drvui->sym_op_signed[j]);	/* add point to cell if inside box */
	    }
	    for (ix[0] = (int) drvui->frames[drvui->frame_no].cryst_lim[0] - 1; ix[0] <= (int) drvui->frames[drvui->frame_no].cryst_lim[3] + 1; ++ix[0]) {	/* step through cells parallel a */
		for (ix[1] = (int) drvui->frames[drvui->frame_no].cryst_lim[1] - 1; ix[1] <= (int) drvui->frames[drvui->frame_no].cryst_lim[4] + 1; ++ix[1]) {	/* parallel b */
		    for (ix[2] = (int) drvui->frames[drvui->frame_no].cryst_lim[2] - 1; ix[2] <= (int) drvui->frames[drvui->frame_no].cryst_lim[5] + 1; ++ix[2]) {	/* parallel c */
			if (abs (ix[0]) + abs (ix[1]) + abs (ix[2]) != 0) {
			    for (j = 0; j < ncell; ++j) {	/*atoms in cell */
				for (k = 0; k <= 2; ++k)
				    vert[k] = ix[k] + drvui->cell_xyz[j][k];
				add_vert (vert, 1000 * atoms + drvui->sym_op_no[j], 0, 1, drvui->sym_op_signed[j]);	/* add point to cell if inside box */
			    }
			}
		    }
		}
	    }
	}
	if ((domolcomp == 0 && drvui->npoly == 1) || saved_nvert == NvertM) {
	    for (atoms = 0; atoms < natom; atoms++) {	/* loop through atoms */
		if (drvui->atoms[atoms].atom_fn != drvui->frame_no)
		    continue;
		for (j = 0; j < 3; j++) {	// restore base atom position
		    drvui->atoms[atoms].atom_xyz[j] = drvui->atoms[atoms].saved_xyz[j];
		}
	    }
	    return;
	}
	if (loop == 0)
	    fprintf (drvui->fcns, "\n");
	loop++;
	fprintf (drvui->fcns,
		 "Recursive Atom Addition, Loop No. %3d, %4d Atoms in List.\n", loop,
		 NvertM - 1);
    }
}				/* end of build_box_contents */


/* ************************************************************** */
/* ************************************************************** */

void
Calc_Rot (float v1[3], float v2[3])

/* routine to calculate rotation axes needed to put v1 parallel to Z
   and the projection of v2 parallel to Y */
{
    float u1[3], u2[3], u3[3], temp;

    double U[3][3];

    int i, j;

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    U[i][j] = 0.0f;	/* clear U */
    temp = 0.0f;
    for (i = 0; i < 3; i++) {	/* convert v1 to Cartesian */
	u1[i] = 0.0f;
	for (j = 0; j < 3; j++)
	    u1[i] += (float) drvui->b_mat[i][j] * v1[j];
	temp += u1[i] * u1[i];
    }
    if (temp <= 0.0) {
	u1[2] = 1.0f;		/* if null input, make u1 = 0,0,1 */
	U[2][2] = 1.0f;
    } else {
	temp = (float) sqrt (temp);
	for (i = 0; i < 3; i++) {	/* make a unit vector */
	    u1[i] /= temp;
	    U[i][2] = u1[i];	/* put column in U */
	}
    }
    for (i = 0; i < 3; i++) {	/* convert v2 to Cartesian */
	u2[i] = 0.0f;
	for (j = 0; j < 3; j++)
	    u2[i] += (float) drvui->b_mat[i][j] * v2[j];
    }
    u3[0] = u2[1] * u1[2] - u1[1] * u2[2];	/* calculate u2 X u1 */
    u3[1] = u2[2] * u1[0] - u2[0] * u1[2];
    u3[2] = u2[0] * u1[1] - u1[0] * u2[1];
    temp = 0.0f;
    for (i = 0; i < 3; i++)
	temp += u3[i] * u3[i];
    if (temp <= 0.0) {
	u3[0] = 1.0f;		/* if null input, make u3 = 1,0,0 */
	U[0][0] = 1.0f;
    } else {
	temp = (float) sqrt (temp);
	for (i = 0; i < 3; i++) {	/* make a unit vector */
	    u3[i] /= temp;
	    U[i][0] = u3[i];	/* add column to U */
	}
    }
    u2[0] = u1[1] * u3[2] - u3[1] * u1[2];	/* calculate u2 = u1 X u3 */
    u2[1] = u1[2] * u3[0] - u1[0] * u3[2];
    u2[2] = u1[0] * u3[1] - u3[0] * u1[1];
    for (i = 0; i < 3; i++)
	U[i][1] = u2[i];
/* make determinant of U > 0 */
    if (determinant (U) < 0.0) {
	for (i = 0; i < 3; i++)
	    U[i][0] = -U[i][0];
    }

/* The rotation required to get the desired orientation, which is described
    by G_Rot, has been set up so G_Rot * U = I. Thus G_Rot = U(inv),
    U is unitary and U(inv) = U(transpose). */

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    G_Rot[i][j] = (float) U[j][i];
/* calculate rotations */
    temp = (float) (G_Rot[1][2] * G_Rot[1][2] + G_Rot[2][2] * G_Rot[2][2]);
    if (temp > 0.001) {
	float temp2;

	Yrot = (float) (atan (G_Rot[0][2] / sqrt (temp)) * RAD);
	temp2 = (float) (G_Rot[1][2] / cos (Yrot / RAD));
	if (temp2 > 1.0)
	    temp2 = 1.0f;
	if (temp2 < -1.0)
	    temp2 = -1.0f;
	Xrot = -(float) asin (temp2) * (float) RAD;
	temp2 = -(float) (G_Rot[0][1] / cos (Yrot / RAD));
	if (temp2 > 1.0)
	    temp2 = 1.0f;
	if (temp2 < -1.0)
	    temp2 = -1.0f;
	Zrot = (float) asin (temp2) * (float) RAD;
    } else {
	Yrot = -(float) asin (G_Rot[0][2]) * (float) RAD;
	Xrot = 0.0f;		/* XXXXXXXXXXX  - May not be right */
	Zrot = 0.0f;		/* XXXXXXXXXXX  - May not be right */
    }
}


/* ************************************************************** */
/* ************************************************************** */

void
convert_ellipsoid (void)

/* routine to convert ellipsoid coefficients Bij or Uij to betaij form */
{
    int i, j, save;

    for (i = 1; i < drvui->n_ellips; ++i) {
	if (drvui->ellips[i].ellips_ismod != 0)
	    continue;
	save = drvui->ellips[i].ell_type / 1000;
	if (drvui->ellips[i].ell_type % 1000 != 0) {	// multiply Uij and Bij by a*(i)*a*(j)/4
	    for (j = 0; j < 3; ++j)
		drvui->ellips[i].ellips[j] *=
		    drvui->rec_lat_con[j] * drvui->rec_lat_con[j] * 0.25f;
	    drvui->ellips[i].ellips[3] *=
		drvui->rec_lat_con[0] * drvui->rec_lat_con[1] * 0.25f;
	    drvui->ellips[i].ellips[4] *=
		drvui->rec_lat_con[0] * drvui->rec_lat_con[2] * 0.25f;
	    drvui->ellips[i].ellips[5] *=
		drvui->rec_lat_con[1] * drvui->rec_lat_con[2] * 0.25f;
	    if (drvui->ellips[i].ell_type % 1000 == 1)	// multiply Uij by 8pi^2
		for (j = 0; j < 6; ++j)
		    drvui->ellips[i].ellips[j] *= 78.9568f;
	}
	drvui->ellips[i].ell_type = 1000 * save;
    }
}				/* end of convert_ellipsoid */

/* ************************************************************** */
/* ************************************************************** */

void
Conv_Sym_Mat (void)

/* routine to convert the rotational part of a symmetry matrix to a Cartesian rotation */
{
    int i, j, k, l;

    double mat[3][3];

    float B[3][3], Binv[3][3];

    float tmp[3][3];

    if (drvui->sys == 5) {	/* hexagonal systems are special */
	for (i = 0; i < 3; i++)
	    for (j = 0; j < 3; j++) {
		B[i][j] = (float) drvui->b_mat[i][j];
		Binv[i][j] = B[i][j];
	    }
	matinv (Binv);
    }
    for (l = 0; l < drvui->ng; l++) {
	for (j = 0; j < 3; j++) {
	    for (k = 0; k < 3; k++) {
		mat[j][k] = drvui->rss[l][j][k];
	    }
	}
	if (determinant (mat) < 0.0) {	/* if improper rotation, negate all elements */
	    for (j = 0; j < 3; j++) {
		for (k = 0; k < 3; k++) {
		    mat[j][k] *= -1.0f;
		}
	    }
	}
	for (j = 0; j < 3; j++)
	    for (k = 0; k < 3; k++)
		drvui->rssC[l][j][k] = (float) mat[j][k];	/* save a copy of the original */

	if (drvui->sys == 5) {	/* hexagonal systems are special */
	    float R[3][3];

	    for (j = 0; j < 3; j++)
		for (k = 0; k < 3; k++)
		    R[j][k] = (float) mat[j][k];
	    matmul (B, R, tmp);
	    matmul (tmp, Binv, R);	/* The Cartesian matrix is B R Binv */
	    for (j = 0; j < 3; j++)
		for (k = 0; k < 3; k++)
		    mat[j][k] = R[j][k];
	}
	for (j = 0; j < 3; j++) {	/* restore matrix */
	    for (k = 0; k < 3; k++) {
		drvui->rss[l][j][k] = (float) mat[j][k];
	    }
	}
    }
}

/* ************************************************************** */
/* ************************************************************** */

float
dist (int j, int k)

/* calculates distance between vertices at s_vert[j] and s_vert[k] */
{
    float d, t;

    int i;

    d = 0.0f;
    for (i = 0; i <= 2; ++i) {
	t = s_vert[3 * j + i] - s_vert[3 * k + i];
	d += t * t;
    }
    return ((float) sqrt (d));
}

/* ************************************************************** */
/* ************************************************************** */

float
dot0_3d (float x0, float y0, float z0, float x1, float y1, float z1,
	 float x2, float y2, float z2)
{
// 3D dot product of (P1-P0) and (P2-P0), John Burkardt

    return (x1 - x0) * (x2 - x0) + (y1 - y0) * (y2 - y0) + (z1 - z0) * (z2 - z0);
}

/* ************************************************************** */
/* ************************************************************** */

int
eigen (float *biso, float beta[6], float valu[3], float vect[3][3])

/* routine to calculate eigenvalues 'valu' and vectors 'v' for
 * array of anisotropic thermal parameters 'beta', 'biso' is the equivalent
 * isotropic B
 *
 * some parts of this routine were copied from ORTEP
 */
{
    static double errnd = 1.0E-9;

    double b0, b1, b2;

    double a1, a2, phi;

    double a[3][3], b[3][3], w[3][3], u[3];

    float v1[3], v2[3], v3[3];

    double tem, sigma, smax;

    double Vect[3][3];

    int i, j, l, ii, iii, i1, imax = 0;

    int neq;

    static double tpi2 = 2.0 * M_PI * M_PI;	/* 2 pi^2 */

/* put beta into square symmetric matrix b (called M in ORTEP) */

    for (i = 0; i < 3; i++)
	b[i][i] = beta[i];
    b[0][1] = b[1][0] = beta[3];
    b[0][2] = b[2][0] = beta[4];
    b[1][2] = b[2][1] = beta[5];
/* multiply b * ginv to get w */
    for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++) {
	    w[i][j] = 0.0;
	    for (l = 0; l < 3; l++) {
		w[i][j] += b[l][i] * drvui->ginv[j][l];
	    }
	}
    }
    sigma = 0.0;
    for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++) {
	    a[i][j] = w[i][j];
	    sigma += a[i][j] * a[i][j];
	}
    }
    if (sigma <= 0.0)
	return (0);		/* error on null 'w' */
    sigma = sqrt (sigma);
/* get coefficients of third-order characteristic equation 
 *        Equation is -L^3 + b2 L^2 + b1 L + b0 = 0
 */
    b2 = a[0][0] + a[1][1] + a[2][2];
    b1 = -a[0][0] * a[1][1] - a[0][0] * a[2][2] - a[1][1] * a[2][2]
	+ a[0][2] * a[2][0] + a[0][1] * a[1][0] + a[1][2] * a[2][1];
    b0 = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2]
	+ a[2][0] * a[1][2] * a[0][1] - a[2][0] * a[0][2] * a[1][1]
	- a[0][0] * a[2][1] * a[1][2] - a[1][0] * a[0][1] * a[2][2];

    *biso = (float) (4.0 * b2 / 3.0);
/* transform to get to form x^3 + a1 x + a2  == 0 */
    a1 = (-3.0 * b1 - b2 * b2) / 3.0;
    a2 = (-2.0 * b2 * b2 * b2 - 9.0 * b2 * b1 - 27.0 * b0) / 27.0;
    tem = (a2 * a2) / 4.0 + (a1 * a1 * a1) / 27.0;
    if (tem > errnd) {
	/* two complex roots */
	fprintf (drvui->flout, "Complex Roots! tem = %f\n", tem);
	return (0);		/* error if complex roots */
    } else if (tem < -errnd) {
	/* 3 real and unequal roots */
	phi = acos (-(a2 / 2.0) / sqrt (-(a1 * a1 * a1) / 27.0));
	u[0] = 2.0 * sqrt (-a1 / 3.0) * cos (phi / 3.0) + b2 / 3.0;
	u[1] = 2.0 * sqrt (-a1 / 3.0) * cos (phi / 3.0 + 120.0 / RAD) + b2 / 3.0;
	u[2] = 2.0 * sqrt (-a1 / 3.0) * cos (phi / 3.0 + 240.0 / RAD) + b2 / 3.0;
    } else {
	double sign_a2 = 1.0;
	/* Workaround brain-dead cbrt routine in Windows */
	if (a2 < 0.0) {
	    sign_a2 = -1.0;
	    a2 = fabs(a2);
	}
	/* 3 real roots with at least 2 being equal */
	u[0] = -2.0 * sign_a2 * cbrt (a2 / 2.0) + b2 / 3.0;
	u[1] = sign_a2 * cbrt (a2 / 2.0) + b2 / 3.0;
	u[2] = u[1];
    }
    for (j = 0; j < 2; j++) {	/* sort roots in increasing order */
	if (u[j] > u[j + 1]) {
	    tem = u[j + 1];
	    u[j + 1] = u[j];
	    u[j] = tem;
	}
    }
    for (i = 0; i < 3; i++)
	if (valu[i] < 0.0)
	    valu[i] = 0.001f;
	else
	    valu[i] = (float) sqrt (u[i] / tpi2);

/* count multiple roots */
    neq = 1;
    for (j = 0; j < 2; j++)
	if (fabs (u[j] - u[j + 1]) < errnd)
	    neq++;
    if (neq > 2) {
/* 3 equal roots */
	for (ii = 0; ii < 3; ii++) {
	    for (i = 0; i < 3; i++)
		vect[i][ii] = 0.0;
	    vect[ii][ii] = 1.0;
	}
	return (1);
    } else if (neq == 2) {
	for (ii = 0; ii < 3; ii++) {
	    v1[ii] = v2[ii] = v3[ii] = 0.0;
	}
	if (drvui->sys == 5) {	// hexagonal 
/* With 2 equal eigenvalues, we have degenerate vectors */
	    if (fabs (u[0] - u[1]) < errnd) {
		v3[2] = 1.0f;	// 3rd vector is unique
		v1[0] = 1.0f;
	    } else {
		v3[0] = 1.0f;	// 1st one is
		v1[2] = 1.0f;
	    }
	} else if (drvui->sys == 6) {		// cubic with atom at xxx
	    if (fabs (u[0] - u[1]) < errnd) {
		v3[0] = v3[1] = v3[2] = (float) sqrt (3.0) / 3.0f;
		v1[0] = -(float) sqrt (2.0) / 2.0f;
		v1[1] = -v1[0];
	    } else {
		v1[0] = v1[1] = v1[2] = (float) sqrt (3.0) / 3.0f;
		v3[0] = -(float) sqrt (2.0) / 2.0f;
		v3[1] = -v3[0];
	    }
	} else {		// tetragonal
	    if (fabs (u[0] - u[1]) < errnd) {
		v3[2] = 1.0f;	// 3rd vector is unique
		v2[1] = 1.0f;
		v1[0] = 1.0f;
	    } else {
		v3[0] = 1.0f;	// 1st one is
		v1[1] = 1.0f;
		v2[2] = 1.0f;
	    }
	}
	vcross (v3, v1, v2);	// v1, v2 and v3 are mutually perp. unit vectors
	vnormalize (v2);
	for (j = 0; j < 3; j++) {
	    vect[0][j] = v1[j];
	    vect[1][j] = v2[j];
	    vect[2][j] = v3[j];
	}
	return (1);
    }
    for (ii = 0; ii < 3; ii++) {
	for (iii = 0; iii < 3; iii++) {
	    for (i = 0; i < 3; i++)
		a[i][iii] = w[i][iii];
	    a[iii][iii] = w[iii][iii] - u[ii];
	}
	smax = 0.0;
	for (i = 0; i < 3; i++) {
	    i1 = (i < 2) ? i + 1 : 0;
	    b[0][i] = a[1][i] * a[2][i1] - a[2][i] * a[1][i1];
	    b[1][i] = a[2][i] * a[0][i1] - a[0][i] * a[2][i1];
	    b[2][i] = a[0][i] * a[1][i1] - a[1][i] * a[0][i1];
	    tem = b[0][i] * b[0][i] + b[1][i] * b[1][i] + b[2][i] * b[2][i];
	    if (tem > smax) {
		smax = tem;
		imax = i;
	    }
	}
	if (smax <= errnd)
	    return (0);
	smax = (float) sqrt (smax);
	for (i = 0; i < 3; i++)
	    v2[i] = (float) (b[i][imax] / smax);
	for (i = 0; i < 3; i++) {	/* Convert v2 to cartesian system */
	    v1[i] = 0.0;
	    for (j = 0; j < 3; j++)
		v1[i] += (float) drvui->b_mat[i][j] * v2[j];
	}
	vnormalize (v1);
	for (i = 0; i < 3; i++)
	    Vect[ii][i] = v1[i];	/* unit vector */
    }
    if (determinant (Vect) < 0.97) {
	for (i = 0; i < 3; i++) {
	    v2[i] = (float) Vect[1][i];
	    v3[i] = (float) Vect[2][i];
	}
	vcross (v2, v3, v1);	/* Set is not orthogonal -the smallest root is probably wrong */
	for (i = 0; i < 3; i++)
	    Vect[0][i] = v1[i];
    }
    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    vect[i][j] = (float) Vect[i][j];
    return (1);
}

/* ************************************************************** */
/* ************************************************************** */

void
find_all_in_box (int i)

/* routine to locate all atoms of type 'i' in the display box */
/* i - atom number */
{
    int n, j, k, l;		/* loop indices */

    int sym_no;			/* symmetry number */

#ifndef WIN32
    int* itmp;

    float* ftmp;
#endif

    for (n = 1; n < NvertM; n++) {	/* loop through the Master list */
	j = (int) drvui->atom_no[n] / 1000;	/* extract atom number */
	if (j == i) {
	    drvui->orig_atom_no[nvert] = i;
	    sym_no = drvui->atom_no[n] - 1000 * j;	/* symmetry operator number */
	    for (k = 0; k < 3; k++) {
		o_vert[3 * nvert + k] = xypos[3 * n + k];
		o_vert_nm[3 * nvert + k] = xypos_nm[3 * n + k];
	    }
	    for (k = 0; k < 3; k++) {
		s_vert[3 * nvert + k] = 0.0f;
		for (l = 0; l < 3; l++) {	/* calculate cartesian coordinates */
		    s_vert[3 * nvert + k] +=
			(float) drvui->b_mat[k][l] * (o_vert[3 * nvert + l] - origin[l]);
		}
	    }
	    vert_sym_nos[nvert] = drvui->atom_so[n];
	    vert_sym_no[nvert++] = sym_no;
	    if (nvert > drvui->verts_alloc) {
		Error_Box
		    ("Too many vertices for dimensions - please report to developers.");
		return;
	    }
	    if (nvert > 2 * NvertM) {
#ifdef WIN32
		Error_Box ("Overrun of vert_sym_no. Please send 'str' file\n"
			   "to Larry.Finger@@lwfinger.net.");
		exit (0);
#else
		itmp = (int *) realloc (vert_sym_no, (nvert + 2) * sizeof (int));
		if (itmp) 
		    vert_sym_no = itmp;
		else
		    Error_Box ("Unable to expand storage for vert_sym_no. (Out of memory)\n");
		itmp = (int *) realloc (vert_sym_nos, (nvert + 2) * sizeof (int));
		if (itmp) 
		    vert_sym_nos = itmp;
		else
		    Error_Box ("Unable to expand storage for vert_sym_nos. (Out of memory)\n");
		itmp = (int *) realloc (drvui->orig_atom_no, (nvert + 2) * sizeof (int));
		if (itmp) 
		    drvui->orig_atom_no = itmp;
		else
		    Error_Box ("Unable to expand storage for orig_atom_no. (Out of memory)\n");
		ftmp = (float *) realloc (o_vert, (6 * (nvert + 2) * sizeof (float)));
		if (ftmp) 
		    o_vert = ftmp;
		else
		    Error_Box ("Unable to expand storage for o_vert. (Out of memory)\n");
		ftmp = (float *) realloc (o_vert_nm, (6 * (nvert + 2) * sizeof (float)));
		if (ftmp) 
		    o_vert_nm = ftmp;
		else
		    Error_Box ("Unable to expand storage for o_vert_nm. (Out of memory)\n");
#endif
	    }
	}
    }
}				/* end of find_all_in_box */

/* ************************************************************** */
/* ************************************************************** */

void
generate_arrows (void)
// routine to generate arrow objects for magnetic moments
{
    int Arrow_Count;		// Counter for number output

    int i, k, l;		// loop counters

    float glr, glb, glg;

    GLUquadricObj *glu_quadObj;

    float base;

    float x[3], xp[3], yp[3];

    int xi, xj, xk;

    float phi, chi;

    int savedNvertM = NvertM;	// save the end of the vertex list

    char col_arrow_v[40];

    char col_arrow_p[40];

    Arrow_Count = 0;

    glu_quadObj = gluNewQuadric ();

    for (i = 0; i < drvui->nmag; ++i) {	// loop through magnetic arrows
	if (drvui->arrows[i].arrow_fn != drvui->frame_no)
	    continue;
	for (k = 0; k < 3; k++) {
	    x[k] = drvui->arrows[i].mag_xp[k] - drvui->xyzoff[k];
	}
	for (xi = -3; xi <= 3; xi++) {
	    xp[0] = (float) xi;
	    for (xj = -3; xj <= 3; xj++) {
		xp[1] = (float) xj;
		for (xk = -3; xk <= 3; xk++) {
		    xp[2] = (float) xk;	// xp is a magnetic cell translation
		    for (k = 0; k < 3; k++) {
			yp[k] = x[k];	// yp will be position in nuclear cell
			for (l = 0; l < 3; l++) {
			    yp[k] += drvui->mag_matrix[k][l] * xp[l];
			}
		    }
		    add_vert (yp, 1000 * (i + drvui->atom_alloc), 0, 0, 0);	// add vector to vertex list - NOT modulated??
		}
	    }
	}
    }				// end of loop through magnetic arrows
    for (k = savedNvertM; k < NvertM; k++) {
	float tmp, tmp1;

	i = drvui->atom_no[k] / 1000 - drvui->atom_alloc;	// get number of magnetic item
	strcpy (col_arrow_v, drvui->arrows[i].col_arrow);
	strcpy (col_arrow_p, col_arrow_v);
	Transform_VRML_Color (col_arrow_v);
	Transform_POV_Color (col_arrow_p);
	for (l = 0; l < 3; l++) {
	    x[l] = drvui->arrows[i].mag_xc[l];
	}
	double t = sqrt (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);

	if (t < 0.002)
	    continue;		// skip this arrow if 0 length
	for (l = 0; l < 3; l++) {
	    x[l] /= (float) t;
	}
	t = sqrt (x[0] * x[0] + x[1] * x[1]);
	chi = (float) atan2 (t, (double) x[2]);
	if (t > 0.001)
	    phi = (float) atan2 (x[0] / t, x[1] / t);
	else
	    phi = 0.0f;
	glLoadName (200000 + i);
	glPushName (k);
	glPushMatrix ();
	(void) sscanf (col_arrow_v, "%f %f %f", &glr, &glg, &glb);
	glColor3f (glr, glg, glb);
	glTranslatef (s_vert[3 * k], s_vert[3 * k + 1], s_vert[3 * k + 2]);
	glRotatef (-phi * (float) RAD, 0.0, 0.0, 1.0f);
	glRotatef (-chi * (float) RAD, 1.0f, 0., 0.);
	glTranslatef (0.0, 0.0, -drvui->arrows[i].arrow_length / 2.0f);
	gluQuadricDrawStyle (glu_quadObj, GLU_FILL);
	gluCylinder (glu_quadObj, drvui->arrows[i].arrow_diam,
		     drvui->arrows[i].arrow_diam, 0.9 * drvui->arrows[i].arrow_length, 10,
		     1);
	glPopMatrix ();
	glPushMatrix ();
	glTranslatef (s_vert[3 * k], s_vert[3 * k + 1], s_vert[3 * k + 2]);
	glRotatef (-phi * (float) RAD, 0.0, 0.0, 1.0f);
	glRotatef (-chi * (float) RAD, 1.0f, 0., 0.);
	glTranslatef (0.0f, 0.0f, drvui->arrows[i].arrow_length * 0.27f);

	gluQuadricDrawStyle (glu_quadObj, GLU_FILL);
	gluCylinder (glu_quadObj, drvui->arrows[i].arrow_diam * 2.0,
		     0.01, 0.4 * drvui->arrows[i].arrow_length, 10, 1);
	glPopMatrix ();
	glPopName ();

// arrow in POV

	if (doPOV) {
	    tmp = 0.5f * drvui->arrows[i].arrow_length;
	    base = 2.0f * drvui->arrows[i].arrow_diam;
	    fprintf (drvui->fpoutp, " object{\n");
	    fprintf (drvui->fpoutp, "   union{\n");
	    fprintf (drvui->fpoutp, "     object{\n");
	    fprintf (drvui->fpoutp,
		     "      cylinder{ < 0,0,%8.5f > , < 0,0,%8.5f >,%8.5f\n",
		     -0.5f * drvui->arrows[i].arrow_length,
		     0.27f * drvui->arrows[i].arrow_length, drvui->arrows[i].arrow_diam);
	    fprintf (drvui->fpoutp, "         texture{pigment{color %s}}\n", col_arrow_p);
	    fprintf (drvui->fpoutp,
		     "     finish{phong %5.2f phong_size %5.2f}}\n     }\n",
		     drvui->Phong_Value, drvui->Phong_Size);
	    tmp1 = 0.27f * drvui->arrows[i].arrow_length;
	    tmp = 0.77f * drvui->arrows[i].arrow_length;
	    fprintf (drvui->fpoutp, "     object{\n");
	    fprintf (drvui->fpoutp,
		     "      cone{<0, 0, %8.5f>, %8.5f, <0, 0, %8.5f>,  0 \n", tmp1, base,
		     tmp);
	    fprintf (drvui->fpoutp, "         texture{pigment{color %s}}\n", col_arrow_p);
	    fprintf (drvui->fpoutp, "    finish{phong %5.2f phong_size %5.2f}}\n     }\n",
		     drvui->Phong_Value, drvui->Phong_Size);
	    fprintf (drvui->fpoutp, "   rotate <%8.5f,0,0>\n", -chi * RAD);
	    fprintf (drvui->fpoutp, "   rotate <0,0,%8.5f>\n", -phi * RAD);
	    fprintf (drvui->fpoutp, "   translate <%8.5f, %8.5f, %8.5f>\n }}\n",
		     s_vert[3 * k], s_vert[3 * k + 1], s_vert[3 * k + 2]);
	}
	if (doAsy) {
	    fprintf (drvui->fpouta, "draw(pic, (%5.2f,%5.2f,%5.2f)--(%5.2f,%5.2f,%5.2f),\n",
		     s_vert[3 * k]-0.8f*drvui->arrows[i].mag_xc[0], 
		     s_vert[3 * k + 1]-0.8f*drvui->arrows[i].mag_xc[1], 
		     s_vert[3 * k + 2]-0.8f*drvui->arrows[i].mag_xc[2],
		     s_vert[3 * k]+1.2f*drvui->arrows[i].mag_xc[0], 
		     s_vert[3 * k + 1]+1.2f*drvui->arrows[i].mag_xc[1], 
		     s_vert[3 * k + 2]+1.2f*drvui->arrows[i].mag_xc[2]);
	    fprintf (drvui->fpouta, "\tlinewidth(%.2f)+rgb(%.2f,%.2f,%.2f),EndArrow3(DefaultHead3,%.2f),currentlight);\n",
		     20.f*drvui->arrows[i].arrow_diam,glr,glg,glb,4.f*20.f*drvui->arrows[i].arrow_diam); 
	}
	// arrow in VRML
	float d1, d2, d3, alpha, temp;

	float cosalp;

	d1 = x[1];
	d2 = -x[0];
	d3 = 0.0f;
	temp = (float) sqrt (d1 * d1 + d2 * d2 + d3 * d3);
	if (temp < 0.001f) {
	    d1 = 1.0;
	    temp = 1.0;
	}
	d1 /= temp;
	d2 /= temp;
	d3 /= temp;
	cosalp = x[2] / (float) sqrt (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);;
	alpha = -(float) acos (cosalp);
	if (doVrml) {
	    if (Vrml2) {
// cylinder for shaft of arrow
		fprintf (drvui->fpoutv, "  Transform{\n");
		fprintf (drvui->fpoutv, "    rotation %6.3f %6.3f %6.3f %7.3f\n",
			 d1, d2, d3, alpha);
		fprintf (drvui->fpoutv, "    translation %5.3f %5.3f %5.3f\n",
			 s_vert[3 * k], s_vert[3 * k + 1], s_vert[3 * k + 2]);
		fprintf (drvui->fpoutv, "    children [\n");
		fprintf (drvui->fpoutv, "      Transform{\n");
		fprintf (drvui->fpoutv, "        rotation 1 0 0 1.5708\n");
		fprintf (drvui->fpoutv, "        translation 0 0 %5.3f\n",
			 -0.05f * drvui->arrows[i].arrow_length);
		fprintf (drvui->fpoutv, "        children [\n");
		fprintf (drvui->fpoutv, "          Shape {\n");
		fprintf (drvui->fpoutv, "            geometry Cylinder {\n");
		fprintf (drvui->fpoutv, "              radius %f\n",
			 drvui->arrows[i].arrow_diam);
		fprintf (drvui->fpoutv, "              height %f\n",
			 0.9f * drvui->arrows[i].arrow_length);
		fprintf (drvui->fpoutv, "            }\n");
		fprintf (drvui->fpoutv,
			 "            appearance Appearance { material Material "
			 "{ diffuseColor %s}}\n", col_arrow_v);
		fprintf (drvui->fpoutv, "          }\n");
		fprintf (drvui->fpoutv, "        ]\n");
		fprintf (drvui->fpoutv, "      }\n");
		fprintf (drvui->fpoutv, "    ]\n");
		fprintf (drvui->fpoutv, "  }\n");
// Cone at end of arrow
		fprintf (drvui->fpoutv, "  Transform{\n");
		fprintf (drvui->fpoutv, "    rotation %6.3f %6.3f %6.3f %7.3f\n",
			 d1, d2, d3, alpha);
		fprintf (drvui->fpoutv, "    translation %5.3f %5.3f %5.3f\n",
			 s_vert[3 * k], s_vert[3 * k + 1], s_vert[3 * k + 2]);
		fprintf (drvui->fpoutv, "    children Transform {\n");
		fprintf (drvui->fpoutv, "      children [\n");
		fprintf (drvui->fpoutv, "        Transform{\n");
		fprintf (drvui->fpoutv, "          rotation 1 0 0 1.5708\n");
		fprintf (drvui->fpoutv, "          translation 0 0 %5.3f\n",
			 0.5f * drvui->arrows[i].arrow_length);
		fprintf (drvui->fpoutv, "          children [\n");
		fprintf (drvui->fpoutv, "            Shape {\n");
		fprintf (drvui->fpoutv, "              geometry Cone {\n");
		fprintf (drvui->fpoutv, "                bottomRadius %f\n",
			 2.0 * drvui->arrows[i].arrow_diam);
		fprintf (drvui->fpoutv, "                height %f\n",
			 drvui->arrows[i].arrow_length * 0.4);
		fprintf (drvui->fpoutv, "              }\n");
		fprintf (drvui->fpoutv,
			 "              appearance Appearance { material Material "
			 "{   diffuseColor %s}}\n", col_arrow_v);
		fprintf (drvui->fpoutv, "            }\n");
		fprintf (drvui->fpoutv, "          ]\n");
		fprintf (drvui->fpoutv, "        }\n");
		fprintf (drvui->fpoutv, "      ]\n");
		fprintf (drvui->fpoutv, "    }\n");
		fprintf (drvui->fpoutv, "  }\n");
	    } else {		// VRML 1
		fprintf (drvui->fpoutv, " Separator{\n");
		fprintf (drvui->fpoutv, "  Transform{\n");
		fprintf (drvui->fpoutv, "   translation %5.3f %5.3f %5.3f\n",
			 s_vert[3 * k], s_vert[3 * k + 1], s_vert[3 * k + 2]);
		fprintf (drvui->fpoutv, "   rotation %6.3f %6.3f %6.3f %7.3f\n  }\n",
			 d1, d2, d3, alpha);
		fprintf (drvui->fpoutv, "  Material {diffuseColor %s}\n", col_arrow_v);
		fprintf (drvui->fpoutv, "  Rotation{rotation 1 0 0 1.5708}\n");
		fprintf (drvui->fpoutv, "  Cylinder {\n");
		fprintf (drvui->fpoutv, "   parts SIDES\n");
		fprintf (drvui->fpoutv, "   radius %f\n", drvui->arrows[i].arrow_diam);
		fprintf (drvui->fpoutv, "   height %f\n",
			 0.9f * drvui->arrows[i].arrow_length);
		fprintf (drvui->fpoutv, "  }\n }\n");
		fprintf (drvui->fpoutv, " Separator{\n");
		fprintf (drvui->fpoutv, "  Transform{\n");
		fprintf (drvui->fpoutv, "   translation %5.3f %5.3f %5.3f\n",
			 s_vert[3 * k], s_vert[3 * k + 1], s_vert[3 * k + 2]);
		fprintf (drvui->fpoutv, "   rotation %6.3f %6.3f %6.3f %7.3f\n  }\n",
			 d1, d2, d3, alpha);
		fprintf (drvui->fpoutv, "  Material {diffuseColor %s}\n", col_arrow_v);
		fprintf (drvui->fpoutv, "  Transform{\n");
		fprintf (drvui->fpoutv, "   rotation 1 0 0 1.5708\n");
		fprintf (drvui->fpoutv, "   translation 0 0 %5.3f\n  }\n",
			 0.5f * drvui->arrows[i].arrow_length);
		fprintf (drvui->fpoutv, "  Cone {\n");
		fprintf (drvui->fpoutv, "   parts SIDES\n");
		fprintf (drvui->fpoutv, "   bottomRadius %f\n",
			 2.0 * drvui->arrows[i].arrow_diam);
		fprintf (drvui->fpoutv, "   height %f\n",
			 drvui->arrows[i].arrow_length * 0.4);
		fprintf (drvui->fpoutv, "  }\n }\n");
	    }
	}
	Arrow_Count++;
    }

    gluDeleteQuadric (glu_quadObj);

    if (Arrow_Count > 0) {
	fprintf (drvui->fcns, "%4d magnetic moment arrows.\n", Arrow_Count);
	fprintf (drvui->flout, "Generated %4d arrows.\n", Arrow_Count);
    }
    NvertM = savedNvertM;	// restore vertex count
}				// end of generate_arrows

/* ************************************************************** */
/* ************************************************************** */

void
generate_bonds (void)

/* routine to generate bond descriptions and add then to edit list */
{
    int start, same_type;

    int i, j, k;		/* loop counters */

    int N_Bond;			/* Counter for number of bonds of a type */

    int l1;			/* loop counter for matrix operations */

    int nvert1;			/* vertex counter */

    int Bond_Count;		/* Counter for bonds output */

    char bnd1[4], bnd2[4];	/* temporary storage of bond labels */

    float d;			/* bond distance */

    float df[3];		/* vector difference */

    float beta;			/* rotation angle */

    float gamma;		/* another rotation angle */

    float at1[3], at2[3];	/*initial coordinates of bond endpoints */

    float clip0, clip1;		/* clipping scale factors */

    int m, l2, l3, out = 0;
    int ell[2] = { 0, 0 };	/* save ellipsoid numbers */

    int is_ellipsoid[2] = { 0, 0 };	/* ellipsoid markers - for moving bond endpoint */
    int o1, p, p1, q, r;	/* loop variables for rotation matrix calculation */

    float elrot0[3][3];		/* rotation matrices of ellipsoids */

    float Z0[3], Z1[3];		/* components of rms displacement along bond */

    int *ellips_id;

    int *ellips_num;

    int jj;

    float factor, rms0, rms1;

    float glr, glg, glb;

    int numdashes;

    float dashes = 5.0f;

    float mf;

    GLUquadricObj *glu_quadObj;

    int omit;

    char col_bond_p[40];

    char col_bond_v[40];

    if (drvui->nbond == 1)
	return;			/* exit if no work to be done */

    Bond_Count = 0;
    if (!
	(ellips_id =
	 (int *) zalloc ((unsigned) ((long) (drvui->verts_alloc * sizeof (int)))))) {
	Error_Box ("Unable to get ellips_id allocation!");
	return;
    }
    if (!
	(ellips_num =
	 (int *) zalloc ((unsigned) ((long) (drvui->verts_alloc * sizeof (int)))))) {
	Error_Box ("Unable to get ellips_num allocation!");
	free (ellips_id);
	return;
    }

    glu_quadObj = gluNewQuadric ();

    for (i = 1; i < drvui->nbond; ++i) {	/* loop through bond types */
	if (drvui->bonds[i].bond_fn != drvui->frame_no)
	    continue;		// skip if not in this frame
	strcpy (col_bond_p, drvui->bonds[i].col_bond);
	strcpy (col_bond_v, col_bond_p);
	Transform_VRML_Color (col_bond_v);
	Transform_POV_Color (col_bond_p);
	nvert = 0;		/* empty the vertex list */
	for (j = 0; j < natom; j++) {	/* loop through atom types */
	    if (drvui->atoms[j].atom_fn != drvui->frame_no)
		continue;
	    int o = drvui->atoms[j].atom_n;

	    if (((o & 255) > 0) || (((o >> 24) & 255) > 0)) {	// sphere or ellipsoid

		o = (drvui->atoms[j].atom_n >> 24) & 255;
		if (o > 0) {	// ellipsoid
		    if (check_atom_name
			(drvui->bonds[i].bond_l1, drvui->ellips[o].ellips_l)) {
			if (drvui->do_ellipsoids == 1 && drvui->El_Cutout == 1 && drvui->ellips[o].ell_type > 100) {	/* skip if B iso */
			    is_ellipsoid[0] = 1;	/* or no cutout */
			    ell[0] = o;
			}
			o = nvert;
			find_all_in_box (j);	// get all atoms of type j
			for (p = o; p < nvert; p++)
			    ellips_id[p] = ell[0];	// save ellipsoid id for these entries
			ellips_num[p] = j;
		    }
		} else {	// sphere
		    o = drvui->atoms[j].atom_n & 255;	// number of sphere
		    if (check_atom_name
			(drvui->bonds[i].bond_l1, drvui->spheres[o].sphere_l)) {
			find_all_in_box (j);	/* get all atoms of type j */
		        is_ellipsoid[0] = 0;
		    }
		}		/* o > 0 */
	    }
	}			/* for (j=0;j<natom;++j) */

	bnd1[0] = drvui->bonds[i].bond_l1[0];
	bnd1[1] = drvui->bonds[i].bond_l1[1];
	bnd1[2] = drvui->bonds[i].bond_l1[2];
	bnd1[3] = drvui->bonds[i].bond_l1[3];
	nvert1 = nvert;
	if (check_atom_name (drvui->bonds[i].bond_l1, drvui->bonds[i].bond_l2)) {
	    bnd2[0] = drvui->bonds[i].bond_l1[0];	/* atoms are same type */
	    bnd2[1] = drvui->bonds[i].bond_l1[1];
	    bnd2[2] = drvui->bonds[i].bond_l1[2];
	    bnd2[3] = drvui->bonds[i].bond_l1[3];
	    is_ellipsoid[1] = is_ellipsoid[0];
	    ell[1] = ell[0];
	    same_type = 1;
	} else {
	    bnd2[0] = drvui->bonds[i].bond_l2[0];
	    bnd2[1] = drvui->bonds[i].bond_l2[1];
	    bnd2[2] = drvui->bonds[i].bond_l2[2];
	    bnd2[3] = drvui->bonds[i].bond_l2[3];
	    same_type = 0;
/* atoms i and j are not the same type */
	    for (j = 0; j < natom; ++j) {	/* loop through atom types */
		if (drvui->atoms[j].atom_fn != drvui->frame_no)
		    continue;
		if (((drvui->atoms[j].atom_n & 255) > 0) || (((drvui->atoms[j].atom_n >> 24) & 255) > 0)) {	// sphere or ellipsoid
		    int o;

		    o = (drvui->atoms[j].atom_n >> 24) & 255;
		    if (o > 0) {	// ellipsoid
			if (check_atom_name
			    (drvui->bonds[i].bond_l2, drvui->ellips[o].ellips_l)) {
			    if (drvui->do_ellipsoids == 1 && drvui->El_Cutout == 1 && drvui->ellips[o].ell_type > 100) {	/* skip if B iso */
				is_ellipsoid[1] = 1;	/* or no cutout */
				ell[1] = o;
			    }
			    o = nvert;
			    find_all_in_box (j);	/* get all atoms of type j */
			    for (p = o; p < nvert; p++)
				ellips_id[p] = ell[1];	/* save ellipsoid type */
			    ellips_num[p] = j;
			}
		    } else {
			o = drvui->atoms[j].atom_n & 255;	// sphere number
			if (check_atom_name
			    (drvui->bonds[i].bond_l2, drvui->spheres[o].sphere_l)) {
			    find_all_in_box (j);	/* get all atoms of type j */
			    is_ellipsoid[1] = 0;
			}
		    }		// o > 0
		}
	    }			/* for (j=0;j<natom;++j) */
	}			/* if atoms are same type */

	if (((nvert1 > 0) && (nvert > nvert1)) || (same_type == 1)) {
	    N_Bond = 0;
	    for (j = 0; j < nvert1; ++j) {	/* loop through atoms of type 1 */
		start = nvert1;
		if (same_type == 1)
		    start = j;
		for (k = start; k < nvert; ++k) {	/* loop on type 2 atoms  */
		    d = dist (j, k);	/* get atomic distance */
		    if ((d >= drvui->bonds[i].bond_min) && (d <= drvui->bonds[i].bond_max)) {	/* here if distance within range */
			df[0] = df[1] = drvui->Bond_Mult * drvui->bonds[i].bond_size;	/* scaling */
			df[2] = 0.5f * d;

			Bond_Count++;	/* update count */
			for (p = 0; p < 3; p++) {	/* initialize endpoints of bond */
			    Z0[p] = s_vert[3 * j + p];
			    Z1[p] = s_vert[3 * k + p];
			}

			if (drvui->do_ellipsoids == 1 && drvui->El_Cutout == 1) {

/* Elipsoid has octant cutout, calculate rms displacement along bond  */

			    /* calculate multiplier - the 0.9 forces the end of the bond closer to the ellipsoid */
			    factor =
				(float) (drvui->Ellipsoid_Scale * 0.9 /
					 (PI * d * 1.41421356));

/* Mean-square displacement of a scaled thermal ellipsoid along the direction of a bond is given by: */
/*                                                                                      */
/*      D = factor^2 * (X'G'RBR'GX)                                                     */
/*                                                                                      */
/*  where X is the interatomic vector given in triclinic fractional coordinates,        */
/*        G is the lattice metric with the ij'th element given by a[i].a[j],            */
/*        R is the rotational part of the symmetry operation that transform atom i,     */
/*        B is a square matrix of the betaij thermal factors,                           */
/*    and a prime indicates the transpose of the matrix.                                */
/*                                                                                      */
			    if (is_ellipsoid[0] == 1) {	/* at one end of bond */
				float B[3][3], BP[3][3];

				o1 = vert_sym_no[j];
				p1 = ellips_id[j];	/* have original ellipsoid number */
				jj = ellips_num[j];
				if (drvui->ellips[p1].ellips_ismod != 0) {	/* apply modulation */
				    float uij[6];

				    float vert[3];

				    vert[0] = o_vert_nm[3 * j];
				    vert[1] = o_vert_nm[3 * j + 1];
				    vert[2] = o_vert_nm[3 * j + 2];
				    if (modulate_uij (vert, p1, jj, o1, uij) == 1) {	/* if n.p.d, keep average value */
					for (p = 0; p < 3; p++)
					    B[p][p] = drvui->ellips[p1].ellips[p];	// extract betaij's into square matrix
					B[0][1] = B[1][0] = drvui->ellips[p1].ellips[3];
					B[0][2] = B[2][0] = drvui->ellips[p1].ellips[4];
					B[1][2] = B[2][1] = drvui->ellips[p1].ellips[5];
				    } else {	/* use modulated values */
					for (p = 0; p < 3; p++)
					    B[p][p] = uij[p];	// extract betaij's into square matrix
					B[0][1] = B[1][0] = uij[3];
					B[0][2] = B[2][0] = uij[4];
					B[1][2] = B[2][1] = uij[5];
				    }
				} else {	/* not modulated */
				    for (p = 0; p < 3; p++)
					B[p][p] = drvui->ellips[p1].ellips[p];	// extract betaij's into square matrix
				    B[0][1] = B[1][0] = drvui->ellips[p1].ellips[3];
				    B[0][2] = B[2][0] = drvui->ellips[p1].ellips[4];
				    B[1][2] = B[2][1] = drvui->ellips[p1].ellips[5];
				}
				for (p = 0; p < 3; p++) {
				    df[p] = 0.0;
				    for (q = 0; q < 3; ++q) {
					elrot0[p][q] = 0.0;
					df[p] +=
					    drvui->ginv[q][p] * (o_vert[3 * k + q] -
								 o_vert[3 * j + q]);
					for (r = 0; r < 3; r++)
					    elrot0[p][q] +=
						B[p][r] * drvui->rssC[o1][q][r];
				    }
				}
				for (p = 0; p < 3; p++) {
				    for (q = 0; q < 3; q++) {
					BP[p][q] = 0.0;
					for (r = 0; r < 3; r++)
					    BP[p][q] += (float) drvui->rssC[o1][p][r] * elrot0[r][q];	/* BP is rotated betij matrix */
				    }
				}
				for (p = 0; p < 3; p++) {
				    Z0[p] = 0.0;
				    for (q = 0; q < 3; q++)
					Z0[p] += df[q] * BP[q][p];
				}
/* rms0 is the scaled rms displacement value of the ellipsoid in the direction of the bond */
				rms0 = factor * (float) sqrt (Z0[0] * df[0] +
							      Z0[1] * df[1] +
							      Z0[2] * df[2]);
/* calculate components of end point at intersection of ellipsoid and bond */
				for (q = 0; q < 3; q++)
				    Z0[q] =
					s_vert[3 * j + q] + (s_vert[3 * k + q] -
							     s_vert[3 * j +
								    q]) * rms0 / d;
			    }
			    /* if is_ellipsoid[0] */
			    if (is_ellipsoid[1] == 1) {	/* ellipsoid at other end of bond */
				float B[3][3], BP[3][3];

				o1 = vert_sym_no[k];
				p1 = ellips_id[k];	/* have original ellipsoid number */
				jj = ellips_num[k];
				if (drvui->ellips[p1].ellips_ismod != 0) {	/* apply modulation if needed */
				    float uij[6];

				    float vert[3];

				    vert[0] = o_vert_nm[3 * k];
				    vert[1] = o_vert_nm[3 * k + 1];
				    vert[2] = o_vert_nm[3 * k + 2];
				    if (modulate_uij (vert, p1, jj, o1, uij) == 1) {	/* if n.p.d, use average value */
					for (p = 0; p < 3; p++)
					    B[p][p] = drvui->ellips[p1].ellips[p];	/* extract betaij's into square matrix */
					B[0][1] = B[1][0] = drvui->ellips[p1].ellips[3];
					B[0][2] = B[2][0] = drvui->ellips[p1].ellips[4];
					B[1][2] = B[2][1] = drvui->ellips[p1].ellips[5];
				    } else {
					for (p = 0; p < 3; p++)
					    B[p][p] = uij[p];	// extract betaij's into square matrix
					B[0][1] = B[1][0] = uij[3];
					B[0][2] = B[2][0] = uij[4];
					B[1][2] = B[2][1] = uij[5];
				    }
				} else {	/* not modulated */
				    for (p = 0; p < 3; p++)
					B[p][p] = drvui->ellips[p1].ellips[p];	/* extract betaij's into square matrix */
				    B[0][1] = B[1][0] = drvui->ellips[p1].ellips[3];
				    B[0][2] = B[2][0] = drvui->ellips[p1].ellips[4];
				    B[1][2] = B[2][1] = drvui->ellips[p1].ellips[5];
				}
				for (p = 0; p < 3; p++) {
				    df[p] = 0.0;
				    for (q = 0; q < 3; q++) {
					elrot0[p][q] = 0.0;
					df[p] +=
					    drvui->ginv[q][p] * (o_vert[3 * k + q] -
								 o_vert[3 * j + q]);
					for (r = 0; r < 3; ++r)
					    elrot0[p][q] +=
						B[p][r] * drvui->rssC[o1][q][r];
				    }
				}
				for (p = 0; p < 3; p++) {
				    for (q = 0; q < 3; q++) {
					BP[p][q] = 0.0;
					for (r = 0; r < 3; r++)
					    BP[p][q] += (float) drvui->rssC[o1][p][r] * elrot0[r][q];	/* BP is rotated betij matrix */
				    }
				}
				for (p = 0; p < 3; p++) {
				    Z1[p] = 0.0;
				    for (q = 0; q < 3; q++)
					Z1[p] += df[q] * BP[q][p];
				}
/* rms1 is the scaled rms displacement value of the ellipsoid in the direction of the bond */

				rms1 = factor * (float) sqrt (Z1[0] * df[0] +
							      Z1[1] * df[1] +
							      Z1[2] * df[2]);

/* calculate components of end point at intersection of ellipsoid and bond */
				for (q = 0; q < 3; q++)
				    Z1[q] =
					s_vert[3 * k + q] - (s_vert[3 * k + q] -
							     s_vert[3 * j +
								    q]) * rms1 / d;
			    }	/* if is_ellipsoid[1] */
			}
			/* if any segmented ellipsoids present */
			for (l1 = 0; l1 < 3; l1++) {
			    at1[l1] = (float) Z0[l1];
			    at2[l1] = (float) Z1[l1];
			    df[l1] = (float) (Z1[l1] - Z0[l1]);
			}

			if (clipflag != 0) {
			    out = 0;
			    l2 = 0;
			    l3 = 0;
			    clip0 = 0.0f;
			    clip1 = 1.0f;
			    for (l1 = 0; l1 < 3; l1++) {
				if (o_vert[3 * j + l1] <
				    drvui->frames[drvui->frame_no].clip_lim[l1] - 0.01
				    || o_vert[3 * j + l1] >
				    drvui->frames[drvui->frame_no].clip_lim[l1 + 3] +
				    0.01)
				    l2 = 1;
				if (o_vert[3 * k + l1] <
				    drvui->frames[drvui->frame_no].clip_lim[l1] - 0.01
				    || o_vert[3 * k + l1] >
				    drvui->frames[drvui->frame_no].clip_lim[l1 + 3] +
				    0.01)
				    l3 = 1;
			    }
			    if (l2 && l3)
				out = 1;	/* both atoms outside */
			    if (l3)
				clip1 = 0.5f;
			    if (l2)
				clip0 = 0.5f;

			    at2[0] = at1[0] + clip1 * df[0];
			    at2[1] = at1[1] + clip1 * df[1];
			    at2[2] = at1[2] + clip1 * df[2];

			    at1[0] = at1[0] + clip0 * df[0];
			    at1[1] = at1[1] + clip0 * df[1];
			    at1[2] = at1[2] + clip0 * df[2];
			}
			/* if clipflag */
			omit = 0;
			for (l1 = 0; l1 < Omit->nomits; l1++) {
			    if (Omit->omit1[l1] == i * 1000
				&& Omit->omit2[l1] == Bond_Count)
				omit = 1;
			}

			if ((clipflag == 0 || out == 0) && omit == 0) {
			    if (doPOV) {
				fprintf (drvui->fpoutp,
					 " /* Bond: %c%c %8.5f %8.5f %8.5f TO %c%c %8.5f %8.5f %8.5f*/\n",
					 bnd1[0], bnd1[1], o_vert[3 * j],
					 o_vert[3 * j + 1], o_vert[3 * j + 2], bnd2[0],
					 bnd2[1], o_vert[3 * k], o_vert[3 * k + 1],
					 o_vert[3 * k + 2]);
				if (drvui->bonds[i].bond_style == 0) {	/* solid bonds */
				    fprintf (drvui->fpoutp,
					     " cylinder{<%8.5f, %8.5f, %8.5f>, <%8.5f, %8.5f, %8.5f>, %8.5f \n",
					     at1[0], at1[1], at1[2], at2[0], at2[1],
					     at2[2],
					     drvui->Bond_Mult *
					     drvui->bonds[i].bond_size);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n",
					     col_bond_p);
				    fprintf (drvui->fpoutp, " }\n");
				} else {	/* dashed bonds */
				    numdashes = drvui->bonds[i].bond_style;
				    dashes = (float) numdashes;

				    for (m = 0; m < 3; m++)
					df[m] = (float) (at2[m] - at1[m]);
				    fprintf (drvui->fpoutp,
					     " cylinder{<%8.5f, %8.5f, %8.5f>, <%8.5f, %8.5f, %8.5f>, %8.5f \n",
					     at1[0], at1[1], at1[2], at1[0] + df[0] / dashes,
					     at1[1] + df[1] / dashes, at1[2] + df[2] / dashes,
					     drvui->Bond_Mult *
					     drvui->bonds[i].bond_size);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n",
					     col_bond_p);
				    fprintf (drvui->fpoutp, " }\n");
				    for (m = 2; m < numdashes - 1; m += 2) {
					mf = (float) m;
					fprintf (drvui->fpoutp,
						 " cylinder{<%8.5f, %8.5f, %8.5f>, <%8.5f, %8.5f, %8.5f>, %8.5f \n",
						 at1[0] + df[0] * mf / dashes,
						 at1[1] + df[1] * mf / dashes,
						 at1[2] + df[2] * mf / dashes,
						 at1[0] + df[0] * (mf + 1.) / dashes,
						 at1[1] + df[1] * (mf + 1.) / dashes,
						 at1[2] + df[2] * (mf + 1.) / dashes,
						 drvui->Bond_Mult *
						 drvui->bonds[i].bond_size);
					fprintf (drvui->fpoutp,
						 "  texture{pigment{color %s  }}\n",
						 col_bond_p);
					fprintf (drvui->fpoutp, " }\n");
				    }
				    fprintf (drvui->fpoutp,
					     " cylinder{<%8.5f, %8.5f, %8.5f>, <%8.5f, %8.5f, %8.5f>, %8.5f \n",
					     at1[0] + df[0] * (dashes - 1.) / dashes,
					     at1[1] + df[1] * (dashes - 1.) / dashes,
					     at1[2] + df[2] * (dashes - 1.) / dashes, 
					     at2[0], at2[1], at2[2],
					     drvui->Bond_Mult *
					     drvui->bonds[i].bond_size);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n",
					     col_bond_p);
				    fprintf (drvui->fpoutp, " }\n");
				}
			    }
			    if (doAsy) {
				(void) sscanf (col_bond_v, "%f %f %f", &glr, &glg, &glb);
				fprintf (drvui->fpouta,
					" // Bond: %c%c %8.5f %8.5f %8.5f TO %c%c %8.5f %8.5f %8.5f\n",
					bnd1[0], bnd1[1], o_vert[3 * j],
					o_vert[3 * j + 1], o_vert[3 * j + 2], bnd2[0],
					bnd2[1], o_vert[3 * k], o_vert[3 * k + 1],
					o_vert[3 * k + 2]);
				if (drvui->bonds[i].bond_style == 0) {	/* solid bonds */
				    fprintf (drvui->fpouta," draw(pic, (%8.5f,%8.5f,%8.5f)--(%8.5f,%8.5f,%8.5f),",
					    at1[0], at1[1], at1[2], at2[0], at2[1], at2[2]);
				    fprintf (drvui->fpouta, "rgb(%4.2f,%4.2f,%4.2f)+linewidth(%5.2f) );\n",
					    glr,glg,glb, 25. * drvui->Bond_Mult *
					    drvui->bonds[i].bond_size);
				} else {	/* dashed bonds */
				    numdashes = drvui->bonds[i].bond_style;
				    dashes = (float) numdashes;

				    for (m = 0; m < 3; m++)
					df[m] = (float) (at2[m] - at1[m]);
				    fprintf (drvui->fpouta," draw(pic, (%8.5f,%8.5f,%8.5f)--(%8.5f,%8.5f,%8.5f),",
					    at1[0], at1[1], at1[2], at1[0] + df[0] / dashes, 
					    at1[1] + df[1] / dashes, at1[2] + df[2] / dashes);
				    fprintf (drvui->fpouta, "rgb(%4.2f,%4.2f,%4.2f)+linewidth(%5.2f)+squarecap );\n",
					    glr,glg,glb, 25. * drvui->Bond_Mult *
					    drvui->bonds[i].bond_size);
				    for (m = 2; m < numdashes - 1; m += 2) {
					mf = (float) m;
					fprintf (drvui->fpouta," draw(pic, (%8.5f,%8.5f,%8.5f)--(%8.5f,%8.5f,%8.5f),",
						at1[0] + df[0] * mf / dashes,
						at1[1] + df[1] * mf / dashes,
						at1[2] + df[2] * mf / dashes,
						at1[0] + df[0] * (mf + 1.) / dashes,
						at1[1] + df[1] * (mf + 1.) / dashes,
						at1[2] + df[2] * (mf + 1.) / dashes);
					fprintf (drvui->fpouta, "rgb(%4.2f,%4.2f,%4.2f)+linewidth(%5.2f)+squarecap );\n",
						glr,glg,glb, 25. * drvui->Bond_Mult *
						drvui->bonds[i].bond_size);
				    }
				    fprintf (drvui->fpouta," draw(pic, (%8.5f,%8.5f,%8.5f)--(%8.5f,%8.5f,%8.5f),",
					     at1[0] + df[0] * (dashes - 1.) / dashes,
					     at1[1] + df[1] * (dashes - 1.) / dashes,
					     at1[2] + df[2] * (dashes - 1.) / dashes, 
					     at2[0], at2[1], at2[2]);
				    fprintf (drvui->fpouta, "rgb(%4.2f,%4.2f,%4.2f)+linewidth(%5.2f)+squarecap );\n",
					    glr,glg,glb, 25. * drvui->Bond_Mult *
					    drvui->bonds[i].bond_size);
				}
			    }
			    N_Bond++;
			    if (doVrml && no_comment == 0)
				fprintf (drvui->fpoutv,
					 "# Bond: %c%c %8.5f %8.5f %8.5f TO %c%c %8.5f %8.5f %8.5f\n",
					 bnd1[0], bnd1[1], o_vert[3 * j],
					 o_vert[3 * j + 1], o_vert[3 * j + 2], bnd2[0],
					 bnd2[1], o_vert[3 * k], o_vert[3 * k + 1],
					 o_vert[3 * k + 2]);
			    for (m = 0; m < 3; m++)
				df[m] = (float) (at2[m] - at1[m]);
			    d = (float) sqrt (df[0] * df[0] + df[1] * df[1] +
					      df[2] * df[2]);
/*
   calculate rotation matrices needed to rotate bond about Z onto Y-Z plane and then about X so that
   bond is parallel to Z.  Finally, apply inverse transformation to cylinder, which started parallel
   to X.  Resulting cylinder will be parallel to vector between atoms.
*/
			    beta = (float) atan2 (df[0], df[1] + 0.0000001f);	/* rotation angle about Z (in radians) */
			    gamma = (float) sqrt (df[0] * df[0] + df[1] * df[1]);
			    gamma = (float) atan2 (gamma, df[2]);	/* Rotation angle about X */
			    glPushMatrix ();
			    glLoadName (i * 1000);
			    glPushName (Bond_Count);
			    (void) sscanf (col_bond_v, "%f %f %f", &glr, &glg, &glb);
			    glColor3f (glr, glg, glb);
			    if (drvui->bonds[i].bond_style == 0) {
				glPushMatrix ();
				glTranslatef (at1[0], at1[1], at1[2]);
				glRotatef (-beta * (float) RAD, 0.0f, 0.0f, 1.);
				glRotatef (-gamma * (float) RAD, 1., 0.0f, 0.0f);
				gluQuadricDrawStyle (glu_quadObj, GLU_FILL);
				gluCylinder (glu_quadObj,
					     drvui->Bond_Mult * drvui->bonds[i].bond_size,
					     drvui->Bond_Mult * drvui->bonds[i].bond_size,
					     d, 10, 1);
				glPopMatrix ();
			    } else {
				numdashes = drvui->bonds[i].bond_style;
				dashes = (float) numdashes;
				glPushMatrix ();
				glTranslatef (at1[0], at1[1], at1[2]);
				glRotatef (-beta * (float) RAD, 0.0f, 0.0f, 1.);
				glRotatef (-gamma * (float) RAD, 1., 0.0f, 0.0f);
				gluQuadricDrawStyle (glu_quadObj, GLU_FILL);
				gluCylinder (glu_quadObj,
					     drvui->Bond_Mult * drvui->bonds[i].bond_size,
					     drvui->Bond_Mult * drvui->bonds[i].bond_size,
					     d / dashes, 10, 1);
				glPopMatrix ();
				for (m = 2; m < numdashes; m += 2) {
				    mf = (float) m;
				    glPushMatrix ();
				    glTranslatef (at1[0] + df[0] * mf / dashes,
						  at1[1] + df[1] * mf / dashes,
						  at1[2] + df[2] * mf / dashes);
				    glRotatef (-beta * (float) RAD, 0.0f, 0.0f, 1.);
				    glRotatef (-gamma * (float) RAD, 1., 0.0f, 0.0f);
				    gluQuadricDrawStyle (glu_quadObj, GLU_FILL);
				    gluCylinder (glu_quadObj,
						 drvui->Bond_Mult *
						 drvui->bonds[i].bond_size,
						 drvui->Bond_Mult *
						 drvui->bonds[i].bond_size, d / dashes,
						 10, 1);
				    glPopMatrix ();
				}
			    }
			    glPopName ();
			    glPopMatrix ();
			    if (doVrml) {
				if (Vrml2) {
				    if (drvui->bonds[i].bond_style == 0) {	/* solid bond */
					fprintf (drvui->fpoutv, "  Transform{\n ");
					fprintf (drvui->fpoutv,
						 "   translation %5.3f %5.3f %5.3f\n",
						 (at1[0] + at2[0]) / 2.,
						 (at1[1] + at2[1]) / 2.,
						 (at1[2] + at2[2]) / 2.);
					fprintf (drvui->fpoutv, "   rotation 0 0 1 %f\n",
						 -beta);
					fprintf (drvui->fpoutv,
						 "   children Transform { rotation 1 0 0 %f\n",
						 1.5708 - gamma);
					fprintf (drvui->fpoutv, " children [ Shape {\n");
					fprintf (drvui->fpoutv,
						 " appearance Appearance { material Material {   diffuseColor %s} }\n",
						 col_bond_v);
					fprintf (drvui->fpoutv, "  geometry Cylinder {");
					fprintf (drvui->fpoutv, "radius %f",
						 drvui->Bond_Mult *
						 drvui->bonds[i].bond_size);
					fprintf (drvui->fpoutv, " height %f", d);
					fprintf (drvui->fpoutv, "}}]}}\n");
				    } else {	/* dashed bonds */
					numdashes = drvui->bonds[i].bond_style;
					dashes = (float) numdashes;
					fprintf (drvui->fpoutv, "  Transform{\n ");
					fprintf (drvui->fpoutv,
						 "   translation %5.3f %5.3f %5.3f\n",
						 at1[0] + df[0] / dashes,
						 at1[1] + df[1] / dashes,
						 at1[2] + df[2] / dashes);
					fprintf (drvui->fpoutv, "   rotation 0 0 1 %f\n",
						 -beta);
					fprintf (drvui->fpoutv,
						 "   children Transform { rotation 1 0 0 %f\n",
						 1.5708 - gamma);
					fprintf (drvui->fpoutv, " children [ Shape {\n");
					fprintf (drvui->fpoutv,
						 " appearance Appearance { material Material {   diffuseColor %s} }\n",
						 col_bond_v);
					fprintf (drvui->fpoutv, "  geometry Cylinder {");
					fprintf (drvui->fpoutv, "radius %f",
						 drvui->Bond_Mult *
						 drvui->bonds[i].bond_size);
					fprintf (drvui->fpoutv, " height %f",
						 d / dashes);
					fprintf (drvui->fpoutv, "}}]}}\n");
					for (m = 2; m < dashes - 1; m += 2) {
					    mf = (float) m;
					    fprintf (drvui->fpoutv, "  Transform{\n ");
					    fprintf (drvui->fpoutv,
						     "   translation %5.3f %5.3f %5.3f\n",
						     at1[0] +
						     df[0] * mf / dashes,
						     at1[1] +
						     df[1] * mf / dashes,
						     at1[2] +
						     df[2] * mf / dashes);
					    fprintf (drvui->fpoutv,
						     "   rotation 0 0 1 %f\n", -beta);
					    fprintf (drvui->fpoutv,
						     "   children Transform { rotation 1 0 0 %f\n",
						     1.5708 - gamma);
					    fprintf (drvui->fpoutv,
						     " children [ Shape {\n");
					    fprintf (drvui->fpoutv,
						     " appearance Appearance { material Material {   diffuseColor %s} }\n",
						     col_bond_v);
					    fprintf (drvui->fpoutv,
						     "  geometry Cylinder {");
					    fprintf (drvui->fpoutv, "radius %f",
						     drvui->Bond_Mult *
						     drvui->bonds[i].bond_size);
					    fprintf (drvui->fpoutv, " height %f",
						     d / dashes);
					    fprintf (drvui->fpoutv, "}}]}}\n");
					}
				    }
				} else {	/* VRML1 */
				    if (drvui->bonds[i].bond_style == 0) {	/* solid bonds */
					fprintf (drvui->fpoutv, " Separator{\n");
					fprintf (drvui->fpoutv, "  Transform{\n ");
					fprintf (drvui->fpoutv,
						 "   translation %5.3f %5.3f %5.3f\n  }\n",
						 (at1[0] + at2[0]) / 2.,
						 (at1[1] + at2[1]) / 2.,
						 (at1[2] + at2[2]) / 2.);
					fprintf (drvui->fpoutv, "  Rotation {\n");
					fprintf (drvui->fpoutv,
						 "   rotation 0 0 1 %f\n  }\n", -beta);
					fprintf (drvui->fpoutv, "  Rotation {\n");
					fprintf (drvui->fpoutv,
						 "   rotation 1 0 0 %f\n  }\n",
						 1.5708 - gamma);
					fprintf (drvui->fpoutv,
						 "  Material {\n   diffuseColor %s\n",
						 col_bond_v);
					fprintf (drvui->fpoutv, "  }\n");
					fprintf (drvui->fpoutv, "  Cylinder {\n");
					fprintf (drvui->fpoutv, "   parts SIDES\n");
					fprintf (drvui->fpoutv, "   radius %f\n",
						 drvui->Bond_Mult *
						 drvui->bonds[i].bond_size);
					fprintf (drvui->fpoutv, "   height %f\n", d);
					fprintf (drvui->fpoutv, "  }\n }\n");
				    } else {	/* dashed bonds */
					fprintf (drvui->fpoutv, " Separator{\n");
					fprintf (drvui->fpoutv, "  Transform{\n ");
					fprintf (drvui->fpoutv,
						 "   translation %5.3f %5.3f %5.3f\n  }\n",
						 at1[0] + df[0] / dashes, 
						 at1[1] + df[1] / dashes,
						 at1[2] + df[2] / dashes);
					fprintf (drvui->fpoutv, "  Rotation {\n");
					fprintf (drvui->fpoutv,
						 "   rotation 0 0 1 %f\n  }\n", -beta);
					fprintf (drvui->fpoutv, "  Rotation {\n");
					fprintf (drvui->fpoutv,
						 "   rotation 1 0 0 %f\n  }\n",
						 1.5708 - gamma);
					fprintf (drvui->fpoutv,
						 "  Material {\n   diffuseColor %s\n",
						 col_bond_v);
					fprintf (drvui->fpoutv, "  }\n");
					fprintf (drvui->fpoutv, "  Cylinder {\n");
					fprintf (drvui->fpoutv, "   parts SIDES\n");
					fprintf (drvui->fpoutv, "   radius %f\n",
						 drvui->Bond_Mult *
						 drvui->bonds[i].bond_size);
					fprintf (drvui->fpoutv, "   height %f\n", d / dashes);
					fprintf (drvui->fpoutv, "  }\n }\n");
					for (m = 2; m < dashes - 1; m += 2) {
					    mf = (float) m;
					    fprintf (drvui->fpoutv, " Separator{\n");
					    fprintf (drvui->fpoutv, "  Transform{\n ");
					    fprintf (drvui->fpoutv,
						     "   translation %5.3f %5.3f %5.3f\n  }\n",
						     at1[0] +
						     df[0] * mf / dashes,
						     at1[1] +
						     df[1] * mf / dashes,
						     at1[2] +
						     df[2] * mf / dashes);
					    fprintf (drvui->fpoutv, "  Rotation {\n");
					    fprintf (drvui->fpoutv,
					 	     "   rotation 0 0 1 %f\n  }\n", -beta);
					    fprintf (drvui->fpoutv, "  Rotation {\n");
					    fprintf (drvui->fpoutv,
						     "   rotation 1 0 0 %f\n  }\n",
						     1.5708 - gamma);
					    fprintf (drvui->fpoutv,
						     "  Material {\n   diffuseColor %s\n",
						     col_bond_v);
					    fprintf (drvui->fpoutv, "  }\n");
					    fprintf (drvui->fpoutv, "  Cylinder {\n");
					    fprintf (drvui->fpoutv, "   parts SIDES\n");
					    fprintf (drvui->fpoutv, "   radius %f\n",
						     drvui->Bond_Mult *
						     drvui->bonds[i].bond_size);
					    fprintf (drvui->fpoutv, "   height %f\n", d / dashes);
					    fprintf (drvui->fpoutv, "  }\n }\n");
					}
				    }
				}
			    }
			}
		    }		/* if distance OK */
		}		/* for k */
	    }			/* for j */
	}			/* if nvert */
    }				/* for i (bond types) */
    free (ellips_id);
    free (ellips_num);
    gluDeleteQuadric (glu_quadObj);
    fprintf (drvui->fcns, "%4d bonds.\n", Bond_Count);
    fprintf (drvui->flout, "\nGenerated %4d bonds.\n\n", Bond_Count);
}				/* end of generate_bonds */

/* ************************************************************** */
/* ************************************************************** */

void
generate_cones (void)

/* routine to generate lone pair cone descriptions and add them to edit list */
{
    int start, same_type;

    int i, j, k, l;		/* loop counters */

    int N_Bond;			/* Counter for number of bonds of a type */

    int l1;			/* loop counter for matrix operations */

    int nvert1;			/* vertex counter */

    int Bond_Count;		/* Counter for bonds output */

    char bnd1[4];		/* temporary storage of bond labels */

    float d = 0;		/* bond distance */

    float df[3];		/* vector difference */

    float beta;			/* rotation angle */

    float gamma;		/* another rotation angle */

    float at1[3], at2[3];	/*initial coordinates of bond endpoints */

    float clip0, clip1;		/* clipping scale factors */

    int m, l2, l3, out = 0;

    int p;			/* loop variables for rotation matrix calculation */

    float Z0[3], Z1[3];		/* components of rms displacement along bond */

    float mybonds[20][20];	/*temporary storage for atoms bonded to lonepair site */

    int nbnds = 0;

    float vecsum[3], vecpro[3];

    float glr, glg, glb;

    GLUquadricObj *glu_quadObj;

    char col_cone_p[40];

    char col_cone_v[40];

    if (drvui->ncone == 1)
	return;			/* exit if no work to be done */

    Bond_Count = 0;

    glu_quadObj = gluNewQuadric ();

    for (i = 1; i < drvui->ncone; ++i) {	// loop through cone types
	if (drvui->cones[i].cone_fn != drvui->frame_no)
	    continue;		// skip if wrong frame
	strcpy (col_cone_p, drvui->cones[i].col_cone);
	strcpy (col_cone_v, col_cone_p);
	Transform_VRML_Color (col_cone_v);
	Transform_POV_Color (col_cone_p);
	nvert = 0;		// empty the vertex list
	for (j = 0; j < natom; ++j) {	// loop through atom types
	    if (drvui->atoms[j].atom_fn != drvui->frame_no)
		continue;
	    if (((drvui->atoms[j].atom_n & 255) > 0) || (((drvui->atoms[j].atom_n >> 24) & 255) > 0)) {	// sphere or ellipsoid
		int o;

		o = drvui->atoms[j].atom_n >> 24;
		if (o > 0) {	// ellipsoid
		    if (check_atom_name
			(drvui->cones[i].cone_l1, drvui->ellips[o].ellips_l))
			find_all_in_box (j);	/* get all atoms of type j */
		} else {
		    o = drvui->atoms[j].atom_n & 255;
		    if (check_atom_name
			(drvui->cones[i].cone_l1, drvui->spheres[o].sphere_l)) {
			find_all_in_box (j);	/* get all atoms of type j */
		    }
		}		// o > 0
	    }
	}			/* for (j=0;j<natom;++j) */

	bnd1[0] = drvui->cones[i].cone_l1[0];
	bnd1[1] = drvui->cones[i].cone_l1[1];
	bnd1[2] = drvui->cones[i].cone_l1[2];
	bnd1[3] = drvui->cones[i].cone_l1[3];
	nvert1 = nvert;
	same_type = 0;

/* atoms i and j are not the same type */
	for (j = 0; j < natom; ++j) {	/* loop through atom types */
	    if (drvui->atoms[j].atom_fn != drvui->frame_no)
		continue;
	    if (((drvui->atoms[j].atom_n & 255) > 0) || (((drvui->atoms[j].atom_n >> 24) & 255) > 0)) {	// sphere or ellipsoid
		find_all_in_box (j);	/* get all atoms of type j */
	    }
	}			/* for (j=0;j<natom;++j) */

	if (((nvert1 > 0) && (nvert > nvert1)) || (same_type == 1)) {
	    N_Bond = 0;
	    nbnds = 0;
	    for (j = 0; j < nvert1; ++j) {	/* loop through atoms of type 1 */
		nbnds = 0;
		start = nvert1;
		if (same_type)
		    start = j;
		start = 0;
		for (k = start; k < nvert; ++k) {	/* loop on type 2 atoms  */
		    d = dist (j, k);	/* get atomic distance */
//          if (d > 0.2 && d <= drvui->cones[i].cone_height) {     /* here if distance within range */
		    if (d > 0.2 && d <= printdist) {	/* here if distance within range */
			mybonds[0][nbnds] = s_vert[3 * k];
			mybonds[1][nbnds] = s_vert[3 * k + 1];
			mybonds[2][nbnds] = s_vert[3 * k + 2];
			nbnds++;
		    }		/* if distance ok */
		}		/* for k */
		/*# */
		if (nbnds < 4 - abs (drvui->cones[i].numlonepairs)) {
		    continue;
		}
		df[0] = df[1] = drvui->Bond_Mult * drvui->bonds[i].bond_size;	/* scaling */
		df[2] = 0.5f * d;

		Bond_Count++;	/* update count */
		for (p = 0; p < 3; p++) {	/* initialize endpoints of bond */
		    Z0[p] = s_vert[3 * j + p];
		}

		for (l = 0; l < abs (drvui->cones[i].numlonepairs); l++) {
		    for (p = 0; p < 3; p++) {	/* initialize endpoints of bond */
			Z1[p] = 0;

			if (drvui->cones[i].numlonepairs == 1) {
			    for (k = 0; k < nbnds; k++)
				Z1[p] -= mybonds[p][k] - s_vert[3 * j + p];
			    Z1[p] = Z0[p] + Z1[p];
			} else if (drvui->cones[i].numlonepairs < 0) {
			    vecpro[p] = 0.;
			    vecpro[0] = (s_vert[3 * j + 1] - mybonds[1][0])
				* (s_vert[3 * j + 2] - mybonds[2][1])
				- (s_vert[3 * j + 1] - mybonds[1][1])
				* (s_vert[3 * j + 2] - mybonds[2][0]);
			    vecpro[1] = (s_vert[3 * j + 2] - mybonds[2][0])
				* (s_vert[3 * j] - mybonds[0][1])
				- (s_vert[3 * j + 2] - mybonds[2][1])
				* (s_vert[3 * j] - mybonds[0][0]);
			    vecpro[2] = (s_vert[3 * j] - mybonds[0][0])
				* (s_vert[3 * j + 1] - mybonds[1][1])
				- (s_vert[3 * j] - mybonds[0][1])
				* (s_vert[3 * j + 1] - mybonds[1][0]);
			    if (l > 0)
				Z1[p] = Z0[p] + vecpro[p];
			    else
				Z1[p] = Z0[p] - vecpro[p];
			} else if (drvui->cones[i].numlonepairs == 2) {
			    vecsum[p] = 0.;
			    vecpro[p] = 0.;
			    for (k = 0; k < nbnds; k++) {
				vecsum[p] += mybonds[p][k] - s_vert[3 * j + p];
			    }

			    vecpro[0] = (s_vert[3 * j + 1] - mybonds[1][0])
				* (s_vert[3 * j + 2] - mybonds[2][1])
				- (s_vert[3 * j + 1] - mybonds[1][1])
				* (s_vert[3 * j + 2] - mybonds[2][0]);
			    vecpro[1] = (s_vert[3 * j + 2] - mybonds[2][0])
				* (s_vert[3 * j] - mybonds[0][1])
				- (s_vert[3 * j + 2] - mybonds[2][1])
				* (s_vert[3 * j] - mybonds[0][0]);
			    vecpro[2] = (s_vert[3 * j] - mybonds[0][0])
				* (s_vert[3 * j + 1] - mybonds[1][1])
				- (s_vert[3 * j] - mybonds[0][1])
				* (s_vert[3 * j + 1] - mybonds[1][0]);
			    switch (l) {
			    case 0:
				Z1[p] = -vecsum[p] + vecpro[p];
				break;
			    case 1:
				Z1[p] = -vecsum[p] - vecpro[p];
				break;
			    }
			    Z1[p] = Z0[p] + Z1[p];
			}
		    }		/* for p */

		    d = (float) sqrt ((Z0[0] - Z1[0]) * (Z0[0] - Z1[0])
				      + (Z0[1] - Z1[1]) * (Z0[1] - Z1[1])
				      + (Z0[2] - Z1[2]) * (Z0[2] - Z1[2]));
		    d = drvui->cones[i].cone_height / d;
		    Z1[0] = Z0[0] + d * (Z1[0] - Z0[0]);
		    Z1[1] = Z0[1] + d * (Z1[1] - Z0[1]);
		    Z1[2] = Z0[2] + d * (Z1[2] - Z0[2]);

		    for (l1 = 0; l1 < 3; l1++) {
			at1[l1] = (float) Z0[l1];
			at2[l1] = (float) Z1[l1];
			df[l1] = (float) (Z1[l1] - Z0[l1]);
		    }

		    if (clipflag != 0) {
			out = 0;
			l2 = 0;
			l3 = 0;
			clip0 = 0.0f;
			clip1 = 1.0f;
			for (l1 = 0; l1 < 3; l1++) {
			    if (o_vert[3 * j + l1] <
				drvui->frames[drvui->frame_no].clip_lim[l1] - 0.01
				|| o_vert[3 * j + l1] >
				drvui->frames[drvui->frame_no].clip_lim[l1 + 3] + 0.01)
				l2 = 1;
			    if (o_vert[3 * k + l1] <
				drvui->frames[drvui->frame_no].clip_lim[l1] - 0.01
				|| o_vert[3 * k + l1] >
				drvui->frames[drvui->frame_no].clip_lim[l1 + 3] + 0.01)
				l3 = 1;
			}
			if (l2 && l3)
			    out = 1;	/* both atoms outside */
			if (l2)
			    out = 1;
			if (l3)
			    clip1 = 0.5f;
			if (l2)
			    clip0 = 0.5f;

			at1[0] = at1[0] + clip0 * df[0];
			at1[1] = at1[1] + clip0 * df[1];
			at1[2] = at1[2] + clip0 * df[2];

			at2[0] = at1[0] + clip1 * df[0];
			at2[1] = at1[1] + clip1 * df[1];
			at2[2] = at1[2] + clip1 * df[2];

		    }
		    /* if clipflag */
		    if (clipflag == 0 || out == 0) {
			if (doPOV) {
			    fprintf (drvui->fpoutp,
				     " /* LonePair on %c%c %8.5f %8.5f %8.5f */\n",
				     bnd1[0], bnd1[1], o_vert[3 * j], o_vert[3 * j + 1],
				     o_vert[3 * j + 2]);
			    fprintf (drvui->fpoutp,
				     " cone{<%8.5f, %8.5f, %8.5f>, %8.5f, <%8.5f, %8.5f, %8.5f>,  %8.5f \n",
				     at1[0], at1[1], at1[2],
				     drvui->Bond_Mult * drvui->cones[i].cone_min, at2[0],
				     at2[1], at2[2],
				     drvui->Bond_Mult * drvui->cones[i].cone_max);
			    fprintf (drvui->fpoutp, "  texture{pigment{color %s  }}\n",
				     col_cone_p);
			    fprintf (drvui->fpoutp, " }\n");
			    fprintf (drvui->fpoutp,
				     " sphere{<%8.5f, %8.5f, %8.5f> %8.5f \n", at2[0],
				     at2[1], at2[2],
				     drvui->Bond_Mult * drvui->cones[i].cone_max);
			    fprintf (drvui->fpoutp, "  texture{pigment{color %s  }}\n",
				     col_cone_p);
			    fprintf (drvui->fpoutp, " }\n");
			}
			N_Bond++;
			if (no_comment == 0 && doVrml)
			    fprintf (drvui->fpoutv,
				     "# LonePair on: %c%c %8.5f %8.5f %8.5f\n", bnd1[0],
				     bnd1[1], o_vert[3 * j], o_vert[3 * j + 1],
				     o_vert[3 * j + 2]);
			for (m = 0; m < 3; m++)
			    df[m] = (float) (at2[m] - at1[m]);
			d = (float) sqrt (df[0] * df[0] + df[1] * df[1] + df[2] * df[2]);
/*
   calculate rotation matrices needed to rotate bond about Z onto Y-Z plane and then about X so that
   bond is parallel to Z.  Finally, apply inverse transformation to cylinder, which started parallel
   to X.  Resulting cylinder will be parallel to vector between atoms.
*/
			beta = (float) atan2 (df[0], df[1] + 0.0000001f);	/* rotation angle about Z (in radians) */
			gamma = (float) sqrt (df[0] * df[0] + df[1] * df[1]);
			gamma = (float) atan2 (gamma, df[2]);	/* Rotation angle about X */
			glPushMatrix ();
			sscanf (col_cone_v, "%f %f %f", &glr, &glg, &glb);
			glColor3f (glr, glg, glb);
			glTranslatef (at1[0], at1[1], at1[2]);
			glRotatef (-beta * (float) RAD, 0.0f, 0.0f, 1.0f);
			glRotatef (-gamma * (float) RAD, 1.0f, 0.0f, 0.0f);
			gluQuadricDrawStyle (glu_quadObj, GLU_FILL);
			gluCylinder (glu_quadObj,
				     drvui->Bond_Mult * drvui->cones[i].cone_min,
				     drvui->Bond_Mult * drvui->cones[i].cone_max, d, 10,
				     1);
			glPopMatrix ();
			glPushMatrix ();
			glTranslatef (at2[0], at2[1], at2[2]);
			glutSolidSphere (drvui->cones[i].cone_max, 10, 10);
			glPopMatrix ();
			if (doVrml) {
			    if (Vrml2) {
				fprintf (drvui->fpoutv, "  Transform{\n ");
				fprintf (drvui->fpoutv,
					 "   translation %5.3f %5.3f %5.3f\n",
					 (at1[0] + at2[0]) / 2., (at1[1] + at2[1]) / 2.,
					 (at1[2] + at2[2]) / 2.);
				fprintf (drvui->fpoutv, "   rotation 0 0 1 %f\n", -beta);
				fprintf (drvui->fpoutv,
					 "   children Transform { rotation 0 1 0 %f\n",
					 3.14);
				fprintf (drvui->fpoutv,
					 "   children Transform { rotation 1 0 0 %f\n",
					 1.5708 + gamma);
				fprintf (drvui->fpoutv, " children [ Shape {\n");
				fprintf (drvui->fpoutv,
					 " appearance Appearance { material Material {   diffuseColor %s} }\n",
					 col_cone_v);
				fprintf (drvui->fpoutv, "  geometry Cone {");
				fprintf (drvui->fpoutv, "bottomRadius %f",
					 drvui->Bond_Mult * drvui->cones[i].cone_max);
				fprintf (drvui->fpoutv, " height %f", d);
				fprintf (drvui->fpoutv, "}}]}}}\n");
				fprintf (drvui->fpoutv, "  Transform {");
				fprintf (drvui->fpoutv,
					 " translation %8.5f %8.5f %8.5f\n", at2[0],
					 at2[1], at2[2]);
				fprintf (drvui->fpoutv, " children [\n Shape {");
				fprintf (drvui->fpoutv,
					 "appearance Appearance {material Material { diffuseColor %s} }\n",
					 col_cone_v);
				fprintf (drvui->fpoutv,
					 "  geometry Sphere { radius %8.5f}\n",
					 drvui->cones[i].cone_max);
				fprintf (drvui->fpoutv, "}]}\n");
			    } else {
				fprintf (drvui->fpoutv, " Separator{\n");
				fprintf (drvui->fpoutv, "  Transform{\n ");
				fprintf (drvui->fpoutv,
					 "   translation %5.3f %5.3f %5.3f\n  }\n",
					 (at1[0] + at2[0]) / 2., (at1[1] + at2[1]) / 2.,
					 (at1[2] + at2[2]) / 2.);
				fprintf (drvui->fpoutv, "  Rotation {\n");
				fprintf (drvui->fpoutv, "   rotation 0 0 1 %f\n  }\n",
					 -beta);
				fprintf (drvui->fpoutv, "  Rotation {\n");
				fprintf (drvui->fpoutv, "   rotation 0 1 0  %f\n  }\n",
					 3.14);
				fprintf (drvui->fpoutv, "  Rotation {\n");
				fprintf (drvui->fpoutv, "   rotation 1 0 0 %f\n  }\n",
					 1.5708 + gamma);
				fprintf (drvui->fpoutv,
					 "  Material {\n   diffuseColor %s\n",
					 col_cone_v);
				fprintf (drvui->fpoutv, "  }\n");
				fprintf (drvui->fpoutv, "  Cone {\n");
				fprintf (drvui->fpoutv, "   parts SIDES\n");
				fprintf (drvui->fpoutv, "   bottomRadius %f\n",
					 drvui->Bond_Mult * drvui->cones[i].cone_max);
				fprintf (drvui->fpoutv, "   height %f\n", d);
				fprintf (drvui->fpoutv, "  }\n }\n");
				fprintf (drvui->fpoutv, " Separator {");
				fprintf (drvui->fpoutv, "  Transform {");
				fprintf (drvui->fpoutv,
					 " translation %8.5f %8.5f %8.5f}\n", at2[0],
					 at2[1], at2[2]);
				fprintf (drvui->fpoutv, "  Material {");
				fprintf (drvui->fpoutv, " shininess 0.3\n");
				fprintf (drvui->fpoutv, " diffuseColor %s}\n",
					 col_cone_v);
				fprintf (drvui->fpoutv, "  Sphere { radius %8.5f}\n",
					 drvui->cones[i].cone_max);
				fprintf (drvui->fpoutv, " }\n");
			    }	/* if (Vrml2) */
			}
			if (doAsy) {
			    fprintf (drvui->fpouta, "draw (pic, shift(%5.2f,%5.2f,%5.2f)*rotate(%.2f,Z)*rotate(%.2f,X)*\n",
				     at2[0],at2[1],at2[2],-beta*RAD,-gamma*RAD);
			    fprintf (drvui->fpouta, "\tscale(%.2f,%.2f,%.2f)*unitcone,rgb(%.2f,%.2f,%.2f));\n",
				     drvui->Bond_Mult * drvui->cones[i].cone_max,drvui->Bond_Mult * drvui->cones[i].cone_max,-d,glr,glg,glb);
			    fprintf (drvui->fpouta, "draw (pic, shift(%5.2f,%5.2f,%5.2f)*scale3(%.2f)*unitsphere,rgb(%.2f,%.2f,%.2f));\n",
				     at2[0],at2[1],at2[2],drvui->cones[i].cone_max,glr,glg,glb);
			}
		    }
		    /* if not clipped */
		}		/* for l */
	    }			/* for j */
	}			/* if nvert */
    }				/* for i */
    gluDeleteQuadric (glu_quadObj);
    fprintf (drvui->fcns, "%4d cones.\n", Bond_Count);
    fprintf (drvui->flout, "\nGenerated %4d cones.\n\n", Bond_Count);
}				/* end of generate_cones */

/* ************************************************************** */
/* ************************************************************** */

void
get_atom_id (void)

/* routine to rework atom numbers to point at polyhedron, sphere
   plane, or ellipsoid characteristics list */
{
    int i, j;			/* local variables */

/* lines commented out with a question mark would be needed for correctness, but
   would break the stishovite example that relies on polyhedral completion carrying
   over from the first frame into the second */

    for (i = 0; i < natom; ++i) {
	if (drvui->atoms[i].atom_fn != drvui->frame_no)
	    continue;
	drvui->atoms[i].sv_atom_n = drvui->atoms[i].atom_n;	/* save original atom number */
	drvui->atoms[i].atom_n = 0;	// clear encoded location
	for (j = 1; j < drvui->nsphere; ++j) {
//?      if (drvui->sphere_fn[j] != drvui->frame_no) continue; 
	    if (check_atom_name (drvui->atoms[i].atom_l, drvui->spheres[j].sphere_l)
		&& (drvui->atoms[i].sv_atom_n == drvui->spheres[j].sphere_n
		    || drvui->spheres[j].sphere_n == -1)) {
		drvui->atoms[i].atom_n = j;	// atom is in sphere list
	    }
	}
	for (j = 1; j < drvui->npoly; ++j) {
//?      if (drvui->poly_fn[j] != drvui->frame_no) continue; 
	    if (check_atom_name (drvui->atoms[i].atom_l, drvui->polyhedra[j].poly_l)) {
		drvui->atoms[i].atom_n += (j << 8);	// atom is in poly. list
	    }
	}
	for (j = 1; j < drvui->nplane; ++j) {
//?      if (drvui->plane_fn[j] != drvui->frame_no) continue; 
	    if (check_atom_name (drvui->atoms[i].atom_l, drvui->planes[j].plane_l)) {
		drvui->atoms[i].atom_n += (j << 16);	// atom is in plane list
	    }
	}
	if (drvui->do_ellipsoids) {
	    for (j = 1; j < drvui->n_ellips; ++j) {
		if (check_atom_name (drvui->atoms[i].atom_l, drvui->ellips[j].ellips_l)
		    && (drvui->atoms[i].sv_atom_n == drvui->ellips[j].ellips_n)) {
		    drvui->atoms[i].atom_n += (j << 24);	// atom is in ellipsoid list
		}
	    }
	}
    }
}				// end of get_atom_id

/* ************************************************************** */
/* ************************************************************** */

void
get_input (int Quick)

/* Routine to process input file, make lattice metric and do preliminary
   ellipsoid processing */
{
    int i, j, k;

    float biso;

    if (!Quick)
	fprintf (drvui->flout,
		 "\n********************* Input: ************************\n");
    read_inp (Quick);
    if (Quick)
	return;
    fprintf (drvui->flout, "*****************************************************\n\n");
    if (!Labels)
	Display_axes = 0;	/* If no labels, no vectors */
    if (drvui->frame_no == 1 && doVrml) {
	if (Vrml2)
	    fprintf (drvui->fcns, "\nGenerating VRML97 Output.\n");
	else
	    fprintf (drvui->fcns, "\nGenerating VRML1.0 Output.\n");
    }
    print_sym ();		/* print symmetry information */
    if (drvui->xyzoff_read)
	drvui->origin1_flag = 0;
    (void) fclose (drvui->fpin);	/* close file */
    make_bmat (drvui->sys, drvui->lat_con, drvui->b_mat, drvui->ginv, drvui->rec_lat_con);	/* create the lattice metric */
    fprintf (drvui->flout, "Lattice metrics\n");
    fprintf (drvui->flout, "%8.3f %8.3f %8.3f\n", drvui->b_mat[0][0], drvui->b_mat[0][1],
	     drvui->b_mat[0][2]);
    fprintf (drvui->flout, "%8.3f %8.3f %8.3f\n", drvui->b_mat[1][0], drvui->b_mat[1][1],
	     drvui->b_mat[1][2]);
    fprintf (drvui->flout, "%8.3f %8.3f %8.3f\n", drvui->b_mat[2][0], drvui->b_mat[2][1],
	     drvui->b_mat[2][2]);
/* set limits of search regions */
    if (packflag == 0 && boxflag == 0 && drvui->frame_no == 1) {	/* neither given and frame 1 */
	for (i = 0; i < 3; i++) {
	    drvui->frames[drvui->frame_no].cryst_lim[i] = origin[i] - 0.5f;	/* set to generate one cell */
	    drvui->frames[drvui->frame_no].cryst_lim[i + 3] = origin[i] + 0.5f;
	}
	packflag = 1;
    }
    if (packflag != 0) {	/* pack command given - calculate approximate bounds */
	float temp, temp1;

	for (i = 0; i < 3; i++) {
	    temp = origin[i] - drvui->frames[drvui->frame_no].cryst_lim[i];
	    temp1 = drvui->frames[drvui->frame_no].cryst_lim[i + 3] - origin[i];
	    if (temp1 > temp)
		temp = temp1;
	    boxlim[i] = (temp + 0.2f) * drvui->lat_con[i] / 0.5f;
	}
    } else {			/* bounds command given - calculate approximate pack limits */
	for (i = 0; i < 3; i++) {
	    drvui->frames[drvui->frame_no].cryst_lim[i] =
		origin[i] - boxlim[i] / drvui->lat_con[i] - 0.2f;
	    drvui->frames[drvui->frame_no].cryst_lim[i + 3] =
		boxlim[i] / drvui->lat_con[i] - origin[i] + 0.2f;
	}
    }
    if ((Options & L_OPT) != 0)
	Calc_Rot (drvui->lookat_v1, drvui->lookat_v2);	// Get rotation angles for "lookat"
    if (drvui->n_ellips > 1) {
	convert_ellipsoid ();	/* change Bij and Uij forms to betaij */
	fprintf (drvui->flout, "\nAnalysis of Thermal Ellipsoids:\n");
	for (i = 1; i < drvui->n_ellips; i++) {
	    if (drvui->ellips[i].ellips_ismod != 0) {
		fprintf (drvui->flout, " (deferred for %c%c%c%c %d due to modulation)\n",
			 drvui->ellips[i].ellips_l[0], drvui->ellips[i].ellips_l[1],
			 drvui->ellips[i].ellips_l[2], drvui->ellips[i].ellips_l[3],
			 drvui->ellips[i].ellips_n);
		continue;	// modulation must be applied to unchanged Uij values
	    }
	    if (eigen
		(&biso, drvui->ellips[i].ellips, drvui->ellips[i].ellips_RMS,
		 drvui->ellips[i].ellips_EV) == 0) {
		fprintf (drvui->flout,
			 "Ellipsoid for %c%c%c%c %d is non-positive definite and has "
			 "been removed from the active list.\n",
			 drvui->ellips[i].ellips_l[0], drvui->ellips[i].ellips_l[1],
			 drvui->ellips[i].ellips_l[2], drvui->ellips[i].ellips_l[3],
			 drvui->ellips[i].ellips_n);
		drvui->ellips[i].ell_type = -100;	// disable this one
		drvui->ellips[i].ellips_RMS[0] = (float) sqrt (biso / 78.957);	/* convert Biso to mu */
	    } else {
		drvui->ellips[i].ell_type += 1;
		fprintf (drvui->flout,
			 "\n%c%c%c%c %3d Equivalent Isotropic B: %6.2f, U: %6.3f\n",
			 drvui->ellips[i].ellips_l[0], drvui->ellips[i].ellips_l[1],
			 drvui->ellips[i].ellips_l[2], drvui->ellips[i].ellips_l[3],
			 drvui->ellips[i].ellips_n, biso, biso / 78.957);
		fprintf (drvui->flout, "      RMS Amplitudes: %6.3f %6.3f %6.3f\n",
			 drvui->ellips[i].ellips_RMS[0], drvui->ellips[i].ellips_RMS[1],
			 drvui->ellips[i].ellips_RMS[2]);
		fprintf (drvui->flout, "      Eigen Vectors (by row):\n");
		for (j = 0; j < 3; j++) {
		    fprintf (drvui->flout, "     ");
		    for (k = 0; k < 3; k++)
			fprintf (drvui->flout, "     %8.5f",
				 drvui->ellips[i].ellips_EV[j][k]);
		    fprintf (drvui->flout, "\n");
		}
	    }
	}
	drvui->Ellipsoid_Scale = P_to_C (drvui->Ellipsoid_Prob);
    }
}				/* end of get_input */

/* ************************************************************** */
/* ************************************************************** */

void
Locate_Triple (void)

/* routine to try to place the vector triple where it will be seen in the
POV diagram */
{
    int i, j;

    float temp, pos[3], cpos[3];

    for (i = 0; i < 3; ++i) {	// convert origin to Cartesian
	pos[i] = 0.0f;
	for (j = 0; j < 3; ++j)
	    pos[i] -= (float) drvui->b_mat[i][j] * origin[j];
    }
    for (i = 0; i <= 2; ++i) {	// calculate position of point after POV rotation
	cpos[i] = 0.0f;
	for (j = 0; j <= 2; ++j)
	    cpos[i] += (float) G_Rot[j][i] * pos[j];
    }
    for (i = 0; i < 3; i++) {
	temp = 0.04f * (POV_Max[i] - POV_Min[i]);

	if (cpos[i] > 0.)
	    pos[i] = POV_Max[i] + temp;
	else
	    pos[i] = POV_Min[i] - temp;
    }
// calculate coordinates of pos before rotation
    for (i = 0; i < 3; i++) {
	offset[i] = 0.0f;
	for (j = 0; j < 3; j++) {
	    offset[i] += (float) G_Rot[i][j] * pos[j];
	}
    }
}

/* ************************************************************** */
/* ************************************************************** */

void
make_bmat (int sys, float lat_con[6], double b_mat[3][3], float ginv[3][3],
	   float rec_lat_con[6])
{
    float snal, snbe, snga, csal, csbe, csga, vol;

    int i, j;

    for (i = 3; i <= 5; ++i)
	if (lat_con[i] <= 0)
	    lat_con[i] = 90.0f;
    if (sys == 5) {
	lat_con[1] = lat_con[0];
	lat_con[5] = 120.0f;
    }
    if (sys == 4)
	lat_con[1] = lat_con[0];
    if (sys == 6) {
	lat_con[1] = lat_con[0];
	lat_con[2] = lat_con[0];
    }
    for (i = 0; i <= 2; ++i)
	for (j = 0; j <= 2; ++j)
	    b_mat[j][i] = 0.0f;	/* initialize matrix */
    b_mat[0][0] = 1.0f;
    csal = (float) cos (lat_con[3] / RAD);
    snal = (float) sin (lat_con[3] / RAD);
    csbe = (float) cos (lat_con[4] / RAD);
    snbe = (float) sin (lat_con[4] / RAD);
    csga = (float) cos (lat_con[5] / RAD);
    snga = (float) sin (lat_con[5] / RAD);
    for (i = 0; i < 3; i++)
	ginv[i][i] = lat_con[i] * lat_con[i];
    ginv[0][1] = ginv[1][0] = lat_con[0] * lat_con[1] * csga;
    ginv[0][2] = ginv[2][0] = lat_con[0] * lat_con[2] * csbe;
    ginv[1][2] = ginv[2][1] = lat_con[1] * lat_con[2] * csal;
    b_mat[0][1] = csga;
    b_mat[1][1] = snga;
    b_mat[0][2] = csbe;
    b_mat[1][2] = (csal - csbe * csga) / snga;
    b_mat[2][2] = (float) sqrt (1.0 - b_mat[0][2] * b_mat[0][2] -
				b_mat[1][2] * b_mat[1][2]);
    for (i = 0; i <= 2; ++i)
	for (j = 0; j <= 2; ++j)
	    b_mat[j][i] *= lat_con[i];
/* calculate cell volume and reciprocal lengths */
    vol = (float) (b_mat[0][0] * (b_mat[1][1] * b_mat[2][2] - b_mat[1][2] * b_mat[2][1])
		   - b_mat[0][1] * (b_mat[1][0] * b_mat[2][2] - b_mat[1][2] * b_mat[2][0])
		   + b_mat[0][2] * (b_mat[1][0] * b_mat[2][1] -
				    b_mat[1][1] * b_mat[2][0]));
    rec_lat_con[0] = lat_con[1] * lat_con[2] * snal / vol;
    rec_lat_con[1] = lat_con[0] * lat_con[2] * snbe / vol;
    rec_lat_con[2] = lat_con[1] * lat_con[0] * snga / vol;
    rec_lat_con[3] = (csbe * csga - csal) / (snbe * snga);
    rec_lat_con[4] = (csal * csga - csbe) / (snal * snga);
    rec_lat_con[5] = (csal * csbe - csga) / (snal * snbe);
    drvui->subsys_ref_volume = 1.0f / vol;	/* save the reciprocal volume */

}				/* end of make_bmat */

/* ************************************************************** */
/* ************************************************************** */

void
modulate_parameters (float vert[3], double *occ, int sym_no, int atom_no)
{
/* routine to generate the occupancy and positional parameters for an atom with
   "average" positional parameters in 'vert' 
   This routine relies very heavily on S. van Smaalen, Z. Kristallogr. 219 (2004) 681-691.
   A second reference used for 2 and 3D modulation is Jakubowicz et al., Phys. Rev. B 63 (2001). */

    double x4[3] = { 0.0, 0.0, 0.0 }, rx4[3] = {
    0.0, 0.0, 0.0};
    double arg;
    float xadd[3] = { 0.0, 0.0, 0.0 };
    float rxadd[3];

    float eps_inv[3][3];

    int id, i1, i2, i3, j, k, theatom;

    int axis;

    *occ = 1.0;
    if (drvui->modulated <= 0)
	return;
    for (j = 0; j < 3; j++)
	for (k = 0; k < 3; k++)
	    eps_inv[j][k] = (float) drvui->ss_m[sym_no][j][k];
    (void) matinv (eps_inv);	/* the 'magical' epsilon^(-1) (van Smaalen eq 14) */
    for (k = 0; k < drvui->modulated; k++) {
	x4[k] = drvui->phaseshift[k];
	for (j = 0; j < 3; j++)
	    x4[k] += drvui->cell_vec[k][j] * vert[j];	/* contains x4=q1.r, x5=q2.r, x6=q3.r */
    }
    for (k = 0; k < 3; k++) {
	rx4[k] = -drvui->ts_m[sym_no][k];
	for (j = 0; j < 3; j++) {
	    rx4[k] += eps_inv[k][j] * x4[j];	/* rx4 is x4 transformed by epsilon^(-1)  */
	}
	rx4[k] = fmod (rx4[k], 1.0);
	if (rx4[k] < 0)
	    rx4[k] += 1.;
    }
    if (drvui->atoms[atom_no].occ_ismod & 2) {	/* check for crenel occupancy modulation */
	float lower = drvui->modulate_x[atom_no].atom_occ_crenel[0];

	float upper = drvui->modulate_x[atom_no].atom_occ_crenel[1];

	if ((rx4[0] < lower) || (rx4[0] > upper))
	    *occ = 0.0;
	if (((upper - 1.0f) > rx4[0]) || ((lower + 1.0f) < rx4[0]))
	    *occ = 1.0;		/* handle the cases when lower < 0 or upper > 1 */
	if (*occ == 0.0)
	    return;		/* done if crenel occ is zero */
    }

    if (drvui->atoms[atom_no].occ_ismod & 1) {	/* check for Fourier occupancy modulation */
	*occ = drvui->atoms[atom_no].occupancy;
	for (j = 0; j < 3; j++)
	    rx4[j] *= 2.0 * PI;
	for (k = 0; k < drvui->no_site_occ; k++) {	/* loop through the fourier mod vectors */
	    id = drvui->modulate_x[k].atom_occpar_id - 1;
	    theatom = drvui->modulate_x[k].atom_occpar_atom;
	    if (theatom != atom_no)
		continue;
	    i1 = drvui->modulate_gbl[id].vector_mult[0];	/* multiplier for first cell mod vector */
	    i2 = drvui->modulate_gbl[id].vector_mult[1];	/* multiplier for second cell mod vector */
	    i3 = drvui->modulate_gbl[id].vector_mult[2];	/* multiplier for third cell mod vector */
	    arg = i1 * rx4[0] + i2 * rx4[1] + i3 * rx4[2];
	    *occ += (float) (drvui->modulate_x[k].atom_occpar[0] * cos (arg)
			     + drvui->modulate_x[k].atom_occpar[1] * sin (arg));
	}
	if (*occ < drvui->atoms[atom_no].min_occ) {
	    *occ = 0.;
	    return;
	}
    }

    if (drvui->atoms[atom_no].atom_ismod & 1) {	/* Fourier position modulation */
	if (!(drvui->atoms[atom_no].occ_ismod & 1)) {	//only if not already done in previous step
	    for (j = 0; j < 3; j++)
		rx4[j] *= 2.0 * PI;
	}
	for (j = 0; j < drvui->no_site_displace; j++) {	/* loop through the fourier mod vectors */
	    id = drvui->modulate_3x[j].atom_modpar_id - 1;
	    axis = drvui->modulate_3x[j].atom_modpar_axis;
	    theatom = drvui->modulate_3x[j].atom_modpar_atom;
	    if (theatom != atom_no)
		continue;
	    i1 = drvui->modulate_gbl[id].vector_mult[0];	/* multiplier for first cell mod vector */
	    i2 = drvui->modulate_gbl[id].vector_mult[1];	/* multiplier for second cell mod vector */
	    i3 = drvui->modulate_gbl[id].vector_mult[2];	/* multiplier for third cell mod vector */
	    arg = i1 * rx4[0] + i2 * rx4[1] + i3 * rx4[2];
	    xadd[axis] += (float) (drvui->modulate_3x[j].atom_modpar[0] * cos (arg)
				   + drvui->modulate_3x[j].atom_modpar[1] * sin (arg));
	}
	for (k = 0; k < 3; k++) {
	    rxadd[k] = 0.0;
	    for (i1 = 0; i1 < 3; i1++) {
		rxadd[k] += (float) drvui->ss[sym_no][k][i1] * xadd[i1];
	    }
	    vert[k] += rxadd[k];
	}
    }

    if (drvui->atoms[atom_no].atom_ismod & 2) {	/* Position modulation by a sawtooth function */
	if ((drvui->atoms[atom_no].occ_ismod & 1) ^ (drvui->atoms[atom_no].atom_ismod & 1)) {	// undo any previous multiplication by 2pi
	    for (j = 0; j < 3; j++)
		rx4[j] /= 2.0 * PI;
	}
	for (k = 0; k < 3; k++)
	    xadd[k] = 2.0f * drvui->modulate_x[atom_no].atom_mod_sawtooth[k]
		* (((float) rx4[k] - drvui->modulate_x[atom_no].atom_mod_sawtooth[3])
		   / drvui->modulate_x[atom_no].atom_mod_sawtooth[4]);

	for (k = 0; k < 3; k++) {
	    xadd[k] = 0.0;
	    for (i1 = 0; i1 < 3; i1++) {
		rxadd[k] += (float) drvui->ss[sym_no][k][i1] * xadd[i1];
	    }
	    vert[k] += rxadd[k];
	}
    }
}

/* ************************************************************** */
/* ************************************************************** */

int
modulate_uij (float vert[3], int ellips_no, int atom_no, int sym_no, float uij[])
{
    double x4[3] = { 0.0, 0.0, 0.0 }, rx4[3] = {
    0.0, 0.0, 0.0};
    double arg;

    float eps_inv[3][3];

    int id, i1, i2, i3, j, k, theatom;

    int term;

    float biso;

//    float uij[6];

    for (j = 0; j < 6; j++)
	uij[j] = drvui->ellips[ellips_no].ellips[j];	/* copy unmodulated values */
    for (k = 0; k < drvui->modulated; k++) {
	x4[k] = drvui->phaseshift[k];
	for (j = 0; j < 3; j++)
	    x4[k] += drvui->cell_vec[k][j] * vert[j];	/* contains x4=q1.r, x5=q2.r, x6=q3.r */
    }
    for (k = 0; k < 3; k++) {
	rx4[k] = drvui->ts_m[sym_no][k];
	for (j = 0; j < 3; j++) {
	    eps_inv[j][k] = (float) drvui->ss_m[sym_no][j][k];
	    rx4[k] += drvui->ss_m[sym_no][k][j] * x4[j];	/* rx4 is x4 transformed by the _mm sym op */
	}
	rx4[k] = fmod (rx4[k], 1.0);
    }
    (void) matinv (eps_inv);	/* the 'magical' epsilon^(-1) (van Smaalen eq 14) */

    for (j = 0; j < 3; j++)
	rx4[j] *= 2.0 * PI;
    for (k = 0; k < drvui->no_site_U_terms; k++) {	/* loop through the fourier mod vectors */
	id = drvui->modulate_3t[k].ellips_modpar_id - 1;
	theatom = drvui->modulate_x[k].ellips_modpar_atom;
	if (theatom != atom_no)
	    continue;
	term = drvui->modulate_3t[k].ellips_modpar_term;
	i1 = drvui->modulate_gbl[id].vector_mult[0];	/* multiplier for first cell mod vector */
	i2 = drvui->modulate_gbl[id].vector_mult[1];	/* multiplier for second cell mod vector */
	i3 = drvui->modulate_gbl[id].vector_mult[2];	/* multiplier for third cell mod vector */
	arg = i1 * rx4[0] + i2 * rx4[1] + i3 * rx4[2];
	uij[term] += (float) (drvui->modulate_3t[k].ellips_modpar[0] * cos (arg)
			      + drvui->modulate_3t[k].ellips_modpar[1] * sin (arg));
    }

// FIXME? Does rotation by symmetry need to be done before making the modulation correction?

// convert from Uij to betaij

    for (j = 0; j < 3; ++j)
	uij[j] *= drvui->rec_lat_con[j] * drvui->rec_lat_con[j] * 0.25f;
    uij[3] *= drvui->rec_lat_con[0] * drvui->rec_lat_con[1] * 0.25f;
    uij[4] *= drvui->rec_lat_con[0] * drvui->rec_lat_con[2] * 0.25f;
    uij[5] *= drvui->rec_lat_con[1] * drvui->rec_lat_con[2] * 0.25f;
// multiply Uij by 8pi^2
    for (j = 0; j < 6; ++j)
	uij[j] *= 78.9568f;

// calculate eigenvalues and eigenvectors

    if (eigen
	(&biso, uij, drvui->ellips[ellips_no].ellips_RMS,
	 drvui->ellips[ellips_no].ellips_EV) == 0) {
	fprintf (drvui->flout,
		 "Ellipsoid for %c%c%c%c %d is non-positive definite and has "
		 "been removed from the active list.\n",
		 drvui->ellips[ellips_no].ellips_l[0],
		 drvui->ellips[ellips_no].ellips_l[1],
		 drvui->ellips[ellips_no].ellips_l[2],
		 drvui->ellips[ellips_no].ellips_l[3], drvui->ellips[ellips_no].ellips_n);
	return (1);
    } else {
	fprintf (drvui->flout, "\n%c%c%c%c %3d Equivalent Isotropic B: %6.2f, U: %6.3f\n",
		 drvui->ellips[ellips_no].ellips_l[0],
		 drvui->ellips[ellips_no].ellips_l[1],
		 drvui->ellips[ellips_no].ellips_l[2],
		 drvui->ellips[ellips_no].ellips_l[3], drvui->ellips[ellips_no].ellips_n,
		 biso, biso / 78.957);
	fprintf (drvui->flout, "      RMS Amplitudes: %6.3f %6.3f %6.3f\n",
		 drvui->ellips[ellips_no].ellips_RMS[0],
		 drvui->ellips[ellips_no].ellips_RMS[1],
		 drvui->ellips[ellips_no].ellips_RMS[2]);
	fprintf (drvui->flout, "      Eigen Vectors (by row):\n");
	for (j = 0; j < 3; j++) {
	    fprintf (drvui->flout, "     ");
	    for (k = 0; k < 3; k++)
		fprintf (drvui->flout, "     %8.5f",
			 drvui->ellips[ellips_no].ellips_EV[j][k]);
	    fprintf (drvui->flout, "\n");
	}
    }
    return 0;

}

/* ************************************************************** */
/* ************************************************************** */

float *
polygon_normal_3d (int n, float v[])
/* helper function for not_in_slab(), simplified from John Burkardt's geometry.c */
{
    int i;

    int j;

    float *normal;

    float normal_norm;

    float *p;

    float *v1;

    float *v2;

    normal = new float[3];

    v1 = new float[3];

    v2 = new float[3];

    p = new float[3];

    normal[0] = normal[1] = normal[2] = 0.;

    for (i = 0; i < 3; i++) {
	v1[i] = v[i + 1 * 3] - v[i + 0 * 3];
    }

    for (j = 2; j < n; j++) {
	for (i = 0; i < 3; i++) {
	    v2[i] = v[i + j * 3] - v[i + 0 * 3];
	}

	vcross (v1, v2, p);

	for (i = 0; i < 3; i++) {
	    normal[i] = normal[i] + p[i];
	}

	for (i = 0; i < 3; i++) {
	    v1[i] = v2[i];
	}

    }
    normal_norm =
	(float) sqrt (normal[0] * normal[0] + normal[1] * normal[1] +
		      normal[2] * normal[2]);

    if (normal_norm != 0.0) {
	for (i = 0; i < 3; i++) {
	    normal[i] = normal[i] / normal_norm;
	}
    }

    delete[]v1;
    delete[]v2;

    return normal;
}

/* ************************************************************** */
/* ************************************************************** */

float
polygon_solid_angle_3d (int n, float v[], float p[3])
/* helper function for not_in_slab(), simplified from John Burkardt's geometry.c */
{
    float a[3];

    double angle;

    float area = 0.0;

    float b[3];

    float c[3];

    int j;

    int jp1;

    float normal1[3];

    float normal1_norm;

    float normal2[3];

    float normal2_norm;

    float *plane;

    float r1[3];

    double s;

    float value;

    if (n < 3)
	return 0.0;

    plane = polygon_normal_3d (n, v);

    a[0] = v[0 + (n - 1) * 3] - v[0 + 0 * 3];
    a[1] = v[1 + (n - 1) * 3] - v[1 + 0 * 3];
    a[2] = v[2 + (n - 1) * 3] - v[2 + 0 * 3];

    for (j = 0; j < n; j++) {
	r1[0] = v[0 + j * 3] - p[0];
	r1[1] = v[1 + j * 3] - p[1];
	r1[2] = v[2 + j * 3] - p[2];

	jp1 = j + 1;
	if (jp1 > n - 1)
	    jp1 = 0;

	b[0] = v[0 + jp1 * 3] - v[0 + j * 3];
	b[1] = v[1 + jp1 * 3] - v[1 + j * 3];
	b[2] = v[2 + jp1 * 3] - v[2 + j * 3];

	vcross (a, r1, normal1);

	normal1_norm =
	    (float) sqrt (normal1[0] * normal1[0] + normal1[1] * normal1[1] +
			  normal1[2] * normal1[2]);

	vcross (r1, b, normal2);

	normal2_norm =
	    (float) sqrt (normal2[0] * normal2[0] + normal2[1] * normal2[1] +
			  normal2[2] * normal2[2]);

	s = vdot (normal1, normal2)
	    / (normal1_norm * normal2_norm);

	if (s <= -1.) {
	    angle = PI;
	} else if (s >= 1.) {
	    angle = 0.;
	} else
	    angle = acos (s);

	vcross (a, plane, c);
	s = vdot (b, c);

	if (0.0 < s)
	    area += (float) (PI - angle);
	else
	    area += (float) (PI + angle);

	a[0] = -b[0];
	a[1] = -b[1];
	a[2] = -b[2];

    }
    area = area - (float) (PI * (double) (n - 2));

    if (0.0 < vdot (plane, r1))
	value = -area;
    else
	value = area;

    delete[]plane;

    return value;
}

/* ************************************************************** */
/* ************************************************************** */

int
not_in_slab (float x, float y, float z)
/* adapted from polyhedron_contains_point_3d (Carvalho et al. Graphics Gems V)
   as contained in John Burkardt's geometry.c on http://www.csit.fsu.edu/~burkardt/ */
{
    float area;

    int face;

    int i, k;

    int node;

    int node_num_face = 4;

    float v_face[3 * 4];

    float p[3];
    static const int face_point[4 * 6] =
	{ 1, 2, 4, 3, 1, 3, 6, 5, 1, 5, 7, 2, 5, 6, 8, 7, 2, 7, 8, 4, 3, 4, 8, 6 };


    p[0] = x;
    p[1] = y;
    p[2] = z;

    area = 0.0;
    for (face = 0; face < 6; face++) {

	for (k = 0; k < 4; k++) {
	    node = face_point[k + face * 4];
	    for (i = 0; i < 3; i++) {
		v_face[i + k * 3] = slabv[i + (node - 1) * 3];
	    }
	}
	area = area + polygon_solid_angle_3d (node_num_face, v_face, p);
    }
    if (area < -2. * PI || 2. * PI < area)
	return 0;
    else
	return 1;
}

/* ************************************************************** */
/* ************************************************************** */

void
Output_Spheres (float *radii, int i)

/* routine to output the spheres with vertices defined in s_vert and o_vert
   The radius of the sphere is in radii */
{
    float glr, glg, glb;

    int k;

    int outside = 0;

    int l1, l2;

    int omit;

    char sphere_col_p[40];

    char sphere_col_v[40];

    strcpy (sphere_col_p, drvui->spheres[i].sphere_col);
    strcpy (sphere_col_v, sphere_col_p);
    Transform_VRML_Color (sphere_col_v);
    Transform_POV_Color (sphere_col_p);

    for (k = 0; k < nvert; ++k) {
	if (clipflag != 0) {
	    outside = 0;
	    for (l1 = 0; l1 < 3; l1++)
		if (o_vert[3 * k + l1] <
		    drvui->frames[drvui->frame_no].clip_lim[l1] - 0.01
		    || o_vert[3 * k + l1] >
		    drvui->frames[drvui->frame_no].clip_lim[l1 + 3] + 0.01)
		    outside = 1;
	}

	omit = 0;
	for (l1 = 0; l1 < Omit->nomits; l1++) {
	    if (Omit->omit1[l1] == i && Omit->omit2[l1] == k)
		omit = 1;
	}

	if (omit == 0 && (clipflag == 0 || outside == 0)) {
	    float vert1[3], vert2[3];

	    for (l1 = 0; l1 < 3; l1++)
		vert1[l1] = o_vert[3 * k + l1];
	    for (l1 = 0; l1 < 3; ++l1) {	/* convert vertex coordinates to Cartesian */
		vert2[l1] = 0.0f;
		for (l2 = 0; l2 < 3; ++l2)
		    vert2[l1] += (float) drvui->b_mat[l1][l2] * (vert1[l2] - origin[l2]);
	    }
	    if (doPOV) {
		fprintf (drvui->fpoutp,
			 " /* Sphere for %c%c%c%c at %8.5f %8.5f %8.5f */ \n",
			 drvui->spheres[i].sphere_l[0], drvui->spheres[i].sphere_l[1],
			 drvui->spheres[i].sphere_l[2], drvui->spheres[i].sphere_l[3],
			 vert1[0], vert1[1], vert1[2]);
		fprintf (drvui->fpoutp, " sphere{<%8.5f, %8.5f, %8.5f>, %8.5f\n",
			 vert2[0], vert2[1], vert2[2], radii[k]);
		fprintf (drvui->fpoutp, "  texture{pigment{color %s  }}\n", sphere_col_p);
		fprintf (drvui->fpoutp, "   finish{phong %5.2f phong_size %5.2f}\n",
			 drvui->Phong_Value, drvui->Phong_Size);
		fprintf (drvui->fpoutp, " }\n");
	    }
	    glPushMatrix ();
	    glLoadName (i);
	    glPushName (k);
	    glTranslatef (vert2[0], vert2[1], vert2[2]);
	    (void) sscanf (sphere_col_v, "%f %f %f", &glr, &glg, &glb);
	    glColor3f (glr, glg, glb);
	    glutSolidSphere (radii[k], 10, 10);
	    glPopName ();
	    glPopMatrix ();
	    if (doAsy) {
		fprintf (drvui->fpouta,
			 " // Sphere for %c%c%c%c at %8.5f %8.5f %8.5f  \n",
			 drvui->spheres[i].sphere_l[0], drvui->spheres[i].sphere_l[1],
			 drvui->spheres[i].sphere_l[2], drvui->spheres[i].sphere_l[3],
			 vert1[0], vert1[1], vert1[2]);
		fprintf (drvui->fpouta, "draw(pic, shift (%8.5f, %8.5f, %8.5f)*scale3(%.2f)*unitsphere,rgb(%4.2f,%4.2f,%4.2f));\n",
			 vert2[0], vert2[1], vert2[2], radii[k],glr,glg,glb);
	    }
	    if (doVrml) {
		if (no_comment == 0)
		    fprintf (drvui->fpoutv,
			     "# Sphere for %c%c%c%c at %8.5f %8.5f %8.5f \n",
			     drvui->spheres[i].sphere_l[0], drvui->spheres[i].sphere_l[1],
			     drvui->spheres[i].sphere_l[2], drvui->spheres[i].sphere_l[3],
			     vert1[0], vert1[1], vert1[2]);
		if (Vrml2) {
		    fprintf (drvui->fpoutv, "  Transform {\n");
		    fprintf (drvui->fpoutv, "   translation %8.5f %8.5f %8.5f\n",
			     vert2[0], vert2[1], vert2[2]);
		    fprintf (drvui->fpoutv, "   children [\n    Shape {\n");
		    fprintf (drvui->fpoutv,
			     "     appearance Appearance {material Material { diffuseColor %s} }\n",
			     sphere_col_v);
		    fprintf (drvui->fpoutv, "     geometry Sphere { radius %8.5f}\n",
			     radii[k]);
		    fprintf (drvui->fpoutv, "    }\n   ]\n  }\n");
		} else {	// VRML1
		    fprintf (drvui->fpoutv, " Separator {\n");
		    fprintf (drvui->fpoutv, "  Transform {");
		    fprintf (drvui->fpoutv, " translation %8.5f %8.5f %8.5f}\n", vert2[0],
			     vert2[1], vert2[2]);
		    fprintf (drvui->fpoutv, "  Material {\n");
		    fprintf (drvui->fpoutv, "   shininess 0.3\n");
		    fprintf (drvui->fpoutv, "   diffuseColor %s\n  }\n", sphere_col_v);
		    fprintf (drvui->fpoutv, "  Sphere { radius %8.5f}\n }\n", radii[k]);
		}
	    }
	}
    }
}

/* ************************************************************** */
/* ************************************************************** */

void
plot_vrml_poly (int polyno)
{
    int i, j, k, l, m, n;

    float *out;

    int *point, equal;

    char edgecolor[40];

/* get space for unique vertices and a list to point to them */

    if (!(out = (float *) zalloc ((unsigned) ((draw_list + 10) * 3 * sizeof (float))))) {
	Error_Box ("Unable to allocate space in plot_vrml_poly.");
	return;
    }
    if (!(point = (int *) zalloc ((unsigned) ((draw_list + 10) * sizeof (int))))) {
	free (out);
	Error_Box ("Unable to allocate space in plot_vrml_poly.");
	return;
    }
    if (!Vrml2)
	fprintf (drvui->fpoutv, "Coordinate3 { point [\n");
    n = 0;
    for (i = 0; i < draw_list; i++) {	/* build list of unique points */
	j = poly_list[i];
/*    fprintf(fpoutv,"# %d-sided polygon:\n",j); */
	for (k = i + 1; k <= i + j; k++) {
	    l = poly_list[k] - 1;
	    for (m = 0; m < 3; m++)
		out[3 * n + m] = s_vert[3 * l + m];
	    equal = 0;
	    for (m = 0; m < n; m++) {	/* check if new */
		if ((fabs (out[3 * m] - out[3 * n]) < 0.0002)
		    && (fabs (out[3 * m + 1] - out[3 * n + 1]) < 0.0002)
		    && (fabs (out[3 * m + 2] - out[3 * n + 2]) < 0.0002)) {
		    point[k] = m;	/* point previously found */
		    equal = 1;
		    break;
		}
	    }
	    if (equal == 0)
		point[k] = n++;	/* new point - save it by incrementing n */
	}
	i = i + j;
    }
    --n;
    if (Vrml2)
	fprintf (drvui->fpoutv,
		 "}\n geometry IndexedFaceSet { coord Coordinate{ point [\n");
    for (l = 0; l < n; l++)	/* output the unique set of coordinates */
	fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n", out[3 * l], out[3 * l + 1],
		 out[3 * l + 2]);
    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f]\n", out[3 * n], out[3 * n + 1],
	     out[3 * n + 2]);
    if (Vrml2) {
	fprintf (drvui->fpoutv, " }\n");
	fprintf (drvui->fpoutv, " coordIndex [\n");
    } else {
	fprintf (drvui->fpoutv, "}\n IndexedFaceSet { coordIndex [\n");
    }
    l = 0;
    for (i = 0; i < draw_list - 1; i++) {
	j = poly_list[i];
	l++;
	for (k = 0; k < j; k++)
	    fprintf (drvui->fpoutv, "%d,", point[l++]);
	if (i + j < draw_list - 2) {
	    fprintf (drvui->fpoutv, "-1,\n");
	} else {
	    if (Vrml2) {
		fprintf (drvui->fpoutv, "-1]\n");
	    } else {
		fprintf (drvui->fpoutv, "-1]\n  }\n");
	    }
	}
	i = i + j;
    }
    if (Vrml2) {
	fprintf (drvui->fpoutv, " solid FALSE\n");
	fprintf (drvui->fpoutv, " convex TRUE\n");
	fprintf (drvui->fpoutv, " creaseAngle 1.5708\n");
	fprintf (drvui->fpoutv, "  }\n }\n");
    }
    if (edges > 0) {
	if (drvui->polyhedra[polyno].poly_rad_edge != 0.)
	    strncpy (edgecolor, drvui->polyhedra[polyno].poly_col_edge, 39);
	else
	    strncpy (edgecolor, drvui->col_edge, 39);
	edgecolor[39] = 0;
	Transform_VRML_Color (edgecolor);
	if (Vrml2) {
	    fprintf (drvui->fpoutv, " Shape {\n");
	    fprintf (drvui->fpoutv,
		     "  geometry IndexedLineSet {\n   coord Coordinate{\n    point [\n");
	    for (l = 0; l < n; l++)	/* output the unique set of coordinates */
		fprintf (drvui->fpoutv, "     %5.3f %5.3f %5.3f,\n", out[3 * l],
			 out[3 * l + 1], out[3 * l + 2]);
	    fprintf (drvui->fpoutv, "     %5.3f %5.3f %5.3f]\n", out[3 * n],
		     out[3 * n + 1], out[3 * n + 2]);
	    fprintf (drvui->fpoutv, "  }\n   coordIndex [\n    ");
	} else {
	    fprintf (drvui->fpoutv, " Material {\n   diffuseColor %s \n }\n", edgecolor);
	    fprintf (drvui->fpoutv, "IndexedLineSet { coordIndex [\n    ");
	}
	l = 0;
	for (i = 0; i < draw_list - 1; i++) {
	    j = poly_list[i];
	    l++;
	    for (k = 0; k < j; k++)
		fprintf (drvui->fpoutv, "%d,", point[l++]);
	    if (i + j < draw_list - 2) {
		fprintf (drvui->fpoutv, "%d,-1,\n    ", point[l - j]);
	    } else {
		if (Vrml2) {
		    fprintf (drvui->fpoutv, "%d,-1]\n  ", point[l - j]);
		} else {
		    fprintf (drvui->fpoutv, "%d,-1]\n  }\n }\n", point[l - j]);
		}
	    }
	    i = i + j;
	}
	if (Vrml2) {
	    fprintf (drvui->fpoutv,
		     "  color Color { color [%s]}\n colorIndex[0]\n colorPerVertex FALSE}\n",
		     edgecolor);
	    fprintf (drvui->fpoutv, "  }\n");
	}
    } else {
	if (!Vrml2)
	    fprintf (drvui->fpoutv, " }\n");
    }
    free (out);			/* release the space */
    free (point);
    return;
}

/* ************************************************************** */
/* ************************************************************** */

void
print_sym (void)

/* print symmetry information */
{
    int i, j, k;

    switch (drvui->sys) {
    case 1:
	fprintf (drvui->flout, "Triclinic - ");
	break;
    case 2:
	fprintf (drvui->flout, "Monoclinic - ");
	break;
    case 3:
	fprintf (drvui->flout, "Orthorhombic - ");
	break;
    case 4:
	fprintf (drvui->flout, "Tetragonal - ");
	break;
    case 5:
	fprintf (drvui->flout, "Hexagonal - ");
	break;
    case 6:
	fprintf (drvui->flout, "Cubic - ");
    }
    if (!drvui->acentric) {
	fprintf (drvui->flout, "Centric - ");
    } else {
	fprintf (drvui->flout, "Acentric - ");
    }
    switch (drvui->nbr) {
    case 1:
	fprintf (drvui->flout, "P lattice");
	break;
    case 2:
	fprintf (drvui->flout, "A-centered");
	break;
    case 3:
	fprintf (drvui->flout, "B-centered");
	break;
    case 4:
	fprintf (drvui->flout, "C-centered");
	break;
    case 5:
	fprintf (drvui->flout, "F-centered");
	break;
    case 6:
	fprintf (drvui->flout, "I-centered");
	break;
    case 7:
	fprintf (drvui->flout, "R-centered");
	break;
    default:
	fprintf (drvui->flout, "**** Unknown Lattice ***");
    }
    fprintf (drvui->flout, " with %d Symmetry Operators\n\n", drvui->ng);
    for (i = 1; i <= drvui->ng; ++i) {

	for (j = 1; j <= 3; ++j) {
	    for (k = 1; k <= 3; ++k) {
		drvui->rss[i - 1][j - 1][k - 1] = (float) drvui->ss[i - 1][j - 1][k - 1];
		if (drvui->ss[i - 1][j - 1][k - 1] != 0) {
		    if (drvui->ss[i - 1][j - 1][k - 1] > 0) {
			fprintf (drvui->flout, "+");
		    } else {
			fprintf (drvui->flout, "-");
		    }
		    switch (k) {
		    case 1:
			fprintf (drvui->flout, "x");
			break;
		    case 2:
			fprintf (drvui->flout, "y");
			break;
		    case 3:
			fprintf (drvui->flout, "z");
		    }
		}
	    }
	    if (drvui->ts[i - 1][j - 1] != 0.0) {
		if (drvui->ts[i - 1][j - 1] == 0.5)
		    fprintf (drvui->flout, "+1/2");
		if (drvui->ts[i - 1][j - 1] == 0.25)
		    fprintf (drvui->flout, "+1/4");
		if (drvui->ts[i - 1][j - 1] == 0.75)
		    fprintf (drvui->flout, "+3/4");
		if (fabs (drvui->ts[i - 1][j - 1] - 0.33333333) < 0.00001)
		    fprintf (drvui->flout, "+1/3");
		if (fabs (drvui->ts[i - 1][j - 1] - 0.66666667) < 0.00001)
		    fprintf (drvui->flout, "+2/3");
		if (fabs (drvui->ts[i - 1][j - 1] - 0.16666667) < 0.00001)
		    fprintf (drvui->flout, "+1/6");
		if (fabs (drvui->ts[i - 1][j - 1] - 0.83333333) < 0.00001)
		    fprintf (drvui->flout, "+5/6");
	    }
	    if (j != 3)
		fprintf (drvui->flout, ",");
	}
	if (i != drvui->ng) {
	    if (i % 4 != 0) {
		fprintf (drvui->flout, " ; ");
	    } else {
		fprintf (drvui->flout, "\n");
	    }
	}
    }
    Conv_Sym_Mat ();
    fprintf (drvui->flout, "\n");
}				/* end of print_sym */
