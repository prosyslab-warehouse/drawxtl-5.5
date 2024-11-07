// $Id: DRAWxtl2.cxx 1114 2011-02-23 20:29:18Z martin $
//
// module drawxtl2.cxx - part of DRAWxtl V5.5
// Coded using the FLTK 1.1.6 widget set
//
//     Larry W. Finger, Martin Kroeker and Brian Toby
//
//
// This module contains the following routines:
//
//  check_dynamic_storage - adds to a dynamic array if needed
//  check_vert_alloc - checks to see if space for vertices
//  determinant - only for 3x3 matrices
//  generate_ellipsoids - builds the ellipsoid descriptions for the output lists
//  generate_gl_texts - adds text labels to openGL rendering
//  generate_planes - builds the plane descriptions for the output lists
//  generate_poly - builds the polyhedra descriptions for the output lists
//  generate_texts - adds text labels to the output lists (not openGL)
//  init_dynamic_storage - initializes the dynamic arrays
//  matinv - invert a 3 x 3 matrix
//
#include "drawxtl.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <FL/glut.H>

#include "draw_ext.h"

#include "DRAWxtl_proto.h"

/* ************************************************************** */
/* ************************************************************** */

void
check_dynamic_storage (void)
{

//arrows

    if (drvui->nmag >= drvui->nmag_alloc) {
	drvui->arrows =
	    (arrow_struct *) realloc (drvui->arrows,
				      (drvui->nmag_alloc +
				       20) * sizeof (struct arrow_struct));
	if (!drvui->arrows) {
	    Error_Box ("Unable to allocate storage for arrows.");
	    exit (0);
	}
	drvui->nmag_alloc += 20;
    }
// atoms

    if (natom >= drvui->atom_alloc) {
	drvui->atoms =
	    (atom_struct *) realloc (drvui->atoms,
				     (drvui->atom_alloc +
				      20) * sizeof (struct atom_struct));
	if (!drvui->atoms) {
	    Error_Box ("Unable to allocate storage for atoms.");
	    exit (0);
	}
	drvui->atom_alloc += 20;
	drvui->modulate_x =
	    (mod_x_struct *) realloc (drvui->modulate_x,
				      drvui->atom_alloc * drvui->mod_x_alloc *
				      sizeof (struct mod_x_struct));
	if (!drvui->modulate_x) {
	    Error_Box ("Unable to allocate storage for modulation vectors.");
	    exit (0);
	}
	drvui->modulate_3x =
	    (mod_3x_struct *) realloc (drvui->modulate_3x,
				       drvui->atom_alloc * 3 * drvui->mod_3x_alloc *
				       sizeof (struct mod_3x_struct));
	if (!drvui->modulate_3x) {
	    Error_Box ("Unable to allocate storage for modulation vectors.");
	    exit (0);
	}
    }
//bonds

    if (drvui->nbond >= drvui->nbond_alloc) {
	drvui->bonds = (bond_struct *) realloc (drvui->bonds, (drvui->nbond_alloc + 20) *
						sizeof (struct bond_struct));
	if (!drvui->bonds) {
	    Error_Box ("Unable to allocate storage for bonds.");
	    exit (0);
	}
	drvui->nbond_alloc += 20;
    }
//cones

    if (drvui->ncone >= drvui->ncone_alloc) {
	drvui->cones = (cone_struct *) realloc (drvui->cones, (drvui->ncone_alloc + 20) *
						sizeof (struct cone_struct));
	if (!drvui->cones) {
	    Error_Box ("Unable to allocate storage for cones.");
	    exit (0);
	}
	drvui->ncone_alloc += 20;
    }
// ellipsoids

    if (drvui->n_ellips >= drvui->ellips_alloc) {
	drvui->ellips =
	    (ellips_struct *) realloc (drvui->ellips,
				       (drvui->ellips_alloc +
					20) * sizeof (struct ellips_struct));
	if (!drvui->ellips) {
	    Error_Box ("Unable to allocate storage for ellipsoids.\n");
	    exit (0);
	}
	drvui->ellips_alloc += 20;
	drvui->modulate_3t =
	    (mod_3t_struct *) realloc (drvui->modulate_3t,
				       drvui->ellips_alloc * 3 * drvui->mod_3t_alloc *
				       sizeof (struct mod_3t_struct));
	if (!drvui->modulate_3t) {
	    Error_Box ("Unable to allocate storage for modulation vectors.");
	    exit (0);
	}
    }
// Fourier contours

    if (drvui->numOfFourierContours >= drvui->num_Fourier_alloc) {
	drvui->fourier =
	    (map_struct *) realloc (drvui->fourier,
				    (drvui->num_Fourier_alloc +
				     20) * sizeof (struct map_struct));
	if (!drvui->fourier) {
	    Error_Box ("Unable to allocate storage for Fourier contours.");
	    exit (0);
	}
	drvui->num_Fourier_alloc += 20;
    }
// frames

    if (drvui->max_frame >= drvui->frame_alloc) {
	drvui->frames = (frame_struct *) realloc (drvui->frames, (drvui->max_frame + 20) *
						  sizeof (struct frame_struct));
	if (!drvui->frames) {
	    Error_Box ("Unable to allocate storage for frames.");
	    exit (0);
	}
	drvui->max_frame += 20;
    }
//labels

    if (drvui->nlabel >= drvui->nlabel_alloc) {
	drvui->labels =
	    (label_struct *) realloc (drvui->labels,
				      (drvui->nlabel_alloc +
				       100) * sizeof (struct label_struct));
	if (!drvui->labels) {
	    Error_Box ("Unable to allocate storage for labels.");
	    exit (0);
	}
	drvui->nlabel_alloc += 100;
    }
//least-squares planes

    if (drvui->nbplane >= drvui->nbplane_alloc) {
	drvui->bplanes =
	    (bplane_struct *) realloc (drvui->bplanes,
				       (drvui->nbplane_alloc +
					20) * sizeof (struct bplane_struct));
	if (!drvui->bplanes) {
	    Error_Box ("Unable to allocate storage for least-squares planes");
	    exit (0);
	}
	drvui->nbplane_alloc += 20;
    }
// modulation

    if (drvui->no_mod_vectors >= drvui->mod_gbl_alloc) {
	drvui->modulate_gbl =
	    (mod_gbl_struct *) realloc (drvui->modulate_gbl,
					(drvui->mod_gbl_alloc +
					 20) * sizeof (struct mod_gbl_struct));
	if (!drvui->modulate_gbl) {
	    Error_Box ("Unable to allocate storage for modulation vectors.");
	    exit (0);
	}
	drvui->mod_gbl_alloc += 20;
    }
    if (drvui->no_site_occ >= drvui->mod_x_alloc) {
	drvui->modulate_x =
	    (mod_x_struct *) realloc (drvui->modulate_x,
				      (drvui->mod_x_alloc +
				       20) * drvui->atom_alloc *
				      sizeof (struct mod_x_struct));
	if (!drvui->modulate_x) {
	    Error_Box ("Unable to allocate storage for modulation vectors.");
	    exit (0);
	}
	drvui->mod_x_alloc += 20;
    }
    if (drvui->no_site_displace >= drvui->mod_3x_alloc) {
	drvui->modulate_3x =
	    (mod_3x_struct *) realloc (drvui->modulate_3x,
				       (drvui->mod_3x_alloc +
					20) * 3 * drvui->atom_alloc *
				       sizeof (struct mod_3x_struct));
	if (!drvui->modulate_3x) {
	    Error_Box ("Unable to allocate storage for modulation vectors.");
	    exit (0);
	}
	drvui->mod_3x_alloc += 20;
    }
    if (drvui->no_site_U_terms >= drvui->mod_3t_alloc) {
	drvui->modulate_3t =
	    (mod_3t_struct *) realloc (drvui->modulate_3t,
				       (drvui->mod_3t_alloc +
					20) * 3 * drvui->ellips_alloc *
				       sizeof (struct mod_3t_struct));
	if (!drvui->modulate_3t) {
	    Error_Box ("Unable to allocate storage for modulation vectors.");
	    exit (0);
	}
	drvui->mod_3t_alloc += 20;
    }
//planes

    if (drvui->nplane >= drvui->nplane_alloc) {
	drvui->planes =
	    (plane_struct *) realloc (drvui->planes,
				      (drvui->nplane_alloc +
				       20) * sizeof (struct plane_struct));
	if (!drvui->planes) {
	    Error_Box ("Unable to allocate storage for planes.");
	    exit (0);
	}
	drvui->nplane_alloc += 20;
    }
//polyhedra

    if (drvui->npoly >= drvui->npoly_alloc) {
	drvui->polyhedra =
	    (poly_struct *) realloc (drvui->polyhedra,
				     (drvui->npoly_alloc +
				      20) * sizeof (struct poly_struct));
	if (!drvui->polyhedra) {
	    Error_Box ("Unable to allocate storage for polyhedra.");
	    exit (0);
	}
	drvui->npoly_alloc += 20;
    }
//spheres

    if (drvui->nsphere >= drvui->nsphere_alloc) {
	drvui->spheres =
	    (sphere_struct *) realloc (drvui->spheres,
				       (drvui->nsphere_alloc +
					20) * sizeof (struct sphere_struct));
	if (!drvui->spheres) {
	    Error_Box ("Unable to allocate storage for spheres.");
	    exit (0);
	}
	drvui->nsphere_alloc += 20;
    }

    if (drvui->natprop >= drvui->natprop_alloc) {
	drvui->atprops =
	    (atprop_struct *) realloc (drvui->atprops,
				       (drvui->natprop_alloc +
					20) * sizeof (struct atprop_struct));
	if (!drvui->atprops) {
	    Error_Box ("Unable to allocate storage for atomic properties.");
	    exit (0);
	}
	drvui->natprop_alloc += 20;
    }
}

/* ************************************************************** */
/* ************************************************************** */

int
check_vert_alloc (int number, int extend_ok)
{
    if (number < drvui->verts_alloc)
	return 1;		/* return success */
    if (!extend_ok)
	return 0;		/* return failure for no space and no extend */
    return 0;
}

/* ************************************************************** */
/* ************************************************************** */

float
determinant (double rot[3][3])
{
    double temp;

    temp = rot[0][0] * (rot[1][1] * rot[2][2] - rot[1][2] * rot[2][1])
	- rot[1][0] * (rot[0][1] * rot[2][2] - rot[0][2] * rot[2][1])
	+ rot[2][0] * (rot[0][1] * rot[1][2] - rot[0][2] * rot[1][1]);
    return (float) temp;
}

/* ************************************************************** */
/* ************************************************************** */

void
generate_ellipsoids (void)

/* routine to generate ellipsoid objects and add then to edit list */
{
    float ev[3];		/* scaled eigenvalue */

    double rot[3][3];		/* rotation matrix for ellipsoid */

    int i, j, k, l, m, n, o, oo;	/* loop counters */
    int Vo[3] = { 0, 0, 0 };
    int n_ellipsoids;

    float EV[3][3];		/* array for eigenvectors */

    float alpha1, alpha2, alpha3;	/* terms needed to calculate rotations */

    float d1, d2, d3, alpha, sinalp, radius, temp;

    float glr, glg, glb;
    GLdouble clipx[4] = { 0., 0., 0., 0. };
    GLdouble clipy[4] = { 0., 0., 0., 0. };
    GLdouble clipz[4] = { 0., 0., 0., 0. };

    int outside = 0;
    int ntest[8][3] =
	{ {1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1}, {-1, -1, 1}, {1, -1, -1},
    {-1, 1, -1}, {-1, -1, -1}
    };				/* list of octant coordinates to check */
    double eye[3], vmin, dist;	/* used to find octant for cutout */

    int best;			/* octant list pointer */

    int omit, l1;

    char string[40];

    float rm[16];

    float cosal2;

    char ellips_col_p[40];

    char ellips_col_v[40];

    char ellipaxis_col_v[40];

    char ellipaxis_col_p[40];

    float uij[6];

    n_ellipsoids = 0;
    radius = 0.00015f * Scale;

    if (drvui->Ellipaxis_width != 0.)
	radius = drvui->Ellipaxis_width;
    fprintf (drvui->flout, "\nGenerated %6.1f percent probability ellipsoids for:\n\n",
	     100. * drvui->Ellipsoid_Prob);

    if (drvui->El_Cutout == 1) {	// if cutouts, update G_Rot
	glPushMatrix ();
	glLoadIdentity ();
	crystal->calculate (rm);	// get rotation matrix
	glPopMatrix ();
	for (m = 0; m < 3; m++) {	//  and copy to G_Rot
	    G_Rot[0][m] = rm[m];
	    G_Rot[1][m] = rm[m + 4];
	    G_Rot[2][m] = rm[m + 8];
	}
    }
    for (i = 1; i < drvui->n_ellips; ++i) {	// loop through ellipsoids
	if (drvui->ellips[i].ell_type < -1)
	    continue;		// This one is non-positive definite
	if (drvui->ellips[i].ell_type < 1000)
	    continue;		// ellipsoid OK but not displayed
	strcpy (ellips_col_p, drvui->ellips[i].ellips_col);
	strcpy (ellips_col_v, drvui->ellips[i].ellips_col);
	Transform_VRML_Color (ellips_col_v);
	Transform_POV_Color (ellips_col_p);
	strcpy (ellipaxis_col_p, drvui->Ellipaxis_color);
	strcpy (ellipaxis_col_v, drvui->Ellipaxis_color);
	Transform_VRML_Color (ellipaxis_col_v);
	Transform_POV_Color (ellipaxis_col_p);

	if (drvui->ellips[i].ellips_ismod == 0) {
	    for (j = 0; j < 3; j++) {
		ev[j] = drvui->ellips[i].ellips_RMS[j] * drvui->Ellipsoid_Scale;
		for (k = 0; k < 3; k++)
		    EV[j][k] = drvui->ellips[i].ellips_EV[j][k];
	    }
	}

	for (j = 0; j < natom; ++j) {	/* loop through atoms */
	    if (drvui->atoms[j].atom_fn != drvui->frame_no)
		continue;
	    nvert = 0;		/* clear the vertex list */
	    if (((drvui->atoms[j].atom_n >> 24) & 255) == i)	/* this ellipsoid */
		find_all_in_box (j);	/* generate all locations for this atom */
	    if (nvert > 0) {	/* if some present */
		if (drvui->do_ellipsoids == 0) {	/* Biso sphere */
		    float *radii;	/* allocate storage for radii of spheres */

		    radii = (float *) zalloc ((unsigned) (nvert * sizeof (float)));
		    for (k = 0; k < nvert; ++k)	/* copy radii info to storage */
			radii[k] =
			    drvui->ellips[i].ellips_RMS[0] * drvui->Ellipsoid_Scale;
		    Output_Spheres (radii, i);
		    free (radii);
		} else {
		    for (k = 0; k < nvert; ++k) {	/* loop through the atom positions */
			if (clipflag != 0) {

			    outside = 0;
			    for (l1 = 0; l1 < 3; l1++)
				if (o_vert[3 * k + l1] <
				    drvui->frames[drvui->frame_no].clip_lim[l1]
				    || o_vert[3 * k + l1] >
				    drvui->frames[drvui->frame_no].clip_lim[l1 + 3])
				    outside = 1;
			}

			omit = 0;
			for (l1 = 0; l1 < Omit->nomits; l1++) {
			    if (Omit->omit1[l1] == i * 100000 && Omit->omit2[l1] == k)
				omit = 1;
			}
			if ((!clipflag || !outside) && !omit) {
			    o = vert_sym_no[k];
			    fprintf (drvui->flout,
				     "  %c%c%c%c %1d at %8.5f %8.5f %8.5f with sym op %d\n",
				     drvui->atoms[j].atom_l[0], drvui->atoms[j].atom_l[1],
				     drvui->atoms[j].atom_l[2], drvui->atoms[j].atom_l[3],
				     drvui->atoms[j].sv_atom_n, o_vert[3 * k],
				     o_vert[3 * k + 1], o_vert[3 * k + 2], o);
			    if (no_comment == 0 && doVrml)
				fprintf (drvui->fpoutv,
					 "# Ellipsoid around %c%c%c%c %1d at %8.5f %8.5f %8.5f \n",
					 drvui->atoms[j].atom_l[0],
					 drvui->atoms[j].atom_l[1],
					 drvui->atoms[j].atom_l[2],
					 drvui->atoms[j].atom_l[3],
					 drvui->atoms[j].sv_atom_n, o_vert[3 * k],
					 o_vert[3 * k + 1], o_vert[3 * k + 2]);
			    if (doPOV)
				fprintf (drvui->fpoutp,
					 "/* Ellipsoid around %c%c%c%c %1d at %8.5f %8.5f %8.5f  */\n",
					 drvui->atoms[j].atom_l[0],
					 drvui->atoms[j].atom_l[1],
					 drvui->atoms[j].atom_l[2],
					 drvui->atoms[j].atom_l[3],
					 drvui->atoms[j].sv_atom_n, o_vert[3 * k],
					 o_vert[3 * k + 1], o_vert[3 * k + 2]);
			    if (doAsy)
				fprintf (drvui->fpouta,
					 " // Ellipsoid around %c%c%c%c %1d at %8.5f %8.5f %8.5f\n",
					 drvui->atoms[j].atom_l[0],
					 drvui->atoms[j].atom_l[1],
					 drvui->atoms[j].atom_l[2],
					 drvui->atoms[j].atom_l[3],
					 drvui->atoms[j].sv_atom_n, o_vert[3 * k],
					 o_vert[3 * k + 1], o_vert[3 * k + 2]);
				
			    if (drvui->ellips[i].ellips_ismod != 0) {
				int jj, kk;

				float vert[3];

				vert[0] = o_vert_nm[3 * k];
				vert[1] = o_vert_nm[3 * k + 1];
				vert[2] = o_vert_nm[3 * k + 2];
				if (modulate_uij (vert, i, j, o, uij) == 1)
				    continue;	// This one is non-positive definite
				for (jj = 0; jj < 3; jj++) {
				    ev[jj] =
					drvui->ellips[i].ellips_RMS[jj] *
					drvui->Ellipsoid_Scale;
				    for (kk = 0; kk < 3; kk++)
					EV[jj][kk] = drvui->ellips[i].ellips_EV[jj][kk];
				}
			    }
			    for (l = 0; l < 3; ++l) {
				for (m = 0; m < 3; ++m) {
				    rot[l][m] = 0.0;
				    for (n = 0; n < 3; ++n)
					rot[l][m] += drvui->rss[o][l][n] * EV[n][m];
				}	// for on m
			    }	// for on l

// make rotation right-handed (det > 0)

			    if (determinant (rot) < 0.0) {
				for (l = 0; l < 3; l++) {	// make rot right handed
				    for (m = 0; m < 3; m++) {
					rot[l][m] *= -1.0;
				    }
				}
			    }
			    if (fabs ((temp = determinant (rot)) - 1.0f) > 1.0e-3)
				fprintf (drvui->flout, "Determinant of rot = %8.5f\n",
					 temp);

// calculate ellipsoid rotation angles for POV and openGL

			    alpha2 = (float) asin (rot[2][0]);
			    cosal2 = (float) cos (alpha2);
			    if (fabs (cosal2) > 0.00001) {
				alpha1 =
				    (float) atan2 (-rot[2][1] / cosal2,
						   rot[2][2] / cosal2);
				alpha3 =
				    (float) atan2 (-rot[1][0] / cosal2,
						   rot[0][0] / cosal2);
			    } else {
				alpha3 = (float) atan2 (rot[0][1], -rot[0][2]);
				alpha1 = 0.0f;
			    }
			    if (drvui->El_Cutout == 1) {	// select octant for cutout
				double fullrot[3][3];

				vmin = -1000000.;	// initialize to impossible value
				best = 0;
				for (m = 0; m < 3; m++) {	// calculate combined rotation matrix
				    for (n = 0; n < 3; n++) {
					fullrot[m][n] = 0.0;
					for (oo = 0; oo < 3; oo++) {
					    fullrot[m][n] += rot[oo][m] * G_Rot[oo][n];
					}
				    }
				}
				if (fabs ((temp = determinant (fullrot)) - 1.0f) > 1.0e-3)
				    fprintf (drvui->flout,
					     "Determinant of fullrot = %8.5f\n", temp);

				eye[0] = fullrot[0][2];	// get the rotated position of the viewer (z axis)
				eye[1] = fullrot[1][2];
				eye[2] = fullrot[2][2];

				for (n = 0; n < 8; n++) {
				    for (m = 0; m < 3; m++)	// take next triple from list
					Vo[m] = ntest[n][m];

				    dist =
					eye[0] * Vo[0] + eye[1] * Vo[1] + eye[2] * Vo[2];
				    if (dist > vmin) {	// maximum is the one we want
					vmin = dist;
					best = n;
				    }
				}

				for (m = 0; m < 3; m++)	// set the octant to display
				    Vo[m] = ntest[best][m];

			    }	// ... if cutout

// output this ellipsoid to the POV file
			    if (doPOV) {
				fprintf (drvui->fpoutp, " object{\n");
				fprintf (drvui->fpoutp, "   union{\n");
				fprintf (drvui->fpoutp, "     object{\n");
				fprintf (drvui->fpoutp, "      cylinder{ < 0,0,-0.5 > , < 0,0,0.5 >, 1. open\n");	// first principal ellipse
				fprintf (drvui->fpoutp,
					 "       texture{pigment{color %s}}\n",
					 ellipaxis_col_p);
				fprintf (drvui->fpoutp,
					 "       scale<%8.5f,%8.5f,%8.5f>\n      }\n     }\n",
					 ev[0], ev[1], radius);
				fprintf (drvui->fpoutp, "     object{\n      cylinder{ < 0,-0.5,0 > , < 0,0.5,0 >, 1. open\n");	// second principal ellipse
				fprintf (drvui->fpoutp,
					 "       texture{pigment{color %s}}\n",
					 ellipaxis_col_p);
				fprintf (drvui->fpoutp,
					 "       scale<%8.5f,%8.5f,%8.5f>\n      }\n     }\n",
					 ev[0], radius, ev[2]);
				fprintf (drvui->fpoutp, "     object{\n      cylinder{ < -0.5,0,0 > , < 0.5,0,0 >, 1. open\n");	// third principal ellipse
				fprintf (drvui->fpoutp,
					 "       texture{pigment{color %s}}\n",
					 ellipaxis_col_p);
				fprintf (drvui->fpoutp,
					 "       scale<%8.5f,%8.5f,%8.5f>\n      }\n     }\n",
					 radius, ev[1], ev[2]);
				fprintf (drvui->fpoutp, "     object{\n");
				if (drvui->El_Cutout)
				    fprintf (drvui->fpoutp, "      difference {\n");
				fprintf (drvui->fpoutp, "       sphere{< 0,0,0 > 1}\n");
				if (drvui->El_Cutout) {
				    fprintf (drvui->fpoutp,
					     "       box{<0,0,0>,<%d,%d,%d>}\n", Vo[0],
					     Vo[1], Vo[2]);
				}
				fprintf (drvui->fpoutp,
					 "       texture{pigment{color %s  }}\n",
					 ellips_col_p);
				fprintf (drvui->fpoutp,
					 "       finish{phong %5.2f phong_size %5.2f reflection 0.1}\n",
					 drvui->Phong_Value, drvui->Phong_Size);
				fprintf (drvui->fpoutp,
					 "       scale<%8.5f,%8.5f,%8.5f>\n", ev[0],
					 ev[1], ev[2]);
				if (drvui->El_Cutout) {
				    fprintf (drvui->fpoutp, "      }\n");
				    fprintf (drvui->fpoutp, "     }\n");
				    fprintf (drvui->fpoutp, "     object{\n");
				    fprintf (drvui->fpoutp, "      cylinder{ < 0,0,-0.5 > , < 0,0,0.5 >, 0.99\n");	// first principal ellipse
				    fprintf (drvui->fpoutp,
					     "       texture{pigment{color %s}}\n",
					     drvui->Cutout_color);
				    fprintf (drvui->fpoutp,
					     "       scale<%8.5f,%8.5f,%8.5f>\n      }\n     }\n",
					     ev[0], ev[1], radius);
				    fprintf (drvui->fpoutp, "     object{\n      cylinder{ < 0,-0.5,0 > , < 0,0.5,0 >, 0.99\n");	// second principal ellipse
				    fprintf (drvui->fpoutp,
					     "       texture{pigment{color %s}}\n",
					     drvui->Cutout_color);
				    fprintf (drvui->fpoutp,
					     "       scale<%8.5f,%8.5f,%8.5f>\n      }\n     }\n",
					     ev[0], radius, ev[2]);
				    fprintf (drvui->fpoutp, "     object{\n      cylinder{ < -0.5,0,0 > , < 0.5,0,0 >, 0.99\n");	// third principal ellipse
				    fprintf (drvui->fpoutp,
					     "       texture{pigment{color %s}}\n",
					     drvui->Cutout_color);
				    fprintf (drvui->fpoutp,
					     "       scale<%8.5f,%8.5f,%8.5f>\n      }\n     }\n",
					     radius, ev[1], ev[2]);
				}
				fprintf (drvui->fpoutp, "      }\n");
				if (!drvui->El_Cutout)
				    fprintf (drvui->fpoutp, "     }\n");
				fprintf (drvui->fpoutp, "   rotate <%8.5f,0,0>\n",
					 -alpha1 * RAD);
				fprintf (drvui->fpoutp, "   rotate <0,%8.5f,0>\n",
					 -alpha2 * RAD);
				fprintf (drvui->fpoutp, "   rotate <0,0,%8.5f>\n",
					 -alpha3 * RAD);
				fprintf (drvui->fpoutp,
					 "   translate <%8.5f, %8.5f, %8.5f>\n }\n",
					 s_vert[3 * k], s_vert[3 * k + 1],
					 s_vert[3 * k + 2]);
			    }
// Output the ellipsoid to the VRML file

			    d1 = (float) (0.5 * (rot[2][1] - rot[1][2]));
			    d2 = (float) (0.5 * (rot[0][2] - rot[2][0]));
			    d3 = (float) (0.5 * (rot[1][0] - rot[0][1]));

// d's contain direction cosines * sin(alpha)

			    alpha = (float) (0.5 * (rot[0][0] + rot[1][1] +
						    rot[2][2] - 1.0f));
			    if (alpha < -1.0f) alpha = -1.0f;
			    if (alpha > 1.0f) alpha = 1.0f;

			    alpha = (float) acos (alpha);
			    sinalp = (float) sin (alpha);
			    if (fabs (sinalp) > 0.0001) {
				d1 = d1 / sinalp;
				d2 = d2 / sinalp;
				d3 = d3 / sinalp;

// Make sure alpha is in correct quadrant

				if (fabs (d1) > 0.0001) {
				    temp = (float) (1.0 - cos (alpha)) * d2 * d3 -
					d1 * (float) sin (alpha);
				    if (fabs (rot[1][2] - temp) > 0.0001)
					alpha = -alpha;
				} else {
				    if (fabs (d2) > 0.0001) {
					temp = (float) (1.0 - cos (alpha)) * d1 *
					    d3 - d2 * (float) sin (alpha);
					if (fabs (rot[2][0] - temp) > 0.0001)
					    alpha = -alpha;
				    } else {
					if (fabs (d3) > 0.0001) {
					    temp = (float) (1.0 - cos (alpha)) *
						d1 * d2 - d3 * (float) sin (alpha);
					    if (fabs (rot[0][1] - temp) > 0.0001)
						alpha = -alpha;
					}
				    }
				}
			    } else {
				d1 = (float) sqrt (fabs (rot[0][0]));
				d2 = (float) sqrt (fabs (rot[1][1]));
				d3 = (float) sqrt (fabs (rot[2][2]));
			    }	//sinalp != 0

			    temp = (float) sqrt (d1 * d1 + d2 * d2 + d3 * d3);
			    d1 = d1 / temp;
			    d2 = d2 / temp;
			    d3 = d3 / temp;

			    if (doVrml) {
				if (!Vrml2)
				    fprintf (drvui->fpoutv, " Separator {\n");
				fprintf (drvui->fpoutv, "  Transform {\n");
				fprintf (drvui->fpoutv,
					 "   rotation %6.3f %6.3f %6.3f %7.3f\n",
					 d1, d2, d3, -alpha);
				fprintf (drvui->fpoutv,
					 "   translation %8.5f %8.5f %8.5f\n",
					 s_vert[3 * k], s_vert[3 * k + 1],
					 s_vert[3 * k + 2]);
				if (strlen (ellips_col_v) == 0)
				    strcpy (ellips_col_v, "1 0 0");
				if (Vrml2) {
				    fprintf (drvui->fpoutv,
					     "   scale %8.5f %8.5f %8.5f\n", ev[0], ev[1],
					     ev[2]);
				    fprintf (drvui->fpoutv, "  children [ Shape {");
				    fprintf (drvui->fpoutv, " appearance Appearance {\n");
				    fprintf (drvui->fpoutv,
					     "  material Material {diffuseColor  %s}\n",
					     ellips_col_v);
				    fprintf (drvui->fpoutv, "  }\n");
				    fprintf (drvui->fpoutv, "  geometry Sphere {");
				    fprintf (drvui->fpoutv, "   radius     1.}}]}\n");
				} else {
				    fprintf (drvui->fpoutv,
					     "   scaleFactor %8.5f %8.5f %8.5f\n", ev[0],
					     ev[1], ev[2]);
				    fprintf (drvui->fpoutv, "  }\n");
				    fprintf (drvui->fpoutv, "  Material {\n");
				    fprintf (drvui->fpoutv, "   shininess 0.3\n");
				    fprintf (drvui->fpoutv, "   diffuseColor  %s\n",
					     ellips_col_v);
				    fprintf (drvui->fpoutv, "  }\n");
				    fprintf (drvui->fpoutv, "  Sphere {\n");
				    fprintf (drvui->fpoutv, "   radius     1\n  }\n");
				    fprintf (drvui->fpoutv, " }\n");
				}
			    }
			    (void) sscanf (ellips_col_v, "%f %f %f", &glr, &glg, &glb);
			    if (doAsy) {
				fprintf (drvui->fpouta," draw(pic, shift(%f,%f,%f)*rotate(%f,Z)*rotate(%f,Y)*rotate(%f,X)*\n",
				s_vert[3*k],s_vert[3*k+1],s_vert[3*k+2],-alpha3,-alpha2,-alpha1);
				fprintf (drvui->fpouta,"\tscale(%f,%f,%f)*unitsphere,rgb(%.2f,%.2f,%.2f));\n",
				    ev[0],ev[1],ev[2],glr,glg,glb);
			    }
// output ellipsoid to screen with openGL
			    glPushMatrix ();
			    glLoadName (i * 100000);
			    glPushName (k);
			    glColor3f (glr, glg, glb);
			    glTranslatef (s_vert[3 * k], s_vert[3 * k + 1],
					  s_vert[3 * k + 2]);
			    glRotatef (-alpha3 * (float) RAD, 0.0f, 0.0f, 1.0f);
			    glRotatef (-alpha2 * (float) RAD, 0.0f, 1.0f, 0.0f);
			    glRotatef (-alpha1 * (float) RAD, 1.0f, 0.0f, 0.0f);
			    glScalef (ev[0], ev[1], ev[2]);
			    if (drvui->El_Cutout) {
				clipx[0] = -(GLdouble) Vo[0];
				clipy[1] = -(GLdouble) Vo[1];
				clipz[2] = -(GLdouble) Vo[2];
				glClipPlane (GL_CLIP_PLANE0, clipx);
				glClipPlane (GL_CLIP_PLANE1, clipy);
				glClipPlane (GL_CLIP_PLANE2, clipz);
				glEnable (GL_CLIP_PLANE0);
				glutSolidSphere (1., 10, 10);
				glDisable (GL_CLIP_PLANE0);
				glEnable (GL_CLIP_PLANE1);
				glutSolidSphere (1., 10, 10);
				glDisable (GL_CLIP_PLANE1);
				glEnable (GL_CLIP_PLANE2);
				glutSolidSphere (1., 10, 10);
				glDisable (GL_CLIP_PLANE2);
				strncpy (string, drvui->Ellipaxis_color, 39);
				trim_string (string, 40);
				Transform_VRML_Color (string);
				(void) sscanf (string, "%f %f %f", &glr, &glg, &glb);
				glColor3f (glr, glg, glb);
				glDisable (GL_LIGHTING);
				glBegin (GL_LINES);
				glVertex3f (0.0f, 0.0f, 0.0f);
				glVertex3f ((float) Vo[0], 0.0f, 0.0f);
				glVertex3f (0.0f, 0.0f, 0.0f);
				glVertex3f (0.0f, (float) Vo[1], 0.0f);
				glVertex3f (0.0f, 0.0f, 0.0f);
				glVertex3f (0.0f, 0.0f, (float) Vo[2]);
				glEnd ();
				glEnable (GL_LIGHTING);
			    } else {
				glutSolidSphere (1., 10, 10);
			    }	// if drvui->El_Cutout
			    glPopName ();
			    glPopMatrix ();
			    n_ellipsoids++;
			}	// if clipflag == 0
		    }		// for on k
		}		// if ellipsoid 
	    }			// if nvert > 0
	}			// for all atoms j of this type 
    }				// for on i - loop through ellipsoids
    fprintf (drvui->fcns, "%4d ellipsoids\n", n_ellipsoids);
    fprintf (drvui->flout, "%4d ellipsoids generated\n", n_ellipsoids);
}				// end of generate_ellipsoids

/* ************************************************************** */
/* ************************************************************** */

void
generate_planes (void)

/* procedure to generate corners of planes */
{
    int i, j, k;		/* loop counters, etc */

    int nvert_start;		/* starting number for cations */

    int planeno;		/* number of current polyhedron */

    int Plane_Count;		/* number of planes output */

    int ii, jj, ll, l;

    int *no, *ns;

    float xm, ym, zm, vx, vy, vz, bm, vxs, vys, vzs, bms, phi, phi0, cosphi, radius;

    int outside;

    int omit, l1;

    float glr, glg, glb;

    static int warn = 0;

    char plane_col_p[40];

    char plane_col_v[40];

    char col_edge_p[40];

    char col_edge_v[40];

    if (drvui->nplane == 1)
	return;			/* exit NOW if no work */
    radius = drvui->rad_edge * (float) (0.005 * Scale);
    if (edges) {
	strcpy (col_edge_p, drvui->col_edge);
	strcpy (col_edge_v, drvui->col_edge);
	Transform_VRML_Color (col_edge_v);
	Transform_POV_Color (col_edge_p);
    }
    Plane_Count = 0;
    for (i = 0; i < 3; ++i)
	boxlim[i] += 3.0f;	/* expand box for cation generation */
    nvert = 0;			/* make list empty */
    for (i = 0; i < natom; ++i) {	/* loop through atoms */
	if (drvui->atoms[i].atom_fn != drvui->frame_no)
	    continue;
	find_all_in_box (i);	/* generate all potential cations */
    }
    for (i = 0; i < 3; ++i)
	boxlim[i] -= 3.0f;	/* shrink box to original size */
    if (nvert == 0)
	return;			/* if no anions, no polyhedra */
    nvert_start = nvert;	/* save start of cation list */
    for (i = 0; i < natom; ++i) {	/* start loop for cations */
	if (drvui->atoms[i].atom_fn != drvui->frame_no)
	    continue;
	nvert = nvert_start;	/* get rid of previous cations */
	planeno = (drvui->atoms[i].atom_n >> 16) & 255;	/* save plane number */
	if (planeno > 0 && drvui->planes[planeno].plane_fn == drvui->frame_no) {	/* greater than zero for polyhedral centers */
	    strcpy (plane_col_p, drvui->planes[planeno].plane_col);
	    strcpy (plane_col_v, drvui->planes[planeno].plane_col);
	    Transform_VRML_Color (plane_col_v);
	    Transform_POV_Color (plane_col_p);
	    find_all_in_box (i);	/* add all polyhedral centers in box */
	    if (nvert_start < nvert) {	/* check for some found */
		draw_list = 0;	/* initialize the draw list */
		for (j = nvert_start; j < nvert; ++j) {	/* loop through cations */
		    outside = 0;
		    if (clipflag != 0) {
			for (ii = 0; ii < 3; ii++) {
			    if (o_vert[3 * j + ii] <
				drvui->frames[drvui->frame_no].clip_lim[ii] - 0.01
				|| o_vert[3 * j + ii] >
				drvui->frames[drvui->frame_no].clip_lim[ii + 3] + 0.01)
				outside = 1;
			}
		    }
		    if (domolcomp && packflag) {
			for (ii = 0; ii < 3; ii++) {
			    if (o_vert[3 * j + ii] <
				drvui->frames[drvui->frame_no].cryst_lim[ii] - 0.01
				|| o_vert[3 * j + ii] >
				drvui->frames[drvui->frame_no].cryst_lim[ii + 3] + 0.01)
				outside = 1;
			}
		    }
		    numb_list = 0;	/* list clear */
		    for (k = 0; k < nvert_start; ++k) {	/* loop through anions */
			if (dist (j, k) <= 0.)
			    continue;
			if (dist (j, k) <= drvui->planes[planeno].plane_size + 0.0001 && !outside)	/* add vertex k to poly list */
			    vertex_list[numb_list++] = k;
		    }

		    omit = 0;
		    for (l1 = 0; l1 < Omit->nomits; l1++) {
			if (Omit->omit1[l1] == (i + 1) * 100000 + (i + 1) * 1000
			    && Omit->omit2[l1] == j)
			    omit = 1;
		    }

		    if (numb_list > 2 && omit == 0) {	/* only do if at least a triangle */
			fprintf (drvui->flout,
				 "Polygon with %d vertices around %c%c%c%c%3d at %8.5f %8.5f %8.5f\n",
				 numb_list, drvui->atoms[i].atom_l[0],
				 drvui->atoms[i].atom_l[1], drvui->atoms[i].atom_l[2],
				 drvui->atoms[i].atom_l[3], drvui->atoms[i].sv_atom_n,
				 s_vert[3 * j], s_vert[3 * j + 1], s_vert[3 * j + 2]);
			if (doVrml) {
			    if (no_comment == 0)
				fprintf (drvui->fpoutv,
					 "# Polygon (%d) around %c%c%c%c%3d at %8.5f %8.5f %8.5f \n",
					 numb_list, drvui->atoms[i].atom_l[0],
					 drvui->atoms[i].atom_l[1],
					 drvui->atoms[i].atom_l[2],
					 drvui->atoms[i].atom_l[3],
					 drvui->atoms[i].sv_atom_n, o_vert[3 * j],
					 o_vert[3 * j + 1], o_vert[3 * j + 2]);
			    if (Vrml2) {
				fprintf (drvui->fpoutv, " Shape {");
				fprintf (drvui->fpoutv, "appearance Appearance {\n");
				fprintf (drvui->fpoutv,
					 "   material Material {diffuseColor %s} \n \n",
					 plane_col_v);
			    } else {
				fprintf (drvui->fpoutv, " Separator {\n        ");
				fprintf (drvui->fpoutv, "  Material {\n");
				fprintf (drvui->fpoutv, "   diffuseColor %s \n  }\n",
					 plane_col_v);
			    }
			}
			Plane_Count++;	/* increment number output */
			poly_list[0] = numb_list;
/* Make vertices consecutive - find center of polygon */
			xm = ym = zm = 0.0f;
			if (!
			    (no =
			     (int *)
			     zalloc ((unsigned) ((4 * numb_list) * sizeof (int))))) {
			    Error_Box ("Unable to get no allocation");
			    return;
			}
			if (!
			    (ns =
			     (int *)
			     zalloc ((unsigned) ((4 * numb_list) * sizeof (int))))) {
			    Error_Box ("Unable to get ns allocation");
			    free (no);
			    return;
			}
			for (ii = 0; ii < numb_list; ii++) {
			    no[ii] = vertex_list[ii];
			    xm = xm + s_vert[3 * no[ii]];
			    ym = ym + s_vert[3 * no[ii] + 1];
			    zm = zm + s_vert[3 * no[ii] + 2];
			}
			xm = xm / (float) numb_list;
			ym = ym / (float) numb_list;
			zm = zm / (float) numb_list;

/* sort vertices in consecutive order based on minimal angles between vectors */

			for (ii = 0; ii < numb_list; ii++)
			    ns[ii] = 0;
			ns[0] = no[0];
			for (ii = 0; ii < numb_list - 1; ii++) {
			    phi0 = 1000.0f;
			    vxs = s_vert[3 * ns[ii]] - xm;
			    vys = s_vert[3 * ns[ii] + 1] - ym;
			    vzs = s_vert[3 * ns[ii] + 2] - zm;
			    bms = (float) sqrt (vxs * vxs + vys * vys + vzs * vzs);
			    for (jj = 0; jj < numb_list; jj++) {
				for (ll = 0; ll < ii; ll++)
				    if (ns[ll] == no[jj])
					jj++;
				if (ns[ii] != no[jj] && jj < numb_list) {
				    vx = s_vert[3 * no[jj]] - xm;
				    vy = s_vert[3 * no[jj] + 1] - ym;
				    vz = s_vert[3 * no[jj] + 2] - zm;
				    bm = (float) sqrt (vx * vx + vy * vy + vz * vz);
				    cosphi =
					(vxs * vx + vys * vy + vzs * vz) / (bms * bm);
				    if (cosphi < -1.)
					cosphi = -1.0f;
				    if (cosphi > 1.)
					cosphi = 1.0f;
				    phi = (float) acos (cosphi) * (float) RAD;
				    if (phi < phi0) {
					phi0 = phi;
					ns[ii + 1] = no[jj];
				    }
				}
			    }
			}
			draw_list = 1;
			poly_list[0] = numb_list;
			for (l = 0; l < numb_list; ++l)
			    poly_list[draw_list++] = ns[l] + 1;
			poly_list[draw_list++] = 0;	/* terminate polygon list */
			if (draw_list > 4 * NvertM) {
			    poly_list =
				(int *)
				zalloc ((unsigned) (4 * (nvert + 2) * sizeof (int)));
			    if (!warn) {
				fprintf (drvui->flout, "Overrun of poly_list.\n");
				fprintf (drvui->flout, "Overrun of poly_list.\n");
				Error_Box
				    ("Overrun of poly_list. Please send 'str' file\n"
				     "to Larry.Finger@@lwfinger.net.");
				warn = 1;
			    }
			}
			if (doVrml)
			    plot_vrml_poly (planeno);
			if (doPOV)
			    fprintf (drvui->fpoutp,
				     " /* Polygon (%d) around %c%c%c%c at %8.5f %8.5f %8.5f */ \n",
				     numb_list, drvui->atoms[i].atom_l[0],
				     drvui->atoms[i].atom_l[1], drvui->atoms[i].atom_l[2],
				     drvui->atoms[i].atom_l[3], o_vert[3 * j],
				     o_vert[3 * j + 1], o_vert[3 * j + 2]);
			if (doAsy)
			    fprintf (drvui->fpouta,
				     " // Polygon (%d) around %c%c%c%c at %8.5f %8.5f %8.5f  \n",
				     numb_list, drvui->atoms[i].atom_l[0],
				     drvui->atoms[i].atom_l[1], drvui->atoms[i].atom_l[2],
				     drvui->atoms[i].atom_l[3], o_vert[3 * j],
				     o_vert[3 * j + 1], o_vert[3 * j + 2]);
			glPushMatrix ();
			glDisable (GL_LIGHTING);
			(void) sscanf (plane_col_v, "%f %f %f", &glr, &glb, &glg);
			glColor3f (glr, glb, glg);
			glLoadName (100000 * (i + 1) + 1000 * (i + 1));
			glPushName (j);
			glBegin (GL_TRIANGLES);
			for (l = 0; l < numb_list - 2; ++l) {
			    if (doPOV) {
				fprintf (drvui->fpoutp,
					 " triangle {< %8.5f, %8.5f, %8.5f>,\n",
					 s_vert[3 * ns[0]], s_vert[3 * ns[0] + 1],
					 s_vert[3 * ns[0] + 2]);
				fprintf (drvui->fpoutp,
					 "<%8.5f, %8.5f, %8.5f>, <%8.5f, %8.5f, %8.5f>\n",
					 s_vert[3 * ns[l + 1]], s_vert[3 * ns[l + 1] + 1],
					 s_vert[3 * ns[l + 1] + 2], s_vert[3 * ns[l + 2]],
					 s_vert[3 * ns[l + 2] + 1],
					 s_vert[3 * ns[l + 2] + 2]);
				fprintf (drvui->fpoutp,
					 "  texture{pigment{color %s  }}\n", plane_col_p);
				fprintf (drvui->fpoutp, " }\n");
			    }
			    if (doAsy) {
				fprintf (drvui->fpouta,
					 " draw(pic, surface ( ( %8.5f, %8.5f, %8.5f)--",
					 s_vert[3 * ns[0]], s_vert[3 * ns[0] + 1],
					 s_vert[3 * ns[0] + 2]);
				fprintf (drvui->fpouta,
					 "( %8.5f, %8.5f, %8.5f)--(%8.5f, %8.5f, %8.5f)--cycle),",
					 s_vert[3 * ns[l + 1]], s_vert[3 * ns[l + 1] + 1],
					 s_vert[3 * ns[l + 1] + 2], s_vert[3 * ns[l + 2]],
					 s_vert[3 * ns[l + 2] + 1],
					 s_vert[3 * ns[l + 2] + 2]);
				fprintf (drvui->fpouta,"rgb( %4.2f,%4.2f,%4.2f) );\n",glr,glg,glb);
			    }
			    glVertex3f (s_vert[3 * ns[0]], s_vert[3 * ns[0] + 1],
					s_vert[3 * ns[0] + 2]);
			    glVertex3f (s_vert[3 * ns[l + 1]], s_vert[3 * ns[l + 1] + 1],
					s_vert[3 * ns[l + 1] + 2]);
			    glVertex3f (s_vert[3 * ns[l + 2]], s_vert[3 * ns[l + 2] + 1],
					s_vert[3 * ns[l + 2] + 2]);
			}
			glEnd ();
			glPopName ();
			glEnable (GL_LIGHTING);
			glPopMatrix ();
			if (edges) {
			    float df[3], at[3];

			    int m;

			    for (l = 0; l < numb_list - 1; ++l) {
				if (doPOV) {
				    fprintf (drvui->fpoutp,
					     " cylinder { <%8.5f, %8.5f, %8.5f>, <%8.5f, %8.5f, %8.5f> , %8.5f \n",
					     s_vert[3 * ns[l]], s_vert[3 * ns[l] + 1],
					     s_vert[3 * ns[l] + 2], s_vert[3 * ns[l + 1]],
					     s_vert[3 * ns[l + 1] + 1],
					     s_vert[3 * ns[l + 1] + 2], radius);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n",
					     col_edge_p);
				    fprintf (drvui->fpoutp, " }\n");
				}
				for (m = 0; m < 3; m++) {
				    at[m] = s_vert[3 * ns[l] + m];
				    df[m] = s_vert[3 * ns[l + 1] + m] - at[m];
				}
				push_cylinder (df, at, radius, col_edge_v);	// do the openGL edge
			    }
			    if (doPOV) {
				fprintf (drvui->fpoutp,
					 " cylinder { <%8.5f, %8.5f, %8.5f>, <%8.5f, %8.5f, %8.5f> , %8.5f \n",
					 s_vert[3 * ns[0]], s_vert[3 * ns[0] + 1],
					 s_vert[3 * ns[0] + 2],
					 s_vert[3 * ns[numb_list - 1]],
					 s_vert[3 * ns[numb_list - 1] + 1],
					 s_vert[3 * ns[numb_list - 1] + 2], radius);
				fprintf (drvui->fpoutp, " texture{pigment{color %s  }}\n",
					 col_edge_p);
				fprintf (drvui->fpoutp, " }\n");
			    }
			    for (m = 0; m < 3; m++) {
				at[m] = s_vert[3 * ns[0] + m];
				df[m] = s_vert[3 * ns[numb_list - 1] + m] - at[m];
			    }
			    push_cylinder (df, at, radius, col_edge_v);	// do the openGL edge
			}
			draw_list = 0;
			free (no);
			free (ns);
		    }		/* numb_list > 2 */
		}		/* loop on j through cations */
	    }			/* nvert_start < nvert */
	}
    }				/* end of loop through cations */
    fprintf (drvui->fcns, "%4d planes.\n", Plane_Count);
    fprintf (drvui->flout, "Generated %4d planes.\n", Plane_Count);
}

/* ************************************************************** */
/* ************************************************************** */

void
generate_gl_texts (void)
/* procedure to render all textual labels in the OpenGL view - this has
   to be done outside the display list so that the labels always face
   toward the viewer during rotations */
{
    int i, n;

    float vert[3];

    char *p;

    float scale;

    // if (nvert == 0) return; // happens when the mapread progress bar forces a 
    // (re)draw before the structure has been read

    display_cursor_text ();
    if (drvui->nlabel == 1)
	return;

    if (Display_axes) {
	nvert = 0;
	generate_triple ();
	draw_GL_triple ();
    }

    for (n = 1; n < drvui->nlabel; n++) {
	if (drvui->labels[n].label_fn != drvui->frame_no)
	    continue;
	nvert = 0;
	if (!Labels && strlen (drvui->labels[n].label_label) == 1) {
	    if (strncmp (drvui->labels[n].label_label, "a", 1))
		nvert = 1;
	    if (strncmp (drvui->labels[n].label_label, "b", 1))
		nvert = 1;
	    if (strncmp (drvui->labels[n].label_label, "c", 1))
		nvert = 1;
	    if (strncmp (drvui->labels[n].label_label, "o", 1))
		nvert = 1;
	    if (nvert == 1)
		continue;
	}
	nvert = 0;		/* initialize vertex list */
	for (i = 0; i < 3; i++) {
	    vert[i] = drvui->labels[n].label_x[i];
	    if (Display_axes && ((n == drvui->triple[0]) || (n == drvui->triple[1])
				 || (n == drvui->triple[2]) || (n == drvui->triple[3])))
		vert[i] += origin[i];
	}
	add_vert_nc (vert);

	if (Display_axes && ((n == drvui->triple[0]) || (n == drvui->triple[1])
			     || (n == drvui->triple[2]) || (n == drvui->triple[3])))
	    for (i = 0; i < 3; i++)
		s_vert[i] += offset[i];

	glPushMatrix ();
	glLoadName (n);
	glTranslatef (s_vert[0], s_vert[1], s_vert[2]);
	scale = 0.005f;
	if (!strcmp (drvui->labels[n].label_label, "triple_vect"))
	    scale *= 0.00001f;	// make "triple_vect" really small
	glRotatef (-(float) atof (drvui->X_Rot->value ()), 1.0f, 0.0f, 0.0f);
	glRotatef (-(float) atof (drvui->Y_Rot->value ()), 0.0f, 1.0f, 0.0f);
	glRotatef (-(float) atof (drvui->Z_Rot->value ()), 0.0f, 0.0f, 1.0f);
	glColor3f (0.0f, 0.0f, 0.0f);
	glDisable (GL_LIGHTING);
	glLineWidth (2.0);
	glScalef (drvui->label_scale * scale, drvui->label_scale * scale,
		  drvui->label_scale * scale);
	for (p = drvui->labels[n].label_label; *p; p++)
	    glutStrokeCharacter (GLUT_STROKE_ROMAN, *p);
	glLineWidth (1.0);
	glEnable (GL_LIGHTING);
	glPopMatrix ();
    }				/* end of loop on n   */
}

/* ************************************************************** */
/* ************************************************************** */

void
generate_poly (void)

/* procedure to generate vertices of polyhedra */
{
    int i, j, k;		/* loop counters, etc */

    int nvert_start;		/* starting number for cations */

    int polyno;			/* number of current polyhedron */

    int Poly_Count;		/* Number of polyhedra output */

    int Max_Vertices;		/* Maximum number of vertices */

    int ii, outside;

    int type_check;		/* if atom type checking needs to be done */

    int omit, l1;

    static int warn = 0;

    float c[3];

    char poly_col_v[40];

    float vert[3];

    if (drvui->npoly == 1)
	return;			/* exit NOW if no work */
    Poly_Count = 0;
    fprintf (drvui->flout,
	     "Polyhedra: DV = Distance Variation, QE = Quadratic Elongation\n"
	     " AV = Angle Variance (see Robinson et al. (1971), Science 172, 567.)\n"
	     "\n Atom   NoVert        Coords              Vol      DV       QE       AV\n");
    for (i = 0; i < 3; ++i)
	boxlim[i] += 3.0f;	/* expand box for anion generation */
    nvert = 0;			/* make list empty */
    for (i = 0; i < natom; ++i) {	/* loop through atoms */
	if (drvui->atoms[i].atom_fn != drvui->frame_no)
	    continue;
	find_all_in_box (i);	/* generate all potential anions */
    }
    for (i = 0; i < 3; ++i)
	boxlim[i] -= 3.0f;	/*shrink box to original size */
    if (nvert == 0)
	return;			/* if no anions, no polyhedra */
    nvert_start = nvert;	/* save start of cation list */
    for (i = 0; i < natom; ++i) {	/* start loop for cations */
	if (drvui->atoms[i].atom_fn != drvui->frame_no)
	    continue;
	nvert = nvert_start;	/* get rid of previous cations */
	type_check = 0;
	polyno = (drvui->atoms[i].atom_n >> 8) & 255;	// save polyhedron number
	if (polyno > 0 && drvui->polyhedra[polyno].poly_fn == drvui->frame_no) {	// greater than zero for polyhedral centers
	    strcpy (poly_col_v, drvui->polyhedra[polyno].poly_col);
	    Transform_VRML_Color (poly_col_v);
	    if (drvui->polyhedra[polyno].poly_t[0] != '\0')
		type_check = 1;
	    find_all_in_box (i);	/* add all polyhedral centers in box */
	    if (nvert_start < nvert) {	/* check for some found */
		draw_list = 0;	/* initialize the draw list */
		Max_Vertices = 0;
		for (j = nvert_start; j < nvert; ++j) {	/* first loop through cations */
		    numb_list = 0;	/* list clear */
		    outside = 0;
		    if (clipflag != 0) {
			for (ii = 0; ii < 3; ii++) {
			    if (o_vert[3 * j + ii] <
				drvui->frames[drvui->frame_no].clip_lim[ii] - 0.01
				|| o_vert[3 * j + ii] >
				drvui->frames[drvui->frame_no].clip_lim[ii + 3] + 0.01)
				outside = 1;
			}
		    }
		    if (outside == 1)
			continue;
		    for (k = 0; k < nvert_start; ++k) {	/* loop through anions */
			float d;

			int kk;

			kk = drvui->orig_atom_no[k];
			if ((d = dist (j, k)) <= drvui->polyhedra[polyno].poly_size) {	/* add vertex k to poly list */
			    if ((!type_check)
				||
				(check_atom_name
				 (drvui->atoms[kk].atom_l,
				  drvui->polyhedra[polyno].poly_t))) {
				if (d >= drvui->polyhedra[polyno].poly_min) {
				    numb_list++;
				}
			    }
			}
		    }
		    if (numb_list > Max_Vertices)
			Max_Vertices = numb_list;
		}		/* end of first pass through list */
		draw_list = 0;	/* initialize the draw list */
		for (j = nvert_start; j < nvert; ++j) {	/* second loop through cations */
		    numb_list = 0;	/* list clear */
		    for (k = 0; k < nvert_start; ++k) {	/* loop through anions */
			float d;

			int kk;

			kk = drvui->orig_atom_no[k];
			outside = 0;
			if (domolcomp && packflag) {
			    for (ii = 0; ii < 3; ii++) {
				if (o_vert[3 * k + ii] <
				    drvui->frames[drvui->frame_no].cryst_lim[ii] - 0.01
				    || o_vert[3 * k + ii] >
				    drvui->frames[drvui->frame_no].cryst_lim[ii + 3] +
				    0.01)
				    outside = 1;
			    }
			}
			if (clipflag != 0) {
			    for (ii = 0; ii < 3; ii++) {
				if (o_vert[3 * k + ii] <
				    drvui->frames[drvui->frame_no].clip_lim[ii] - 0.01
				    || o_vert[3 * k + ii] >
				    drvui->frames[drvui->frame_no].clip_lim[ii + 3] +
				    0.01)
				    outside = 1;
			    }
			}
			if ((d = dist (j, k)) <= drvui->polyhedra[polyno].poly_size + 0.000001 && !outside) {	/* add vertex k to poly list */
			    if (d >= drvui->polyhedra[polyno].poly_min) {
				if ((type_check == 0)
				    || check_atom_name (drvui->atoms[kk].atom_l,
							drvui->polyhedra[polyno].poly_t))
				    vertex_list[numb_list++] = k;
				if (numb_list > 4 * NvertM) {
				    vertex_list =
					(int *)
					zalloc ((unsigned)
						(4 * (nvert + 2) * sizeof (int)));
				    if (!warn) {
					fprintf (drvui->flout,
						 "Overrun of vertex_list.\n");
					fprintf (drvui->flout,
						 "Overrun of vertex_list.\n");
					Error_Box
					    ("Overrun of vertex_list.\nPlease send 'str' file "
					     "to Larry.Finger@@lwfinger.net.");
					warn = 1;
				    }
				}
			    }
			}
		    }

		    omit = 0;
		    for (l1 = 0; l1 < Omit->nomits; l1++) {
			if (Omit->omit1[l1] == (i + 1) * 100000 + (i + 1)
			    && Omit->omit2[l1] == j)
			    omit = 1;
		    }
		    for (l1 = 0; l1 < 3; l1++) {
			vert[l1] = o_vert[3 * j + l1];
			c[l1] = s_vert[3 * j + l1];
		    }
		    if (numb_list > 3 && numb_list == Max_Vertices && omit == 0) {	/* only do if polyhedron complete */
			fprintf (drvui->flout, " %c%c%c%c%2d%4d %8.5f %8.5f %8.5f",
				 drvui->atoms[i].atom_l[0], drvui->atoms[i].atom_l[1],
				 drvui->atoms[i].atom_l[2], drvui->atoms[i].atom_l[3],
				 drvui->atoms[i].sv_atom_n, numb_list, vert[0], vert[1],
				 vert[2]);
			if (doPOV)
			    fprintf (drvui->fpoutp,
				     " /* Polyhedron (%d) around %c%c%c%c%3d at %8.5f %8.5f %8.5f */ \n",
				     numb_list, drvui->atoms[i].atom_l[0],
				     drvui->atoms[i].atom_l[1], drvui->atoms[i].atom_l[2],
				     drvui->atoms[i].atom_l[3], drvui->atoms[i].sv_atom_n,
				     vert[0], vert[1], vert[2]);
			if (doAsy)
			    fprintf (drvui->fpouta,
				     " // Polyhedron (%d) around %c%c%c%c%3d at %8.5f %8.5f %8.5f  \n",
				     numb_list, drvui->atoms[i].atom_l[0],
				     drvui->atoms[i].atom_l[1], drvui->atoms[i].atom_l[2],
				     drvui->atoms[i].atom_l[3], drvui->atoms[i].sv_atom_n,
				     vert[0], vert[1], vert[2]);
			if (doVrml) {
			    if (no_comment == 0)
				fprintf (drvui->fpoutv,
					 "# Polyhedron (%d) around %c%c%c%c%3d at %8.5f %8.5f %8.5f \n",
					 numb_list, drvui->atoms[i].atom_l[0],
					 drvui->atoms[i].atom_l[1],
					 drvui->atoms[i].atom_l[2],
					 drvui->atoms[i].atom_l[3],
					 drvui->atoms[i].sv_atom_n, vert[0], vert[1],
					 vert[2]);
			    if (Vrml2) {
				fprintf (drvui->fpoutv, " Shape {");
				fprintf (drvui->fpoutv,
					 "appearance Appearance {\n  material Material { diffuseColor  %s} \n \n",
					 poly_col_v);
			    } else {
				fprintf (drvui->fpoutv, " Separator {\n    ");
				fprintf (drvui->fpoutv,
					 "  Material {\n   diffuseColor  %s \n  }\n",
					 poly_col_v);
			    }
			}
			build_poly_list (numb_list, polyno, i + 1, j, c);	/* select polygons */
			poly_list[draw_list++] = 0;	/* terminate polygon list */
			if (doVrml)
			    plot_vrml_poly (polyno);	/* pov and openGL part already done in build_poly_list */
			draw_list = 0;
			Poly_Count++;
		    }
		}
	    }			/* nvert_start < nvert */
	}
    }				/* end of loop through cations */
    fprintf (drvui->fcns, "%4d polyhedra.\n", Poly_Count);
    fprintf (drvui->flout, "\nGenerated %4d polyhedra.\n", Poly_Count);
}

/* ************************************************************** */
/* ************************************************************** */

void
generate_texts (void)

/* procedure to create POV and VRML command sequences for rendering of
   user-defined labels (labeltext command) - the GL rendering is done
   in generate_gl_texts */
{
    int i, n;

    int Label_Count = 0;

    float vert[3];

    if (drvui->nlabel == 1)
	return;

    for (n = 1; n < drvui->nlabel; n++) {

	float size;

	if (drvui->labels[n].label_fn != drvui->frame_no)
	    continue;
	if (drvui->triple[0] == n)
	    continue;		// skip "vector_trip"

	nvert = 0;
	if (!Labels && strlen (drvui->labels[n].label_label) == 1) {
	    if (strncmp (drvui->labels[n].label_label, "a", 1))
		nvert = 1;
	    if (strncmp (drvui->labels[n].label_label, "b", 1))
		nvert = 1;
	    if (strncmp (drvui->labels[n].label_label, "c", 1))
		nvert = 1;
	    if (strncmp (drvui->labels[n].label_label, "o", 1))
		nvert = 1;
	    if (nvert == 1)
		continue;
	}
	nvert = 0;		/* initialize vertex list */
	for (i = 0; i < 3; i++) {
	    vert[i] = drvui->labels[n].label_x[i];
	    if (Display_axes && ((n == drvui->triple[0]) || (n == drvui->triple[1])
				 || (n == drvui->triple[2]) || (n == drvui->triple[3])))
		vert[i] += origin[i];
	}
	add_vert_nc (vert);

	if (Display_axes && ((n == drvui->triple[0]) || (n == drvui->triple[1])
			     || (n == drvui->triple[2]) || (n == drvui->triple[3])))
	    for (i = 0; i < 3; i++)
		s_vert[i] += offset[i];
	Label_Count++;
	fprintf (drvui->flout, "Label[%d] %g %g %g %s\n", n, vert[0], vert[1], vert[2],
		 drvui->labels[n].label_label);
	if (doPOV)
	    fprintf (drvui->fpoutp, "/* Labels */\n");
	if (doVrml)
	    fprintf (drvui->fpoutv, "# Labels\n");
	size = drvui->label_scale * 8.0f * Text_Size;
	if (doPOV) {
	    fprintf (drvui->fpoutp, "  text { ttf \"crystal.ttf\",\"%s\" 0.15,0\n",
		     drvui->labels[n].label_label);
	    fprintf (drvui->fpoutp, "    scale <%4.2f,%4.2f,%4.2f>\n", size, size, size);
	    fprintf (drvui->fpoutp, "    rotate <0, 0, -zrot>\n");
	    fprintf (drvui->fpoutp, "    rotate <0, -yrot, 0>\n");
	    fprintf (drvui->fpoutp, "    rotate <-xrot,0, 0>\n");
	    fprintf (drvui->fpoutp,
		     "    translate <%8.5f,%8.5f,%8.5f> pigment{color Black}\n",
		     s_vert[0], s_vert[1], s_vert[2]);
	    fprintf (drvui->fpoutp, "  }\n");
	}
	if (doAsy) {
	    fprintf (drvui->fpouta, " label(pic, \"%s\",(%8.5f,%8.5f,%8.5f));\n",
		drvui->labels[n].label_label, s_vert[0], s_vert[1], s_vert[2]);
	}
	if (doVrml) {
	    if (Vrml2) {
		fprintf (drvui->fpoutv, "   Transform {");
		fprintf (drvui->fpoutv, " translation %8.5f %8.5f %8.5f \n",
			 s_vert[0], s_vert[1], s_vert[2]);
		fprintf (drvui->fpoutv, " rotation 1 0 0 %f\n", -xrot / RAD);
		fprintf (drvui->fpoutv, "children Transform { rotation 0 1 0 %f\n",
			 -yrot / RAD);
		fprintf (drvui->fpoutv, "children Transform { rotation 0 0 1 %f\n",
			 -zrot / RAD);
		fprintf (drvui->fpoutv, "  children [ \n");
		fprintf (drvui->fpoutv, "   Billboard {\n    axisOfRotation 0 1 0\n");
		fprintf (drvui->fpoutv, "    children [ \n");
		fprintf (drvui->fpoutv, "     Shape {\n");
		fprintf (drvui->fpoutv, "      appearance Appearance{\n");
		fprintf (drvui->fpoutv,
			 "       material Material {diffuseColor 0 0 0}\n");
		fprintf (drvui->fpoutv, "      }\n");
		fprintf (drvui->fpoutv, "      geometry Text{string [\"%s\"]}\n",
			 drvui->labels[n].label_label);
		fprintf (drvui->fpoutv, "     }\n    ]\n");
		fprintf (drvui->fpoutv, "   }\n  ]\n }}}\n");
	    } else {
		fprintf (drvui->fpoutv, " Separator {\n  Translation {");
		fprintf (drvui->fpoutv, " translation %8.5f %8.5f %8.5f }\n",
			 s_vert[0], s_vert[1], s_vert[2]);
		fprintf (drvui->fpoutv, "  Rotation{ rotation 0 0 1 %8.5f}\n",
			 -zrot / RAD);
		fprintf (drvui->fpoutv, "  Rotation{ rotation 0 1 0 %8.5f}\n",
			 -yrot / RAD);
		fprintf (drvui->fpoutv, "  Rotation{ rotation 1 0 0 %8.5f}\n",
			 -xrot / RAD);
		fprintf (drvui->fpoutv, "  Scale{ scaleFactor %4.2f %4.2f %4.2f }\n",
			 Text_Size, Text_Size, Text_Size);
		fprintf (drvui->fpoutv, "  Material{ diffuseColor 0. 0. 0. }\n");
		fprintf (drvui->fpoutv, "  AsciiText { string  \"%s\" }\n",
			 drvui->labels[n].label_label);
		fprintf (drvui->fpoutv, " }\n");
	    }
	}
    }				/* end of loop on n   */
    fprintf (drvui->fcns, "%4d labels.\n", Label_Count);
    fprintf (drvui->flout, "Generated %4d labels.\n", Label_Count);
}

/* ************************************************************** */
/* ************************************************************** */

void
init_dynamic_storage ()
{

// atoms

    drvui->verts_alloc = MAX_VERTS;
    drvui->atoms = (atom_struct *) zalloc (1000 * sizeof (struct atom_struct));
    if (!drvui->atoms) {
	printf ("\n*** Unable to allocate storage for atoms.\n");
	exit (0);
    }
    drvui->atom_alloc = 1000;

//arrows

    drvui->arrows = (arrow_struct *) zalloc (20 * sizeof (struct arrow_struct));
    if (!drvui->arrows) {
	printf ("\n*** Unable to allocate storage for arrows.\n");
	exit (0);
    }
    drvui->nmag_alloc = 20;
    drvui->nmag = 0;

//bonds

    drvui->bonds = (bond_struct *) zalloc (20 * sizeof (struct bond_struct));
    if (!drvui->bonds) {
	printf ("\n*** Unable to allocate storage for bonds.\n");
	exit (0);
    }
    drvui->nbond_alloc = 20;
    drvui->nbond = 0;

//cones

    drvui->cones = (cone_struct *) zalloc (20 * sizeof (struct cone_struct));
    if (!drvui->cones) {
	printf ("\n*** Unable to allocate storage for cones.\n");
	exit (0);
    }
    drvui->ncone_alloc = 20;
    drvui->ncone = 0;

// ellipsoids

    drvui->ellips = (ellips_struct *) zalloc (200 * sizeof (struct ellips_struct));
    if (!drvui->ellips) {
	printf ("\n*** Unable to allocate storage for ellipsoids.\n");
	exit (0);
    }
    drvui->ellips_alloc = 200;
    drvui->n_ellips = 0;

// Fourier map contours

    drvui->fourier = (struct map_struct *) zalloc (20 * sizeof (struct map_struct));
    if (!drvui->fourier) {
	printf ("\n*** Unable to allocate storage for Fourier contours.\n");
	exit (0);
    }
    drvui->num_Fourier_alloc = 20;

// frame data

    drvui->frames = (struct frame_struct *) zalloc (20 * sizeof (struct frame_struct));
    if (!drvui->frames) {
	printf ("\n*** Unable to allocate storage for frames.\n");
	exit (0);
    }
    drvui->frame_alloc = 20;


//labels

    drvui->labels = (struct label_struct *) zalloc (100 * sizeof (struct label_struct));
    if (!drvui->labels) {
	printf ("\n*** Unable to allocate storage for labels.\n");
	exit (0);
    }
    drvui->nlabel_alloc = 100;

//least-squares planes

    drvui->bplanes = (bplane_struct *) zalloc (20 * sizeof (struct bplane_struct));
    if (!drvui->bplanes) {
	printf ("\n*** Unable to allocate storage for least-squares planes.\n");
	exit (0);
    }
    drvui->nbplane_alloc = 20;

//planes

    drvui->planes = (plane_struct *) zalloc (20 * sizeof (struct plane_struct));
    if (!drvui->planes) {
	printf ("\n*** Unable to allocate storage for planes.\n");
	return;
    }
    drvui->nplane_alloc = 20;

//polyhedra

    drvui->polyhedra = (poly_struct *) zalloc (20 * sizeof (struct poly_struct));
    if (!drvui->polyhedra) {
	printf ("\n*** Unable to allocate storage for polyhedra.\n");
	exit (0);
    }
    drvui->npoly_alloc = 20;

    drvui->polyedges = (edge_struct *) zalloc (20 * sizeof (struct edge_struct));
    if (!drvui->polyedges) {
	printf ("\n*** Unable to allocate storage for polyhedron edge parameters.\n");
	exit (0);
    }
    drvui->nedge_alloc = 20;

//spheres

    drvui->spheres = (sphere_struct *) zalloc (20 * sizeof (struct sphere_struct));
    if (!drvui->spheres) {
	printf ("\n*** Unable to allocate storage for spheres.\n");
	exit (0);
    }
    drvui->nsphere_alloc = 20;

// global modulation parameters

    drvui->modulate_gbl =
	(struct mod_gbl_struct *) zalloc (20 * sizeof (struct mod_gbl_struct));
    if (!drvui->modulate_gbl) {
	printf ("\n*** Unable to allocate storage for global modulation parameters.\n");
	exit (0);
    }
    drvui->mod_gbl_alloc = 20;
    drvui->no_mod_vectors = 0;

// X modulation parameters

    drvui->modulate_x =
	(struct mod_x_struct *) zalloc (drvui->mod_gbl_alloc * drvui->atom_alloc *
					sizeof (struct mod_x_struct));
    if (!drvui->modulate_x) {
	printf ("\n*** Unable to allocate storage for X modulation parameters.\n");
	exit (0);
    }
    drvui->mod_x_alloc = 20;
    drvui->no_site_occ = 0;

// 3X modulation parameters

    drvui->modulate_3x =
	(struct mod_3x_struct *) zalloc (drvui->mod_gbl_alloc * drvui->atom_alloc *
					 sizeof (struct mod_3x_struct));
    if (!drvui->modulate_3x) {
	printf ("\n*** Unable to allocate storage for 3X modulation parameters.\n");
	exit (0);
    }
    drvui->mod_3x_alloc = 20;
    drvui->no_site_displace = 0;

// 3T modulation parameters

    drvui->modulate_3t =
	(struct mod_3t_struct *) zalloc (drvui->mod_gbl_alloc * drvui->ellips_alloc *
					 sizeof (struct mod_3t_struct));
    if (!drvui->modulate_3t) {
	printf ("\n*** Unable to allocate storage for 3T modulation parameters.\n");
	exit (0);
    }
    drvui->mod_3t_alloc = 20;
    drvui->no_site_U_terms = 0;

    drvui->atprops = (atprop_struct *) zalloc (20 * sizeof (struct atprop_struct));
    if (!drvui->atprops) {
	printf ("\n*** Unable to allocate storage for atomic parameters.\n");
	exit (0);
    }
    drvui->natprop_alloc = 20;

}

/* ************************************************************** */
/* ************************************************************** */

float
matinv (float a[3][3])
{
/* invert a 3 x 3 float matrix - a^(-1) replaces a,
   returned value of determinant is zero if matrix is singular */
    float det;

    double ad[3][3];

    float ainv[3][3];

    int i, j;

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    ad[i][j] = a[i][j];
    if (!(det = determinant (ad)))
	return 0.0f;		/* a is singular */
    ainv[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
    ainv[0][1] = -(a[0][1] * a[2][2] - a[0][2] * a[2][1]) / det;
    ainv[0][2] = (a[0][1] * a[1][2] - a[1][1] * a[0][2]) / det;
    ainv[1][0] = -(a[1][0] * a[2][2] - a[1][2] * a[2][0]) / det;
    ainv[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det;
    ainv[1][2] = -(a[0][0] * a[1][2] - a[1][0] * a[0][2]) / det;
    ainv[2][0] = -(a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
    ainv[2][1] = -(a[0][0] * a[2][1] - a[2][0] * a[0][1]) / det;
    ainv[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    a[i][j] = ainv[i][j];
    return det;
}

/* ************************************************************** */
/* ************************************************************** */

void
matmul (float a[3][3], float b[3][3], float c[3][3])
{
/* Multiply matrix a * b and return the result in c. Matrix c may be
 * the same as a or b */
    float tmp[3][3];

    int i, j, k;

    for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++) {
	    tmp[i][j] = 0.0;
	    for (k = 0; k < 3; k++) {
		tmp[i][j] += a[i][k] * b[k][j];
	    }
	}
    }
/* copy tmp to c */
    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    c[i][j] = tmp[i][j];
}
