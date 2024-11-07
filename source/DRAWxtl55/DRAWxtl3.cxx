// $Id: DRAWxtl3.cxx 1114 2011-02-23 20:29:18Z martin $
//
// module drawxtl3.cxx - part of DRAWxtl V5.5
// Coded using the FLTK 1.1.6 widget set
//
//     Larry W. Finger, Martin Kroeker and Brian Toby
//
//
// This module contains the following routines:
//
//  build_poly_list - generate the list of polygon vertices that describe the polyhedron
//  draw_cell - draw the unit-cell boundaries
//  draw_GL_triple - draw the axes for the vector triple in openGL
//  generate_spheres - output the list of spheres for the display lists 
//  generate_triple - generate the vertex locations for the triple vector
//  push_cylinder - routine to push an openGL cylinder that represents a unit-cell edge
//  generate_lsq_planes - generate best fitting planes 
//  generate_aimsurf - generate surfaces for AIM datasets
//  calculate_voids - calculate cavity volumes for rendering
//  generate_voids - draw cavity volumes
//  calc_simplevoids - calculate cavity volume by iterative voxel tests
//  generate_simplevoids - draw cavity voxels
//  calculate_sas - calculate and draw solvent accessible surfaces using random sampling
//  calculate_msms - calculate accessible surface using Sanners' msms program 
//  generate_msms - draw solvent accessible surface from Sanners' msms program 

#include "drawxtl.h"
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
#define rint(x) (float)(int)(x < 0 ? x - 0.5 : x + 0.5)
#endif
#include "draw_ext.h"

#include "DRAWxtl_proto.h"

/* ************************************************************** */
/* ************************************************************** */

void
build_poly_list (int numb_list, int polyno, int i0, int i1, float cx[3])

/* routine to take list of 'numb_list' indices in 's_vert'
   and find each plane on surface of convex polyhedron.  It then builds
   a list of pointers to these vertices for drawing program.  The algorithm
   is as follows:
   
   1.  If there are 4 vertices, the polyhedron is a tetrahedron and all
       combinations of 3 vertices constitute a face.

   2.  For more than 4 vertices, some planes will be interior.  These are
       found by choosing three points (indicated by i,j,k) and forming two
       vectors defined by V1 = P(j) - P(i) and V2 = P(k) - P(i).  For each
       additional point (l), a third vector V3 = P(l) - P(i) is formed and
       the volume enclosed by these three vectors V = V3 . V1 x V2 is
       calculated.  If this volume is zero, point l lies in the plane.  If
       the plane is on the surface of the polyhedron, all non-zero volumes
       must have the same sign, which may be positive or negative.
*/
/*  numb_list - number of candidate vertices */
/*  polyno - number of polygon - for polygon color */
/*  i0 - pointer to central atom */
/*  i1 - pointer to an atom in the polyhedron */
/*  cx - orthogonal coordinates of center of polyhedron */
{
    int i, j, k, l, m, n;	/* local loop counters */

    float v0[3], v1[3], v2[3], v3[3], v4[3];	/* vectors used in calculation */

    float p[3];			/* coordinates of point i */

    float vol, sign_z;

    int ind;			/* indicator of surface plane */

    int *no;			/* place to keep vertex id's */

    int imk, jmk, kmk;

    int ii, jj, ll;

    int jplane;

    double phi, cosphi;

    float xm, ym, zm, vx, vy, vz, bm, phi0, vxs, vys, vzs, bms;

    float d1, d2, d3, radius;

    int *ns;

    char edgecolor[40];

    static int warn = 0;

    float glr = 0.0f, glg = 0.0f, glb = 0.0f;

    float Vol = 0.0f;		/* volume of polyhedron */

    float Area, Z0, Z;

    int ih[3];

    float sumth = 0.0f, sumth2 = 0.0f, sumd2 = 0.0f, sumd = 0.0f;;
    int nangle = 0;

    int nm, mn;

    float va[3], vb[3], temp, c1, c2 = 0.0f, vlo, qe = 0.0f, siga = 0.0f, sigd = 0.0f;

    char poly_col_p[40];

    char poly_col_v[40];

    strcpy (poly_col_p, drvui->polyhedra[polyno].poly_col);
    strcpy (poly_col_v, poly_col_p);
    Transform_VRML_Color (poly_col_v);
    Transform_POV_Color (poly_col_p);

    if (!(no = (int *) zalloc ((unsigned) ((4 * numb_list) * sizeof (int))))) {
	Error_Box ("Unable to get no allocation!");
	return;
    }
    if (!(ns = (int *) zalloc ((unsigned) ((4 * numb_list) * sizeof (int))))) {
	Error_Box ("Unable to get ns allocation!");
	free (no);
	return;
    }
    radius = drvui->rad_edge * 0.005f * Scale;
    memset (no, 0, sizeof (no));
    memset (ns, 0, sizeof (ns));

    glLoadName (i0 * 100000 + i0);
    glPushName (i1);

    if (drvui->polyhedra[polyno].poly_rad_edge != 0.)
	radius = drvui->polyhedra[polyno].poly_rad_edge * 0.005f * Scale;

    Vol = 0.0f;
    for (i = 0; i < numb_list - 2; i++) {
	ih[0] = i;
	for (j = i + 1; j < numb_list - 1; j++) {
	    ih[1] = j;
	    for (k = j + 1; k < numb_list; k++) {
		ih[2] = k;
		for (m = 0; m < 3; m++) {
		    v1[m] =
			s_vert[3 * vertex_list[j] + m] - s_vert[3 * vertex_list[i] + m];
		    v2[m] =
			s_vert[3 * vertex_list[k] + m] - s_vert[3 * vertex_list[i] + m];
		    v4[m] = s_vert[3 * vertex_list[i] + m] - cx[m];
		}
		vcross (v1, v2, v3);
		Area = 0.5f * vlength (v3);
		Z0 = 0.5f * vdot (v3, v4) / Area;
		if (fabs (Z0) < 0.001)
		    continue;	/* 3 points are in a line */
		for (l = 0; l < numb_list; l++) {
		    if (l == i || l == j || l == k)
			continue;
		    for (m = 0; m < 3; m++) {
			v0[m] =
			    s_vert[3 * vertex_list[i] + m] - s_vert[3 * vertex_list[l] +
								    m];
		    }
		    Z = 0.5f * vdot (v0, v3);
		    if (Z * Z0 < 0.0f)
			goto next_ver;	/* not all points on same side */
		    for (m = 0; m < 3; m++) {
			v3[m] *= 0.5f / Area;	/* directions cosines of plane normal */
		    }
		}
		Vol += Area * (float) (fabs (Z0) / 3.0);
		for (l = 0; l < 2; l++) {
		    nm = ih[l];
		    for (m = l + 1; m < 3; m++) {
			mn = ih[m];
			for (n = 0; n < 3; n++) {
			    va[n] = s_vert[3 * vertex_list[nm] + n] - cx[n];
			    vb[n] = s_vert[3 * vertex_list[mn] + n] - cx[n];
			}
			temp = vdot (va, vb) / (vlength (va) * vlength (vb));
			temp = min (temp, 1.0f);
			temp = max (-1.0f, temp);
			temp =
			    (float) (atan2 (sqrt (1.0f - temp * temp), temp) * 180.0 /
				     PI);
			sumth += temp;
			sumth2 += temp * temp;
			nangle++;
		    }
		}
	      next_ver:;
	    }
	}
    }
    for (i = 0; i < numb_list; i++) {
	for (n = 0; n < 3; n++) {
	    va[n] = s_vert[3 * vertex_list[i] + n] - cx[n];
	}
	temp = vlength (va);
	sumd2 += temp * temp;
	sumd += temp;

    }
    if (numb_list == 4) {
	c1 = 9.0f * (float) (sqrt (3.0) / 8.0);
	c2 = 109.471f;
    } else if (numb_list == 6) {
	c1 = 0.75f;
	c2 = 90.0f;
    } else {
	c1 = 0.0f;
    }
    if (c1) {
	vlo = (float) (exp (2.0 * log (c1 * Vol) / 3.0));
	qe = sumd2 / ((float) numb_list * vlo);
	nangle = (nangle + 1) / 2;
	sumth2 *= 0.5f;
	sumth *= 0.5f;
	siga = (sumth2 - 2.0f * c2 * sumth + nangle * c2 * c2) / (float) (nangle - 1);
	sigd = (sumd2 - sumd * sumd / (float) numb_list) / (float) (nangle - 1);
    }
    for (i = 0; i < numb_list - 2; ++i) {
	no[0] = vertex_list[i];	/* first point */
	if (DepthCue != 0.) {
/*calculate z coordinate after rotation, and scale edge width accordingly*/
	    radius = 0.0f;
	    for (l = 0; l <= 2; l++) {
		radius += (float) G_Rot[l][2] * s_vert[3 * no[0] + l];
	    }
	    radius = drvui->rad_edge * 0.005f * Scale + DepthCue * radius;
	}
	for (m = 0; m <= 2; ++m)
	    p[m] = s_vert[3 * no[0] + m];
	for (j = i + 1; j < numb_list - 1; ++j) {
	    no[1] = vertex_list[j];	/* second point */
	    for (k = 0; k <= 2; ++k)
		v1[k] = s_vert[3 * no[1] + k] - p[k];
	    for (k = j + 1; k < numb_list; ++k) {
		if (numb_list == 4) {	/* tetrahedron */
		    poly_list[draw_list++] = 3;	/* number of vertices for triangle */
		    poly_list[draw_list++] = vertex_list[i] + 1;
		    poly_list[draw_list++] = vertex_list[j] + 1;
		    poly_list[draw_list++] = vertex_list[k] + 1;
		    if (draw_list > 4 * NvertM) {
			poly_list =
			    (int *) zalloc ((unsigned) (4 * (nvert + 2) * sizeof (int)));
			if (!warn) {
			    fprintf (drvui->flout, "Overrun of poly_list.\n");
			    fprintf (drvui->flout, "Overrun of poly_list.\n");
			    Error_Box ("Overrun of poly_list. Please send 'str' file\n"
				       "to Larry.Finger@@lwfinger.net.");
			    warn = 1;
			}
		    }

		    imk = vertex_list[i];
		    jmk = vertex_list[j];
		    kmk = vertex_list[k];
		    for (m = 0; m < 3; m++) {
			v1[m] = s_vert[3 * imk + m] - s_vert[3 * jmk + m];
			v2[m] = s_vert[3 * imk + m] - s_vert[3 * kmk + m];
			v3[m] = s_vert[3 * kmk + m] - s_vert[3 * jmk + m];
		    }
		    d1 = vlength (v1);
		    d2 = vlength (v2);
		    d3 = vlength (v3);
		    if ((d1 > 0.01) && (d2 > 0.01) && (d3 > 0.01)) {
			if (doPOV)
			    fprintf (drvui->fpoutp,
				     " triangle{< %8.5f, %8.5f, %8.5f>,\n",
				     s_vert[3 * imk], s_vert[3 * imk + 1],
				     s_vert[3 * imk + 2]);
			if (doPOV)
			    fprintf (drvui->fpoutp,
				     "<%8.5f, %8.5f, %8.5f>, <%8.5f, %8.5f, %8.5f>  texture{pigment{color %s  }}\n",
				     s_vert[3 * jmk], s_vert[3 * jmk + 1],
				     s_vert[3 * jmk + 2], s_vert[3 * kmk],
				     s_vert[3 * kmk + 1], s_vert[3 * kmk + 2],
				     poly_col_p);
			if (doPOV)
			    fprintf (drvui->fpoutp, " }\n");
			glPushMatrix ();
			(void) sscanf (poly_col_v, "%f %f %f", &glr, &glg, &glb);
			if (doAsy) {
			    fprintf (drvui->fpouta,
				     " draw(pic, surface ( ( %8.5f, %8.5f, %8.5f)--",
				     s_vert[3 * imk], s_vert[3 * imk + 1],
				     s_vert[3 * imk + 2]);
			    fprintf (drvui->fpouta,
				     "(%8.5f, %8.5f, %8.5f)--(%8.5f, %8.5f, %8.5f)--cycle),",
				     s_vert[3 * jmk], s_vert[3 * jmk + 1],
				     s_vert[3 * jmk + 2], s_vert[3 * kmk],
				     s_vert[3 * kmk + 1], s_vert[3 * kmk + 2]);
			    fprintf (drvui->fpouta,"rgb(%4.2f,%4.2f,%4.2f) );\n",glr,glg,glb);
			}
			glColor3f (glr, glg, glb);
			glDisable (GL_LIGHTING);
			glBegin (GL_TRIANGLES);
			glVertex3f (s_vert[3 * imk], s_vert[3 * imk + 1],
				    s_vert[3 * imk + 2]);
			glVertex3f (s_vert[3 * jmk], s_vert[3 * jmk + 1],
				    s_vert[3 * jmk + 2]);
			glVertex3f (s_vert[3 * kmk], s_vert[3 * kmk + 1],
				    s_vert[3 * kmk + 2]);
			glEnd ();
			glEnable (GL_LIGHTING);
			glPopMatrix ();
			if (edges) {
			    float df[3], at[3];

			    int m;

			    char color[128];

			    if (drvui->polyhedra[polyno].poly_rad_edge > 0.005f)
				strncpy (edgecolor,
					 drvui->polyhedra[polyno].poly_col_edge, 39);
			    else
				strncpy (edgecolor, drvui->col_edge, 39);
			    trim_string (edgecolor, 40);
			    if (!strlen (edgecolor))
				strcpy (edgecolor, "White");
			    strcpy (color, edgecolor);
			    Transform_VRML_Color (color);
			    Transform_POV_Color (edgecolor);
			    if (doPOV)
				fprintf (drvui->fpoutp,
					 " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
					 s_vert[3 * imk], s_vert[3 * imk + 1],
					 s_vert[3 * imk + 2], s_vert[3 * jmk],
					 s_vert[3 * jmk + 1], s_vert[3 * jmk + 2],
					 radius);
			    if (doPOV)
				fprintf (drvui->fpoutp,
					 "  texture{pigment{color %s  }}\n", edgecolor);
			    if (doPOV)
				fprintf (drvui->fpoutp, " }\n");
			    for (m = 0; m < 3; m++) {
				at[m] = s_vert[3 * imk + m];
				df[m] = s_vert[3 * jmk + m] - s_vert[3 * imk + m];
			    }
			    push_cylinder (df, at, radius, color);
			    if (doPOV) {
				fprintf (drvui->fpoutp,
					 " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
					 s_vert[3 * jmk], s_vert[3 * jmk + 1],
					 s_vert[3 * jmk + 2], s_vert[3 * kmk],
					 s_vert[3 * kmk + 1], s_vert[3 * kmk + 2],
					 radius);
				fprintf (drvui->fpoutp,
					 "  texture{pigment{color %s  }}\n", edgecolor);
				fprintf (drvui->fpoutp, " }\n");
			    }
			    for (m = 0; m < 3; m++) {
				at[m] = s_vert[3 * jmk + m];
				df[m] = s_vert[3 * kmk + m] - s_vert[3 * jmk + m];
			    }
			    push_cylinder (df, at, radius, color);
			    if (doPOV) {
				fprintf (drvui->fpoutp,
					 " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f \n",
					 s_vert[3 * kmk], s_vert[3 * kmk + 1],
					 s_vert[3 * kmk + 2], s_vert[3 * imk],
					 s_vert[3 * imk + 1], s_vert[3 * imk + 2],
					 radius);
				fprintf (drvui->fpoutp,
					 "  texture{pigment{color %s  }}\n", edgecolor);
				fprintf (drvui->fpoutp, " }\n");
			    }
			    for (m = 0; m < 3; m++) {
				at[m] = s_vert[3 * kmk + m];
				df[m] = s_vert[3 * imk + m] - s_vert[3 * kmk + m];
			    }
			    push_cylinder (df, at, radius, color);
			}
		    }
		} else {
		    no[2] = vertex_list[k];	/* third point */
		    jplane = 3;
		    for (l = 0; l <= 2; ++l)
			v2[l] = s_vert[3 * no[2] + l] - p[l];	/* get v2 */
		    sign_z = 0.0f;
		    ind = 0;
		    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
		    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
		    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
		    for (l = 0; l < numb_list; ++l) {
			vol = 0.0f;
			for (m = 0; m <= 2; ++m)
			    vol = vol + (s_vert[3 * vertex_list[l] + m] - p[m]) * v3[m];
			if ((float) fabs (vol) < drvui->polylimit) {
			    if (l > k) {
				no[jplane++] = vertex_list[l];	/* add to plane list */
			    } else {
				if ((l != i) && (l != j) && (l != k)) {
				    ind = 1;	/* we already have this plane */
				    l = numb_list;	/* skip rest of loop */
				}
			    }
			} else {
/* point is not in plane - check if first out-of-plane point found */
			    if (sign_z == 0.0) {
				sign_z = vol;	/* yes - set sign */
			    } else {
/* if here, not first out of plane - check that vol and sign_z have same sign. 
   If not, set ind=1 to indicate not a plane */
				if (vol * sign_z < 0.0) {
				    ind = 1;
				    l = numb_list;	/* skip rest of inner loop */
				}
			    }
			}
		    }		/* end of loop on l */
		    if (ind == 0) {	/* 0 means a surface plane */
			poly_list[draw_list++] = jplane;	/* number of vertices in plane */
			for (l = 0; l < jplane; ++l)
			    poly_list[draw_list++] = no[l] + 1;
			if (doPOV)
			    fprintf (drvui->fpoutp, "/* %d-sided polygon */\n", jplane);
/* Make vertices consecutive - find center of polygon */
			xm = ym = zm = 0.0f;
			for (ii = 0; ii < jplane; ii++) {
			    xm = xm + s_vert[3 * no[ii]];
			    ym = ym + s_vert[3 * no[ii] + 1];
			    zm = zm + s_vert[3 * no[ii] + 2];
			}
			xm = xm / (float) jplane;
			ym = ym / (float) jplane;
			zm = zm / (float) jplane;

/* sort vertices in consecutive order based on minimal angles between vectors */

			for (ii = 0; ii < jplane; ii++)
			    ns[ii] = 0;
			ns[0] = no[0];
			for (ii = 0; ii < jplane - 1; ii++) {
			    phi0 = 1000.0f;
			    vxs = s_vert[3 * ns[ii]] - xm;
			    vys = s_vert[3 * ns[ii] + 1] - ym;
			    vzs = s_vert[3 * ns[ii] + 2] - zm;
			    bms = (float) sqrt (vxs * vxs + vys * vys + vzs * vzs);
			    for (jj = 0; jj < jplane; jj++) {
				for (ll = 0; ll < ii; ll++) {
				    if (ns[ll] == no[jj])
					jj++;
				    if (jj >= jplane)
					break;
				}
				if (jj < jplane) {
				    if (ns[ii] != no[jj]) {
					vx = s_vert[3 * no[jj]] - xm;
					vy = s_vert[3 * no[jj] + 1] - ym;
					vz = s_vert[3 * no[jj] + 2] - zm;
					bm = (float) sqrt (vx * vx + vy * vy + vz * vz);
					cosphi =
					    (vxs * vx + vys * vy + vzs * vz) / (bms * bm);
					if (cosphi < -1.)
					    cosphi = -1.;
					if (cosphi > 1.)
					    cosphi = 1.;
					phi = acos (cosphi) * RAD;
					if (phi < phi0) {
					    phi0 = (float) phi;
					    ns[ii + 1] = no[jj];
					}
				    }
				}
			    }
			}
			if (numb_list != 4) {
			    draw_list = draw_list - jplane;
			    for (l = 0; l < jplane; ++l)
				poly_list[draw_list++] = ns[l] + 1;
			}

			glPushMatrix ();
			(void) sscanf (poly_col_v, "%f %f %f", &glr, &glg, &glb);
			glDisable (GL_LIGHTING);
			glColor3f (glr, glg, glb);
			glBegin (GL_TRIANGLES);

			for (l = 0; l < jplane - 2; ++l) {
			    d1 = (float) fabs (s_vert[3 * ns[0]] -
					       s_vert[3 * ns[l + 1]]) +
				(float) fabs (s_vert[3 * ns[0] + 1] -
					      s_vert[3 * ns[l + 1] + 1]) +
				(float) fabs (s_vert[3 * ns[0] + 2] -
					      s_vert[3 * ns[l + 1] + 2]);
			    d2 = (float) fabs (s_vert[3 * ns[0]] -
					       s_vert[3 * ns[l + 2]]) +
				(float) fabs (s_vert[3 * ns[0] + 1] -
					      s_vert[3 * ns[l + 2] + 1]) +
				(float) fabs (s_vert[3 * ns[0] + 2] -
					      s_vert[3 * ns[l + 2] + 2]);
			    d3 = (float) fabs (s_vert[3 * ns[l + 2]] -
					       s_vert[3 * ns[l + 1]]) +
				(float) fabs (s_vert[3 * ns[l + 2] + 1] -
					      s_vert[3 * ns[l + 1] + 1]) +
				(float) fabs (s_vert[3 * ns[l + 2] + 2] -
					      s_vert[3 * ns[l + 1] + 2]);
			    if ((d1 > 0.01) && (d2 > 0.01) && (d3 > 0.01)) {
				if (doPOV) {
				    fprintf (drvui->fpoutp,
					     " triangle {< %8.5f, %8.5f, %8.5f>,\n",
					     s_vert[3 * ns[0]], s_vert[3 * ns[0] + 1],
					     s_vert[3 * ns[0] + 2]);
				    fprintf (drvui->fpoutp,
					     "<%8.5f, %8.5f, %8.5f>, <%8.5f, %8.5f, %8.5f>\n",
					     s_vert[3 * ns[l + 1]],
					     s_vert[3 * ns[l + 1] + 1],
					     s_vert[3 * ns[l + 1] + 2],
					     s_vert[3 * ns[l + 2]],
					     s_vert[3 * ns[l + 2] + 1],
					     s_vert[3 * ns[l + 2] + 2]);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n",
					     poly_col_p);
				    fprintf (drvui->fpoutp, " }\n");
				}
				if (doAsy) {
			    	    fprintf (drvui->fpouta,
					" draw(pic, surface ( ( %8.5f, %8.5f, %8.5f)--",
					     s_vert[3 * ns[0]], s_vert[3 * ns[0] + 1],
					     s_vert[3 * ns[0] + 2]);
				    fprintf (drvui->fpouta,
					"(%8.5f, %8.5f, %8.5f)--(%8.5f, %8.5f, %8.5f)--cycle),",
					     s_vert[3 * ns[l + 1]],
					     s_vert[3 * ns[l + 1] + 1],
					     s_vert[3 * ns[l + 1] + 2],
					     s_vert[3 * ns[l + 2]],
					     s_vert[3 * ns[l + 2] + 1],
					     s_vert[3 * ns[l + 2] + 2]);
				    fprintf (drvui->fpouta,"rgb(%4.2f,%4.2f,%4.2f) );\n",glr,glg,glb);
				}
/* 
      need to add something like
      glNormal3f(xm-xc,ym-yc,zm-zc ); 
      for lighting -  xc,yc,zc being the center of the polyhedron
      (which is not normally available to build_poly_list)
*/
				glVertex3f (s_vert[3 * ns[0]], s_vert[3 * ns[0] + 1],
					    s_vert[3 * ns[0] + 2]);
				glVertex3f (s_vert[3 * ns[l + 1]],
					    s_vert[3 * ns[l + 1] + 1],
					    s_vert[3 * ns[l + 1] + 2]);
				glVertex3f (s_vert[3 * ns[l + 2]],
					    s_vert[3 * ns[l + 2] + 1],
					    s_vert[3 * ns[l + 2] + 2]);
			    }
			}
			glEnd ();
			glEnable (GL_LIGHTING);
			glPopMatrix ();
			if (edges) {
			    float at[3], df[3];

			    int m;

			    char color[128];

			    if (drvui->polyhedra[polyno].poly_rad_edge >= 0.001f)
				strncpy (edgecolor,
					 drvui->polyhedra[polyno].poly_col_edge, 39);
			    else
				strncpy (edgecolor, drvui->col_edge, 39);
			    trim_string (edgecolor, 40);
			    if (!strlen (edgecolor))
				strcpy (edgecolor, "White");
			    strcpy (color, edgecolor);
			    Transform_VRML_Color (color);
			    Transform_POV_Color (edgecolor);
			    for (l = 0; l < jplane - 1; ++l) {
				d1 = (float) fabs (s_vert[3 * ns[l]] -
						   s_vert[3 * ns[l + 1]]) +
				    (float) fabs (s_vert[3 * ns[l] + 1] -
						  s_vert[3 * ns[l + 1] + 1]) +
				    (float) fabs (s_vert[3 * ns[l] + 2] -
						  s_vert[3 * ns[l + 1] + 2]);

				if (d1 > 0.01) {
				    if (doPOV) {
					fprintf (drvui->fpoutp,
						 " cylinder { <%8.5f, %8.5f, %8.5f>, <%8.5f, %8.5f, %8.5f> , %8.5f \n",
						 s_vert[3 * ns[l]], s_vert[3 * ns[l] + 1],
						 s_vert[3 * ns[l] + 2],
						 s_vert[3 * ns[l + 1]],
						 s_vert[3 * ns[l + 1] + 1],
						 s_vert[3 * ns[l + 1] + 2], radius);
					fprintf (drvui->fpoutp,
						 "  texture{pigment{color %s  }}\n",
						 edgecolor);
					fprintf (drvui->fpoutp, " }\n");
				    }
				    for (m = 0; m < 3; m++) {
					at[m] = s_vert[3 * ns[l] + m];
					df[m] = s_vert[3 * ns[l + 1] + m] - at[m];
				    }
				    push_cylinder (df, at, radius, color);
				}
			    }
			    d1 = (float) fabs (s_vert[3 * ns[0]] -
					       s_vert[3 * ns[jplane - 1]]) +
				(float) fabs (s_vert[3 * ns[0] + 1] -
					      s_vert[3 * ns[jplane - 1] + 1]) +
				(float) fabs (s_vert[3 * ns[0] + 2] -
					      s_vert[3 * ns[jplane - 1] + 2]);
			    if (d1 > 0.01) {
				if (doPOV) {
				    fprintf (drvui->fpoutp,
					     " cylinder { <%8.5f, %8.5f, %8.5f>, <%8.5f, %8.5f, %8.5f> , %8.5f \n",
					     s_vert[3 * ns[0]], s_vert[3 * ns[0] + 1],
					     s_vert[3 * ns[0] + 2],
					     s_vert[3 * ns[jplane - 1]],
					     s_vert[3 * ns[jplane - 1] + 1],
					     s_vert[3 * ns[jplane - 1] + 2], radius);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n",
					     edgecolor);
				    fprintf (drvui->fpoutp, " }\n");
				}
				for (m = 0; m < 3; m++) {
				    at[m] = s_vert[3 * ns[0] + m];
				    df[m] = s_vert[3 * ns[jplane - 1] + m] - at[m];
				}
				push_cylinder (df, at, radius, color);
			    }
			}
		    }
		}		/* end of numb_list != 4 */
	    }			/* end of for (k.... */
	}
    }
    if (c1) {
	fprintf (drvui->flout, " %7.3f %8.5f %8.5f %8.3f\n", Vol, sigd, qe, siga);
    } else {
	fprintf (drvui->flout, " %7.3f %8.5f\n", Vol, sigd);
    }

    glPopName ();
    free (no);
    free (ns);
    return;
}				/* end of build_poly_list */

/* ************************************************************** */
/* ************************************************************** */

void
draw_cell (int docell)

/* procedure to add unit-cell outline to geometry file */
/* docell - non-zero if unit-cell to be drawn */
{
    int i, j, k;

    float vert[3];

    float radius;

//  static char abc[] = "cbao";
    float df[3];

    float at[3];

    char col_cell_v[40];

    char col_cell_p[40];

    int nocell;

    radius = rad_cell * 0.01f * Scale;
    if (!drvui->no_subsys)
	drvui->no_subsys = 1;

    for (nocell = 0; nocell < drvui->no_subsys; nocell++) {
	nvert = 0;		/* initialize vertex list */
	for (i = 0; i <= 2; ++i)
	    vert[i] = 0.0f;
	sub_add_vert_nc (vert, nocell);	// place 0,0,0 in list
	vert[2] = 1.0f;
	sub_add_vert_nc (vert, nocell);	// place 0,0,1 in list
	vert[1] = 1.0f;
	sub_add_vert_nc (vert, nocell);	// place 0,1,1 in list
	vert[2] = 0.0f;
	sub_add_vert_nc (vert, nocell);	// place 0,1,0 in list
	vert[1] = 0.0f;
	sub_add_vert_nc (vert, nocell);	// place 0,0,0 in list
	vert[0] = 1.0f;
	sub_add_vert_nc (vert, nocell);	// place 1,0,0 in list
	vert[1] = 1.0f;
	sub_add_vert_nc (vert, nocell);	// place 1,1,0 in list
	vert[2] = 1.0f;
	sub_add_vert_nc (vert, nocell);	// place 1,1,1 in list
	vert[1] = 0.0f;
	sub_add_vert_nc (vert, nocell);	/* place 1,0,1 in list */
	vert[2] = 0.0f;
	sub_add_vert_nc (vert, nocell);	/* place 1,0,0 in list */
	vert[1] = 1.0f;
	sub_add_vert_nc (vert, nocell);	/* place 1,1,0 in list */
	vert[0] = 0.0f;
	sub_add_vert_nc (vert, nocell);	/* place 0,1,0 in list */
	vert[2] = 1.0f;
	sub_add_vert_nc (vert, nocell);	/* place 0,1,1 in list */
	vert[0] = 1.0f;
	sub_add_vert_nc (vert, nocell);	/* place 1,1,1 in list */
	vert[1] = 0.0f;
	sub_add_vert_nc (vert, nocell);	/* place 1,0,1 in list */
	vert[0] = 0.0f;
	sub_add_vert_nc (vert, nocell);	/* place 0,0,1 in list */

	if (docell) {
	    strcpy (col_cell_v, drvui->col_cell);
	    Transform_VRML_Color (col_cell_v);
	    strcpy (col_cell_p, drvui->col_cell);
	    Transform_POV_Color (col_cell_p);
	    if (doVrml) {
		if (Vrml2) {
		    fprintf (drvui->fpoutv, " Shape {\n");
		    fprintf (drvui->fpoutv,
			     "  geometry IndexedLineSet{ coord Coordinate{ point[\n");
		} else {
		    fprintf (drvui->fpoutv, " Separator{\n");
		    fprintf (drvui->fpoutv, "  Material{\n");
		    fprintf (drvui->fpoutv, "   diffuseColor %s\n", col_cell_v);
		    fprintf (drvui->fpoutv, "  }\n");
		    fprintf (drvui->fpoutv, "  Coordinate3{ point[\n");
		}
		fprintf (drvui->fpoutv, "   %8.5f %8.5f %8.5f,\n", s_vert[0], s_vert[1],
			 s_vert[2]);
		fprintf (drvui->fpoutv, "   %8.5f %8.5f %8.5f,\n", s_vert[3], s_vert[4],
			 s_vert[5]);
		fprintf (drvui->fpoutv, "   %8.5f %8.5f %8.5f,\n", s_vert[6], s_vert[7],
			 s_vert[8]);
		fprintf (drvui->fpoutv, "   %8.5f %8.5f %8.5f,\n", s_vert[9], s_vert[10],
			 s_vert[11]);
		fprintf (drvui->fpoutv, "   %8.5f %8.5f %8.5f,\n", s_vert[15], s_vert[16],
			 s_vert[17]);
		fprintf (drvui->fpoutv, "   %8.5f %8.5f %8.5f,\n", s_vert[18], s_vert[19],
			 s_vert[20]);
		fprintf (drvui->fpoutv, "   %8.5f %8.5f %8.5f,\n", s_vert[21], s_vert[22],
			 s_vert[23]);
		fprintf (drvui->fpoutv, "   %8.5f %8.5f %8.5f]\n  }\n", s_vert[24],
			 s_vert[25], s_vert[26]);
		if (Vrml2) {
		    fprintf (drvui->fpoutv, "  coordIndex[");
		} else {
		    fprintf (drvui->fpoutv, "  IndexedLineSet{coordIndex[");
		}
		fprintf (drvui->fpoutv, "0,1,2,3,0,4,5,6,7,1,2,6,5,3,0,4,7,-1]\n");
		if (Vrml2) {
		    fprintf (drvui->fpoutv,
			     "  color Color { color [%s]}\n  colorIndex[0]\n"
			     "  colorPerVertex FALSE\n", col_cell_v);
		}
		fprintf (drvui->fpoutv, " }\n");
		fprintf (drvui->fpoutv, "}\n");
	    }
	    for (i = 0; i < 12; i++) {	// generate unit-cell edges for openGL
		j = i;
		if (i == 9)
		    j = 10;
		if (i == 10)
		    j = 12;
		k = j + 1;
		if (i == 11) {
		    j = 8;
		    k = 1;
		}
		at[0] = s_vert[3 * j];	// beginning vertex
		at[1] = s_vert[3 * j + 1];
		at[2] = s_vert[3 * j + 2];
		df[0] = s_vert[3 * k] - at[0];	// difference from start to end
		df[1] = s_vert[3 * k + 1] - at[1];
		df[2] = s_vert[3 * k + 2] - at[2];
		push_cylinder (df, at, radius, col_cell_v);	// draw the cylinder
	    }
	    if (doPOV) {
		fprintf (drvui->fpoutp, "  cylinder{<%f,%f,%f>,<%f,%f,%f>,%f \n",
			 s_vert[0], s_vert[1], s_vert[2], s_vert[3], s_vert[4],
			 s_vert[5], radius);
		fprintf (drvui->fpoutp, "   texture{pigment{color %s }}\n", col_cell_p);
		fprintf (drvui->fpoutp, "   no_shadow no_reflection\n  }\n");
		fprintf (drvui->fpoutp, "  cylinder{<%f,%f,%f>,<%f,%f,%f>,%f\n",
			 s_vert[3], s_vert[4], s_vert[5], s_vert[6], s_vert[7],
			 s_vert[8], radius);
		fprintf (drvui->fpoutp, "   texture{pigment{color %s }}\n", col_cell_p);
		fprintf (drvui->fpoutp, "   no_shadow no_reflection\n  }\n");
		fprintf (drvui->fpoutp, "  cylinder{<%f,%f,%f>,<%f,%f,%f>,%f\n",
			 s_vert[6], s_vert[7], s_vert[8], s_vert[9], s_vert[10],
			 s_vert[11], radius);
		fprintf (drvui->fpoutp, "   texture{pigment{color %s }}\n", col_cell_p);
		fprintf (drvui->fpoutp, "   no_shadow no_reflection\n  }\n");
		fprintf (drvui->fpoutp, "  cylinder{<%f,%f,%f>,<%f,%f,%f>, %f\n",
			 s_vert[9], s_vert[10], s_vert[11], s_vert[12], s_vert[13],
			 s_vert[14], radius);
		fprintf (drvui->fpoutp, "   texture{pigment{color %s }}\n", col_cell_p);
		fprintf (drvui->fpoutp, "   no_shadow no_reflection\n  }\n");
		fprintf (drvui->fpoutp, "  cylinder{<%f,%f,%f>,<%f,%f,%f>,%f\n",
			 s_vert[12], s_vert[13], s_vert[14], s_vert[15], s_vert[16],
			 s_vert[17], radius);
		fprintf (drvui->fpoutp, "   texture{pigment{color %s }}\n", col_cell_p);
		fprintf (drvui->fpoutp, "   no_shadow no_reflection\n  }\n");
		fprintf (drvui->fpoutp, "  cylinder{<%f,%f,%f>,<%f,%f,%f>,%f\n",
			 s_vert[15], s_vert[16], s_vert[17], s_vert[18], s_vert[19],
			 s_vert[20], radius);
		fprintf (drvui->fpoutp, "   texture{pigment{color %s }}\n", col_cell_p);
		fprintf (drvui->fpoutp, "   no_shadow no_reflection\n  }\n");
		fprintf (drvui->fpoutp, "  cylinder{<%f,%f,%f>,<%f,%f,%f>,%f\n",
			 s_vert[18], s_vert[19], s_vert[20], s_vert[21], s_vert[22],
			 s_vert[23], radius);
		fprintf (drvui->fpoutp, "   texture{pigment{color %s }}\n", col_cell_p);
		fprintf (drvui->fpoutp, "   no_shadow no_reflection\n  }\n");
		fprintf (drvui->fpoutp, "  cylinder{<%f,%f,%f>,<%f,%f,%f>,%f\n",
			 s_vert[21], s_vert[22], s_vert[23], s_vert[24], s_vert[25],
			 s_vert[26], radius);
		fprintf (drvui->fpoutp, "   texture{pigment{color %s }}\n", col_cell_p);
		fprintf (drvui->fpoutp, "   no_shadow no_reflection\n  }\n");
		fprintf (drvui->fpoutp, "  cylinder{<%f,%f,%f>,<%f,%f,%f>,%f\n",
			 s_vert[24], s_vert[25], s_vert[26], s_vert[27], s_vert[28],
			 s_vert[29], radius);
		fprintf (drvui->fpoutp, "   texture{pigment{color %s }}\n", col_cell_p);
		fprintf (drvui->fpoutp, "   no_shadow no_reflection\n  }\n");
		fprintf (drvui->fpoutp, "  cylinder{<%f,%f,%f>,<%f,%f,%f>,%f\n",
			 s_vert[30], s_vert[31], s_vert[32], s_vert[33], s_vert[34],
			 s_vert[35], radius);
		fprintf (drvui->fpoutp, "   texture{pigment{color %s }}\n", col_cell_p);
		fprintf (drvui->fpoutp, "   no_shadow no_reflection\n  }\n");
		fprintf (drvui->fpoutp, "  cylinder{<%f,%f,%f>,<%f,%f,%f>,%f\n",
			 s_vert[36], s_vert[37], s_vert[38], s_vert[39], s_vert[40],
			 s_vert[41], radius);
		fprintf (drvui->fpoutp, "   texture{pigment{color %s }}\n", col_cell_p);
		fprintf (drvui->fpoutp, "   no_shadow no_reflection\n  }\n");
		fprintf (drvui->fpoutp, "  cylinder{<%f,%f,%f>,<%f,%f,%f>,%f\n",
			 s_vert[24], s_vert[25], s_vert[26], s_vert[3], s_vert[4],
			 s_vert[5], radius);
		fprintf (drvui->fpoutp, "   texture{pigment{color %s }}\n", col_cell_p);
		fprintf (drvui->fpoutp, "   no_shadow no_reflection\n  }\n");
	    }
	    if (doAsy) {
		float r,g,b;
		sscanf (col_cell_v,"%f %f %f",&r,&g,&b);
		fprintf (drvui->fpouta, " draw(pic, (%8.5f,%8.5f,%8.5f)--", s_vert[0], s_vert[1],
			 s_vert[2]);
		fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--", s_vert[3], s_vert[4],
			 s_vert[5]);
		fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--", s_vert[6], s_vert[7],
			 s_vert[8]);
		fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--cycle,", s_vert[9], s_vert[10],
			 s_vert[11]);
		fprintf (drvui->fpouta,"rgb(%4.2f,%4.2f,%4.2f) );\n",r,g,b);
		fprintf (drvui->fpouta, " draw(pic, (%8.5f,%8.5f,%8.5f)--", s_vert[15], s_vert[16],
			 s_vert[17]);
		fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--", s_vert[18], s_vert[19],
			 s_vert[20]);
		fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--", s_vert[21], s_vert[22],
			 s_vert[23]);
		fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--cycle,", s_vert[24], s_vert[25],
			 s_vert[26]);
		fprintf (drvui->fpouta,"rgb(%4.2f,%4.2f,%4.2f) );\n",r,g,b);
		fprintf (drvui->fpouta, " draw(pic, (%8.5f,%8.5f,%8.5f)--", s_vert[0], s_vert[1],
			 s_vert[2]);
		fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f),", s_vert[15], s_vert[16],
			 s_vert[17]);
		fprintf (drvui->fpouta,"rgb(%4.2f,%4.2f,%4.2f) );\n",r,g,b);
		fprintf (drvui->fpouta, " draw(pic, (%8.5f,%8.5f,%8.5f)--", s_vert[3], s_vert[4],
			 s_vert[5]);
		fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f),", s_vert[24], s_vert[25],
			 s_vert[26]);
		fprintf (drvui->fpouta,"rgb(%4.2f,%4.2f,%4.2f) );\n",r,g,b);
		fprintf (drvui->fpouta, " draw(pic, (%8.5f,%8.5f,%8.5f)--", s_vert[6], s_vert[7],
			 s_vert[8]);
		fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f),", s_vert[21], s_vert[22],
			 s_vert[23]);
		fprintf (drvui->fpouta,"rgb(%4.2f,%4.2f,%4.2f) );\n",r,g,b);
		fprintf (drvui->fpouta, " draw(pic, (%8.5f,%8.5f,%8.5f)--", s_vert[9], s_vert[10],
			 s_vert[11]);
		fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f),", s_vert[18], s_vert[19],
			 s_vert[20]);
		fprintf (drvui->fpouta,"rgb(%4.2f,%4.2f,%4.2f) );\n",r,g,b);
	    }
	}			/* if (docell) */
    }				/* for(nocell  */

    if (Display_axes) {

	nvert = 0;		/* initialize vertex list */
	generate_triple ();
	if (radius == 0.)
	    radius = 0.0002f * Scale;
	strcpy (col_cell_v, drvui->col_cell);
	Transform_VRML_Color (col_cell_v);
	strcpy (col_cell_p, drvui->col_cell);
	Transform_POV_Color (col_cell_p);

/* draw vector triple in POV */
	if (doPOV) {
	    fprintf (drvui->fpoutp, "/* drawing vector triple */");
	    fprintf (drvui->fpoutp, "  cylinder{<%f,%f,%f>,<%f,%f,%f>,%f \n", s_vert[0],
		     s_vert[1], s_vert[2], s_vert[3], s_vert[4], s_vert[5], radius);
	    fprintf (drvui->fpoutp, "   texture{pigment{color %s }}\n", col_cell_p);
	    fprintf (drvui->fpoutp, "   no_shadow no_reflection\n  }\n");
	    fprintf (drvui->fpoutp, "  cylinder{<%f,%f,%f>,<%f,%f,%f>,%f \n", s_vert[0],
		     s_vert[1], s_vert[2], s_vert[6], s_vert[7], s_vert[8], radius);
	    fprintf (drvui->fpoutp, "   texture{pigment{color %s }}\n", col_cell_p);
	    fprintf (drvui->fpoutp, "   no_shadow no_reflection\n  }\n");
	    fprintf (drvui->fpoutp, "  cylinder{<%f,%f,%f>,<%f,%f,%f>,%f \n", s_vert[0],
		     s_vert[1], s_vert[2], s_vert[9], s_vert[10], s_vert[11], radius);
	    fprintf (drvui->fpoutp, "   texture{pigment{color %s }}\n", col_cell_p);
	    fprintf (drvui->fpoutp, "   no_shadow no_reflection\n  }\n");
	}
/* draw vector triple in Asymptote */
	if (doAsy) {
	    float r,g,b;
	    sscanf (col_cell_v,"%f %f %f",&r,&g,&b);
	    fprintf (drvui->fpouta, "// vector triple\n");
	    fprintf (drvui->fpouta, "  draw(pic, (%8.5f,%8.5f,%8.5f)--(%8.5f,%8.5f,%8.5f), \n", 
		     s_vert[0], s_vert[1], s_vert[2], s_vert[3], s_vert[4], s_vert[5]);
	    fprintf (drvui->fpouta,"rgb(%4.2f,%4.2f,%4.2f) );\n",r,g,b);
	    fprintf (drvui->fpouta, "  draw(pic, (%8.5f,%8.5f,%8.5f)--(%8.5f,%8.5f,%8.5f), \n", 
		     s_vert[0], s_vert[1], s_vert[2], s_vert[6], s_vert[7], s_vert[8]);
	    fprintf (drvui->fpouta,"rgb(%4.2f,%4.2f,%4.2f) );\n",r,g,b);
	    fprintf (drvui->fpouta, "  draw(pic, (%8.5f,%8.5f,%8.5f)--(%8.5f,%8.5f,%8.5f), \n", 
		     s_vert[0], s_vert[1], s_vert[2], s_vert[9], s_vert[10], s_vert[11]);
	    fprintf (drvui->fpouta,"rgb(%4.2f,%4.2f,%4.2f) );\n",r,g,b);
	}
/* draw vector triple in VRML */
	if (doVrml) {
	    if (Vrml2) {
		fprintf (drvui->fpoutv, " Shape{\n");
		fprintf (drvui->fpoutv,
			 "  geometry IndexedLineSet { coord Coordinate{ point [\n");
		fprintf (drvui->fpoutv, "   %8.5f %8.5f %8.5f,\n", s_vert[0], s_vert[1],
			 s_vert[2]);
		fprintf (drvui->fpoutv, "   %8.5f %8.5f %8.5f,\n", s_vert[3], s_vert[4],
			 s_vert[5]);
		fprintf (drvui->fpoutv, "   %8.5f %8.5f %8.5f,\n", s_vert[6], s_vert[7],
			 s_vert[8]);
		fprintf (drvui->fpoutv, "   %8.5f %8.5f %8.5f]\n  }\n", s_vert[9],
			 s_vert[10], s_vert[11]);
		fprintf (drvui->fpoutv, "  coordIndex [");
		fprintf (drvui->fpoutv, "0,3,0,2,0,1,-1]\n");
		fprintf (drvui->fpoutv, "   color Color { color [%s]}\n   colorIndex[0]\n"
			 "   colorPerVertex FALSE\n", col_cell_v);
		fprintf (drvui->fpoutv, "  }\n");
		fprintf (drvui->fpoutv, " }\n");
	    } else {
		fprintf (drvui->fpoutv, " Separator{\n");
		fprintf (drvui->fpoutv, "  Material{\n");
		fprintf (drvui->fpoutv, "   diffuseColor %s\n", col_cell_v);
		fprintf (drvui->fpoutv, "  }\n");
		fprintf (drvui->fpoutv, "  Coordinate3{ point[\n");
		fprintf (drvui->fpoutv, "   %8.5f %8.5f %8.5f,\n", s_vert[0], s_vert[1],
			 s_vert[2]);
		fprintf (drvui->fpoutv, "   %8.5f %8.5f %8.5f,\n", s_vert[3], s_vert[4],
			 s_vert[5]);
		fprintf (drvui->fpoutv, "   %8.5f %8.5f %8.5f,\n", s_vert[6], s_vert[7],
			 s_vert[8]);
		fprintf (drvui->fpoutv, "   %8.5f %8.5f %8.5f]\n  }\n", s_vert[9],
			 s_vert[10], s_vert[11]);
		fprintf (drvui->fpoutv, "  IndexedLineSet{coordIndex[");
		fprintf (drvui->fpoutv, "0,3,0,2,0,1,-1]\n");
		fprintf (drvui->fpoutv, "  }\n");
		fprintf (drvui->fpoutv, " }\n");
	    }
	}
    }

    nvert = 0;			/* empty vertex list */
}				/* end of draw_cell */

/* ************************************************************** */
/* ************************************************************** */

void
draw_GL_triple (void)
{
// draw the lines for the vector triple in openGL
    int i;

    float at[3], df[3];

    float radius = rad_cell * 0.01f * Scale;

    char col_cell_v[40];

    strcpy (col_cell_v, drvui->col_cell);
    Transform_VRML_Color (col_cell_v);

    if (radius == 0.)
	radius = 0.0002f * Scale;

    for (i = 0; i < 3; i++) {
	at[i] = s_vert[i];
	df[i] = s_vert[3 + i] - at[i];
    }
    push_cylinder (df, at, radius, col_cell_v);
    for (i = 0; i < 3; i++) {
	at[i] = s_vert[i];
	df[i] = s_vert[6 + i] - at[i];
    }
    push_cylinder (df, at, radius, col_cell_v);
    for (i = 0; i < 3; i++) {
	at[i] = s_vert[i];
	df[i] = s_vert[9 + i] - at[i];
    }
    push_cylinder (df, at, radius, col_cell_v);
}

/* ************************************************************** */
/* ************************************************************** */

void
generate_spheres (void)

/* routine to generate spherical objects and add then to edit list */
{
    float *radii;		/* dynamic array for sphere radii */

    int Sphere_Count;		/* Counter for number output */

    int i, j, k;		/* loop counters */

    Sphere_Count = 0;
    for (i = 1; i < drvui->nsphere; ++i) {	/* loop through spheres */
	if (drvui->spheres[i].sphere_fn != drvui->frame_no)
	    continue;		// skip if not in this frame
	nvert = 0;		/* clear the vertex list */
	for (j = 0; j < natom; ++j) {	/* loop through atoms */
	    if (drvui->atoms[j].atom_fn != drvui->frame_no)
		continue;
	    if ((drvui->atoms[j].atom_n & 255) == i)
		find_all_in_box (j);
	}
	if (nvert > 0) {	/* if spheres present */
/* allocate storage for radii of spheres */
	    if (!(radii = (float *) zalloc ((unsigned) (nvert * sizeof (float))))) {
		Error_Box ("Unable to allocate memory for storage of spheres.");
		return;
	    };
	    for (k = 0; k < nvert; ++k) {	/* copy radii info to storage */
		radii[k] = drvui->spheres[i].sphere_size * drvui->Sphere_Mult;
		if (drvui->spheres[i].sphere_size < 0.)
		    radii[k] *= -drvui->vert_occ[k];
	    }
	    Output_Spheres (radii, i);
	    free (radii);
	}
	Sphere_Count += nvert;
    }
    if (Sphere_Count > 0) {
	fprintf (drvui->fcns, "%4d spheres.\n", Sphere_Count);
	fprintf (drvui->flout, "Generated %4d spheres.\n", Sphere_Count);
    }
}				/* end of generate_spheres */

/* ************************************************************** */
/* ************************************************************** */

void
generate_triple (void)
{
// generate the vertices for the vector triple

    int i, j;

    float vert[3];

    nvert = 0;
    for (i = 0; i < 3; i++)
	vert[i] = 0.0f;
    add_vert_nc (vert);		/* place origin at 0,0,0 in list */
    vert[2] = 1.0f / drvui->lat_con[2];
    add_vert_nc (vert);		/* place 1 A vector along 001 */
    vert[2] = 0.0f;
    vert[1] = 1.0f / drvui->lat_con[1];
    add_vert_nc (vert);		/* place 1 A vector along 010 */
    vert[1] = 0.0f;
    vert[0] = 1.0f / drvui->lat_con[0];
    add_vert_nc (vert);		/* place 1 A vector along 100 */
    vert[2] = 1.4f / drvui->lat_con[2];
    vert[0] = 0.0f;
    add_vert_nc (vert);		/* place 1.4 A vector along 001 */
    vert[2] = 0.0f;
    vert[1] = 1.4f / drvui->lat_con[1];
    add_vert_nc (vert);		/* place 1.4 A vector along 010 */
    vert[1] = 0.0f;
    vert[0] = 1.4f / drvui->lat_con[0];
    add_vert_nc (vert);		/* place 1.4 A vector along 100 */
    for (i = nvert - 1; i >= 0; --i) {	/* Translate coordinates */
	for (j = 0; j < 3; j++) {
	    s_vert[3 * i + j] += offset[j] - s_vert[j];
	}
    }
}

/* ************************************************************** */
/* ************************************************************** */

void
push_cylinder (float df[3], float at[3], float radius, char *color)
// routine to push an openGL cylinder that represents a unit-cell edge
{
    float beta, gamma;

    float glr, glg, glb;

    float d;			//,radius;

    GLUquadricObj *glu_quadObj;

//    radius = radius * 0.01f *Scale;
    d = (float) sqrt (df[0] * df[0] + df[1] * df[1] + df[2] * df[2]);
    beta = (float) atan2 (df[0], df[1] + 0.0000001f);	/* rotation angle about Z (in radians) */
    gamma = (float) sqrt (df[0] * df[0] + df[1] * df[1]);
    gamma = (float) atan2 (gamma, df[2]);	/* Rotation angle about X */
    glu_quadObj = gluNewQuadric ();
    glPushMatrix ();
    (void) sscanf (color, "%f %f %f", &glr, &glg, &glb);
    glColor3f (glr, glg, glb);
    glTranslatef (at[0], at[1], at[2]);
    glRotatef (-beta * (float) RAD, 0.0f, 0.0f, 1.0f);
    glRotatef (-gamma * (float) RAD, 1.0f, 0.0f, 0.0f);
    gluQuadricDrawStyle (glu_quadObj, GLU_FILL);
    gluCylinder (glu_quadObj, radius, radius, d, 10, 1);
    gluDeleteQuadric (glu_quadObj);
    glPopMatrix ();
}

/* ************************************************************** */
/* ************************************************************** */

int
p_eigen (float beta[6], float valu[3], float vect[3][3])

/* routine to calculate eigenvalues 'valu' and vectors 'v' */
/* copy of eigen without check for n.p.d'ness */
{
    static double errnd = 0.000007;

    double b0, b1, b2, c0, c1;

    double a[3][3], b[3][3], w[3][3], u[3], v[3], z[3];

    double x, y, tem, sigma, vnew, vold, smax;

    double Vect[3][3];

    int i, j, l, ii, iii, i1, imax = 0;

    static double tpi2 = 19.7392088;	/* 1/2 pi^2 */

#define phif(z) ((b2 - z) * z + b1) * z + b0

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
        Equation is -L^3 + b2 L^2 + b1 L + b0 = 0
*/
    b2 = a[0][0] + a[1][1] + a[2][2];
    b1 = -a[0][0] * a[1][1] - a[0][0] * a[2][2] - a[1][1] * a[2][2]
	+ a[0][2] * a[2][0] + a[0][1] * a[1][0] + a[1][2] * a[2][1];
    b0 = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2]
	+ a[2][0] * a[1][2] * a[0][1] - a[2][0] * a[0][2] * a[1][1]
	- a[0][0] * a[2][1] * a[1][2] - a[1][0] * a[0][1] * a[2][2];


/* first root by bisection */

    x = 0.0;
    y = sigma;
    tem = phif (sigma);
    vnew = 0.0;
    if (b0 != 0.0) {
	if ((b0 < 0.0) && (tem <= 0.0))
	    y = -y;
	if (b0 > 0.0) {
	    y = 0.0;
	    x = sigma;
	    if (tem > 0.0)
		x = -x;
	}

/* phif(y) > 0 and phif(x) < 0  */

	for (j = 0; j < 80; ++j) {
	    if ((tem = phif (vnew)) == 0.0)
		break;
	    if (tem < 0.0) {
		x = vnew;
	    } else
		y = vnew;
	    vold = vnew;
	    vnew = (x + y) * 0.5;
	    tem = fabs (vold - vnew);
	    if (tem <= errnd)
		break;
	    if (vold != 0.0) {
		if (fabs (tem / vold) <= errnd)
		    break;
	    }
	}
	if (j == 80) {
	    fprintf (drvui->flout, "No convergence.\n");
	    return (0);
	}
    }
    u[2] = vnew;
    c1 = b2 - vnew;
    c0 = b1 + c1 * vnew;
    tem = c1 * c1 + 4.0 * c0;
    if (tem < -0.0001) {
	fprintf (drvui->flout, "Complex Roots! tem = %f\n", tem);
	return (0);		/* error if complex roots */
    }
    if (tem < 0.0)
	tem = 0.0;
    tem = sqrt (tem);
    u[0] = 0.5 * (c1 - tem);
    u[1] = 0.5 * (c1 + tem);
    for (j = 0; j < 2; j++) {	/* sort roots in increasing order */
	if (u[j] > u[2]) {
	    tem = u[2];
	    u[2] = u[j];
	    u[j] = tem;
	}
    }

/* count multiple roots */
    tem = 100.0 * errnd;
    j = 0;
    if ((u[1] - u[0] - tem) < 0)
	j++;
    if ((u[2] - u[1] - tem) < 0)
	j += 2;
    if (j == 3) {
/* 3 equal roots */
	for (ii = 0; ii < 3; ii++) {
	    for (i = 0; i < 3; i++)
		vect[i][ii] = 0.0;
	    vect[ii][ii] = 1.0;
	    valu[ii] = (float) sqrt (u[ii] / tpi2);
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
	if (smax <= 0.0)
	    return (0);
	smax = (float) sqrt (smax);
	for (i = 0; i < 3; i++)
	    v[i] = b[i][imax] / smax;
	axeqb (a, z, v);
/* convert direct-space vector z into Cartesian system */
	for (i = 0; i < 3; i++) {
	    v[i] = 0.0;
	    for (j = 0; j < 3; j++)
		v[i] += drvui->b_mat[i][j] * z[j];
	}
	tem = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	tem = (float) sqrt (tem);
	for (i = 0; i < 3; i++)
	    Vect[ii][i] = v[i] / tem;	/* unit vector */
    }
    for (i = 0; i < 3; i++)
	valu[i] = (float) sqrt (u[i] / tpi2);
/* look for degenerate vectors */
    if (fabs (Vect[0][0] * Vect[1][0] + Vect[0][1] * Vect[1][1]
	      + Vect[0][2] * Vect[1][2]) > 0.001) {	/* set #1 = 2 x 3 */
	Vect[0][0] = Vect[1][1] * Vect[2][2] - Vect[2][1] * Vect[1][2];
	Vect[0][1] = Vect[2][0] * Vect[1][2] - Vect[1][0] * Vect[2][2];
	Vect[0][2] = Vect[1][0] * Vect[2][1] - Vect[2][0] * Vect[1][1];
    }
    if (fabs (Vect[2][0] * Vect[1][0] + Vect[2][1] * Vect[1][1]
	      + Vect[2][2] * Vect[1][2]) > 0.001) {	/* set #3 = 1 x 2 */
	Vect[2][0] = Vect[0][1] * Vect[1][2] - Vect[1][1] * Vect[0][2];
	Vect[2][1] = Vect[1][0] * Vect[0][2] - Vect[0][0] * Vect[1][2];
	Vect[2][2] = Vect[0][0] * Vect[1][1] - Vect[1][0] * Vect[0][1];
    }
    if (determinant (Vect) < 0.0) {
	for (i = 0; i < 3; i++)
	    for (j = 0; j < 3; j++)
		Vect[i][j] *= -1.0;
    }
    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    vect[i][j] = (float) Vect[i][j];
    return (1);
}

/* ************************************************************** */
/* ************************************************************** */

void
generate_lsq_planes (void)
// routine to calculate and display a least squares plane through a set of atoms
{
    int i, j, jj, k, l, m;

    float d;

    float glr, glg, glb;

    char thecolor[40];

    float c[3], s[6], vects[3][3], vals[3];

    double vx[50][15], vy[50][15], vz[50][15];

    double vsx[50][15], vsy[50][15], vsz[50][15];

    double t0[3], t1[3], t2[3], t3[3], p0[3], p1[3], p2[3], p3[3];

    int nsets = 0;

    double sx, sy, sz;

    int nx, add, nn[50];

    int numplanes = 0;

    double *pn0, *pn1, *pn2;

    double v1n, v2n, dot, angval, cosarg;

    if (drvui->nbplane == 1)
	return;

    pn0 = (double *) zalloc (sizeof (double) * drvui->nplane_alloc * 50);
    pn1 = (double *) zalloc (sizeof (double) * drvui->nplane_alloc * 50);
    pn2 = (double *) zalloc (sizeof (double) * drvui->nplane_alloc * 50);

    for (m = 1; m < drvui->nbplane; m++) {

	nsets = 0;

// first pass through atom list to find maximum number of planes to expect
// (numbers may vary for individual atoms due to view box clipping)

	for (i = 0; i <= natom; ++i) {

	    if (drvui->atoms[i].atom_fn != drvui->frame_no)
		continue;

	    add = 0;
	    for (j = 0; j < drvui->bplanes[m].nbatoms; j++) {	// only if part of l.s. plane
		if ((check_atom_name
		     (drvui->atoms[i].atom_l, drvui->bplanes[m].bplane_t[j]) == 1)
		    && (drvui->atoms[i].sv_atom_n == drvui->bplanes[m].bplane_n[j]
			|| drvui->bplanes[m].bplane_n[j] < 0))
		    add = 1;
	    }
	    if (add == 0)
		continue;

	    nvert = 0;
	    find_all_in_box (i);	// find all copies of this atom in the display box

	    if (nvert == 0) {
		Error_Box
		    ("At least one atom of the l.s. plane is not in the current display box");
		free (pn0);
		free (pn1);
		free (pn2);
		return;
	    }
	    if (nvert > nsets)
		nsets = nvert;
	}

	if (nsets > 50) {
	    free (pn0);
	    free (pn1);
	    free (pn2);
	    Error_Box ("Too many copies (max. 50) of a best fitting plane.");
	    return;
	}
// second pass through master atom list:
// assemble atom lists for each copy of the current plane in the view box

	nx = 0;

	for (i = 0; i <= natom; ++i) {

	    if (drvui->atoms[i].atom_fn != drvui->frame_no)
		continue;

	    add = 0;
	    for (j = 0; j < drvui->bplanes[m].nbatoms; j++) {	// only if part of l.s. plane
		if ((check_atom_name
		     (drvui->atoms[i].atom_l, drvui->bplanes[m].bplane_t[j]) == 1)
		    && (drvui->atoms[i].sv_atom_n == drvui->bplanes[m].bplane_n[j]
			|| drvui->bplanes[m].bplane_n[j] < 0))
		    add = 1;
	    }
	    if (add == 0)
		continue;

	    nvert = 0;
	    find_all_in_box (i);	// find all copies of this atom in the display box

	    if (nx == 0) {	// when lists are empty, add one copy to each
		for (jj = 0; jj < nsets; jj++) {
		    vx[jj][0] = o_vert[3 * jj] - origin[0];
		    vy[jj][0] = o_vert[3 * jj + 1] - origin[1];
		    vz[jj][0] = o_vert[3 * jj + 2] - origin[2];
		    vsx[jj][0] = s_vert[3 * jj];
		    vsy[jj][0] = s_vert[3 * jj + 1];
		    vsz[jj][0] = s_vert[3 * jj + 2];
		    nn[jj] = 1;
		}
		nx = 1;
	    }

	    for (jj = 0; jj < nvert; jj++) {	// find a list that already contains
		add = 0;	// atoms within bonding distance
		for (k = 0; k < nsets; k++) {
		    for (l = 0; l < nn[k]; l++) {	// for all atoms already in this list
			d = (float) ((s_vert[3 * jj] - vsx[k][l]) * (s_vert[3 * jj] -
								     vsx[k][l])
				     + (s_vert[3 * jj + 1] -
					vsy[k][l]) * (s_vert[3 * jj + 1] - vsy[k][l])
				     + (s_vert[3 * jj + 2] -
					vsz[k][l]) * (s_vert[3 * jj + 2] - vsz[k][l]));
			if (d < 0.05)
			    break;	// already in list
			if (d < printdist * printdist) {
			    int nxx = nn[k];

			    vx[k][nxx] = o_vert[3 * jj] - origin[0];
			    vy[k][nxx] = o_vert[3 * jj + 1] - origin[1];
			    vz[k][nxx] = o_vert[3 * jj + 2] - origin[2];
			    vsx[k][nxx] = s_vert[3 * jj];
			    vsy[k][nxx] = s_vert[3 * jj + 1];
			    vsz[k][nxx] = s_vert[3 * jj + 2];
			    nn[k]++;
			    add = 1;
			    break;
			}
		    }
		    if (add == 1)
			continue;	// add each atom to one list only
		}		// for all plane lists
	    }			// for all copies of this atom in box
	    nx++;
	    if (nx > 14) {
		free (pn0);
		free (pn1);
		free (pn2);
		Error_Box ("Cannot handle more than 15 atoms in a best fitting plane.");
		return;
	    }
	}			// loop over all atoms i in master list

// calculate and render the best fitting plane for each set of atoms

	for (i = 0; i < nsets; i++) {
	    nx = nn[i];
	    if (nx < drvui->bplanes[m].nbatoms)
		continue;	// drop incomplete (clipped) planes

	    c[0] = c[1] = c[2] = 0.;
	    sx = sy = sz = 0.;
	    s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = 0.;
	    for (k = 0; k < nx; k++) {
		c[0] += (float) vsx[i][k];
		c[1] += (float) vsy[i][k];
		c[2] += (float) vsz[i][k];
		sx += vx[i][k];
		sy += vy[i][k];
		sz += vz[i][k];
	    }
	    c[0] /= nx;
	    c[1] /= nx;
	    c[2] /= nx;
	    sx /= nx;
	    sy /= nx;
	    sz /= nx;
	    for (k = 0; k < nx; k++) {
		s[0] += (float) ((sx - vx[i][k]) * (sx - vx[i][k]));
		s[1] += (float) ((sy - vy[i][k]) * (sy - vy[i][k]));
		s[2] += (float) ((sz - vz[i][k]) * (sz - vz[i][k]));
		s[3] += (float) ((sx - vx[i][k]) * (sy - vy[i][k]));
		s[4] += (float) ((sx - vx[i][k]) * (sz - vz[i][k]));
		s[5] += (float) ((sy - vy[i][k]) * (sz - vz[i][k]));
	    }

	    if (p_eigen (s, vals, vects) == 0) {
		continue;
	    }
	    d = (float) sqrt (vects[0][0] * vects[0][0] + vects[0][1] * vects[0][1] +
			      vects[0][2] * vects[0][2]);
	    vects[0][0] /= d;
	    vects[0][1] /= d;
	    vects[0][2] /= d;
	    d = (float) sqrt (vects[1][0] * vects[1][0] + vects[1][1] * vects[1][1] +
			      vects[1][2] * vects[1][2]);
	    vects[1][0] /= d;
	    vects[1][1] /= d;
	    vects[1][2] /= d;
	    d = (float) sqrt (vects[2][0] * vects[2][0] + vects[2][1] * vects[2][1] +
			      vects[2][2] * vects[2][2]);
	    vects[2][0] /= d;
	    vects[2][1] /= d;
	    vects[2][2] /= d;
	    fprintf (drvui->flout,
		     "least squares plane %2d(%2d): center at    %f %f %f\n", m, i + 1,
		     c[0], c[1], c[2]);
	    fprintf (drvui->flout, "                            plane normal %f %f %f\n",
		     vects[0][0], vects[0][1], vects[0][2]);
	    fprintf (drvui->flout, "                            1st vector   %f %f %f\n",
		     vects[1][0], vects[1][1], vects[1][2]);
	    fprintf (drvui->flout, "                            2nd vector   %f %f %f\n",
		     vects[2][0], vects[2][1], vects[2][2]);

	    // apply scaling factors to generate initial cornerpoints of polygon 
	    t0[0] = c[0] + drvui->bplanes[m].bplane_d1 * vects[2][0];
	    t0[1] = c[1] + drvui->bplanes[m].bplane_d1 * vects[2][1];
	    t0[2] = c[2] + drvui->bplanes[m].bplane_d1 * vects[2][2];
	    t1[0] = c[0] + drvui->bplanes[m].bplane_d2 * vects[1][0];
	    t1[1] = c[1] + drvui->bplanes[m].bplane_d2 * vects[1][1];
	    t1[2] = c[2] + drvui->bplanes[m].bplane_d2 * vects[1][2];
	    t2[0] = c[0] - drvui->bplanes[m].bplane_d1 * vects[2][0];
	    t2[1] = c[1] - drvui->bplanes[m].bplane_d1 * vects[2][1];
	    t2[2] = c[2] - drvui->bplanes[m].bplane_d1 * vects[2][2];
	    t3[0] = c[0] - drvui->bplanes[m].bplane_d2 * vects[1][0];
	    t3[1] = c[1] - drvui->bplanes[m].bplane_d2 * vects[1][1];
	    t3[2] = c[2] - drvui->bplanes[m].bplane_d2 * vects[1][2];
	    // treat these as midpoints on the edges of the final rectangle
	    // to make width and height factors behave as expected 
	    for (jj = 0; jj < 3; jj++) {
		p0[jj] = c[jj] + (t0[jj] - t1[jj]) / 2.;
		p1[jj] = c[jj] + (t1[jj] - t2[jj]) / 2.;
		p2[jj] = c[jj] + (t2[jj] - t3[jj]) / 2.;
		p3[jj] = c[jj] + (t3[jj] - t0[jj]) / 2.;
	    }

	    strncpy (thecolor, drvui->bplanes[m].bplane_col, 39);
	    Transform_VRML_Color (thecolor);
	    sscanf (thecolor, "%f %f %f", &glr, &glg, &glb);

	    glPushMatrix ();
	    glDisable (GL_LIGHTING);
	    glColor3f (glr, glg, glb);
	    glBegin (GL_QUADS);
	    glNormal3f (c[0] + vects[0][0], c[1] + vects[0][1], c[2] + vects[0][2]);
	    glVertex3f ((float) p0[0], (float) p0[1], (float) p0[2]);
	    glVertex3f ((float) p1[0], (float) p1[1], (float) p1[2]);
	    glVertex3f ((float) p2[0], (float) p2[1], (float) p2[2]);
	    glVertex3f ((float) p3[0], (float) p3[1], (float) p3[2]);
	    glEnd ();
	    glEnable (GL_LIGHTING);
	    glPopMatrix ();

	    if (doPOV) {
		fprintf (drvui->fpoutp, "/* L.S. Plane */\n");
		fprintf (drvui->fpoutp, "triangle {< %8.5f, %8.5f, %8.5f>,\n",
			 c[0], c[1], c[2]);
		fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f>, <%8.5f, %8.5f, %8.5f>\n",
			 p0[0], p0[1], p0[2], p1[0], p1[1], p1[2]);
		fprintf (drvui->fpoutp, "  texture{pigment{color %s  }}\n",
			 drvui->bplanes[m].bplane_col);
		fprintf (drvui->fpoutp, " }\n");
		fprintf (drvui->fpoutp, "triangle {< %8.5f, %8.5f, %8.5f>,\n",
			 c[0], c[1], c[2]);
		fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f>, <%8.5f, %8.5f, %8.5f>\n",
			 p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
		fprintf (drvui->fpoutp, "  texture{pigment{color %s  }}\n",
			 drvui->bplanes[m].bplane_col);
		fprintf (drvui->fpoutp, " }\n");
		fprintf (drvui->fpoutp, "triangle {< %8.5f, %8.5f, %8.5f>,\n",
			 c[0], c[1], c[2]);
		fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f>, <%8.5f, %8.5f, %8.5f>\n",
			 p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]);
		fprintf (drvui->fpoutp, "  texture{pigment{color %s  }}\n",
			 drvui->bplanes[m].bplane_col);
		fprintf (drvui->fpoutp, " }\n");
		fprintf (drvui->fpoutp, "triangle {< %8.5f, %8.5f, %8.5f>,\n",
			 c[0], c[1], c[2]);
		fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f>, <%8.5f, %8.5f, %8.5f>\n",
			 p3[0], p3[1], p3[2], p0[0], p0[1], p0[2]);
		fprintf (drvui->fpoutp, "  texture{pigment{color %s  }}\n",
			 drvui->bplanes[m].bplane_col);
		fprintf (drvui->fpoutp, " }\n");
	    }

	    if (doAsy) {
		fprintf (drvui->fpouta, " // L.S. Plane \n");
		fprintf (drvui->fpouta, " draw(pic, surface ( (%8.5f, %8.5f, %8.5f)--",
			 p0[0], p0[1], p0[2]);
		fprintf (drvui->fpouta, "(%8.5f, %8.5f, %8.5f)--(%8.5f, %8.5f, %8.5f)--",
			 p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
		fprintf (drvui->fpouta, "(%8.5f, %8.5f, %8.5f)--cycle ),",
			 p3[0], p3[1], p3[2]);
		fprintf (drvui->fpouta, "rgb( %4.2f,%4.2f,%4.2f) );\n", glr,glg,glb);
	    }

	    if (doVrml) {
		fprintf (drvui->fpoutv, "# L.S. Plane\n");
		if (Vrml2) {
		    fprintf (drvui->fpoutv, " Shape {");
		    fprintf (drvui->fpoutv, "appearance Appearance {\n");
		    fprintf (drvui->fpoutv,
			     "   material Material {diffuseColor %s} \n \n", thecolor);
		} else {
		    fprintf (drvui->fpoutv, " Separator {\n        ");
		    fprintf (drvui->fpoutv, "  Material {\n");
		    fprintf (drvui->fpoutv, "   diffuseColor %s \n  }\n", thecolor);
		}
		if (Vrml2)
		    fprintf (drvui->fpoutv,
			     "}\n geometry IndexedFaceSet { coord Coordinate{ point [\n");
		else
		    fprintf (drvui->fpoutv, "Coordinate3 { point [\n");

		fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n", p0[0], p0[1], p0[2]);
		fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n", p1[0], p1[1], p1[2]);
		fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n", p2[0], p2[1], p2[2]);
		fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f]\n", p3[0], p3[1], p3[2]);

		if (Vrml2) {
		    fprintf (drvui->fpoutv, " }\n");
		    fprintf (drvui->fpoutv, " coordIndex [0,1,2,3,-1]\n");
		    fprintf (drvui->fpoutv, " solid FALSE\n");
		    fprintf (drvui->fpoutv, " convex TRUE\n");
		    fprintf (drvui->fpoutv, " creaseAngle 1.5708\n");
		    fprintf (drvui->fpoutv, "  }\n }\n");
		} else {
		    fprintf (drvui->fpoutv,
			     "}\n IndexedFaceSet { coordIndex [0,1,2,3,-1]\n");
		    fprintf (drvui->fpoutv, "  }\n }\n");
		}
	    }
	    pn0[50 * (m - 1) + i] = vects[0][0];
	    pn1[50 * (m - 1) + i] = vects[0][1];
	    pn2[50 * (m - 1) + i] = vects[0][2];
	    numplanes++;

	}			// for all l.s. planes
	if (nsets > 1) {
	    int idx = 50 * (m - 1);

	    fprintf (drvui->flout, "Angles between symmetry equivalent planes:\n");
	    for (i = 0; i < nsets; i++) {
		if (nn[i] < drvui->bplanes[m].nbatoms)
		    continue;	// drop incomplete (clipped) planes
		v1n =
		    pn0[idx + i] * pn0[idx + i] + pn1[idx + i] * pn1[idx + i] + pn2[idx +
										    i] *
		    pn2[idx + i];
		v1n = sqrt (v1n);
		for (j = i + 1; j < nsets; j++) {
		    if (nn[j] < drvui->bplanes[m].nbatoms)
			continue;	// drop incomplete (clipped) planes
		    v2n =
			pn0[idx + j] * pn0[idx + j] + pn1[idx + j] * pn1[idx + j] +
			pn2[idx + j] * pn2[idx + j];
		    v2n = sqrt (v2n);
		    dot =
			pn0[idx + i] * pn0[idx + j] + pn1[idx + i] * pn1[idx + j] +
			pn2[idx + i] * pn2[idx + j];
		    if (v1n * v2n < 1.e-6) {
			angval = 0.;
		    } else {
			cosarg = dot / (v1n * v2n);
			if (cosarg > 1.)
			    cosarg = 1.;
			if (cosarg < -1.)
			    cosarg = -1.;
			angval = RAD * acos (cosarg);
		    }
		    fprintf (drvui->flout, "Plane %d - Plane %d: %5.3f\n", i, j, angval);
		}
	    }
	}

    }				// for all bestplane keywords

    if (drvui->nbplane > 2) {
	fprintf (drvui->flout, "Angles between unique planes:\n");
	m = drvui->nbplane;
	for (i = 1; i < m; i++) {
	    v1n =
		pn0[50 * (i - 1)] * pn0[50 * (i - 1)] +
		pn1[50 * (i - 1)] * pn1[50 * (i - 1)] +
		pn2[50 * (i - 1)] * pn2[50 * (i - 1)];
	    v1n = sqrt (v1n);
	    for (j = i + 1; j < m; j++) {
		v2n =
		    pn0[50 * (j - 1)] * pn0[50 * (j - 1)] +
		    pn1[50 * (j - 1)] * pn1[50 * (j - 1)] +
		    pn2[50 * (j - 1)] * pn2[50 * (j - 1)];
		v2n = sqrt (v2n);
		dot =
		    pn0[50 * (i - 1)] * pn0[50 * (j - 1)] +
		    pn1[50 * (i - 1)] * pn1[50 * (j - 1)] +
		    pn2[50 * (i - 1)] * pn2[50 * (j - 1)];
		if (v1n * v2n < 1.e-6) {
		    angval = 0.;
		} else {
		    cosarg = dot / (v1n * v2n);
		    if (cosarg > 1.)
			cosarg = 1.;
		    if (cosarg < -1.)
			cosarg = -1.;
		    angval = RAD * acos (cosarg);
		}
		fprintf (drvui->flout, "Plane %d - Plane %d: %5.3f\n", i, j, angval);
	    }
	}
    }

    fprintf (drvui->flout, "Generated %d best fitting planes.\n", numplanes);
    free (pn0);
    free (pn1);
    free (pn2);
}

/* ************************************************************** */
/* ************************************************************** */

void
generate_aimsurf (void)
// routine to display polyhedra corresponding to atomic AIM basins
{
    int i, j, k, l, m;

    int ii, jj, kk, ku;

    int kk0 = 0, ku0 = 0;

    float glr, glg, glb;

    int numhulls = 0;

    int nv;

    int o;

    float p1[3], p2[3], p3[3], p4[3];

    float sp1[3], sp2[3], sp3[3], sp4[3];

    float rotmat[3][3];

    char surf_col_p[40];

    char surf_col_v[40];

    if (drvui->nsurf == 1)
	return;

    for (i = 1; i < drvui->nsurf; i++) {

	strcpy (surf_col_p, drvui->surfcolor[i]);
	strcpy (surf_col_v, surf_col_p);
	Transform_POV_Color (surf_col_p);
	Transform_VRML_Color (surf_col_v);
	sscanf (surf_col_v, "%f %f %f", &glr, &glg, &glb);


	for (j = 0; j <= natom; ++j) {

	    if (drvui->atoms[j].atom_fn != drvui->frame_no)
		continue;

	    if ((check_atom_name (drvui->atoms[j].atom_l, drvui->surfatom[i]) == 1)
		&& (drvui->atoms[j].sv_atom_n == drvui->surfnum[i])) {

		nvert = 0;
		find_all_in_box (j);	// find all copies of this atom in the display box

		if (nvert == 0) {
		    break;
		}

		nv = drvui->ntet[i] * drvui->nphi[i] - 1;

		if (drvui->surftype[i] == 0)
		    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
		else if (drvui->surftype[i] == 1)
		    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
		glPushMatrix ();
		glColor3f (glr, glg, glb);
		glDisable (GL_LIGHTING);
		if (drvui->surftype[i] < 2)
		    glBegin (GL_QUADS);
		else
		    glBegin (GL_POINTS);

		for (k = 0; k < nvert; k++) {

		    o = vert_sym_no[k];
		    for (l = 0; l < 3; l++)
			for (m = 0; m < 3; m++)
			    rotmat[l][m] = (float) drvui->ss[o][l][m];
		    if (vert_sym_nos[k] < 0) {	// inversion
			rotmat[0][0] *= -1;
			rotmat[1][1] *= -1;
			rotmat[2][2] *= -1;
		    }

		    kk = -2;
		    for (ii = 0; ii < drvui->ntet[i] - 1; ii++) {
			kk++;
			for (jj = 0; jj < drvui->nphi[i] - 1; jj++) {
			    kk++;
			    ku = kk + 1 + drvui->nphi[i];
			    if (jj == 0) {
				kk0 = kk;
				ku0 = ku;
			    }
			    if (ku > nv)
				continue;
			    p1[0] = drvui->surfx[i][kk];
			    p1[1] = drvui->surfy[i][kk];
			    p1[2] = drvui->surfz[i][kk];
			    p2[0] = drvui->surfx[i][kk + 1];
			    p2[1] = drvui->surfy[i][kk + 1];
			    p2[2] = drvui->surfz[i][kk + 1];
			    p3[0] = drvui->surfx[i][ku];
			    p3[1] = drvui->surfy[i][ku];
			    p3[2] = drvui->surfz[i][ku];
			    p4[0] = drvui->surfx[i][kk + drvui->nphi[i]];
			    p4[1] = drvui->surfy[i][kk + drvui->nphi[i]];
			    p4[2] = drvui->surfz[i][kk + drvui->nphi[i]];

			    sp1[0] = sp1[1] = sp1[2] = 0.;
			    sp2[0] = sp2[1] = sp2[2] = 0.;
			    sp3[0] = sp3[1] = sp3[2] = 0.;
			    sp4[0] = sp4[1] = sp4[2] = 0.;

			    for (l = 0; l < 3; ++l) {
				for (m = 0; m < 3; ++m) {
				    sp1[l] += rotmat[l][m] * p1[m];
				    sp2[l] += rotmat[l][m] * p2[m];
				    sp3[l] += rotmat[l][m] * p3[m];
				    sp4[l] += rotmat[l][m] * p4[m];
				}
			    }
			    if (drvui->surftype[i] == 2) {
				glVertex3f (s_vert[3 * k] +
					    (sp1[0] + sp2[0] + sp3[0] + sp3[0]) / 4.0f,
					    s_vert[3 * k + 1] + (sp1[1] + sp2[1] +
								 sp3[1] + sp4[1]) / 4.0f,
					    s_vert[3 * k + 2] + (sp1[2] + sp2[2] +
								 sp3[2] + sp4[2]) / 4.0f);
			    } else {
				glVertex3f (s_vert[3 * k] + sp1[0],
					    s_vert[3 * k + 1] + sp1[1],
					    s_vert[3 * k + 2] + sp1[2]);
				glVertex3f (s_vert[3 * k] + sp2[0],
					    s_vert[3 * k + 1] + sp2[1],
					    s_vert[3 * k + 2] + sp2[2]);
				glVertex3f (s_vert[3 * k] + sp3[0],
					    s_vert[3 * k + 1] + sp3[1],
					    s_vert[3 * k + 2] + sp3[2]);
				glVertex3f (s_vert[3 * k] + sp4[0],
					    s_vert[3 * k + 1] + sp4[1],
					    s_vert[3 * k + 2] + sp4[2]);
			    }
			    if (doPOV) {
				if (drvui->surftype[i] == 0) {
				    fprintf (drvui->fpoutp,
					     " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
					     s_vert[3 * k] + sp1[0],
					     s_vert[3 * k + 1] + sp1[1],
					     s_vert[3 * k + 2] + sp1[2],
					     s_vert[3 * k] + sp2[0],
					     s_vert[3 * k + 1] + sp2[1],
					     s_vert[3 * k + 2] + sp2[2], 0.0001 * Scale);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n }\n",
					     surf_col_p);
				    fprintf (drvui->fpoutp,
					     " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
					     s_vert[3 * k] + sp2[0],
					     s_vert[3 * k + 1] + sp2[1],
					     s_vert[3 * k + 2] + sp2[2],
					     s_vert[3 * k] + sp3[0],
					     s_vert[3 * k + 1] + sp3[1],
					     s_vert[3 * k + 2] + sp3[2], 0.0001 * Scale);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n }\n",
					     surf_col_p);
				    fprintf (drvui->fpoutp,
					     " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
					     s_vert[3 * k] + sp3[0],
					     s_vert[3 * k + 1] + sp3[1],
					     s_vert[3 * k + 2] + sp3[2],
					     s_vert[3 * k] + sp4[0],
					     s_vert[3 * k + 1] + sp4[1],
					     s_vert[3 * k + 2] + sp4[2], 0.0001 * Scale);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n }\n",
					     surf_col_p);
				    fprintf (drvui->fpoutp,
					     " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
					     s_vert[3 * k] + sp4[0],
					     s_vert[3 * k + 1] + sp4[1],
					     s_vert[3 * k + 2] + sp4[2],
					     s_vert[3 * k] + sp1[0],
					     s_vert[3 * k + 1] + sp1[1],
					     s_vert[3 * k + 2] + sp1[2], 0.0001 * Scale);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n }\n",
					     surf_col_p);
				} else if (drvui->surftype[i] == 1) {
				    fprintf (drvui->fpoutp,
					     " triangle { <%f,%f,%f>,<%f,%f,%f>,<%f,%f,%f>\n",
					     s_vert[3 * k] + sp1[0],
					     s_vert[3 * k + 1] + sp1[1],
					     s_vert[3 * k + 2] + sp1[2],
					     s_vert[3 * k] + sp2[0],
					     s_vert[3 * k + 1] + sp2[1],
					     s_vert[3 * k + 2] + sp2[2],
					     s_vert[3 * k] + sp3[0],
					     s_vert[3 * k + 1] + sp3[1],
					     s_vert[3 * k + 2] + sp3[2]);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n }\n",
					     surf_col_p);
				    fprintf (drvui->fpoutp,
					     " triangle { <%f,%f,%f>,<%f,%f,%f>,<%f,%f,%f>\n",
					     s_vert[3 * k] + sp3[0],
					     s_vert[3 * k + 1] + sp3[1],
					     s_vert[3 * k + 2] + sp3[2],
					     s_vert[3 * k] + sp4[0],
					     s_vert[3 * k + 1] + sp4[1],
					     s_vert[3 * k + 2] + sp4[2],
					     s_vert[3 * k] + sp1[0],
					     s_vert[3 * k + 1] + sp1[1],
					     s_vert[3 * k + 2] + sp1[2]);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n }\n",
					     surf_col_p);
				} else {
				    fprintf (drvui->fpoutp, " sphere { <%f,%f,%f>,0.03\n",
					     s_vert[3 * k] + (sp1[0] + sp2[0] + sp3[0] +
							      sp4[0]) / 4.0f,
					     s_vert[3 * k + 1] + (sp1[1] + sp2[1] +
								  sp3[1] + sp4[1]) / 4.0f,
					     s_vert[3 * k + 2] + (sp1[2] + sp2[2] +
								  sp3[2] +
								  sp4[2]) / 4.0f);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n }\n",
					     surf_col_p);
				}
			    }
			    if (doAsy) {
				if (drvui->surftype[i] == 0) {
				    fprintf (drvui->fpouta,
					     " draw(pic, (%f,%f,%f)--(%f,%f,%f)--\n",
					     s_vert[3 * k] + sp1[0],
					     s_vert[3 * k + 1] + sp1[1],
					     s_vert[3 * k + 2] + sp1[2],
					     s_vert[3 * k] + sp2[0],
					     s_vert[3 * k + 1] + sp2[1],
					     s_vert[3 * k + 2] + sp2[2]);
				    fprintf (drvui->fpouta,
					     "(%f,%f,%f)--(%f,%f,%f)--cycle,\n",
					     s_vert[3 * k] + sp3[0],
					     s_vert[3 * k + 1] + sp3[1],
					     s_vert[3 * k + 2] + sp3[2],
					     s_vert[3 * k] + sp4[0],
					     s_vert[3 * k + 1] + sp4[1],
					     s_vert[3 * k + 2] + sp4[2]);
				    fprintf (drvui->fpouta,
					     " rgb(%.2f,%.2f,%.2f)+linewidth(%.2f));\n",
					     glr, glg, glb, 0.0001 * Scale);
				} else if (drvui->surftype[i] == 1) {
				    fprintf (drvui->fpouta,
					     " draw(pic, surface( (%f,%f,%f)--(%f,%f,%f)--\n",
					     s_vert[3 * k] + sp1[0],
					     s_vert[3 * k + 1] + sp1[1],
					     s_vert[3 * k + 2] + sp1[2],
					     s_vert[3 * k] + sp2[0],
					     s_vert[3 * k + 1] + sp2[1],
					     s_vert[3 * k + 2] + sp2[2]);
				    fprintf (drvui->fpouta,
					     "(%f,%f,%f)--(%f,%f,%f)--cycle),\n",
					     s_vert[3 * k] + sp3[0],
					     s_vert[3 * k + 1] + sp3[1],
					     s_vert[3 * k + 2] + sp3[2],
					     s_vert[3 * k] + sp4[0],
					     s_vert[3 * k + 1] + sp4[1],
					     s_vert[3 * k + 2] + sp4[2]);
				    fprintf (drvui->fpouta,
					     " rgb(%.2f,%.2f,%.2f)+linewidth(%.2f));\n",
					     glr, glg, glb, 0.0001 * Scale);
				} else {
				    fprintf (drvui->fpouta, " draw(pic, shift(%8.5f, %8.5f, %8.5f)*scale3(%.2f)*unitsphere,",
					     s_vert[3 * k] + (sp1[0] + sp2[0] + sp3[0] +
							      sp4[0]) / 4.0f,
					     s_vert[3 * k + 1] + (sp1[1] + sp2[1] +
								  sp3[1] + sp4[1]) / 4.0f,
					     s_vert[3 * k + 2] + (sp1[2] + sp2[2] +
								  sp3[2] +
								  sp4[2]) / 4.0f, 0.03);
				    fprintf (drvui->fpouta,
					     " rgb(%.2f, %.2f, %.2f));\n", 
					     glr, glg, glb);
				}
			    }
			    if (doVrml) {
				if (!Vrml2) {
				    fprintf (drvui->fpoutv, " Separator {\n");
				    fprintf (drvui->fpoutv,
					     " Material {  diffuseColor %s }\n",
					     surf_col_v);
				    fprintf (drvui->fpoutv, "Coordinate3 { \n point[\n");
				} else {
				    fprintf (drvui->fpoutv, " Shape {\n");
				    if (drvui->surftype[i] == 0)
					fprintf (drvui->fpoutv,
						 "  geometry IndexedLineSet { coord Coordinate{ point [\n");
				    else if (drvui->surftype[i] == 1)
					fprintf (drvui->fpoutv,
						 "  geometry IndexedFaceSet { coord Coordinate{ point [\n");
				    else
					fprintf (drvui->fpoutv,
						 "  geometry PointSet { coord Coordinate{ point [\n");
				}
				if (drvui->surftype[i] == 2) {
				    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
					     s_vert[3 * k] + (sp1[0] + sp2[0] + sp3[0] +
							      sp4[0]) / 4.0f,
					     s_vert[3 * k + 1] + (sp1[0] + sp1[1] +
								  sp2[1] + sp4[1]) / 4.0f,
					     s_vert[3 * k + 2] + (sp1[2] + sp2[2] +
								  sp3[2] +
								  sp4[2]) / 4.0f);
				} else {
				    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
					     s_vert[3 * k] + sp1[0],
					     s_vert[3 * k + 1] + sp1[1],
					     s_vert[3 * k + 2] + sp1[2]);
				    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
					     s_vert[3 * k] + sp2[0],
					     s_vert[3 * k + 1] + sp2[1],
					     s_vert[3 * k + 2] + sp2[2]);
				    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
					     s_vert[3 * k] + sp3[0],
					     s_vert[3 * k + 1] + sp3[1],
					     s_vert[3 * k + 2] + sp3[2]);
				    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f]\n",
					     s_vert[3 * k] + sp4[0],
					     s_vert[3 * k + 1] + sp4[1],
					     s_vert[3 * k + 2] + sp4[2]);
				}
				if (Vrml2) {
				    if (drvui->surftype[i] < 2)
					fprintf (drvui->fpoutv,
						 "  }\n  coordIndex [ 0,1,2,3,-1]\n color Color { color [%s]"
						 "}\n colorIndex [0]\n colorPerVertex FALSE\n }\n }\n",
						 surf_col_v);
				    else
					fprintf (drvui->fpoutv,
						 "  }\n }\n appearance Appearance {material Material { emissiveColor %s }}\n}\n",
						 surf_col_v);
				} else {
				    if (drvui->surftype[i] == 0)
					fprintf (drvui->fpoutv,
						 "}\n IndexedLineSet { coordIndex [0,1,2,3,-1] }\n}\n");
				    else if (drvui->surftype[i] == 1)
					fprintf (drvui->fpoutv,
						 "}\n IndexedFaceSet { coordIndex [0,1,2,3,-1] }\n}\n");
				    else
					fprintf (drvui->fpoutv,
						 "}\n PointSet { \n startIndex 0\n numPoints -1\n }\n}\n");
				}
			    }
			}	// loop on jj
			if (drvui->surftype[i] < 2) {	//additional step to close the mesh
			    p1[0] = drvui->surfx[i][kk + 1];
			    p1[1] = drvui->surfy[i][kk + 1];
			    p1[2] = drvui->surfz[i][kk + 1];
			    p2[0] = drvui->surfx[i][kk0];
			    p2[1] = drvui->surfy[i][kk0];
			    p2[2] = drvui->surfz[i][kk0];
			    p3[0] = drvui->surfx[i][ku0 - 1];
			    p3[1] = drvui->surfy[i][ku0 - 1];
			    p3[2] = drvui->surfz[i][ku0 - 1];
			    p4[0] = drvui->surfx[i][kk + 1 + drvui->nphi[i]];
			    p4[1] = drvui->surfy[i][kk + 1 + drvui->nphi[i]];
			    p4[2] = drvui->surfz[i][kk + 1 + drvui->nphi[i]];

			    sp1[0] = sp1[1] = sp1[2] = 0.;
			    sp2[0] = sp2[1] = sp2[2] = 0.;
			    sp3[0] = sp3[1] = sp3[2] = 0.;
			    sp4[0] = sp4[1] = sp4[2] = 0.;

			    for (l = 0; l < 3; ++l) {
				for (m = 0; m < 3; ++m) {
				    sp1[l] += rotmat[l][m] * p1[m];
				    sp2[l] += rotmat[l][m] * p2[m];
				    sp3[l] += rotmat[l][m] * p3[m];
				    sp4[l] += rotmat[l][m] * p4[m];
				}
			    }
			    glVertex3f (s_vert[3 * k] + sp1[0],
					s_vert[3 * k + 1] + sp1[1],
					s_vert[3 * k + 2] + sp1[2]);
			    glVertex3f (s_vert[3 * k] + sp2[0],
					s_vert[3 * k + 1] + sp2[1],
					s_vert[3 * k + 2] + sp2[2]);
			    glVertex3f (s_vert[3 * k] + sp3[0],
					s_vert[3 * k + 1] + sp3[1],
					s_vert[3 * k + 2] + sp3[2]);
			    glVertex3f (s_vert[3 * k] + sp4[0],
					s_vert[3 * k + 1] + sp4[1],
					s_vert[3 * k + 2] + sp4[2]);

			    if (doPOV) {
				if (drvui->surftype[i] == 0) {
				    fprintf (drvui->fpoutp,
					     " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
					     s_vert[3 * k] + sp1[0],
					     s_vert[3 * k + 1] + sp1[1],
					     s_vert[3 * k + 2] + sp1[2],
					     s_vert[3 * k] + sp2[0],
					     s_vert[3 * k + 1] + sp2[1],
					     s_vert[3 * k + 2] + sp2[2], 0.0001 * Scale);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n }\n",
					     surf_col_p);
				    fprintf (drvui->fpoutp,
					     " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
					     s_vert[3 * k] + sp2[0],
					     s_vert[3 * k + 1] + sp2[1],
					     s_vert[3 * k + 2] + sp2[2],
					     s_vert[3 * k] + sp3[0],
					     s_vert[3 * k + 1] + sp3[1],
					     s_vert[3 * k + 2] + sp3[2], 0.0001 * Scale);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n }\n",
					     surf_col_p);
				    fprintf (drvui->fpoutp,
					     " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
					     s_vert[3 * k] + sp3[0],
					     s_vert[3 * k + 1] + sp3[1],
					     s_vert[3 * k + 2] + sp3[2],
					     s_vert[3 * k] + sp4[0],
					     s_vert[3 * k + 1] + sp4[1],
					     s_vert[3 * k + 2] + sp4[2], 0.0001 * Scale);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n }\n",
					     surf_col_p);
				    fprintf (drvui->fpoutp,
					     " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
					     s_vert[3 * k] + sp4[0],
					     s_vert[3 * k + 1] + sp4[1],
					     s_vert[3 * k + 2] + sp4[2],
					     s_vert[3 * k] + sp1[0],
					     s_vert[3 * k + 1] + sp1[1],
					     s_vert[3 * k + 2] + sp1[2], 0.0001 * Scale);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n }\n",
					     surf_col_p);
				} else {
				    fprintf (drvui->fpoutp,
					     " triangle { <%f,%f,%f>,<%f,%f,%f>,<%f,%f,%f>\n",
					     s_vert[3 * k] + sp1[0],
					     s_vert[3 * k + 1] + sp1[1],
					     s_vert[3 * k + 2] + sp1[2],
					     s_vert[3 * k] + sp2[0],
					     s_vert[3 * k + 1] + sp2[1],
					     s_vert[3 * k + 2] + sp2[2],
					     s_vert[3 * k] + sp3[0],
					     s_vert[3 * k + 1] + sp3[1],
					     s_vert[3 * k + 2] + sp3[2]);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n }\n",
					     surf_col_p);
				    fprintf (drvui->fpoutp,
					     " triangle { <%f,%f,%f>,<%f,%f,%f>,<%f,%f,%f>\n",
					     s_vert[3 * k] + sp3[0],
					     s_vert[3 * k + 1] + sp3[1],
					     s_vert[3 * k + 2] + sp3[2],
					     s_vert[3 * k] + sp4[0],
					     s_vert[3 * k + 1] + sp4[1],
					     s_vert[3 * k + 2] + sp4[2],
					     s_vert[3 * k] + sp1[0],
					     s_vert[3 * k + 1] + sp1[1],
					     s_vert[3 * k + 2] + sp1[2]);
				    fprintf (drvui->fpoutp,
					     "  texture{pigment{color %s  }}\n }\n",
					     surf_col_p);
				}
			    }
			    if (doAsy) {
				if (drvui->surftype[i] == 0) {
				    fprintf (drvui->fpouta,
					     " draw(pic, (%f,%f,%f)--(%f,%f,%f)--\n",
					     s_vert[3 * k] + sp1[0],
					     s_vert[3 * k + 1] + sp1[1],
					     s_vert[3 * k + 2] + sp1[2],
					     s_vert[3 * k] + sp2[0],
					     s_vert[3 * k + 1] + sp2[1],
					     s_vert[3 * k + 2] + sp2[2]);
				    fprintf (drvui->fpouta,
					     "(%f,%f,%f)--(%f,%f,%f)--cycle,\n",
					     s_vert[3 * k] + sp3[0],
					     s_vert[3 * k + 1] + sp3[1],
					     s_vert[3 * k + 2] + sp3[2],
					     s_vert[3 * k] + sp4[0],
					     s_vert[3 * k + 1] + sp4[1],
					     s_vert[3 * k + 2] + sp4[2]);
				    fprintf (drvui->fpouta,
					     " rgb(%.2f,%.2f,%.2f)+linewidth(%.2f));\n",
					     glr, glg, glb, 0.0001 * Scale);
				} else if (drvui->surftype[i] == 1) {
				    fprintf (drvui->fpouta,
					     " draw(pic, surface( (%f,%f,%f)--(%f,%f,%f)--\n",
					     s_vert[3 * k] + sp1[0],
					     s_vert[3 * k + 1] + sp1[1],
					     s_vert[3 * k + 2] + sp1[2],
					     s_vert[3 * k] + sp2[0],
					     s_vert[3 * k + 1] + sp2[1],
					     s_vert[3 * k + 2] + sp2[2]);
				    fprintf (drvui->fpouta,
					     "(%f,%f,%f)--(%f,%f,%f)--cycle),\n",
					     s_vert[3 * k] + sp3[0],
					     s_vert[3 * k + 1] + sp3[1],
					     s_vert[3 * k + 2] + sp3[2],
					     s_vert[3 * k] + sp4[0],
					     s_vert[3 * k + 1] + sp4[1],
					     s_vert[3 * k + 2] + sp4[2]);
				    fprintf (drvui->fpouta,
					     " rgb(%.2f,%.2f,%.2f)+linewidth(%.2f));\n",
					     glr, glg, glb, 0.0001 * Scale);
				} else {
				    fprintf (drvui->fpouta, " draw(pic, shift(%8.5f, %8.5f, %8.5f)*scale3(%.2f)*unitsphere,",
					     s_vert[3 * k] + (sp1[0] + sp2[0] + sp3[0] +
							      sp4[0]) / 4.0f,
					     s_vert[3 * k + 1] + (sp1[1] + sp2[1] +
								  sp3[1] + sp4[1]) / 4.0f,
					     s_vert[3 * k + 2] + (sp1[2] + sp2[2] +
								  sp3[2] +
								  sp4[2]) / 4.0f, 0.03);
				    fprintf (drvui->fpouta,
					     " rgb(%.2f, %.2f, %.2f));\n", 
					     glr, glg, glb);
				}
			    }
			    if (doVrml) {
				if (!Vrml2) {
				    fprintf (drvui->fpoutv, " Separator {\n");
				    fprintf (drvui->fpoutv,
					     " Material {  diffuseColor %s }\n",
					     surf_col_v);
				    fprintf (drvui->fpoutv, "Coordinate3 { \n point[\n");
				} else {
				    fprintf (drvui->fpoutv, " Shape {\n");
				    if (drvui->surftype[i] == 0)
					fprintf (drvui->fpoutv,
						 "  geometry IndexedLineSet { coord Coordinate{ point [\n");
				    else
					fprintf (drvui->fpoutv,
						 "  geometry IndexedFaceSet { coord Coordinate{ point [\n");
				}

				fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
					 s_vert[3 * k] + sp1[0],
					 s_vert[3 * k + 1] + sp1[1],
					 s_vert[3 * k + 2] + sp1[2]);
				fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
					 s_vert[3 * k] + sp2[0],
					 s_vert[3 * k + 1] + sp2[1],
					 s_vert[3 * k + 2] + sp2[2]);
				fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
					 s_vert[3 * k] + sp3[0],
					 s_vert[3 * k + 1] + sp3[1],
					 s_vert[3 * k + 2] + sp3[2]);
				fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f]\n",
					 s_vert[3 * k] + sp4[0],
					 s_vert[3 * k + 1] + sp4[1],
					 s_vert[3 * k + 2] + sp4[2]);

				if (Vrml2) {
				    fprintf (drvui->fpoutv,
					     "  }\n  coordIndex [ 0,1,2,3,-1]\n color Color { color [%s]"
					     "}\n colorIndex [0]\n colorPerVertex FALSE\n }\n }\n",
					     surf_col_v);
				} else {
				    if (drvui->surftype[i] == 0)
					fprintf (drvui->fpoutv,
						 "}\n IndexedLineSet { coordIndex [0,1,2,3,-1] }\n}\n");
				    else
					fprintf (drvui->fpoutv,
						 "}\n IndexedFaceSet { coordIndex [0,1,2,3,-1] }\n}\n");
				}
			    }
			}	// additional step to close the mesh
		    }		// loop on ii - for all vertices of this hull

		    numhulls++;

		}		//for all equivalent sites k
		glEnd ();
		glEnable (GL_LIGHTING);
		glPopMatrix ();
		glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

	    }			// if matching entry found
	}			// for all atoms j
    }				// for all surfaces i

    fprintf (drvui->flout, "Generated %d AIM surface hulls.\n", numhulls);
}

void
calculate_voids (void)
{

    if (drvui->voidflag <= 0)
	return;
    switch (drvui->voidflag) {
    case 1:
	calc_simplevoids ();
	break;
    case 2:
	calculate_msms ();
	break;
    case 3:
	drvui->voidflag = -3;
	break;
    default:
	break;
    }
    return;
}

void
generate_voids (void)
{
    switch (drvui->voidflag) {
    case -1:
	generate_simplevoids ();
	return;
    case -2:
	generate_msms ();
	return;
    case -3:
	calculate_sas ();
	return;
    default:
	return;
    }
}

// classical kitaigorodskij algorithm - apply grid and check 
// each gridpoint for intersection with an atom of the structure
void
calc_simplevoids (void)
{
    int i, j, k, l;

    int l1, l2;

    int onvert;

    float p[3], p2[3];

    float d2, dlim;

    float vert2[3];

    float saved_boxlim[3];

    float saved_crystlim[6];

    float radius;

    float *fp; 

    int *ip;

    int saved_nvert = NvertM;

    ip = (int *) realloc (drvui->orig_atom_no, (4 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	drvui->orig_atom_no = ip;
    }
    ip = (int *) realloc (vert_sym_no, (4 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	vert_sym_no = ip;
    }
    ip = (int *) realloc (vert_sym_nos, (4 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	vert_sym_nos = ip;
    }
    fp = (float *) realloc (o_vert, (12 * NvertM) * sizeof (float));
    if (!fp) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	o_vert = fp;
    }
    fp = (float *) realloc (o_vert_nm, (12 * NvertM) * sizeof (float));
    if (!fp) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	o_vert_nm = fp;
    }
    nvert = 0;			/* clear the vertex list */
    for (j = 0; j < 3; j++) {
	saved_boxlim[j] = boxlim[j];
	saved_crystlim[j] = drvui->frames[drvui->frame_no].cryst_lim[j];
	saved_crystlim[j + 3] = drvui->frames[drvui->frame_no].cryst_lim[j + 3];
	drvui->frames[drvui->frame_no].cryst_lim[j] = -0.25;
	drvui->frames[drvui->frame_no].cryst_lim[j + 3] = 1.25;
	boxlim[j] += 10.;
    }
    Progress_Window (-1, "Computing Cavity Volumes", 130.0f);
    Fl::flush ();
    NvertM = 1;
    build_box_contents ();

    ip = (int *) realloc (drvui->orig_atom_no, (2 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	drvui->orig_atom_no = ip;
    }
    ip = (int *) realloc (vert_sym_no, (2 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	vert_sym_no = ip;
    }
    ip = (int *) realloc (vert_sym_nos, (2 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	vert_sym_nos = ip;
    }
    fp = (float *) realloc (o_vert, 6 * (NvertM + 20) * sizeof (float));
    if (!fp) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	o_vert = fp;
    }
    fp = (float *) realloc (o_vert_nm, 6 * (NvertM + 20) * sizeof (float));
    if (!fp) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	o_vert_nm = fp;
    }
    float *radi = (float *) malloc (NvertM * sizeof (float));
    if (!radi) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    }
    
    for (j = 0; j < natom; ++j) {	/* loop through atoms */
	onvert = nvert;
	if (drvui->atoms[j].atom_fn != drvui->frame_no)
	    continue;
	radius = drvui->atoms[j].radius;
	find_all_in_box (j);
	for (k = onvert; k < nvert; k++)
	    radi[k] = radius;
    }

    float stepx = 1.f / (float) drvui->voidgrid[0];

    float stepy = 1.f / (float) drvui->voidgrid[1];

    float stepz = 1.f / (float) drvui->voidgrid[2];

    Progress_Window (0, NULL, 1.0f);

    for (i = 0; i < drvui->voidgrid[0]; i++) {
	p[0] = stepx * (float) i;

	Progress_Window (0, NULL, p[0] * 100.0f);

	for (j = 0; j < drvui->voidgrid[1]; j++) {
	    p[1] = stepy * (float) j;
	    for (k = 0; k < drvui->voidgrid[2]; k++) {
		p[2] = stepz * (float) k;
		drvui->voidmap[i][j][k] = (char) 1;

		for (l1 = 0; l1 < 3; ++l1) {	/* convert gridpoint coordinates to Cartesian */
		    vert2[l1] = 0.0f;
		    for (l2 = 0; l2 < 3; ++l2)
			vert2[l1] += (float) drvui->b_mat[l1][l2] * (p[l2] - origin[l2]);
		}
		p2[0] = vert2[0];
		p2[1] = vert2[1];
		p2[2] = vert2[2];
		for (l = 1; l < nvert; l++) {
		    dlim = (drvui->probesize + radi[l]) * (drvui->probesize + radi[l]);

		    d2 = (p2[0] - s_vert[3 * l]) * (p2[0] - s_vert[3 * l])
			+ (p2[1] - s_vert[3 * l + 1]) * (p2[1] - s_vert[3 * l + 1])
			+ (p2[2] - s_vert[3 * l + 2]) * (p2[2] - s_vert[3 * l + 2]);
		    if (d2 <= dlim) {
			drvui->voidmap[i][j][k] = (char) 0;
			break;
		    }
		}
	    }
	}
    }

    for (j = 0; j < 3; j++) {
	//boxlim[j]=saved_boxlim[j];
	drvui->frames[drvui->frame_no].cryst_lim[j] = saved_crystlim[j];
	drvui->frames[drvui->frame_no].cryst_lim[j + 3] = saved_crystlim[j + 3];
    }

    ip = (int *) realloc (drvui->orig_atom_no, (2 * saved_nvert) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	drvui->orig_atom_no = ip;
    }
    ip = (int *) realloc (vert_sym_no, (2 * saved_nvert) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	vert_sym_no = ip;
    }
    ip = (int *) realloc (vert_sym_nos, (2 * saved_nvert) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	vert_sym_nos = ip;
    }

    NvertM = 1;
    build_box_contents ();

    Progress_Window (0, NULL, 110.0f);

    l = 0;
    for (i = 0; i < drvui->voidgrid[0]; i++) {
	for (j = 0; j < drvui->voidgrid[1]; j++) {
	    for (k = 0; k < drvui->voidgrid[2]; k++) {
		if (drvui->voidmap[i][j][k] == 1) {
		    l++;
		}
	    }
	}
    }

    Progress_Window (0, NULL, 120.0f);
    for (i = 1; i < drvui->voidgrid[0] - 1; i++) {
	for (j = 1; j < drvui->voidgrid[1] - 1; j++) {
	    for (k = 1; k < drvui->voidgrid[2] - 1; k++) {
		if (drvui->voidmap[i][j][k] == 0)
		    continue;
		if (drvui->voidmap[i - 1][j][k] == 0)
		    continue;
		if (drvui->voidmap[i + 1][j][k] == 0)
		    continue;
		if (drvui->voidmap[i][j - 1][k] == 0)
		    continue;
		if (drvui->voidmap[i][j + 1][k] == 0)
		    continue;
		if (drvui->voidmap[i][j][k - 1] == 0)
		    continue;
		if (drvui->voidmap[i][j][k + 1] == 0)
		    continue;
		drvui->voidmap[i][j][k] = 2;
	    }
	}
    }
    int gridvol = drvui->voidgrid[0] * drvui->voidgrid[1] * drvui->voidgrid[2];

    if (drvui->voiddata1)
	free (drvui->voiddata1);
    if (drvui->voiddata2)
	free (drvui->voiddata2);
    drvui->voiddata1 = (char *) malloc (255 * sizeof (char));
    drvui->voiddata2 = (char *) malloc (255 * sizeof (char));
    sprintf (drvui->voiddata1, "Void voxels %d (of %d), void percentage %5.2f\n", l,
	     gridvol, 100. * (float) l / (float) gridvol);
    float vol =
	(float) (drvui->b_mat[0][0] *
		 (drvui->b_mat[1][1] * drvui->b_mat[2][2] -
		  drvui->b_mat[1][2] * drvui->b_mat[2][1])
		 - drvui->b_mat[0][1] * (drvui->b_mat[1][0] * drvui->b_mat[2][2] -
					 drvui->b_mat[1][2] * drvui->b_mat[2][0])
		 + drvui->b_mat[0][2] * (drvui->b_mat[1][0] * drvui->b_mat[2][1] -
					 drvui->b_mat[1][1] * drvui->b_mat[2][0]));
    sprintf (drvui->voiddata2,
	     "Unit cell volume %5.2f A^3, voxel volume %5.2f, overall void volume %5.2f A^3\n",
	     vol, vol / (float) (gridvol), (float) l * vol / (float) (gridvol));
    free (radi);
    drvui->voidflag = -1;
    Progress_Window (-2, NULL, 130.0f);
    return;
}



void
generate_simplevoids ()
{
    int i, j, k, l1, l2;

    int ii, jj, kk;

    float p[3], vert2[3];

    float stepx = 1.f / (float) drvui->voidgrid[0];

    float stepy = 1.f / (float) drvui->voidgrid[1];

    float stepz = 1.f / (float) drvui->voidgrid[2];

    float glr, glg, glb;

    char vcolor[40];

    float probe_r = drvui->probesize * drvui->Sphere_Mult;

    float step2 = probe_r / 2.f;

    int gridx = drvui->voidgrid[0];

    int gridy = drvui->voidgrid[1];

    int gridz = drvui->voidgrid[2];

    float x0[3], x1[3], n[3];

    strcpy (vcolor, drvui->voidcolor);
    Transform_VRML_Color (vcolor);

    (void) sscanf (vcolor, "%f %f %f", &glr, &glg, &glb);
    glColor3f (glr, glg, glb);

    int minx = (int) (drvui->frames[drvui->frame_no].cryst_lim[0] * gridx);

    int maxx = (int) (drvui->frames[drvui->frame_no].cryst_lim[3] * gridx);

    int miny = (int) (drvui->frames[drvui->frame_no].cryst_lim[1] * gridy);

    int maxy = (int) (drvui->frames[drvui->frame_no].cryst_lim[4] * gridy);

    int minz = (int) (drvui->frames[drvui->frame_no].cryst_lim[2] * gridz);

    int maxz = (int) (drvui->frames[drvui->frame_no].cryst_lim[5] * gridz);

    glBegin (GL_QUADS);
    //glBegin(GL_POINTS);
    for (i = minx; i < maxx; i++) {
	for (j = miny; j < maxy; j++) {
	    for (k = minz; k < maxz; k++) {
		ii = i % gridx;
		if (ii < 0)
		    ii += gridx;
		jj = j % gridy;
		if (jj < 0)
		    jj += gridy;
		kk = k % gridz;
		if (kk < 0)
		    kk += gridz;
		if (drvui->voidmap[ii][jj][kk] == 1) {
		    p[0] = stepx * (float) i;
		    p[1] = stepy * (float) j;
		    p[2] = stepz * (float) k;
		    if (p[0] < drvui->frames[drvui->frame_no].cryst_lim[0]
			|| p[0] > drvui->frames[drvui->frame_no].cryst_lim[3])
			continue;
		    if (p[1] < drvui->frames[drvui->frame_no].cryst_lim[1]
			|| p[1] > drvui->frames[drvui->frame_no].cryst_lim[4])
			continue;
		    if (p[2] < drvui->frames[drvui->frame_no].cryst_lim[2]
			|| p[2] > drvui->frames[drvui->frame_no].cryst_lim[5])
			continue;

		    for (l1 = 0; l1 < 3; ++l1) {	/* convert vertex coordinates to Cartesian */
			vert2[l1] = 0.0f;
			for (l2 = 0; l2 < 3; ++l2)
			    vert2[l1] +=
				(float) drvui->b_mat[l1][l2] * (p[l2] - origin[l2]);
		    }

		    if (ii==0 || (ii > 0 && drvui->voidmap[ii - 1][jj][kk] == 0)) {
			x0[0] = vert2[0] - step2;
			x0[1] = vert2[1] - step2;
			x0[2] = vert2[2] - step2;
			x1[0] = vert2[0] - step2;
			x1[1] = vert2[1] - step2;
			x1[2] = vert2[2] + step2;
			vcross (x0, x1, n);
			vnormalize (n);
			glNormal3fv (n);
			glVertex3f (vert2[0] - step2, vert2[1] - step2, vert2[2] - step2);
			glVertex3f (vert2[0] - step2, vert2[1] - step2, vert2[2] + step2);
			glVertex3f (vert2[0] - step2, vert2[1] + step2, vert2[2] + step2);
			glVertex3f (vert2[0] - step2, vert2[1] + step2, vert2[2] - step2);
			if (doPOV) {
			    fprintf (drvui->fpoutp, "polygon{5,\n");
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				     vert2[0] - step2, vert2[1] - step2,
				     vert2[2] - step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				     vert2[0] - step2, vert2[1] - step2,
				     vert2[2] + step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				     vert2[0] - step2, vert2[1] + step2,
				     vert2[2] + step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				     vert2[0] - step2, vert2[1] + step2,
				     vert2[2] - step2);
			    fprintf (drvui->fpoutp,
				     "<%5.3f,%5.3f,%5.3f>\ntexture{pigment{color %s}}}\n",
				     vert2[0] - step2, vert2[1] - step2, vert2[2] - step2,
				     drvui->voidcolor);
			}
		        if (doAsy) {
			    fprintf (drvui->fpouta, "draw(pic, surface( (%8.5f, %8.5f, %8.5f)--",
				     vert2[0] - step2, vert2[1] - step2,
				     vert2[2] - step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f>--\n",
				     vert2[0] - step2, vert2[1] - step2,
				     vert2[2] + step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f>--",
				     vert2[0] - step2, vert2[1] + step2,
				     vert2[2] + step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f>--cycle),\n",
				     vert2[0] - step2, vert2[1] + step2,
				     vert2[2] - step2);
			    fprintf (drvui->fpouta, "rgb(%.2f, %.2f, %.2f));\n", glr, glg, glb);
			}
			if (doVrml) {
			    if (!Vrml2) {
				fprintf (drvui->fpoutv, " Separator {\n");
				fprintf (drvui->fpoutv,
					 " Material {  diffuseColor %s }\n", vcolor);
				fprintf (drvui->fpoutv, "Coordinate3 { \n point[\n");
			    } else {
				fprintf (drvui->fpoutv, " Shape {\n");
				fprintf (drvui->fpoutv,
					 "  geometry IndexedFaceSet { coord Coordinate{ point [\n");
			    }
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				     vert2[0] - step2, vert2[1] - step2,
				     vert2[2] - step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				     vert2[0] - step2, vert2[1] - step2,
				     vert2[2] + step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				     vert2[0] - step2, vert2[1] + step2,
				     vert2[2] + step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f]\n",
				     vert2[0] - step2, vert2[1] + step2,
				     vert2[2] - step2);
			    if (Vrml2) {
				fprintf (drvui->fpoutv,
					 "  }\n  coordIndex [ 0,1,2,3,-1]\n color Color { color [%s]"
					 "}\n colorIndex [0]\n colorPerVertex FALSE\n }\n }\n",
					 vcolor);
			    } else {
				fprintf (drvui->fpoutv,
					 "}\n IndexedFaceSet { coordIndex [0,1,2,3,-1] }\n}\n");
			    }
			}
		    }
		    if (ii == gridx-1 || (ii < gridx - 1 && drvui->voidmap[ii + 1][jj][kk] == 0)) {
			x0[0] = vert2[0] + step2;
			x0[1] = vert2[1] - step2;
			x0[2] = vert2[2] - step2;
			x1[0] = vert2[0] + step2;
			x1[1] = vert2[1] - step2;
			x1[2] = vert2[2] + step2;
			vcross (x0, x1, n);
			vnormalize (n);
			glNormal3fv (n);
			glVertex3f (vert2[0] + step2, vert2[1] - step2,
				    vert2[2] - step2);
			glVertex3f (vert2[0] + step2, vert2[1] - step2,
				    vert2[2] + step2);
			glVertex3f (vert2[0] + step2, vert2[1] + step2,
				    vert2[2] + step2);
			glVertex3f (vert2[0] + step2, vert2[1] + step2,
				    vert2[2] - step2);
			if (doPOV) {
			    fprintf (drvui->fpoutp, "polygon{5,\n");
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpoutp,
				    "<%5.3f,%5.3f,%5.3f>\ntexture{pigment{color %s}}}\n",
				    vert2[0] + step2, vert2[1] - step2,
					 vert2[2] - step2, drvui->voidcolor);
			}
			if (doAsy) {
			    fprintf (drvui->fpouta, "draw(pic, surface( (%8.5f, %8.5f, %8.5f)--",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpouta, "<%8.5f,%8.5f,%8.5f)--\n",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--cycle),\n",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpouta, "rgb(%.2f, %.2f, %.2f));\n", glr, glg, glb);
			}
			if (doVrml) {
			    if (!Vrml2) {
				fprintf (drvui->fpoutv, " Separator {\n");
				fprintf (drvui->fpoutv,
					" Material {  diffuseColor %s }\n", vcolor);
				fprintf (drvui->fpoutv, "Coordinate3 { \n point[\n");
			    } else {
				fprintf (drvui->fpoutv, " Shape {\n");
				fprintf (drvui->fpoutv,
					"  geometry IndexedFaceSet { coord Coordinate{ point [\n");
			    }
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f]\n",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] - step2);
			    if (Vrml2) {
				fprintf (drvui->fpoutv,
					"  }\n  coordIndex [ 0,1,2,3,-1]\n color Color { color [%s]"
					"}\n colorIndex [0]\n colorPerVertex FALSE\n }\n }\n",
					vcolor);
			    } else {
				fprintf (drvui->fpoutv,
					"}\n IndexedFaceSet { coordIndex [0,1,2,3,-1] }\n}\n");
			    }
			}
		    }
		    if (jj==0 || (jj > 0 && drvui->voidmap[ii][jj - 1][kk] == 0)) {
			x0[0] = vert2[0] - step2;
			x0[1] = vert2[1] - step2;
			x0[2] = vert2[2] - step2;
			x1[0] = vert2[0] - step2;
			x1[1] = vert2[1] - step2;
			x1[2] = vert2[2] + step2;
			vcross (x0, x1, n);
			vnormalize (n);
			glNormal3fv (n);
			glVertex3f (vert2[0] - step2, vert2[1] - step2,
				    vert2[2] - step2);
			glVertex3f (vert2[0] - step2, vert2[1] - step2,
				    vert2[2] + step2);
			glVertex3f (vert2[0] + step2, vert2[1] - step2,
				    vert2[2] + step2);
			glVertex3f (vert2[0] + step2, vert2[1] - step2,
				    vert2[2] - step2);
			if (doPOV) {
			    fprintf (drvui->fpoutp, "polygon{5,\n");
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] - step2, vert2[1] - step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] - step2, vert2[1] - step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpoutp,
				    "<%5.3f,%5.3f,%5.3f>\ntexture{pigment{color %s}}}\n",
				    vert2[0] - step2, vert2[1] - step2,
				    vert2[2] - step2, drvui->voidcolor);
			}
			if (doAsy) {
			    fprintf (drvui->fpouta, "draw(pic, surface( (%8.5f, %8.5f, %8.5f)--",
				    vert2[0] - step2, vert2[1] - step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--\n",
				    vert2[0] - step2, vert2[1] - step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--cycle)\n",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpouta, "rgb(%.2f, %.2f, %.2f));\n", glr, glg, glb);
			}
			if (doVrml) {
			    if (!Vrml2) {
				fprintf (drvui->fpoutv, " Separator {\n");
				fprintf (drvui->fpoutv,
					" Material {  diffuseColor %s }\n", vcolor);
				fprintf (drvui->fpoutv, "Coordinate3 { \n point[\n");
			    } else {
				fprintf (drvui->fpoutv, " Shape {\n");
				fprintf (drvui->fpoutv,
					"  geometry IndexedFaceSet { coord Coordinate{ point [\n");
			    }
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				    vert2[0] - step2, vert2[1] - step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				    vert2[0] - step2, vert2[1] - step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f]\n",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] - step2);
			    if (Vrml2) {
				fprintf (drvui->fpoutv,
					"  }\n  coordIndex [ 0,1,2,3,-1]\n color Color { color [%s]"
					"}\n colorIndex [0]\n colorPerVertex FALSE\n }\n }\n",
					vcolor);
			    } else {
				fprintf (drvui->fpoutv,
					"}\n IndexedFaceSet { coordIndex [0,1,2,3,-1] }\n}\n");
			    }
			}
		    }
		    if (jj == gridy-1 || (jj < gridy - 1 && drvui->voidmap[ii][jj + 1][kk] == 0)) {
			x0[0] = vert2[0] - step2;
			x0[1] = vert2[1] + step2;
			x0[2] = vert2[2] - step2;
			x1[0] = vert2[0] - step2;
			x1[1] = vert2[1] + step2;
			x1[2] = vert2[2] + step2;
			vcross (x0, x1, n);
			vnormalize (n);
			glNormal3fv (n);
			glVertex3f (vert2[0] - step2, vert2[1] + step2,
				    vert2[2] - step2);
			glVertex3f (vert2[0] - step2, vert2[1] + step2,
				    vert2[2] + step2);
			glVertex3f (vert2[0] + step2, vert2[1] + step2,
				    vert2[2] + step2);
			glVertex3f (vert2[0] + step2, vert2[1] + step2,
				    vert2[2] - step2);
			if (doPOV) {
			    fprintf (drvui->fpoutp, "polygon{5,\n");
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] - step2, vert2[1] + step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] - step2, vert2[1] + step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpoutp,
				    "<%5.3f,%5.3f,%5.3f>\ntexture{pigment{color %s}}}\n",
				    vert2[0] - step2, vert2[1] + step2,
				    vert2[2] - step2, drvui->voidcolor);
			}
			if (doAsy) {
			    fprintf (drvui->fpouta, "draw(pic, surface( (%8.5f, %8.5f, %8.5f)--",
				    vert2[0] - step2, vert2[1] + step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--\n",
				    vert2[0] - step2, vert2[1] + step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--cycle),\n",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpouta, "rgb(%.2f, %.2f, %.2f));\n", glr, glg, glb);
			}
			if (doVrml) {
			    if (!Vrml2) {
				fprintf (drvui->fpoutv, " Separator {\n");
				fprintf (drvui->fpoutv,
					" Material {  diffuseColor %s }\n", vcolor);
				fprintf (drvui->fpoutv, "Coordinate3 { \n point[\n");
			    } else {
				fprintf (drvui->fpoutv, " Shape {\n");
				fprintf (drvui->fpoutv,
					"  geometry IndexedFaceSet { coord Coordinate{ point [\n");
			    }
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				    vert2[0] - step2, vert2[1] + step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				    vert2[0] - step2, vert2[1] + step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f]\n",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] - step2);
			    if (Vrml2) {
				fprintf (drvui->fpoutv,
					"  }\n  coordIndex [ 0,1,2,3,-1]\n color Color { color [%s]"
					"}\n colorIndex [0]\n colorPerVertex FALSE\n }\n }\n",
					vcolor);
			    } else {
				fprintf (drvui->fpoutv,
					"}\n IndexedFaceSet { coordIndex [0,1,2,3,-1] }\n}\n");
			    }
			}
		    }
		    if (kk==0 || (kk > 0 &&  drvui->voidmap[ii][jj][kk - 1] == 0)) {
			x0[0] = vert2[0] - step2;
			x0[1] = vert2[1] - step2;
			x0[2] = vert2[2] - step2;
			x1[0] = vert2[0] - step2;
			x1[1] = vert2[1] + step2;
			x1[2] = vert2[2] - step2;
			vcross (x0, x1, n);
			vnormalize (n);
			glNormal3fv (n);
			glVertex3f (vert2[0] - step2, vert2[1] - step2,
				    vert2[2] - step2);
			glVertex3f (vert2[0] - step2, vert2[1] + step2,
				    vert2[2] - step2);
			glVertex3f (vert2[0] + step2, vert2[1] + step2,
				    vert2[2] - step2);
			glVertex3f (vert2[0] + step2, vert2[1] - step2,
				    vert2[2] - step2);
			if (doPOV) {
			    fprintf (drvui->fpoutp, "polygon{5,\n");
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] - step2, vert2[1] - step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] - step2, vert2[1] + step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpoutp,
				    "<%5.3f,%5.3f,%5.3f>\ntexture{pigment{color %s}}}\n",
				    vert2[0] - step2, vert2[1] - step2,
				    vert2[2] - step2, drvui->voidcolor);
			}
			if (doAsy) {
			    fprintf (drvui->fpouta, "draw(pic, surface( (%8.5f, %8.5f, %8.5f)--",
				    vert2[0] - step2, vert2[1] - step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--\n",
				    vert2[0] - step2, vert2[1] + step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--cycle),\n",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpouta, "rgb(%.2f, %.2f, %.2f));\n", glr, glg, glb);
			}
			if (doVrml) {
			    if (!Vrml2) {
				fprintf (drvui->fpoutv, " Separator {\n");
				fprintf (drvui->fpoutv,
					 " Material {  diffuseColor %s }\n", vcolor);
				fprintf (drvui->fpoutv, "Coordinate3 { \n point[\n");
			    } else {
				fprintf (drvui->fpoutv, " Shape {\n");
				fprintf (drvui->fpoutv,
					 "  geometry IndexedFaceSet { coord Coordinate{ point [\n");
			    }
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				    vert2[0] - step2, vert2[1] - step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				    vert2[0] - step2, vert2[1] + step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] - step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f]\n",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] - step2);
			    if (Vrml2) {
				fprintf (drvui->fpoutv,
					"  }\n  coordIndex [ 0,1,2,3,-1]\n color Color { color [%s]"
					"}\n colorIndex [0]\n colorPerVertex FALSE\n }\n }\n",
					vcolor);
			    } else {
				fprintf (drvui->fpoutv,
					"}\n IndexedFaceSet { coordIndex [0,1,2,3,-1] }\n}\n");
			    }
			}
		    }
		    if (kk==gridz-1 || (kk < gridz - 1 && drvui->voidmap[ii][jj][kk + 1] == 0)) {
			x0[0] = vert2[0] - step2;
			x0[1] = vert2[1] - step2;
			x0[2] = vert2[2] + step2;
			x1[0] = vert2[0] - step2;
			x1[1] = vert2[1] + step2;
			x1[2] = vert2[2] + step2;
			vcross (x0, x1, n);
			vnormalize (n);
			glNormal3fv (n);
			glVertex3f (vert2[0] - step2, vert2[1] - step2,
				    vert2[2] + step2);
			glVertex3f (vert2[0] - step2, vert2[1] + step2,
				    vert2[2] + step2);
			glVertex3f (vert2[0] + step2, vert2[1] + step2,
				    vert2[2] + step2);
			glVertex3f (vert2[0] + step2, vert2[1] - step2,
				    vert2[2] + step2);
			if (doPOV) {
			    fprintf (drvui->fpoutp, "polygon{5,\n");
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] - step2, vert2[1] - step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] - step2, vert2[1] + step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutp, "<%5.3f,%5.3f,%5.3f>,\n",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutp,
				    "<%5.3f,%5.3f,%5.3f>\ntexture{pigment{color %s}}}\n",
				    vert2[0] - step2, vert2[1] - step2,
				    vert2[2] + step2, drvui->voidcolor);
			}
			if (doAsy) {
			    fprintf (drvui->fpouta, "draw(pic, surface( (%8.5f, %8.5f, %8.5f)--",
				    vert2[0] - step2, vert2[1] - step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--\n",
				    vert2[0] - step2, vert2[1] + step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpouta, "(%8.5f,%8.5f,%8.5f)--cycle),\n",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpouta, "rgb(%.2f, %.2f, %.2f));\n", glr, glg, glb);
			}
			if (doVrml) {
			    if (!Vrml2) {
				fprintf (drvui->fpoutv, " Separator {\n");
				fprintf (drvui->fpoutv,
					" Material {  diffuseColor %s }\n", vcolor);
				fprintf (drvui->fpoutv, "Coordinate3 { \n point[\n");
			    } else {
				fprintf (drvui->fpoutv, " Shape {\n");
				fprintf (drvui->fpoutv,
					"  geometry IndexedFaceSet { coord Coordinate{ point [\n");
			    }
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				    vert2[0] - step2, vert2[1] - step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				    vert2[0] - step2, vert2[1] + step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n",
				    vert2[0] + step2, vert2[1] + step2,
				    vert2[2] + step2);
			    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f]\n",
				    vert2[0] + step2, vert2[1] - step2,
				    vert2[2] + step2);
			    if (Vrml2) {
				fprintf (drvui->fpoutv,
					"  }\n  coordIndex [ 0,1,2,3,-1]\n color Color { color [%s]"
					"}\n colorIndex [0]\n colorPerVertex FALSE\n }\n }\n",
					vcolor);
			    } else {
				fprintf (drvui->fpoutv,
					"}\n IndexedFaceSet { coordIndex [0,1,2,3,-1] }\n}\n");
			    }
			}
		    }
		}
	    }
	}
    }
    glEnd ();
    fprintf (drvui->flout, "%s", drvui->voiddata1);
    fprintf (drvui->flout, "%s", drvui->voiddata2);
}

void
calculate_sas (void)
{
    int i, j, k;

    int l1, l2;

    int onvert;

    int ncount;

    int nsample = drvui->voidgrid[0];

    int deny;

    float rprobe = drvui->probesize;

    float radius;

    int iseed = drvui->voidgrid[1];

    float phi, theta, costheta;

    float dist2;

    float dx, dy, dz;

    float xpoint, ypoint, zpoint;

    float vert2[3];

    float saved_boxlim[3];

    float saved_crystlim[6];

    float vol, sfrac, sjreal, stotal;

    float *radi;

    float rad;

    float *fp;

    int *ip;

    char tmpfile[256];

    if (nsample <= 0)
	nsample = 5000;
    if (iseed == 0)
	iseed = -5213150;
    strcpy (tmpfile, drvui->Cur_Root);
    strcat (tmpfile, ".coord");
    FILE *tdf = fopen (tmpfile, "w");

    Progress_Window (-1, "Calculating SAS...", 100.0f);
    Fl::flush ();

    for (l1 = 0; l1 < 3; ++l1) {	/* convert vertex coordinates to Cartesian */
	vert2[l1] = 0.0f;
	for (l2 = 0; l2 < 3; ++l2)
	    vert2[l1] += (float) (drvui->b_mat[l1][l2] * (0. - origin[l2]));
    }
    float cxmin = vert2[0];

    float cymin = vert2[1];

    float czmin = vert2[2];

    for (l1 = 0; l1 < 3; ++l1) {	/* convert vertex coordinates to Cartesian */
	vert2[l1] = 0.0f;
	for (l2 = 0; l2 < 3; ++l2)
	    vert2[l1] += (float) (drvui->b_mat[l1][l2] * (1. - origin[l2]));
    }
    float cxmax = vert2[0];

    float cymax = vert2[1];

    float czmax = vert2[2];

    cxmin += cxmax;
    cxmax *= 2.;
    cymin += cymax;
    cymax *= 2.;
    czmin += czmax;
    czmax *= 2.;
    ip = (int *) realloc (drvui->orig_atom_no, (4 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	fclose (tdf);
	return;
    } else {
	drvui->orig_atom_no = ip;
    }
    ip = (int *) realloc (vert_sym_no, (4 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	fclose (tdf);
	return;
    } else {
	vert_sym_no = ip;
    }
    ip = (int *) realloc (vert_sym_nos, (4 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	fclose (tdf);
	return;
    } else {
	vert_sym_nos = ip;
    }
    fp = (float *) realloc (o_vert, (12 * NvertM) * sizeof (float));
    if (!fp) {
	Error_Box("Unable to allocate memory for volume test");
	fclose (tdf);
	return;
    } else {
	o_vert = fp;
    }
    fp = (float *) realloc (o_vert_nm, (12 * NvertM) * sizeof (float));
    if (!fp) {
	Error_Box("Unable to allocate memory for volume test");
	fclose (tdf);
	return;
    } else {
	o_vert_nm = fp;
    }
    nvert = 0;			/* clear the vertex list */
    for (j = 0; j < 3; j++) {
	saved_boxlim[j] = boxlim[j];
	saved_crystlim[j] = drvui->frames[drvui->frame_no].cryst_lim[j];
	saved_crystlim[j + 3] = drvui->frames[drvui->frame_no].cryst_lim[j + 3];
	drvui->frames[drvui->frame_no].cryst_lim[j] = 0.0;
	drvui->frames[drvui->frame_no].cryst_lim[j + 3] = 1.0;
    }
    NvertM = 1;
    nvert = 0;
    onvert = 0;
    build_box_contents ();
    ip = (int *) realloc (drvui->orig_atom_no, (2 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	fclose (tdf);
	return;
    } else {
	drvui->orig_atom_no = ip;
    }
    ip = (int *) realloc (vert_sym_no, (2 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	fclose (tdf);
	return;
    } else {
	vert_sym_no = ip;
    }
    ip = (int *) realloc (vert_sym_nos, (2 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	fclose (tdf);
	return;
    } else {
	vert_sym_nos = ip;
    }
    fp = (float *) realloc (o_vert, 6 * (NvertM + 20) * sizeof (float));
    if (!fp) {
	Error_Box("Unable to allocate memory for volume test");
	fclose (tdf);
	return;
    } else {
	o_vert = fp;
    }
    fp = (float *) realloc (o_vert_nm, 6 * (NvertM + 20) * sizeof (float));
    if (!fp) {
	Error_Box("Unable to allocate memory for volume test");
	fclose (tdf);
	return;
    } else {
	o_vert_nm = fp;
    }
    radi = (float *) malloc (NvertM * sizeof (float));
    if (!radi) {
	Error_Box("Unable to allocate memory for volume test");
	fclose (tdf);
	return;
    }

    fprintf (tdf, "%10d\n", NvertM - 1);

    for (j = 0; j < natom; ++j) {	/* loop through atoms */
	onvert = nvert;
	if (drvui->atoms[j].atom_fn != drvui->frame_no)
	    continue;
	radius = drvui->atoms[j].radius;
	find_all_in_box (j);
	for (k = onvert; k < nvert; k++) {
	    radi[k] = radius;
	    fprintf (tdf, "%4d %11.5f %11.5f %11.5f     %s       0.00   0   0\n",
		     k + 1, s_vert[3 * k], s_vert[3 * k + 1], s_vert[3 * k + 2],
		     drvui->atoms[j].atom_l);
	}
    }
    fclose (tdf);

    Progress_Window (0, NULL, 10.0f);
//    glPointSize(3);
//    glColor3f (0, 200, 0);
//    glBegin(GL_POINTS);

    sfrac = 0;
    stotal = 0.;
    srand (iseed);
    int ntotal = 0;

    for (i = 0; i < nvert; i++) {
	rad = rprobe + radi[i];
	ncount = 0;
	for (j = 0; j < nsample; j++) {	//samples
	    phi = (float) (rand () / RAND_MAX * 2.0 * M_PI);
	    costheta = (float) (1. - rand () / RAND_MAX * 2.0);
	    theta = acosf (costheta);
	    xpoint = (float) (sin (theta) * cos (phi));
	    ypoint = (float) (sin (theta) * sin (phi));
	    zpoint = costheta;
	    xpoint *= rad;
	    ypoint *= rad;
	    zpoint *= rad;

	    xpoint += (s_vert[3 * i] + cxmax / 2.f);
	    ypoint += (s_vert[3 * i + 1] + cymax / 2.f);
	    zpoint += (s_vert[3 * i + 2] + czmax / 2.f);

// apply periodic boundary conditions
	    if (xpoint < 0.)
		xpoint += cxmax;
	    if (xpoint >= cxmax)
		xpoint -= cxmax;
	    if (ypoint < 0.)
		ypoint += cymax;
	    if (ypoint >= cymax)
		ypoint -= cymax;
	    if (zpoint < 0.)
		zpoint += czmax;
	    if (zpoint >= czmax)
		zpoint -= czmax;

	    deny = 0;
	    for (k = 0; k < nvert; k++) {
		if (i == k)
		    continue;

		dx = xpoint - (s_vert[3 * k] + cxmax / 2.f);
		dx -= cxmax * rint (2. * dx / cxmax);

		dy = ypoint - (s_vert[3 * k + 1] + cymax / 2.f);
		dy -= cymax * rint (2. * dy / cymax);

		dz = zpoint - (s_vert[3 * k + 2] + czmax / 2.f);
		dz -= czmax * rint (2. * dz / czmax);

		dist2 = dx * dx + dy * dy + dz * dz;
		if (sqrt (dist2) < 0.999 * (radi[k] + rprobe)) {
		    deny = 1;
		    break;
		}
	    }

	    if (deny == 0) {
//          glVertex3f(xpoint-cxmax/2.f,ypoint-cymax/2.f,zpoint-czmax/2.f);
		ntotal++;
		if (doPOV) {
		    fprintf (drvui->fpoutp, " sphere{<%8.5f, %8.5f, %8.5f>, %8.5f\n",
			     xpoint - cxmax / 2., ypoint - cymax / 2.,
			     zpoint - czmax / 2., rprobe);
		    fprintf (drvui->fpoutp, "  texture{pigment{color %s  }} }\n",
			     drvui->voidcolor);
		}
		ncount++;
	    }
	}			// for j samples on this atom

	sfrac = (float) ncount / (float) nsample;
	sjreal = (float) (4.0 * M_PI * rad * rad * sfrac);
	stotal += sjreal;
	Progress_Window (0, NULL, (float) i / (float) nvert * 100.0f);
    }				// for i atoms in cell
    vol =
	(float) (drvui->b_mat[0][0] *
		 (drvui->b_mat[1][1] * drvui->b_mat[2][2] -
		  drvui->b_mat[1][2] * drvui->b_mat[2][1])
		 - drvui->b_mat[0][1] * (drvui->b_mat[1][0] * drvui->b_mat[2][2] -
					 drvui->b_mat[1][2] * drvui->b_mat[2][0])
		 + drvui->b_mat[0][2] * (drvui->b_mat[1][0] * drvui->b_mat[2][1] -
					 drvui->b_mat[1][1] * drvui->b_mat[2][0]));

    float stotalreduced = stotal / vol * 1.E4f;

    if (drvui->voiddata1)
	free (drvui->voiddata1);
    if (drvui->voiddata2)
	free (drvui->voiddata2);
    drvui->voiddata1 = (char *) malloc (255 * sizeof (char));
    drvui->voiddata2 = (char *) malloc (255 * sizeof (char));
    sprintf (drvui->voiddata1, " accessible surface area (A^2): %5.2f\n", stotal);
    sprintf (drvui->voiddata2, " accessible surface area per volume: %5.2f\n",
	     stotalreduced);
    fprintf (drvui->flout, "%s", drvui->voiddata1);
    fprintf (drvui->flout, "%s", drvui->voiddata2);
//    glEnd();
    free (radi);
    Progress_Window (-2, NULL, 100.0f);
}


//***************
void
calculate_msms (void)
{
    int j, k;

    int onvert;

    float radius;

    float saved_boxlim[3];

    float saved_crystlim[6];

    char msmsfile[256];

    char msmslog[256];

    char cmd[512];

    char pr[10];

    float sesvol, sesa;

    float *fp;

    int *ip;

    FILE *tdf;

    int saved_nvert = NvertM;

    strcpy (msmsfile, drvui->Cur_Root);
    strcat (msmsfile, ".xyzrn");
    
    strcpy (msmslog, drvui->Cur_Root);
    strcat (msmslog, ".msmslog");

    Progress_Window (-1, "Running MSMS...", 100.0f);
    Fl::flush ();
    ip = (int *) realloc (drvui->orig_atom_no, (4 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	drvui->orig_atom_no = ip;
    }
    ip = (int *) realloc (vert_sym_no, (4 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	vert_sym_no = ip;
    }
    ip = (int *) realloc (vert_sym_nos, (4 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	vert_sym_nos = ip;
    }
    fp = (float *) realloc (o_vert, (12 * NvertM) * sizeof (float));
    if (!fp) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	o_vert = fp;
    }
    fp = (float *) realloc (o_vert_nm, (12 * NvertM) * sizeof (float));
    if (!fp) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	o_vert_nm = fp;
    }
    nvert = 0;			/* clear the vertex list */
    Progress_Window (0, NULL, 10.0f);

    for (j = 0; j < 3; j++) {
	saved_boxlim[j] = boxlim[j];
	saved_crystlim[j] = drvui->frames[drvui->frame_no].cryst_lim[j];
	saved_crystlim[j + 3] = drvui->frames[drvui->frame_no].cryst_lim[j + 3];
	drvui->frames[drvui->frame_no].cryst_lim[j] = 0.0f;
	drvui->frames[drvui->frame_no].cryst_lim[j + 3] = 1.0f;
    }
    NvertM = 1;
    build_box_contents ();
    ip = (int *) realloc (drvui->orig_atom_no, (2 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	drvui->orig_atom_no = ip;
    }
    ip = (int *) realloc (vert_sym_no, (2 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	vert_sym_no = ip;
    }
    ip = (int *) realloc (vert_sym_nos, (2 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	vert_sym_nos = ip;
    }
    fp = (float *) realloc (o_vert, 6 * (NvertM + 20) * sizeof (float));
    if (!fp) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	o_vert = fp;
    }
    fp = (float *) realloc (o_vert_nm, 6 * (NvertM + 20) * sizeof (float));
    if (!fp) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	o_vert_nm = fp;
    }
    onvert = 0;

    tdf = fopen (msmsfile, "w");

    for (j = 0; j < natom; ++j) {	/* loop through atoms */
	onvert = nvert;
	if (drvui->atoms[j].atom_fn != drvui->frame_no)
	    continue;
	radius = drvui->atoms[j].radius;
	find_all_in_box (j);
	for (k = onvert; k < nvert; k++) {
	    fprintf (tdf, "%11.5f %11.5f %11.5f  %5.3f  1  %s\n",
		     s_vert[3 * k], s_vert[3 * k + 1], s_vert[3 * k + 2], radius,
		     drvui->atoms[j].atom_l);
	}
    }
    fclose (tdf);

    Progress_Window (0, NULL, 20.0f);

    nvert = 0;
    for (j = 0; j < 3; j++) {
	boxlim[j]=saved_boxlim[j];
	drvui->frames[drvui->frame_no].cryst_lim[j] = saved_crystlim[j];
	drvui->frames[drvui->frame_no].cryst_lim[j + 3] = saved_crystlim[j + 3];
    }
    ip = (int *) realloc (drvui->orig_atom_no, (2 * saved_nvert) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	drvui->orig_atom_no = ip;
    }
    ip = (int *) realloc (vert_sym_no, (2 * saved_nvert) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	vert_sym_no = ip;
    }
    ip = (int *) realloc (vert_sym_nos, (2 * saved_nvert) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	vert_sym_nos = ip;
    }
    fp = (float *)realloc (o_vert, (6 * (saved_nvert+20) * sizeof (float)));
    if (!fp) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	o_vert = fp;
    }
    fp = (float *)realloc (o_vert_nm, (6 * (saved_nvert+20) * sizeof (float)));
    if (!fp) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	o_vert_nm = fp;
    }

    NvertM = 1;
    build_box_contents ();
    ip = (int *) realloc (drvui->orig_atom_no, (2 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	drvui->orig_atom_no = ip;
    }
    ip = (int *) realloc (vert_sym_no, (2 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	vert_sym_no = ip;
    }
    ip = (int *) realloc (vert_sym_nos, (2 * NvertM) * sizeof (int));
    if (!ip) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	vert_sym_nos = ip;
    }
    fp = (float *)realloc (o_vert, (6 * (NvertM+20) * sizeof (float)));
    if (!fp) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	o_vert = fp;
    }
    fp = (float *)realloc (o_vert_nm, (6 * (NvertM+20) * sizeof (float)));
    if (!fp) {
	Error_Box("Unable to allocate memory for volume test");
	return;
    } else {
	o_vert_nm = fp;
    }

    sprintf (pr, "%5.3f", drvui->probesize);
#ifdef WIN32
    strcpy (cmd, "\"\"");	// build the command string
    strcat (cmd, drvui->MSMS_Path);	//   for Windows
    strcat (cmd, "\"");
    strcat (cmd, " ");
#else
    strcpy (cmd, drvui->MSMS_Path);	// build command string
    strcat (cmd, " ");		//    for Linux
#endif
    strcat (cmd, "-probe_radius ");
    strcat (cmd, pr);
    strcat (cmd, " -if ");
    strcat (cmd, msmsfile);
    strcat (cmd, " -of ");
    strcat (cmd, drvui->Cur_Root);
    strcat (cmd, " >");
    strcat (cmd, msmslog);
#ifdef WIN32
    strcat (cmd, "\"\"");
#endif
    if (system (cmd) != 0) {	// call the MSMS program
	Error_Box ("An error occurred running the MSMS program."
		   "\nThe most probable cause is that the path is not correct."
		   "\nCheck settings.");
	Progress_Window (-2, NULL, 100.0f);
    } else {
	Progress_Window (0, NULL, 90.0f);
	tdf = fopen (msmslog, "r");
	while (!feof (tdf)) {
	    if (fgets (cmd, 200, tdf) == NULL)
		break;
	    if (!strncmp (cmd, "NUMERICAL VOL", 13)) {
		fgets (cmd, 200, tdf);
		fgets (cmd, 200, tdf);
		sscanf (cmd, "%*d %*f %f %f", &sesvol, &sesa);
		break;
	    }
	}
	if (drvui->voiddata1)
	    free (drvui->voiddata1);
	drvui->voiddata1 = (char *) malloc (255 * sizeof (char));
	sprintf (drvui->voiddata1, "Calculated SES volume %5.2f, SES area %5.2f\n",
		 sesvol, sesa);
	Progress_Window (-2, NULL, 100.0f);
	drvui->voidflag = -2;
    }

}

void
generate_msms (void)
{
    char msmsvert[256];

    char msmsface[256];

    char line[100];

    int numverts;

    int numfaces;

    float *tvert;

    float *tnorm;

    int *tname;

    int i, j, k;

    int i1, i2, i3;

    FILE *msmsfile;

    float glr, glg, glb;

    char vcolor[40];

    char atname[5];

    int byatom;

    byatom = 0;
    strcpy (vcolor, drvui->voidcolor);
    if (strstr (vcolor,"byatom")) {
	byatom = 1;
    } else {
	Transform_VRML_Color (vcolor);
        (void) sscanf (vcolor, "%f %f %f", &glr, &glg, &glb);
	glColor3f (glr, glg, glb);
    }

    strcpy (msmsvert, drvui->Cur_Root);
    strcat (msmsvert, ".vert");
    strcpy (msmsface, drvui->Cur_Root);
    strcat (msmsface, ".face");

    msmsfile = fopen (msmsvert, "r");
    if (!msmsfile) {
	Error_Box ("MSMS output (vertex file) not found");
	return;
    }
    (void) fgets (line, 100, msmsfile);
    (void) fgets (line, 100, msmsfile);
    (void) fgets (line, 100, msmsfile);
    sscanf (line, "%d", &numverts);
    numverts += 1;
    tvert = (float *) malloc (3 * numverts * sizeof (float));
    if (!tvert) {
	Error_Box("Unable to allocate memory for volume test");
	fclose (msmsfile);
	return;
    }
    tnorm = (float *) malloc (3 * numverts * sizeof (float));
    if (!tnorm) {
	Error_Box("Unable to allocate memory for volume test");
        free (tvert);
	fclose (msmsfile);
	return;
    }
    tname = (int *) malloc (numverts * sizeof(int));
    if (!tname) {
	Error_Box("Unable to allocate memory for volume test");
        free (tvert);
	free (tnorm);
	fclose (msmsfile);
	return;
    }
    j = 3;
    for (i = 0; i < numverts; i++) {
	if (!fgets (line, 100, msmsfile)) {
	    numverts=i;
	    break;
	}
	sscanf (line, "%f %f %f %f %f %f %*d %*d %*d %s\n", &tvert[j], &tvert[j + 1], &tvert[j + 2],
		&tnorm[j], &tnorm[j + 1], &tnorm[j + 2],atname);

	while(strlen(atname)<4)strcat(atname," ");
	for (k=1;k<drvui->nsphere;k++) {
	    if (check_atom_name (atname, drvui->spheres[k].sphere_l)) {
		tname[i+1]=k;
		break;
	    }
	}
	j += 3;
    }
    fclose (msmsfile);

    msmsfile = fopen (msmsface, "r");
    if (!msmsfile) {
	Error_Box ("MSMS output (face file) not found");
	free (tnorm);
	free (tvert);
	free (tname);
	return;
    }
    (void) fgets (line, 100, msmsfile);
    (void) fgets (line, 100, msmsfile);
    (void) fgets (line, 100, msmsfile);
    sscanf (line, "%d", &numfaces);
//    glShadeModel(GL_SMOOTH);
    glBegin (GL_TRIANGLES);

    if (byatom == 0) {
	for (i = 0; i < numfaces; i++) {
	    fgets (line, 100, msmsfile);
	    sscanf (line, "%d %d %d\n", &i1, &i2, &i3);
 	    glNormal3f (tnorm[3 * i1], tnorm[3 * i1 + 1], tnorm[3 * i1 + 2]);
	    glVertex3f (tvert[3 * i1], tvert[3 * i1 + 1], tvert[3 * i1 + 2]);
	    glNormal3f (tnorm[3 * i2], tnorm[3 * i2 + 1], tnorm[3 * i2 + 2]);
	    glVertex3f (tvert[3 * i2], tvert[3 * i2 + 1], tvert[3 * i2 + 2]);
	    glNormal3f (tnorm[3 * i3], tnorm[3 * i3 + 1], tnorm[3 * i3 + 2]);
	    glVertex3f (tvert[3 * i3], tvert[3 * i3 + 1], tvert[3 * i3 + 2]);

	    if (doVrml) {
		if (!Vrml2) {
		    fprintf (drvui->fpoutv, " Separator {\n");
		    fprintf (drvui->fpoutv, " Material {  diffuseColor %s }\n", vcolor);
		    fprintf (drvui->fpoutv, "Coordinate3 { point [\n");
		} else {
		    fprintf (drvui->fpoutv, " Shape {\n");
		    fprintf (drvui->fpoutv, "  appearance Appearance {\n");
		    fprintf (drvui->fpoutv, "   material Material {diffuseColor %s}\n",
			     vcolor);
		    fprintf (drvui->fpoutv,
		    	     "  }\n  geometry IndexedFaceSet { coord Coordinate{ point [\n");
		}
		fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n", tvert[3 * i1],
			tvert[3 * i1 + 1], tvert[3 * i1 + 2]);
		fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n", tvert[3 * i2],
			tvert[3 * i2 + 1], tvert[3 * i2 + 2]);
		fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f]\n", tvert[3 * i3],
			tvert[3 * i3 + 1], tvert[3 * i3 + 2]);
		if (Vrml2) {
		    fprintf (drvui->fpoutv,
			 "  }\n  coordIndex [ 0,1,2,-1]\n solid FALSE\n convex"
			 " TRUE\n creaseAngle 1.57075\n }\n }\n");
		} else {
		    fprintf (drvui->fpoutv,
			 "}\n IndexedFaceSet { coordIndex [0,1,2,-1] }\n}\n");
		}
	    }
	    if (doPOV) {
		fprintf (drvui->fpoutp, "smooth_triangle {\n");
		fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f> , ", tvert[3 * i1],
			tvert[3 * i1 + 1], tvert[3 * i1 + 2]);
		fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f> , ", tnorm[3 * i1],
			tnorm[3 * i1 + 1], tnorm[3 * i1 + 2]);
		fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f> , ", tvert[3 * i2],
			tvert[3 * i2 + 1], tvert[3 * i2 + 2]);
		fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f> , ", tnorm[3 * i2],
			tnorm[3 * i2 + 1], tnorm[3 * i2 + 2]);
		fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f> , ", tvert[3 * i3],
			tvert[3 * i3 + 1], tvert[3 * i3 + 2]);
		fprintf (drvui->fpoutp,
		 	"<%8.5f, %8.5f, %8.5f>\n texture{pigment{color %s}}\n }\n",
		 	tnorm[3 * i3], tnorm[3 * i3 + 1], tnorm[3 * i3 + 2],
			drvui->voidcolor);
	    }
	    if (doAsy) {
		fprintf (drvui->fpouta, "draw(pic, surface( (%8.5f, %8.5f, %8.5f)--", tvert[3 * i1],
			tvert[3 * i1 + 1], tvert[3 * i1 + 2]);
		fprintf (drvui->fpouta, "(%8.5f, %8.5f, %8.5f)--\n", tvert[3 * i2],
			tvert[3 * i2 + 1], tvert[3 * i2 + 2]);
		fprintf (drvui->fpouta, "(%8.5f, %8.5f, %8.5f)--cycle),", tvert[3 * i3],
			tvert[3 * i3 + 1], tvert[3 * i3 + 2]);
		fprintf (drvui->fpouta, "rgb(%.2f, %.2f, %.2f));\n", glr, glg, glb);
	    }
	}
    } else {
	float r[3], g[3], b[3];
	if (doPOV) {
	    for (k=0;k<drvui->nsphere;k++)
		fprintf(drvui->fpoutp,"#declare T%d=texture{ pigment {color %s}}\n",
			k,drvui->spheres[k].sphere_col);
	    fprintf (drvui->fpoutp, "mesh {\n");
	}
	for (i = 0; i < numfaces; i++) {
	    fgets (line, 100, msmsfile);
	    sscanf (line, "%d %d %d\n", &i1, &i2, &i3);
            k=tname[i1];
            strcpy (vcolor, drvui->spheres[k].sphere_col);
            Transform_VRML_Color(vcolor);
            (void) sscanf (vcolor, "%f %f %f", &r[0], &g[0], &b[0]);
            glColor3f (r[0], g[0], b[0]);
	    glNormal3f (tnorm[3 * i1], tnorm[3 * i1 + 1], tnorm[3 * i1 + 2]);
	    glVertex3f (tvert[3 * i1], tvert[3 * i1 + 1], tvert[3 * i1 + 2]);
            k=tname[i2];
            strcpy (vcolor, drvui->spheres[k].sphere_col);
            Transform_VRML_Color(vcolor);
            (void) sscanf (vcolor, "%f %f %f", &r[1], &g[1], &g[1]);
            glColor3f (r[1], g[1], b[1]);
	    glNormal3f (tnorm[3 * i2], tnorm[3 * i2 + 1], tnorm[3 * i2 + 2]);
	    glVertex3f (tvert[3 * i2], tvert[3 * i2 + 1], tvert[3 * i2 + 2]);
            k=tname[i3];
            strcpy (vcolor, drvui->spheres[k].sphere_col);
            Transform_VRML_Color(vcolor);
            (void) sscanf (vcolor, "%f %f %f", &r[2], &g[2], &b[2]);
            glColor3f (r[2], g[2], b[2]);
	    glNormal3f (tnorm[3 * i3], tnorm[3 * i3 + 1], tnorm[3 * i3 + 2]);
	    glVertex3f (tvert[3 * i3], tvert[3 * i3 + 1], tvert[3 * i3 + 2]);

	    if (doVrml) {
		if (!Vrml2) {
		    fprintf (drvui->fpoutv, " Separator {\n");
		    fprintf (drvui->fpoutv, " Material {  diffuseColor %s }\n", vcolor);
		    fprintf (drvui->fpoutv, "Coordinate3 { point [\n");
		} else {
		    fprintf (drvui->fpoutv, " Shape {\n");
		    fprintf (drvui->fpoutv, "  appearance Appearance {\n");
		    fprintf (drvui->fpoutv, "   material Material {diffuseColor %s}\n",
		  	     vcolor);
		    fprintf (drvui->fpoutv,
		  	     "  }\n  geometry IndexedFaceSet { coord Coordinate{ point [\n");
		}
		fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n", tvert[3 * i1],
			tvert[3 * i1 + 1], tvert[3 * i1 + 2]);
		fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n", tvert[3 * i2],
			tvert[3 * i2 + 1], tvert[3 * i2 + 2]);
		fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f]\n", tvert[3 * i3],
			tvert[3 * i3 + 1], tvert[3 * i3 + 2]);
		if (Vrml2) {
		    fprintf (drvui->fpoutv,
			    "  }\n  coordIndex [ 0,1,2,-1]\n solid FALSE\n convex"
			    " TRUE\n creaseAngle 1.57075\n }\n }\n");
		} else {
		    fprintf (drvui->fpoutv,
			 "}\n IndexedFaceSet { coordIndex [0,1,2,-1] }\n}\n");
		}
	    }
	    if (doPOV) {
		fprintf (drvui->fpoutp, "smooth_triangle {\n");
		fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f> , ", tvert[3 * i1],
			tvert[3 * i1 + 1], tvert[3 * i1 + 2]);
		fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f> , ", tnorm[3 * i1],
			tnorm[3 * i1 + 1], tnorm[3 * i1 + 2]);
		fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f> , ", tvert[3 * i2],
			tvert[3 * i2 + 1], tvert[3 * i2 + 2]);
		fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f> , ", tnorm[3 * i2],
			tnorm[3 * i2 + 1], tnorm[3 * i2 + 2]);
		fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f> , ", tvert[3 * i3],
			tvert[3 * i3 + 1], tvert[3 * i3 + 2]);
		fprintf (drvui->fpoutp,
			"<%8.5f, %8.5f, %8.5f>\n texture_list {T%d T%d T%d}\n }\n",
			tnorm[3 * i3],tnorm[3 * i3 + 1], tnorm[3 * i3 + 2],
			tname[i1],tname[i2],tname[i3]);
	    }
	    if (doAsy) {
		fprintf (drvui->fpouta, "draw(pic, surface( (%8.5f, %8.5f, %8.5f)--", tvert[3 * i1],
			tvert[3 * i1 + 1], tvert[3 * i1 + 2]);
		fprintf (drvui->fpouta, "(%8.5f, %8.5f, %8.5f)--", tvert[3 * i2],
			tvert[3 * i2 + 1], tvert[3 * i2 + 2]);
		fprintf (drvui->fpouta, "(%8.5f, %8.5f, %8.5f)--cycle),\n", tvert[3 * i3],
			tvert[3 * i3 + 1], tvert[3 * i3 + 2]);
		fprintf (drvui->fpouta, "new pen[] {rgb(%.2f, %.2f, %.2f),rgb(%.2f,%.2f,%.2f),",
			r[0], g[0], b[0], r[1], g[1], b[1]);
		fprintf (drvui->fpouta, "rgb(%.2f,%.2f,%.2f)});\n",
			r[2], g[2], b[2]);
	    }
	}
	if (doPOV) fprintf(drvui->fpoutp," }\n"); //close the mesh object
    }
    fclose (msmsfile);
    glEnd ();
//    glShadeModel(GL_FLAT);

    free (tvert);
    free (tnorm);
    free (tname);
    return;
}
