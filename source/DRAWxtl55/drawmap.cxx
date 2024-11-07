// $Id: drawmap.cxx 1064 2010-10-22 19:38:20Z martin $
//
// drawmap.cpp - routine for DRAWxtl V5.5 - the GUI version 
//
// Coded using the FLTK 1.1.6 widget set
//
//     Larry W. Finger, Martin Kroeker and Brian Toby
//
// This module includes Fourier support routines
//
// routines contained within this file:
//
// generate_map - Generate Marching Cubes triangles from a map
// InterpolateMap - interpolate map value
// LookupMap - return the value of the mmap at one of the calculated points
// Maximize_rho - position the cursor at the local maximum (Q & D routine)

#include "drawxtl.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "draw_ext.h"
#include "drawmap.h"
#include "DRAWxtlViewUI.h"
#include "MC.h"

#include "DRAWxtl_proto.h"

/* *************************************************************************************** */
/* Generate Marching Cubes triangles from a map with contour minValue*/
/* *************************************************************************************** */

void
generate_map (float minValue, int Solid, char *Color, char *BackColor)
{
    int numOfTriangles;

    float step[3];

    float Red, Green, Blue;

    int nxMin, nxMax, nyMin, nyMax, nzMin, nzMax;	/* max and min steps along a, b & c */

    int snx, sny, snz;		/* number of steps along a, b & c */

    unsigned int sny_snz;	/* sny x snz */

    int i, j, k, ijk, ijkmax, ii, jj;

    unsigned int ni, nj, nk, si, sj, sk, sni, snj, sind;

    mp4Vector *PointsList;

    float fvert[3], cvert[3];

    TRIANGLE *pTriangles;

    char string[40];

    char Color2[40];

    float d1, d2, d3;

    int iproj = 0;

    int iMax;

    int jMax;

    int kMax;

    if (doVrml) {
	if (!drvui->fpoutv) {
	    Error_Box ("Invalid vrml file path in drawmap.");
	    return;
	}
    }
    if (!FourierPt)
	return;

    drvui->mainWindow->cursor (FL_CURSOR_WAIT);
    strcpy (string, Color);
    strcpy (Color2, Color);
    Transform_VRML_Color (string);
    Transform_POV_Color (Color2);
    (void) sscanf (string, "%f %f %f", &Red, &Green, &Blue);

    step[0] = 1.0f / (float) mapstep_a;
    step[1] = 1.0f / (float) mapstep_b;
    step[2] = 1.0f / (float) mapstep_c;
    nxMin = (int) (max (drvui->frames[drvui->frame_no].map_lim[0], drvui->frames[drvui->frame_no].cryst_lim[0]) * mapstep_a);
    nxMax = (int) (min (drvui->frames[drvui->frame_no].map_lim[3], drvui->frames[drvui->frame_no].cryst_lim[3]) * mapstep_a);
    nyMin = (int) (max (drvui->frames[drvui->frame_no].map_lim[1], drvui->frames[drvui->frame_no].cryst_lim[1]) * mapstep_b);
    nyMax = (int) (min (drvui->frames[drvui->frame_no].map_lim[4], drvui->frames[drvui->frame_no].cryst_lim[4]) * mapstep_b);
    nzMin = (int) (max (drvui->frames[drvui->frame_no].map_lim[2], drvui->frames[drvui->frame_no].cryst_lim[2]) * mapstep_c);
    nzMax = (int) (min (drvui->frames[drvui->frame_no].map_lim[5], drvui->frames[drvui->frame_no].cryst_lim[5]) * mapstep_c);
    snx = abs (nxMax - nxMin + 1);
    sny = abs (nyMax - nyMin + 1);
    snz = abs (nzMax - nzMin + 1);
    iMax = nxMax + 1;
    jMax = nyMax + 1;
    kMax = nzMax + 1;
    if (snx == 1) {
	iproj = 1;
	snx++;
	iMax++;
    }
    if (sny == 1) {
	iproj = 2;
	sny++;
	jMax++;
    }
    if (snz == 1) {
	iproj = 3;
	snz++;
	kMax++;
    }
    sny_snz = sny * snz;

    ijkmax= (int) (drvui->frames[drvui->frame_no].cryst_lim[3] * (mapstep_a-1)) * mapstep_b * mapstep_c
	    + (int) (drvui->frames[drvui->frame_no].cryst_lim[4] * (mapstep_b-1)) * mapstep_c
	    + (int) (drvui->frames[drvui->frame_no].cryst_lim[5] * (mapstep_c-1));

/* create a new set of Fourier points */

    PointsList =
	(mp4Vector *) zalloc (sizeof (mp4Vector) * (snx + 1) * (sny + 1) * (snz + 1));
    if (PointsList == NULL) {
	Error_Box ("Error allocating space for PointsList vectors.");
	return;
    }
    for (i = nxMin, si = 0; i <= iMax; i++, si++) {
	ni = i;
	ni = (5 * mapstep_a + ni) % mapstep_a;	/* this will 'wrap' around any value (negative or positive) */
	sni = si * sny_snz;
	fvert[0] = i * step[0];
	if (iproj == 1) {
	    ni = nxMin;
	    fvert[0] = ni * step[0] + 0.001f;
	}
	for (j = nyMin, sj = 0; j <= jMax; j++, sj++) {
	    nj = j;
	    nj = (5 * mapstep_b + nj) % mapstep_b;
	    snj = sj * snz;
	    fvert[1] = j * step[1];
	    if (iproj == 2) {
		nj = nyMin;
		fvert[1] = nj * step[1] + 0.001f;
	    }
	    for (k = nzMin, sk = 0; k <= kMax; k++, sk++) {
		nk = k;
		nk = (5 * mapstep_c + nk) % mapstep_c;
		sind = sni + snj + sk;
		fvert[2] = k * step[2];
		if (iproj == 3) {
		    nk = nzMin;
		    fvert[2] = nk * step[2] + 0.001f;
		}

/* convert Fractional to Orthonormal Coords */

		for (ii = 0; ii <= 2; ++ii) {	/* convert vertex coordinates to Cartesian */
		    cvert[ii] = 0.0f;
		    for (jj = 0; jj <= 2; ++jj)
			cvert[ii] +=
			    (float) drvui->b_mat[ii][jj] * (fvert[jj] - origin[jj]);
		}
		PointsList[sind].x = cvert[0];
		PointsList[sind].y = cvert[1];
		PointsList[sind].z = cvert[2];
		ijk = ni * mapstep_b * mapstep_c + nj * mapstep_c + nk;
		if (ijk > ijkmax) 
		    ijk-=ijkmax;
		PointsList[sind].val = FourierPt[ijk];
	    }
	}
    }

/* create contours from grid of points with matching cubes */

    pTriangles =
	MC (snx - 1, sny - 1, snz - 1, step[0], step[1], step[2], minValue, PointsList,
	    numOfTriangles);
    free (PointsList);

    if (numOfTriangles <= 0)
	return;

    glPushMatrix ();
    if (!drvui->frames[drvui->frame_no].slice) {
	if (!Solid) {
	    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	    glDisable (GL_LIGHTING);
	} else
	    glMaterialf (GL_FRONT, GL_SHININESS, 0.0);	/* disable shinyness */
	glColor3f (Red, Green, Blue);

	glBegin (GL_TRIANGLES);

	for (i = 0; i < numOfTriangles; i++) {
	    d1 = (pTriangles[i].p[0].x - pTriangles[i].p[1].x) *
		(pTriangles[i].p[0].x - pTriangles[i].p[1].x) +
		(pTriangles[i].p[0].y - pTriangles[i].p[1].y) *
		(pTriangles[i].p[0].y - pTriangles[i].p[1].y) +
		(pTriangles[i].p[0].z - pTriangles[i].p[1].z) *
		(pTriangles[i].p[0].z - pTriangles[i].p[1].z);
	    d2 = (pTriangles[i].p[2].x - pTriangles[i].p[1].x) *
		(pTriangles[i].p[2].x - pTriangles[i].p[1].x) +
		(pTriangles[i].p[2].y - pTriangles[i].p[1].y) *
		(pTriangles[i].p[2].y - pTriangles[i].p[1].y) +
		(pTriangles[i].p[2].z - pTriangles[i].p[1].z) *
		(pTriangles[i].p[2].z - pTriangles[i].p[1].z);
	    d3 = (pTriangles[i].p[2].x - pTriangles[i].p[0].x) *
		(pTriangles[i].p[2].x - pTriangles[i].p[0].x) +
		(pTriangles[i].p[2].y - pTriangles[i].p[0].y) *
		(pTriangles[i].p[2].y - pTriangles[i].p[0].y) +
		(pTriangles[i].p[2].z - pTriangles[i].p[0].z) *
		(pTriangles[i].p[2].z - pTriangles[i].p[0].z);

// skip excursions completely across the cell when a contour reaches the edge
	    if ((d1 > 5.*snx*snx* step[0] ) || (d2 > 5.*sny*sny* step[1]) ||
		(d3 > 5.*snz*snz* step[2]))
		continue;

	    if (minValue > 0) {
		if (doVrml) {
		    if (!Vrml2) {
			fprintf (drvui->fpoutv, " Separator {\n");
			fprintf (drvui->fpoutv, " Material {  diffuseColor %s }\n", string);
			fprintf (drvui->fpoutv, "Coordinate3 { point [\n");
		    } else {
			fprintf (drvui->fpoutv, " Shape {\n");
			if (Solid) {
			    fprintf (drvui->fpoutv, "  appearance Appearance {\n");
			    fprintf (drvui->fpoutv,
				    "   material Material {diffuseColor %s}\n", string);
			    fprintf (drvui->fpoutv,
				    "  }\n  geometry IndexedFaceSet { coord Coordinate{ point [\n");
			} else
			    fprintf (drvui->fpoutv,
				    "  geometry IndexedLineSet { coord Coordinate{ point [\n");
		    }
		}
		if (doPOV) {
		    if (!Solid) {
			if (d1 > 0.00001f) {
			    fprintf (drvui->fpoutp, " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
				    pTriangles[i].p[0].x, pTriangles[i].p[0].y,
				    pTriangles[i].p[0].z, pTriangles[i].p[1].x,
		 		    pTriangles[i].p[1].y, pTriangles[i].p[1].z,
				    0.0001 * Scale);
			    fprintf (drvui->fpoutp, "  texture{pigment{color %s  }}\n }\n",
				    Color2);
			}
			if (d2 > 0.00001f) {
			    fprintf (drvui->fpoutp, " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
				    pTriangles[i].p[1].x, pTriangles[i].p[1].y,
				    pTriangles[i].p[1].z, pTriangles[i].p[2].x,
				    pTriangles[i].p[2].y, pTriangles[i].p[2].z,
				    0.0001 * Scale);
			    fprintf (drvui->fpoutp, "  texture{pigment{color %s  }}\n }\n",
				    Color2);
			}
			if (d3 > 0.00001f) {
			    fprintf (drvui->fpoutp, " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
				    pTriangles[i].p[2].x, pTriangles[i].p[2].y,
				    pTriangles[i].p[2].z, pTriangles[i].p[0].x,
				    pTriangles[i].p[0].y, pTriangles[i].p[0].z,
				    0.0001 * Scale);
			    fprintf (drvui->fpoutp, "  texture{pigment{color %s  }}\n }\n",
				    Color2);
			}
		    }
		    if (Solid)
			fprintf (drvui->fpoutp, "smooth_triangle {\n");
		}
		for (j = 0; j < 3; j++) {
		    if (Solid)
			glNormal3f (pTriangles[i].norm[j].x, pTriangles[i].norm[j].y,
				    pTriangles[i].norm[j].z);
		    glVertex3f (pTriangles[i].p[j].x, pTriangles[i].p[j].y,
				pTriangles[i].p[j].z);
		    if (Solid && doPOV) {
			fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f>,", pTriangles[i].p[j].x,
				pTriangles[i].p[j].y, pTriangles[i].p[j].z);
			fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f>", pTriangles[i].norm[j].x,
				pTriangles[i].norm[j].y, pTriangles[i].norm[j].z);
		    }
		    if (doVrml)
			fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f", pTriangles[i].p[j].x,
				pTriangles[i].p[j].y, pTriangles[i].p[j].z);

		    if (j < 2) {
			if (Solid && doPOV)
			    fprintf (drvui->fpoutp, ",");
			if (doVrml)
			    fprintf (drvui->fpoutv, ",\n");
		    } else {
			if (Solid && doPOV) {
			    fprintf (drvui->fpoutp, "\n texture{pigment{color %s }}",
				     Color2);
			    if (BackColor)
				fprintf (drvui->fpoutp, "\n interior_texture{pigment{color %s }}\n }\n",
					 BackColor);
			    else
				fprintf(drvui->fpoutp, "\n }\n");
                        }
			if (doVrml)
			    fprintf (drvui->fpoutv, "]\n");
		    }
		}		// for (j... 
		if (doAsy) {
		    if (!Solid)
			fprintf (drvui->fpouta, "draw(pic, (%5.2f,%5.2f,%5.2f)--(%5.2f,%5.2f,%5.2f)--(%5.2f,%5.2f,%5.2f)--cycle,rgb(%.2f,%.2f,%.2f)+linewidth(%.2f));\n",
				 pTriangles[i].p[0].x, pTriangles[i].p[0].y,pTriangles[i].p[0].z,
				 pTriangles[i].p[1].x, pTriangles[i].p[1].y,pTriangles[i].p[1].z,
				 pTriangles[i].p[2].x, pTriangles[i].p[2].y,pTriangles[i].p[2].z,
				 Red,Green,Blue,0.001*Scale);
		    else
			fprintf (drvui->fpouta, "draw(pic, surface( (%5.2f,%5.2f,%5.2f)--(%5.2f,%5.2f,%5.2f)--(%5.2f,%5.2f,%5.2f)--cycle),rgb(%.2f,%.2f,%.2f));\n",
				 pTriangles[i].p[0].x, pTriangles[i].p[0].y,pTriangles[i].p[0].z,
				 pTriangles[i].p[1].x, pTriangles[i].p[1].y,pTriangles[i].p[1].z,
				 pTriangles[i].p[2].x, pTriangles[i].p[2].y,pTriangles[i].p[2].z,
				 Red,Green,Blue);
		}
		if (doVrml) {
		    if (Vrml2) {
			if (Solid) {
			    fprintf (drvui->fpoutv,
				    "  }\n  coordIndex [ 0,1,2,-1]\n solid FALSE\n convex"
				    " TRUE\n creaseAngle 1.57075\n }\n }\n");
			} else {
			    fprintf (drvui->fpoutv,
				    "  }\n  coordIndex [ 0,1,2,-1]\n color Color { color [%s]"
				    "}\n colorIndex [0]\n colorPerVertex FALSE\n }\n }\n",
				    string);
			}
		    } else {
			if (Solid) {
			    fprintf (drvui->fpoutv,
				    "}\n IndexedFaceSet { coordIndex [0,1,2,-1] }\n}\n");
			} else {
			    fprintf (drvui->fpoutv,
				    "}\n IndexedLineSet { coordIndex [0,1,2,-1] }\n}\n");
			}
		    }
		}
	    } else {		// if (minValue > 0)
		if (doVrml) {
		    if (!Vrml2) {
			fprintf (drvui->fpoutv, " Separator {\n");
			fprintf (drvui->fpoutv, " Material {  diffuseColor %s }\n", string);
			fprintf (drvui->fpoutv, "Coordinate3 { point [\n");
		    } else {
			fprintf (drvui->fpoutv, " Shape {\n");
			if (Solid) {
			    fprintf (drvui->fpoutv, "  appearance Appearance {\n");
			    fprintf (drvui->fpoutv,
				    "   material Material {diffuseColor %s}\n", string);
			    fprintf (drvui->fpoutv,
				    "  }\n  geometry IndexedFaceSet { coord Coordinate{ point [\n");
			} else
			    fprintf (drvui->fpoutv,
				    "  geometry IndexedLineSet { coord Coordinate{ point [\n");
		    }
		}
		if (doPOV) {
		    if (!Solid) {
			if (d1 > 0.00001f) {	// skip zero length cylinder
			    fprintf (drvui->fpoutp, " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
				    pTriangles[i].p[0].x, pTriangles[i].p[0].y,
				    pTriangles[i].p[0].z, pTriangles[i].p[1].x,
				    pTriangles[i].p[1].y, pTriangles[i].p[1].z,
				    0.0001 * Scale);
			    fprintf (drvui->fpoutp, "  texture{pigment{color %s  }}\n }\n",
				    Color2);
			}
			if (d2 > 0.00001f) {	// skip zero length cylinder
			    fprintf (drvui->fpoutp, " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
				    pTriangles[i].p[1].x, pTriangles[i].p[1].y,
				    pTriangles[i].p[1].z, pTriangles[i].p[2].x,
				    pTriangles[i].p[2].y, pTriangles[i].p[2].z,
				    0.0001 * Scale);
			    fprintf (drvui->fpoutp, "  texture{pigment{color %s  }}\n }\n",
				    Color2);
			}
			if (d3 > 0.00001f) {	// skip zero length cylinder
			    fprintf (drvui->fpoutp, " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
				    pTriangles[i].p[2].x, pTriangles[i].p[2].y,
				    pTriangles[i].p[2].z, pTriangles[i].p[0].x,
				    pTriangles[i].p[0].y, pTriangles[i].p[0].z,
				    0.0001 * Scale);
			    fprintf (drvui->fpoutp, "  texture{pigment{color %s  }}\n }\n",
				    Color2);
			}
		    }
		    if (Solid)
			fprintf (drvui->fpoutp, "smooth_triangle {\n");
		}
		for (j = 2; j >= 0; j--) {
		    if (Solid)
			glNormal3f (pTriangles[i].norm[j].x, pTriangles[i].norm[j].y,
				    pTriangles[i].norm[j].z);
		    glVertex3f (pTriangles[i].p[j].x, pTriangles[i].p[j].y,
				    pTriangles[i].p[j].z);

		    if (Solid && doPOV) {
			fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f>,", pTriangles[i].p[j].x,
				pTriangles[i].p[j].y, pTriangles[i].p[j].z);
			fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f>", pTriangles[i].norm[j].x,
				pTriangles[i].norm[j].y, pTriangles[i].norm[j].z);
		    }
		    if (doVrml)
			fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f", pTriangles[i].p[j].x,
				pTriangles[i].p[j].y, pTriangles[i].p[j].z);
		    if (j > 0) {
			if (Solid && doPOV)
			    fprintf (drvui->fpoutp, ",");
			if (doVrml)
			    fprintf (drvui->fpoutv, ",\n");
		    } else {
			if (Solid && doPOV) {
			    fprintf (drvui->fpoutp, "\n texture{pigment{color %s }}",
				     Color2);
			    if (BackColor)
				fprintf (drvui->fpoutp, "\n interior_texture{pigment{color %s }}\n }\n",
					 BackColor);
			    else
				fprintf(drvui->fpoutp, "\n }\n");
                        }
			if (doVrml)
			    fprintf (drvui->fpoutv, "]\n");
		    }
		}
		if (doAsy) {
		    if (!Solid)
			fprintf (drvui->fpouta, "draw(pic, (%5.2f,%5.2f,%5.2f)--(%5.2f,%5.2f,%5.2f)--(%5.2f,%5.2f,%5.2f)--cycle,rgb(%.2f,%.2f,%.2f)+linewidth(%.2f));\n",
				 pTriangles[i].p[0].x, pTriangles[i].p[0].y,pTriangles[i].p[0].z,
				 pTriangles[i].p[1].x, pTriangles[i].p[1].y,pTriangles[i].p[1].z,
				 pTriangles[i].p[2].x, pTriangles[i].p[2].y,pTriangles[i].p[2].z,
				 Red,Green,Blue,0.001*Scale);
		    else
			fprintf (drvui->fpouta, "draw(pic, surface( (%5.2f,%5.2f,%5.2f)--(%5.2f,%5.2f,%5.2f)--(%5.2f,%5.2f,%5.2f)--cycle),rgb(%.2f,%.2f,%.2f));\n",
				 pTriangles[i].p[0].x, pTriangles[i].p[0].y,pTriangles[i].p[0].z,
				 pTriangles[i].p[1].x, pTriangles[i].p[1].y,pTriangles[i].p[1].z,
				 pTriangles[i].p[2].x, pTriangles[i].p[2].y,pTriangles[i].p[2].z,
				 Red,Green,Blue);
		}
		if (doVrml) {
		    if (Vrml2) {
			if (Solid) {
			    fprintf (drvui->fpoutv,
				    "}\n coordIndex [ 1,2,3,-1]\n solid FALSE\n convex TRUE\n"
				    " creaseAngle 1.57075\n }\n }\n");
			} else {
			    fprintf (drvui->fpoutv,
				    "  }\n  coordIndex [ 0,1,2,-1]\n color Color { color [%s]"
				    "}\n colorIndex [0]\n colorPerVertex FALSE\n }\n }\n",
				    string);
			}
		    } else {
			if (Solid) {
			    fprintf (drvui->fpoutv,
				    "}\n IndexedFaceSet { coordIndex [1,2,3,-1] }\n}\n");
			} else {
			    fprintf (drvui->fpoutv,
				    "}\n IndexedLineSet { coordIndex [1,2,3,-1] }\n}\n");
			}
		    }
		}
	    }			// if (minValue > 0)
	}
	glEnd ();
	if (!Solid) {
	    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	    glEnable (GL_LIGHTING);
	}

    } else { // 2D slices
	mpVector p0,p1;

	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	glDisable (GL_LIGHTING);
	glMaterialf (GL_FRONT, GL_SHININESS, 0.0);	/* disable shinyness */
	glColor3f (Red, Green, Blue);
	glBegin (GL_LINES);
	for (i = 0; i < numOfTriangles; i++) {
	    d1 = (pTriangles[i].p[0].x - pTriangles[i].p[1].x) *
		(pTriangles[i].p[0].x - pTriangles[i].p[1].x) +
		(pTriangles[i].p[0].y - pTriangles[i].p[1].y) *
		(pTriangles[i].p[0].y - pTriangles[i].p[1].y) +
		(pTriangles[i].p[0].z - pTriangles[i].p[1].z) *
		(pTriangles[i].p[0].z - pTriangles[i].p[1].z);
	    d2 = (pTriangles[i].p[2].x - pTriangles[i].p[1].x) *
		(pTriangles[i].p[2].x - pTriangles[i].p[1].x) +
		(pTriangles[i].p[2].y - pTriangles[i].p[1].y) *
		(pTriangles[i].p[2].y - pTriangles[i].p[1].y) +
		(pTriangles[i].p[2].z - pTriangles[i].p[1].z) *
		(pTriangles[i].p[2].z - pTriangles[i].p[1].z);
	    d3 = (pTriangles[i].p[2].x - pTriangles[i].p[0].x) *
		(pTriangles[i].p[2].x - pTriangles[i].p[0].x) +
		(pTriangles[i].p[2].y - pTriangles[i].p[0].y) *
		(pTriangles[i].p[2].y - pTriangles[i].p[0].y) +
		(pTriangles[i].p[2].z - pTriangles[i].p[0].z) *
		(pTriangles[i].p[2].z - pTriangles[i].p[0].z);

// skip excursions completely across the cell when a contour reaches the edge
	    if ((d1 > 5.*snx*snx* step[0] ) || (d2 > 5.*sny*sny* step[1]) ||
		(d3 > 5.*snz*snz* step[2]))
		continue;

	    if (ContourFacet(pTriangles[i],&p0,&p1) == 2) {
		glVertex3f(p0.x,p0.y,p0.z);
		glVertex3f(p1.x,p1.y,p1.z);

		if (doVrml) {
		    if (!Vrml2) {
			fprintf (drvui->fpoutv, " Separator {\n");
			fprintf (drvui->fpoutv, " Material {  diffuseColor %s }\n", string);
			fprintf (drvui->fpoutv, "Coordinate3 { point [\n");
		    } else {
			fprintf (drvui->fpoutv, " Shape {\n");
			fprintf (drvui->fpoutv,
				"  geometry IndexedLineSet { coord Coordinate{ point [\n");
		    }
		    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f,\n", p0.x,p0.y,p0.z);
		    fprintf (drvui->fpoutv, " %5.3f %5.3f %5.3f]\n", p1.x,p1.y,p1.z);
		    if (Vrml2) {
			fprintf (drvui->fpoutv,
				"  }\n  coordIndex [ 0,1,-1]\n color Color { color [%s]"
				"}\n colorIndex [0]\n colorPerVertex FALSE\n }\n }\n",
				string);
		    } else {
			fprintf (drvui->fpoutv,
				"}\n IndexedLineSet { coordIndex [0,1,-1] }\n}\n");
		    }
		}
		if (doPOV) {
		    if (fabs(p0.x-p1.x) > 0.00001f || fabs(p0.y-p1.y) > 0.00001f 
			|| fabs(p0.z-p1.z) >0.00001f) {
			fprintf (drvui->fpoutp, " cylinder { <%f,%f,%f>,<%f,%f,%f>,%f\n",
				p0.x, p0.y,p0.z, p1.x,p1.y,p1.z,
				0.0001 * Scale);
			fprintf (drvui->fpoutp, "  texture{pigment{color %s  }}\n }\n",
				Color2);
		    }
		}
		if (doAsy) {
		    if (fabs(p0.x-p1.x) > 0.00001f || fabs(p0.y-p1.y) > 0.00001f 
			|| fabs(p0.z-p1.z) >0.00001f) {
			fprintf (drvui->fpouta, " draw(pic, (<%8.5f,%8.5f,%8.5f)--(%8.5f,%8.5f,%8.5f),\n",
				p0.x, p0.y,p0.z, p1.x,p1.y,p1.z);
			fprintf (drvui->fpouta,"rgb(%4.2f,%4.2f,%4.2f) );\n",Red,Green,Blue);
		    }
		}
	    }
	}
        glEnd();
        glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
        glEnable (GL_LIGHTING);
    }
    glPopMatrix ();
    free (pTriangles);		/* free the contour array */
    drvui->mainWindow->cursor (FL_CURSOR_DEFAULT);
}

/* interpolate rho on a 3-D grid */

float
InterpolateMap (float x, float y, float z)
{

/* for simplicity, use weighting appropriate for a cubic grid (should be better than no interpolation) */

    int ix, iy, iz;

    int ix1, iy1, iz1;

    float Xmin, Ymin, Zmin, Xmax, Ymax, Zmax, AvgRho;

/* Compute the grid points on either size of x, then map those points onto the 3D array indices */

    ix = ((int) (x * mapstep_a));
    Xmin = ((float) ix) / mapstep_a;
    Xmax = ((float) (ix + 1)) / mapstep_a;
    ix = ((int) (x * mapstep_a)) % mapstep_a;
    ix1 = (ix + 1) % mapstep_a;

/* ditto for y & z */

    iy = ((int) (y * mapstep_b));
    Ymin = ((float) iy) / mapstep_b;
    Ymax = ((float) (iy + 1)) / mapstep_b;
    iy = ((int) (y * mapstep_b)) % mapstep_b;
    iy1 = (iy + 1) % mapstep_b;
    iz = ((int) (z * mapstep_c));
    Zmin = ((float) iz) / mapstep_c;
    Zmax = ((float) (iz + 1)) / mapstep_c;
    iz = ((int) (z * mapstep_c)) % mapstep_c;
    iz1 = (iz + 1) % mapstep_c;

    AvgRho = FourierPt[LookupMap (ix, iy, iz)] * (Xmax - x) * (Ymax - y) * (Zmax - z);
    AvgRho += FourierPt[LookupMap (ix1, iy, iz)] * (x - Xmin) * (Ymax - y) * (Zmax - z);
    AvgRho += FourierPt[LookupMap (ix1, iy1, iz)] * (x - Xmin) * (y - Ymin) * (Zmax - z);
    AvgRho += FourierPt[LookupMap (ix1, iy, iz1)] * (x - Xmin) * (Ymax - y) * (z - Zmin);
    AvgRho += FourierPt[LookupMap (ix, iy1, iz)] * (Xmax - x) * (y - Ymin) * (Zmax - z);
    AvgRho += FourierPt[LookupMap (ix, iy1, iz1)] * (Xmax - x) * (y - Ymin) * (z - Zmin);
    AvgRho += FourierPt[LookupMap (ix, iy, iz1)] * (Xmax - x) * (Ymax - y) * (z - Zmin);
    AvgRho += FourierPt[LookupMap (ix1, iy1, iz1)] * (x - Xmin) * (y - Ymin) * (z - Zmin);

    AvgRho /= (Xmax - Xmin) * (Ymax - Ymin) * (Zmax - Zmin);

    return AvgRho;
}

/* lookup a point in the Fourier map */

int
LookupMap (int ix, int iy, int iz)
{

    return abs (ix * mapstep_b * mapstep_c + iy * mapstep_c +
		((mapstep_c + iz % mapstep_c) % mapstep_c));
}

int
Maximize_rho (int Sense)
{
// find min/max in rho - Sense is -1 if min, +1 if max
    float rhomax = -1.0e15f;

    int i, j, k, l;

    float w0[3][3][3];

    int maxpt;

    float param[7];

    int map[3];

    if (!ReadFourMap)
	return (-1);
// extract 3x3x3 array of points around the current position
    for (i = 0; i < 7; i++)
	param[i] = 0.0f;
    map[0] = (int) ((float) mapstep_a * cur_cen[0] + 0.5);
    map[1] = (int) ((float) mapstep_b * cur_cen[1] + 0.5);
    map[2] = (int) ((float) mapstep_c * cur_cen[2] + 0.5);
    while (1) {
	rhomax = -1.0e15f;
	maxpt = -1;
	for (i = 0; i < 3; i++) {
	    for (j = 0; j < 3; j++) {
		for (k = 0; k < 3; k++) {
		    l = (map[0] + i - 1) * mapstep_b * mapstep_c
			+ (map[1] + j - 1) * mapstep_c + map[2] + k - 1;
		    w0[i][j][k] = (float) Sense *FourierPt[l];

		    if (w0[i][j][k] > rhomax) {
			rhomax = w0[i][j][k];
			maxpt = (i * 100) + (j * 10) + k;
		    }
		}
	    }
	}
	if (maxpt == 111)
	    break;		// maximum in middle
	k = maxpt % 10;		// shift maximum to middle
	j = maxpt / 10 % 10;
	i = maxpt / 100;
	map[2] += k - 1;
	map[1] += j - 1;
	map[0] += i - 1;
	i = 0;
	if (map[0] < 1) {
	    map[0] = 1;
	    i = 1;
	}
	if (map[1] < 1) {
	    map[1] = 1;
	    i = 1;
	}
	if (map[2] < 1) {
	    map[2] = 1;
	    i = 1;
	}
	if (map[0] > mapstep_a - 1) {
	    map[0] = mapstep_a - 1;
	    i = 1;
	}
	if (map[1] > mapstep_b - 1) {
	    map[1] = mapstep_b - 1;
	    i = 1;
	}
	if (map[2] > mapstep_c - 1) {
	    map[2] = mapstep_c - 1;
	    i = 1;
	}
	if (i)
	    break;
	cur_cen[2] = map[2] / (float) mapstep_c;
	cur_cen[1] = map[1] / (float) mapstep_b;
	cur_cen[0] = map[0] / (float) mapstep_a;
    }
    param[3] = -w0[1][1][1] + (w0[0][1][1] + w0[2][1][1]) * 0.5f;
    param[4] = -w0[1][1][1] + (w0[1][0][1] + w0[1][2][1]) * 0.5f;
    param[5] = -w0[1][1][1] + (w0[1][1][0] + w0[1][1][2]) * 0.5f;
    if (param[3] != 0.0f)
	param[0] = (w0[0][1][1] - w0[2][1][1]) * 0.25f / param[3];
    if (param[4] != 0.0f)
	param[1] = (w0[1][0][1] - w0[1][2][1]) * 0.25f / param[4];
    if (param[5] != 0.0f)
	param[2] = (w0[1][1][0] - w0[1][1][2]) * 0.25f / param[5];
    param[6] = w0[1][1][1] - param[3] * param[0] * param[0]
	- param[4] * param[1] * param[1] - param[5] * param[2] * param[2];

    cur_cen[0] += param[0] / mapstep_a;
    cur_cen[1] += param[1] / mapstep_b;
    cur_cen[2] += param[2] / mapstep_c;
    i = 0;
    for (k = 0; k < nvert; k++) {
	if (fabs (o_vert[3 * k] - cur_cen[0]) + fabs (o_vert[3 * k + 1] - cur_cen[1])
	    + fabs (o_vert[3 * k + 2] - cur_cen[2]) < 0.001)
	    i = k;
    }
    if (i)
	return i;		// position already in vertex list
    drvui->orig_atom_no[nvert] = -1;	// new position - add it
    for (k = 0; k < 3; k++) {
	o_vert[3 * nvert + k] = cur_cen[k];
	s_vert[3 * nvert + k] = 0.0f;
	for (l = 0; l < 3; l++) {	/* calculate cartesian coordinates */
	    s_vert[3 * nvert + k] +=
		(float) drvui->b_mat[k][l] * (cur_cen[l] - origin[l]);
	}
    }
    vert_sym_no[nvert++] = 0;
    if (!check_vert_alloc (nvert, 1))
	return (-1);
    return (nvert - 1);
}

// 3D variant of Paul Bourke's conrec taken from
//   http://local.wasp.uwa.edu.au/~pbourke/papers/conrec/
//
//-------------------------------------------------------------------------
//   Create a contour slice through a 3 vertex facet "p"
//   Given the normal of the cutting plane "n" and a point on the plane "p0"
//   Return
//       0 if the contour plane doesn't cut the facet
//       2 if it does cut the facet, the contour line segment is p1->p2
//      -1 for an unexpected occurence
//   If a vertex touches the contour plane nothing need to be drawn!?
//   Note: the following has been written as a "stand along" piece of
//   code that will work but is far from efficient....
//
int ContourFacet (TRIANGLE tri, mpVector *p1, mpVector *p2)
{
#define SIGN(x) (x>0? 1: 0)
    double A,B,C,D;
    double side[3];
  
//      Determine the equation of the plane as
//      Ax + By + Cz + D = 0
    A = drvui->frames[drvui->frame_no].planeeq[0];
    B = drvui->frames[drvui->frame_no].planeeq[1];
    C = drvui->frames[drvui->frame_no].planeeq[2];
    D = drvui->frames[drvui->frame_no].planeeq[3];

//      Evaluate the equation of the plane for each vertex
//      If side < 0 then it is on the side to be retained
//      else it is to be clipped

    side[0] = A*tri.p[0].x + B*tri.p[0].y + C*tri.p[0].z + D;
    side[1] = A*tri.p[1].x + B*tri.p[1].y + C*tri.p[1].z + D;
    side[2] = A*tri.p[2].x + B*tri.p[2].y + C*tri.p[2].z + D;
 
//   Are all the vertices on the same side 
    if (side[0] >= 0 && side[1] >= 0 && side[2] >= 0)
	return(0);
    if (side[0] <= 0 && side[1] <= 0 && side[2] <= 0)
	return(0);

//   Is p0 the only point on a side by itself 
    if ((SIGN(side[0]) != SIGN(side[1])) && (SIGN(side[0]) != SIGN(side[2]))) {
	p1->x = (float)(tri.p[0].x - side[0] * (tri.p[2].x - tri.p[0].x) / (side[2] - side[0]));
	p1->y = (float)(tri.p[0].y - side[0] * (tri.p[2].y - tri.p[0].y) / (side[2] - side[0]));
	p1->z = (float)(tri.p[0].z - side[0] * (tri.p[2].z - tri.p[0].z) / (side[2] - side[0]));
	p2->x = (float)(tri.p[0].x - side[0] * (tri.p[1].x - tri.p[0].x) / (side[1] - side[0]));
	p2->y = (float)(tri.p[0].y - side[0] * (tri.p[1].y - tri.p[0].y) / (side[1] - side[0]));
	p2->z = (float)(tri.p[0].z - side[0] * (tri.p[1].z - tri.p[0].z) / (side[1] - side[0]));
	return(2);
    }

// Is p1 the only point on a side by itself 
    if ((SIGN(side[1]) != SIGN(side[0])) && (SIGN(side[1]) != SIGN(side[2]))) {
	p1->x = (float)(tri.p[1].x - side[1] * (tri.p[2].x - tri.p[1].x) / (side[2] - side[1]));
	p1->y = (float)(tri.p[1].y - side[1] * (tri.p[2].y - tri.p[1].y) / (side[2] - side[1]));
	p1->z = (float)(tri.p[1].z - side[1] * (tri.p[2].z - tri.p[1].z) / (side[2] - side[1]));
	p2->x = (float)(tri.p[1].x - side[1] * (tri.p[0].x - tri.p[1].x) / (side[0] - side[1]));
	p2->y = (float)(tri.p[1].y - side[1] * (tri.p[0].y - tri.p[1].y) / (side[0] - side[1]));
	p2->z = (float)(tri.p[1].z - side[1] * (tri.p[0].z - tri.p[1].z) / (side[0] - side[1]));
	return(2);
    }

// Is p2 the only point on a side by itself 
    if ((SIGN(side[2]) != SIGN(side[0])) && (SIGN(side[2]) != SIGN(side[1]))) {
	p1->x = (float)(tri.p[2].x - side[2] * (tri.p[0].x - tri.p[2].x) / (side[0] - side[2]));
	p1->y = (float)(tri.p[2].y - side[2] * (tri.p[0].y - tri.p[2].y) / (side[0] - side[2]));
	p1->z = (float)(tri.p[2].z - side[2] * (tri.p[0].z - tri.p[2].z) / (side[0] - side[2]));
	p2->x = (float)(tri.p[2].x - side[2] * (tri.p[1].x - tri.p[2].x) / (side[1] - side[2]));
	p2->y = (float)(tri.p[2].y - side[2] * (tri.p[1].y - tri.p[2].y) / (side[1] - side[2]));
	p2->z = (float)(tri.p[2].z - side[2] * (tri.p[1].z - tri.p[2].z) / (side[1] - side[2]));
	return(2);
    }

// Shouldn't get here 
    return(-1);
}

/*                                                                         */
/* BW and "Hot-Cold" Colormaps from Paul Bourke's                          */
/* http://local.wasp.uwa.edu.au/~pbourke/texture_colour/colourramp/        */
/*                                                                         */
void colorramp(float rho,float *r, float *g, float *b)
{
float dv=Map_Info.rhomx-Map_Info.rhomn;

    if (drvui->frames[drvui->frame_no].slice == 3) {
	*r=*g=*b=(rho-Map_Info.rhomn)/(Map_Info.rhomx-Map_Info.rhomn);
	return;
    }
    *r = *g = *b = 1.0f;
    if (rho < (Map_Info.rhomn +0.25f *dv) ) {
	*r = 0.;
	*g = 4.0f * (rho - Map_Info.rhomn)/dv;
    } else if (rho < (Map_Info.rhomn +0.5f * dv)) {
	*r = 0.;
	*b = 1.0f + 4.0f * (Map_Info.rhomn +0.25f * dv - rho) / dv;
    } else if (rho < (Map_Info.rhomn + 0.75f * dv)) {
	*r = 4.0f * (rho -Map_Info.rhomn - 0.5f * dv) / dv;
	*b = 0.0f;
    } else {
	*g = 1.0f + 4.0f * (Map_Info.rhomn + 0.75f * dv - rho) /dv; 
	*b = 0.0f;
    }
}

void
Add_mapslice (int type)
{
// add a 2d map section through the last three atoms 

   float p[3],q[3],pq[3];

    if (!ReadFourMap)
	return;
    if (cur_atom[3]>0) {
	p[0]=o_vert[3*cur_atom[1]]   - cur_cen[0];
	p[1]=o_vert[3*cur_atom[1]+1] - cur_cen[1];
	p[2]=o_vert[3*cur_atom[1]+2] - cur_cen[2];
	q[0]=o_vert[3*cur_atom[2]]   - cur_cen[0];
	q[1]=o_vert[3*cur_atom[2]+1] - cur_cen[1];
	q[2]=o_vert[3*cur_atom[2]+2] - cur_cen[2];
    } else {
	p[0]=o_vert[3*cur_atom[0]]   - cur_cen[0];
	p[1]=o_vert[3*cur_atom[0]+1] - cur_cen[1];
	p[2]=o_vert[3*cur_atom[0]+2] - cur_cen[2];
	q[0]=o_vert[3*cur_atom[1]]   - cur_cen[0];
	q[1]=o_vert[3*cur_atom[1]+1] - cur_cen[1];
	q[2]=o_vert[3*cur_atom[1]+2] - cur_cen[2];
    }
    vcross(p,q,pq);
    vnormalize(pq);
    drvui->frames[drvui->frame_no].mapslice[0]=cur_cen[0];
    drvui->frames[drvui->frame_no].mapslice[1]=cur_cen[1];
    drvui->frames[drvui->frame_no].mapslice[2]=cur_cen[2];
    drvui->frames[drvui->frame_no].mapnorm[0]=pq[0];
    drvui->frames[drvui->frame_no].mapnorm[1]=pq[1];
    drvui->frames[drvui->frame_no].mapnorm[2]=pq[2];
    drvui->frames[drvui->frame_no].slice=type;
    drvui->Str_File_Changed = 1;
    Update_Str (0);
    Generate_Drawing (0);
}

void Mul_Rv(float R[3][3], float inp[3], float oup[3])
{
    /* Multiply R * inp => oup */
    int i, j;

    for (i = 0; i < 3; i++) {
	oup[i] = 0.0f;
	for (j = 0; j < 3; j++)
	    oup[i] += R[i][j] * inp[j];
    }
}

void Mul_Rinvv(float R[3][3], float inp[3], float oup[3])
{
    /* Multiply RT * inp => oup */
    int i, j;

    for (i = 0; i < 3; i++) {
	oup[i] = 0.0f;
	for (j = 0; j < 3; j++)
	    oup[i] += R[j][i] * inp[j];
    }
}

void generate_slice (void) 
{
    /* Calculate the density in the plane described by the planeeq information.
     * It lies a distance planeeq[3] from the origin, and has a normal N
     * given by planeeq[0-2].
     * Three different coordinate systems will be used: (1) The crystallographic
     * system x, (2) the corresponding Cartesian system X, and (3) a rotated
     * system X' where the slice is perpendicular to the Z' axis.
     * Routine Convert_Cryst_Cart() converts from crystal to Cartesian.
     * Routine Convert_Cart_Cryst() converts from Cartesian to crystal.
     * Matrix R converts from X to X' and its transpose converts from X' to X.
     */
    int i,j,k,l;
    int dx,dy,dxy;
    float pp[4][3];
    float q[4][3];
    float rho[4],rhor[4],rhog[4],rhob[4];
    float *prhor = NULL, *prhog = NULL, *prhob = NULL;
    int vertn = 0;
    double phi, chi, cosp, sinp, cosc, sinc;
    float min_x = 999999.0f;
    float min_y = 999999.0f;
    float max_x = -999999.0f;
    float max_y = -999999.0f;
    float xlim[2], ylim[2], zlim[2];
    float R[3][3], X[3], XP[3], x[3];
    float bmat_inv[3][3];
    int offset[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
    float D;

    drvui->mainWindow->cursor (FL_CURSOR_WAIT);
    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    bmat_inv[i][j] = (float)drvui->b_mat[i][j];
    matinv(bmat_inv);
    /* The plane normal should be of unit length, but make sure */
    vnormalize(drvui->frames[drvui->frame_no].planeeq);
    D = drvui->frames[drvui->frame_no].planeeq[3];
    /* Build R */
    phi = atan2(drvui->frames[drvui->frame_no].planeeq[1],
		drvui->frames[drvui->frame_no].planeeq[0]);
    cosp = (float)cos(phi);
    sinp = (float)sin(phi);
    chi = atan2((cosp * drvui->frames[drvui->frame_no].planeeq[0] +
		  sinp * drvui->frames[drvui->frame_no].planeeq[1]),
		  -drvui->frames[drvui->frame_no].planeeq[2]);
    sinc = (float)sin(chi);
    cosc = (float)cos(chi);
    R[0][0] = (float)(cosc * cosp);
    R[0][1] = (float)(cosc * sinp);
    R[0][2] = (float)sinc;
    R[1][0] = -(float)sinp;
    R[1][1] = (float)cosp;
    R[1][2] = 0.0f;
    R[2][0] = -(float)(sinc * cosp);
    R[2][1] = -(float)(sinc * sinp);
    R[2][2] = (float)cosc;
    /* find extreme limits in X' for the display box in x */
    xlim[0] = drvui->frames[drvui->frame_no].cryst_lim[0];
    xlim[1] = drvui->frames[drvui->frame_no].cryst_lim[3];
    ylim[0] = drvui->frames[drvui->frame_no].cryst_lim[1];
    ylim[1] = drvui->frames[drvui->frame_no].cryst_lim[4];
    zlim[0] = drvui->frames[drvui->frame_no].cryst_lim[2];
    zlim[1] = drvui->frames[drvui->frame_no].cryst_lim[5];
    for (i = 0; i < 2; i++) {
	x[0] = xlim[i];
	for (j = 0; j < 2; j++) {
	    x[1] = ylim[j];
	    for (k = 0; k < 2; k++) {
		x[2] = zlim[k];
		Convert_Cryst_Cart(drvui->b_mat, x, X, origin);
		Mul_Rv(R, X, XP);
		if (XP[0] < min_x)
		    min_x = XP[0];
		if (XP[0] > max_x)
		    max_x = XP[0];
		if (XP[1] < min_y)
		    min_y = XP[1];
		if (XP[1] > max_y)
		    max_y = XP[1];
	    }
	}
    }

    dx = 50 * (int)(max_x - min_x);
    dy = 50 * (int)(max_y - min_y);
    if (dx > dy) 
	dy = dx;
    else 
	dx = dy;
    fprintf(drvui->flout, "No. of points in 2D section in each direction %d\n", dx);
    min_x = 3.0f * min_x;
    max_x = 3.0f * max_x;
    min_y = 3.0f * min_y;
    max_y = 3.0f * max_y;

//  dxy = dx * dy * 4; // theoretical maximum if no clipping occurs

    dxy=0;
    for (i = 0; i <= dx; i++) {
	XP[0] = min_x + (float)(i) * (max_x - min_x) / (float)(dx);
	for (j = 0; j <= dy; j++) {
	    XP[1] = min_y + (float)(j) * (max_y - min_y) / (float)(dy);
	    XP[2] = D;	/* Now have slice point in XP */
	    for (k = 0; k < 4; k++) {
		float XPP[3];
		XPP[2] = XP[2];
		XPP[1] = XP[1] + offset[k][1] * (max_y - min_y) / (float)(dx);
		XPP[0] = XP[0] + offset[k][0] * (max_x - min_x) / (float)(dy);
		Mul_Rinvv(R, XPP, X);	/* Convert XPP to X */
		Convert_Cart_Cryst(bmat_inv, X, x, origin); /* and X to x */
		if (x[0] < drvui->frames[drvui->frame_no].cryst_lim[0] ||
		    x[0] < drvui->frames[drvui->frame_no].map_lim[0] ||
		    x[0] > drvui->frames[drvui->frame_no].cryst_lim[3] ||
		    x[0] > drvui->frames[drvui->frame_no].map_lim[3])
		    goto eout;
		if (x[1] < drvui->frames[drvui->frame_no].cryst_lim[1] ||
		    x[1] < drvui->frames[drvui->frame_no].map_lim[1] ||
		    x[1] > drvui->frames[drvui->frame_no].cryst_lim[4] ||
		    x[1] > drvui->frames[drvui->frame_no].map_lim[4])
		    goto eout;
		if (x[2] < drvui->frames[drvui->frame_no].cryst_lim[2] ||
		    x[2] < drvui->frames[drvui->frame_no].map_lim[2] ||
		    x[2] > drvui->frames[drvui->frame_no].cryst_lim[5] ||
		    x[2] > drvui->frames[drvui->frame_no].map_lim[5])
	            goto eout;
	    }
	dxy += 4;
eout:;
	}
    }
    fprintf(drvui->flout, "Allocation size for 2D variables (bytes) is %d\n", (int)(3*dxy*sizeof(float)));

    if (dxy == 0) {
	Error_Box("Error: No part of the mapslice would be visible!\nChoose different parameters.");
	return;
    }


    if (doPOV || doVrml) {

	prhor = (float *)malloc(dxy*sizeof(float));
	if (!prhor) {
	    Error_Box("Error: Unable to allocate memory for map slice!");
	    return;
 	}
	prhog = (float *)malloc(dxy*sizeof(float));
 	if (!prhog) {
	    Error_Box("Error: Unable to allocate memory for map slice!");
	    free(prhor);
	    return;
 	}
	prhob = (float *)malloc(dxy*sizeof(float));
	if (!prhob) {
	    Error_Box("Error: Unable to allocate memory for map slice!");
	    free(prhor);
	    free(prhog);
	    return;
	}

	if (doPOV)
	    fprintf(drvui->fpoutp,"mesh2{\n  vertex_vectors { %d,\n",dxy);

	if (doVrml) {
	    if (!Vrml2) {
		fprintf (drvui->fpoutv, " Separator {\n");
		fprintf (drvui->fpoutv, "Coordinate3 { point [\n");
	    } else {
		fprintf (drvui->fpoutv, " Shape {\n");
		fprintf (drvui->fpoutv,
			"  geometry IndexedFaceSet { coord Coordinate{ point [\n");
	    }
	}
    }
    glBegin(GL_TRIANGLES);
    glNormal3f(drvui->frames[drvui->frame_no].mapnorm[0],
		drvui->frames[drvui->frame_no].mapnorm[1],
		drvui->frames[drvui->frame_no].mapnorm[2]);
	
    for (i = 0; i <= dx; i++) {
	XP[0] = min_x + (float)(i) * (max_x - min_x) / (float)(dx);
	for (j = 0; j <= dy; j++) {
	    XP[1] = min_y + (float)(j) * (max_y - min_y) / (float)(dy);
	    XP[2] = D;	/* Now have slice point in XP */
	    for (k = 0; k < 4; k++) {
		float XPP[3];
		XPP[2] = XP[2];
		XPP[1] = XP[1] + offset[k][1] * (max_y - min_y) / (float)(dx);
		XPP[0] = XP[0] + offset[k][0] * (max_x - min_x) / (float)(dy);
		Mul_Rinvv(R, XPP, X);	/* Convert XPP to X */
		Convert_Cart_Cryst(bmat_inv, X, x, origin); /* and X to x */
		if (x[0] < drvui->frames[drvui->frame_no].cryst_lim[0] ||
		    x[0] < drvui->frames[drvui->frame_no].map_lim[0] ||
		    x[0] > drvui->frames[drvui->frame_no].cryst_lim[3] ||
		    x[0] > drvui->frames[drvui->frame_no].map_lim[3])
			goto endy;
		if (x[1] < drvui->frames[drvui->frame_no].cryst_lim[1] ||
		    x[1] < drvui->frames[drvui->frame_no].map_lim[1] ||
		    x[1] > drvui->frames[drvui->frame_no].cryst_lim[4] ||
		    x[1] > drvui->frames[drvui->frame_no].map_lim[4])
			goto endy;
		if (x[2] < drvui->frames[drvui->frame_no].cryst_lim[2] ||
		    x[2] < drvui->frames[drvui->frame_no].map_lim[2] ||
		    x[2] > drvui->frames[drvui->frame_no].cryst_lim[5] ||
		    x[2] > drvui->frames[drvui->frame_no].map_lim[5])
			goto endy;
		for (l = 0; l < 3; l++) {
		    q[k][l] = x[l];
		    pp[k][l] = X[l];
		}
	    }
	    rho[0] = InterpolateMap(q[0][0], q[0][1], q[0][2]);
	    colorramp(rho[0], &rhor[0], &rhog[0], &rhob[0]);
	    rho[1] = InterpolateMap(q[1][0], q[1][1], q[1][2]);
	    colorramp(rho[1], &rhor[1], &rhog[1], &rhob[1]);
	    rho[2] = InterpolateMap(q[2][0], q[2][1], q[2][2]);
	    colorramp(rho[2], &rhor[2], &rhog[2], &rhob[2]);
	    rho[3] = InterpolateMap(q[3][0], q[3][1], q[3][2]);
	    colorramp(rho[3], &rhor[3], &rhog[3], &rhob[3]);
	    if (doPOV || doVrml) {
		prhor[vertn] = rhor[0];
		prhog[vertn] = rhog[0];
		prhob[vertn++] = rhob[0];
		prhor[vertn] = rhor[1];
		prhog[vertn] = rhog[1];
		prhob[vertn++] = rhob[1];
		prhor[vertn] = rhor[2];
		prhog[vertn] = rhog[2];
		prhob[vertn++] = rhob[2];
		prhor[vertn] = rhor[3];
		prhog[vertn] = rhog[3];
		prhob[vertn++] = rhob[3];
	    }

	    glColor3f(rhor[0], rhog[0], rhob[0]);
	    glVertex3f(pp[0][0], pp[0][1], pp[0][2]);
	    glColor3f(rhor[1], rhog[1], rhob[1]);
	    glVertex3f(pp[1][0],pp[1][1],pp[1][2]);
	    glColor3f(rhor[2], rhog[2], rhob[2]);
	    glVertex3f(pp[2][0], pp[2][1], pp[2][2]);

	    glColor3f(rhor[0], rhog[0], rhob[0]);
	    glVertex3f(pp[0][0], pp[0][1], pp[0][2]);
	    glColor3f(rhor[2], rhog[2], rhob[2]);
	    glVertex3f(pp[2][0], pp[2][1], pp[2][2]);
	    glColor3f(rhor[3], rhog[3], rhob[3]);
	    glVertex3f(pp[3][0], pp[3][1], pp[3][2]);

	    if (doPOV) {
		fprintf (drvui->fpoutp, "<%8.5f, %8.5f, %8.5f> <%8.5f,"
			 "%8.5f, %8.5f> <%8.5f, %8.5f, %8.5f> <%8.5f,%8.5f,%8.5f>\n", 
			 pp[0][0],pp[0][1],pp[0][2],pp[1][0],pp[1][1],pp[1][2],
			 pp[2][0],pp[2][1],pp[2][2],pp[3][0],pp[3][1],pp[3][2]);
	    }	
	    if (doAsy) {
		fprintf (drvui->fpouta, "draw (pic, surface ( (%8.5f, %8.5f, %8.5f)--(%8.5f,"
			 "%8.5f, %8.5f)--(%8.5f, %8.5f, %8.5f)--(%8.5f,%8.5f,%8.5f)--cycle),\n", 
			 pp[0][0],pp[0][1],pp[0][2],pp[1][0],pp[1][1],pp[1][2],
			 pp[2][0],pp[2][1],pp[2][2],pp[3][0],pp[3][1],pp[3][2]);
		fprintf(drvui->fpouta, "new pen[] {rgb(%.2f,%.2f,%.2f),rgb(%.2f,%.2f,%.2f),",
			rhor[0],rhog[0],rhob[0],rhor[1],rhog[1],rhob[1]);
		fprintf(drvui->fpouta, "rgb(%.2f,%.2f,%.2f),rgb(%.2f,%.2f,%.2f)});\n",
			rhor[2],rhog[2],rhob[2],rhor[3],rhog[3],rhob[3]);
	    }	
	    if (doVrml) {
		fprintf (drvui->fpoutv, "%8.5f  %8.5f  %8.5f,   %8.5f  %8.5f"
			"%8.5f,    %8.5f  %8.5f  %8.5f,    %8.5f  %8.5f  %8.5f", 
			pp[0][0],pp[0][1],pp[0][2],pp[1][0],pp[1][1],pp[1][2],
			pp[2][0],pp[2][1],pp[2][2],pp[3][0],pp[3][1],pp[3][2]);
		if (vertn < dxy-1) {
		    fprintf(drvui->fpoutv,",\n");
		} else {
		    fprintf(drvui->fpoutv,"]}\n");
		}
	    }
endy:;
	}
    }
    glEnd();

    if (doPOV) {
	fprintf(drvui->fpoutp," }\n texture_list { %d,\n ",dxy);
	for (k=0;k<dxy-1;k++)
	    fprintf(drvui->fpoutp," texture {pigment { color red %.2f green %.2f blue %.2f }},\n",
		    prhor[k],prhog[k],prhob[k]);
    	fprintf(drvui->fpoutp," texture {pigment { color red %.2f green %.2f blue %.2f }}\n }\n",
		prhor[dxy-1],prhog[dxy-1],prhob[dxy-1]);

	fprintf(drvui->fpoutp," face_indices { %d,\n ",dxy/2);
	int kk=0;
	int k=0;
	while (kk<dxy-6) {   // bracketed triple is coordinate, other texture
	    fprintf(drvui->fpoutp,"<%d,%d,%d>,%d,%d,%d,",kk,kk+1,kk+2,kk,kk+1,kk+2);
	    fprintf(drvui->fpoutp,"<%d,%d,%d>,%d,%d,%d,",kk,kk+2,kk+3,kk,kk+1,kk+2);
	    kk+=4;
	    k++;
	    if (k==2) {
		fprintf(drvui->fpoutp,"\n ");
		k=0;
	    }
	}
	fprintf(drvui->fpoutp,"<%d,%d,%d>,%d,%d,%d,<%d,%d,%d>,%d,%d,%d\n }\n}\n",
		dxy-4,dxy-3,dxy-2,dxy-4,dxy-3,dxy-2,dxy-4,dxy-2,dxy-1,dxy-4,dxy-2,dxy-1);
    }

    if (doVrml) {
	if (!Vrml2) {
	    fprintf (drvui->fpoutv, " MaterialBinding { value PER_VERTEX }\n  Material {  diffuseColor[ \n");
	    for (k=0;k<dxy-1;k++)
		fprintf(drvui->fpoutv," %.2f %.2f %.2f,\n",
			prhor[k],prhog[k],prhob[k]);
    	    fprintf(drvui->fpoutv," %.2f %.2f %.2f] }\n",
		    prhor[dxy-1],prhog[dxy-1],prhob[dxy-1]);
	    fprintf(drvui->fpoutv," IndexedFaceSet{");
	}
	fprintf(drvui->fpoutv," coordIndex[");
	int kk=0;
	while (kk<dxy-6) {
	    for (k=0;k<3;k++){
		fprintf(drvui->fpoutv,"%d,%d,%d,-1,%d,%d,%d,-1,",kk,kk+1,kk+2,kk,kk+2,kk+3);
		kk+=4;
	    }
	    fprintf(drvui->fpoutv,"\n ");
	}
	fprintf(drvui->fpoutv,"%d,%d,%d,-1,%d,%d,%d,-1]\n",
		dxy-4,dxy-3,dxy-2,dxy-4,dxy-2,dxy-1);
	if (Vrml2) {
	    fprintf (drvui->fpoutv, " color Color { color[\n");
	    for (k=0;k<dxy-1;k++)
		fprintf(drvui->fpoutv," %.2f %.2f %.2f,\n",
			prhor[k],prhog[k],prhob[k]);
    	    fprintf(drvui->fpoutv," %.2f %.2f %.2f]}\n",
		    prhor[dxy-1],prhog[dxy-1],prhob[dxy-1]);
	    fprintf(drvui->fpoutv," colorPerVertex TRUE\n");
	    fprintf(drvui->fpoutv," convex TRUE\n");
	    fprintf(drvui->fpoutv," solid FALSE\n");
	}
	fprintf(drvui->fpoutv," }\n }\n");
    }
    free (prhor);
    free (prhog);
    free (prhob);
    drvui->mainWindow->cursor (FL_CURSOR_DEFAULT);
}

void MapLegend (void) 
{
    char label[80]="";
    char *p;
    GLfloat fw=0.0;
    float glr, glg, glb, rho;
    float dv;
    int i;
    
    dv = Map_Info.rhomx - Map_Info.rhomn;  
    glMatrixMode (GL_PROJECTION);
    glPushMatrix ();
    glLoadIdentity ();
    glOrtho (0, 1, 0, 1, -1, 1);
    glMatrixMode (GL_MODELVIEW);
    glPushMatrix ();
    glLoadIdentity ();
    glDisable (GL_LIGHTING);
    glBegin (GL_LINES);
    for (i = 0; i < 256; i++) {
	rho = Map_Info.rhomn + 0.0039f * (float)i * dv;
	colorramp (rho, &glr, &glg, &glb);
	glColor3f (glr, glg, glb);
	glVertex3f (0.05f, 0.55f + .00156f * (float)i, 0.0f); 
	glVertex3f (0.1f, 0.55f + .00156f * (float)i, 0.0f); 
    }
    glEnd ();
    glColor3f (0., 0., 0.);
    glBegin (GL_LINES);
    glVertex3f (0.03f, 0.55f + 0.48f * (0.0f - Map_Info.rhomn) / dv, 0.0f);
    glVertex3f (0.12f, 0.55f + 0.48f * (0.0f - Map_Info.rhomn) / dv, 0.0f);
    glEnd();

    for (i = 0; i < 6; i++) {
	fw = (0.2f * (float)i *dv) + Map_Info.rhomn;
	sprintf (label, "% 5.3f", fw);
	glPushMatrix ();
	glTranslatef (0.1f, 0.55f + 0.08f * (float)i, 0.0f);
	glLineWidth (2.0);
	glScalef (drvui->label_scale * 0.0002f, drvui->label_scale * 0.0002f, drvui->label_scale * 0.0002f);
	for (p = label; *p; p++) 
	    glutStrokeCharacter (GLUT_STROKE_ROMAN, *p);
	glScalef (1.0f, 1.0f, 1.0f);
	glLineWidth (1.0f);
	glPopMatrix ();
    }
    glEnable (GL_LIGHTING);
    glMatrixMode (GL_PROJECTION);
    glPopMatrix ();
    glMatrixMode (GL_MODELVIEW);
    glPopMatrix ();
}
