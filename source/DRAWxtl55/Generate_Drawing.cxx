// $Id: Generate_Drawing.cxx 1117 2011-02-26 14:56:50Z martin $
//
// module Generate_Drawing.cxx - part of DRAWxtl V5.5
// Coded using the FLTK 1.1.6 widget set
//
//     Larry W. Finger, Martin Kroeker and Brian Toby
//
// This module contains the following routines:
//
//  sub_add_nvert - add vertices with subcell mods.
//  sub_add_nvert_nc - add vertices with subcell mods - no checking
//  do_drawing - routine to put together all the elements of the drawing
//  Generate_Drawing - create the drawing for the GUI version of DRAWxtl
//  generate_slab - routine to generate slab parameters from the input values
//  show_slab - draw outline if selected slab if slabmode == 2
//
#include "drawxtl.h"
#include <stdio.h>
#include "draw_ext.h"
#include <stdlib.h>
#include <string.h>
#include <FL/glut.H>
#include <FL/gl.h>
#include <FL/fl_ask.H>
#include <math.h>
#include "DRAWxtlViewUI.h"
#include "CrystalView.h"

#include "DRAWxtl_proto.h"

/* ************************************************************** */
/* ************************************************************** */

void
do_drawing (void)
// routine to put together the elements of the drawing
{
    int i;

    nvert = 0;			// empty vertex list
    if (!drvui->subsys_vol[0])
	drvui->subsys_vol[0] = drvui->subsys_ref_volume;
    if (drvui->frame_no == 1)
	draw_cell (docell);	// generate unit-cell corners and draw it
    nvert = 0;			// empty vertex list
    generate_spheres ();	// generate all spheres to be drawn within box
    if (drvui->do_ellipsoids && drvui->n_ellips > 1)
	generate_ellipsoids ();	// now the ellipsoids
    if (drvui->nmag > 0)
	generate_arrows ();	// magnetic moment arrows, if any
    generate_bonds ();		// generate all bonds to be drawn within box
    generate_cones ();		// generate all lonepairs to be drawn within box
    generate_poly ();		// generate all polyhedra to be drawn
    generate_planes ();		// now the planes
    generate_texts ();		// any textual elements
    generate_lsq_planes ();	// any least squares planes
    generate_aimsurf ();	// any surfaces
    generate_voids ();

    if (slabmode > 1)
	show_slab ();		// show outline of cutout section 
    if (ReadFourMap) {
	if (drvui->frames[drvui->frame_no].slice > 1) {
	    generate_slice();
	} else if (!drvui->Fourier2d) {	// 3d map
	    for (i = 1; i <= drvui->numOfFourierContours; ++i) {
// generate & display a Fourier map
		generate_map (drvui->fourier[i].FourierContourLevel,
			      drvui->fourier[i].FourierContourSolid,
			      drvui->fourier[i].FourierContourColor,
			      drvui->fourier[i].FourierBackColor);
	    }
	} else {		// generate 2d map contours
	    for (i = 1; i <= drvui->numOfFourierContours; ++i) {
		float level = drvui->fourier[i].FourierContourLevel;

		for (;;) {
		    char color[] = "";
		    generate_map (level, 0, drvui->fourier[i].FourierContourColor, color);
		    level += drvui->fourier[i].FourierContourStep;
		    if (level > drvui->fourier[i].FourierContourTop)
			break;
		}
	    }
	}
    }
}

/* ************************************************************** */
/* ************************************************************** */

void
sub_add_vert (float *vert, int no_cell)
{

/* transform from x,y,z in a',b',c' (vert) to x',y',z' in a,b,c (new_vert) */

    float new_vert[3];

    float F1 = drvui->subsys_fact[no_cell][0][0];

    float F2 = drvui->subsys_fact[no_cell][0][1];

    float F3 = drvui->subsys_fact[no_cell][0][2];

    float F4 = drvui->subsys_fact[no_cell][1][0];

    float F5 = drvui->subsys_fact[no_cell][1][1];

    float F6 = drvui->subsys_fact[no_cell][1][2];

    float F7 = drvui->subsys_fact[no_cell][2][0];

    float F8 = drvui->subsys_fact[no_cell][2][1];

    float F9 = drvui->subsys_fact[no_cell][2][2];


    new_vert[0] = (vert[0] * (F5 * F9 - F6 * F8)
		   + vert[1] * (F3 * F8 - F2 * F9)
		   + vert[2] * (F2 * F6 - F3 * F5)) * drvui->subsys_ref_volume
	/ drvui->subsys_vol[no_cell];
    new_vert[1] = (vert[0] * (F6 * F7 - F4 * F9)
		   + vert[1] * (F1 * F9 - F3 * F7)
		   + vert[2] * (F3 * F4 - F1 * F6)) * drvui->subsys_ref_volume
	/ drvui->subsys_vol[no_cell];
    new_vert[2] = (vert[0] * (F4 * F8 - F5 * F7)
		   + vert[1] * (F2 * F7 - F1 * F8)
		   + vert[2] * (F1 * F5 - F2 * F4)) * drvui->subsys_ref_volume
	/ drvui->subsys_vol[no_cell];

    add_vert (new_vert, 0, 1, 0, 0);
}

/* ************************************************************** */
/* ************************************************************** */

void
sub_add_vert_nc (float *vert, int no_cell)
{
    float new_vert[3];

    float F1 = drvui->subsys_fact[no_cell][0][0];

    float F2 = drvui->subsys_fact[no_cell][0][1];

    float F3 = drvui->subsys_fact[no_cell][0][2];

    float F4 = drvui->subsys_fact[no_cell][1][0];

    float F5 = drvui->subsys_fact[no_cell][1][1];

    float F6 = drvui->subsys_fact[no_cell][1][2];

    float F7 = drvui->subsys_fact[no_cell][2][0];

    float F8 = drvui->subsys_fact[no_cell][2][1];

    float F9 = drvui->subsys_fact[no_cell][2][2];

/* transform from x,y,z in a',b',c' (vert) to x',y',z' in a,b,c (new_vert) */

    new_vert[0] = (vert[0] * (F5 * F9 - F6 * F8)
		   + vert[1] * (F3 * F8 - F2 * F9)
		   + vert[2] * (F2 * F6 - F3 * F5)) * drvui->subsys_ref_volume
	/ drvui->subsys_vol[no_cell];
    new_vert[1] = (vert[0] * (F6 * F7 - F4 * F9)
		   + vert[1] * (F1 * F9 - F3 * F7)
		   + vert[2] * (F3 * F4 - F1 * F6)) * drvui->subsys_ref_volume
	/ drvui->subsys_vol[no_cell];
    new_vert[2] = (vert[0] * (F4 * F8 - F5 * F7)
		   + vert[1] * (F2 * F7 - F1 * F8)
		   + vert[2] * (F1 * F5 - F2 * F4)) * drvui->subsys_ref_volume
	/ drvui->subsys_vol[no_cell];

    add_vert_nc (new_vert);
}

/* ************************************************************** */
/* ************************************************************** */

void
Generate_Drawing (int keepmatrix)
{
    char filename[256] = "";

    int i, m;

    float temp;

    float rm[16];

    char col_bg_v[40];

    char col_bg_p[40];

    static int in_progress;

    if (strlen (drvui->Cur_File) == 0) {
	return;
    }

    if (in_progress) {
#ifndef WIN32
	fprintf (stderr, "Generate_Drawing: calculation already in progress\n");
#endif
	return;
    }
    in_progress = 1;
    if (drvui->crystalDL) {
	glDeleteLists (drvui->crystalDL, 1);	// flush display list
	drvui->crystalDL = glGenLists (1);	// create a new one
    }
    if (!(drvui->fpin = fopen (drvui->Cur_Temp, "r"))) {
	Error_Box ("Cannot open structure file");
	goto out;
    }
    strcpy (filename, drvui->Cur_Root);
    strcat (filename, ".out");
    if (!(drvui->flout = fopen (filename, "w"))) {
	Error_Box ("Cannot open listing file");
	goto out;
    }
    if (doVrml) {
	strcpy (filename, drvui->Cur_Root);
	if (X3D)
	    strcat (filename, ".x3dv");
	else
	    strcat (filename, ".wrl");
	if (!(drvui->fpoutv = fopen (filename, "w"))) {
	    Error_Box ("Cannot open VRML output file\n");
	    goto out;
	}
    }
    if (doPOV) {
	if (!drvui->automation) {
	    strcpy (filename, drvui->Cur_Root);
	    strcat (filename, ".pov");
	} else
	    strcpy (filename, drvui->automate_name);
	if (!(drvui->fpoutp = fopen (filename, "w"))) {
	    Error_Box ("Cannot open POV output file\n");
	    goto out;
	}
    }
    if (doAsy) {
	strcpy (filename, drvui->Cur_Root);
	strcat (filename, ".asy");
	if (!(drvui->fpouta = fopen (filename, "w"))) {
	    Error_Box ("Cannot open Asymptote output file\n");
	    goto out;
	}
    }
    strcpy (filename, drvui->Cur_Root);
    strcat (filename, ".cns");
    if (!(drvui->fcns = fopen (filename, "w"))) {
	Error_Box ("Cannot open console listing file\n");
	goto out;
    }

    if (drvui->nsurf > 1) {
	for (i = 1; i < drvui->nsurf; i++) {
	    free (drvui->surfx[i]);
	    free (drvui->surfy[i]);
	    free (drvui->surfz[i]);
	}
    }

    if (drvui->atom_no) {
	free (drvui->atom_no);
	drvui->atom_no = NULL;
    }
    if (!
	(drvui->atom_no =
	 (int *) zalloc ((unsigned) ((long) (drvui->verts_alloc * sizeof (int)))))) {
	Error_Box ("Unable to get initial atom_no allocation\n");
	goto out;
    }

    if (drvui->atom_so) {
	free (drvui->atom_so);
	drvui->atom_so = NULL;
    }
    if (!
	(drvui->atom_so =
	 (int *) zalloc ((unsigned) ((long) (drvui->verts_alloc * sizeof (int)))))) {
	Error_Box ("Unable to get initial atom_so allocation\n");
	goto out;
    }

    if (xypos) {
	free (xypos);
    }
    if (!
	(xypos =
	 (float *)
	 zalloc ((unsigned) (3 * (long) (drvui->verts_alloc * sizeof (float)))))) {
	Error_Box ("Unable to get initial xypos allocation\n");
	goto out;
    }
    if (xypos_nm) {
	free (xypos_nm);
    }
    if (drvui->voidflag == 1 && drvui->voidmap) {
	for ( i = 0; i < drvui->voidgrid[0]; i++) {
	    for (int j = 0; j < drvui->voidgrid[1]; j++) 
		free (drvui->voidmap[i][j]);
	    free (drvui->voidmap[i]);
	}
	free (drvui->voidmap);
	drvui->voidmap = NULL;
    }
    if (!
	(xypos_nm =
	 (float *)
	 zalloc ((unsigned) (3 * (long) (drvui->verts_alloc * sizeof (float)))))) {
	Error_Box ("Unable to get initial xypos_nm allocation\n");
	goto out;
    }
    fprintf (drvui->flout, "\n\n                DRAWxtl\n");
    fprintf (drvui->flout, "  Crystal Structure Drawing Program %s\n\n", VERSION);
    fprintf (drvui->flout,
	     "  Copyright 1996, 1997, 2000-2011 by Larry Finger, Martin Kroeker, and Brian Toby\n");
    fprintf (drvui->flout, "  (all rights reserved).\n\n");
    fprintf (drvui->flout, "     Combined POVxtl and VRMLxtl\n");
    fprintf (drvui->flout,
	     "POVxtl  - Scenes for Persistence of Vision Ray Tracer (POV-Ray) V2.0\n");
    fprintf (drvui->flout, "VRMLxtl - Virtual Reality Modeling Language V1.0/V2.0\n");
    fprintf (drvui->fcns, "\n\n                DRAWxtl\n");
    fprintf (drvui->fcns, "  Crystal Structure Drawing Program %s\n\n", VERSION);
    fprintf (drvui->fcns,
	     "  Copyright 1996, 1997, 2000-2011 by Larry Finger, Martin Kroeker, and Brian Toby\n");
    fprintf (drvui->fcns, "  (all rights reserved).\n\n");
    fprintf (drvui->fcns, "     Combined POVxtl and VRMLxtl\n");
    fprintf (drvui->fcns,
	     "POVxtl  - Scenes for Persistence of Vision Ray Tracer (POV-Ray) V2.0\n");
    fprintf (drvui->fcns, "VRMLxtl - Virtual Reality Modeling Language V1.0/V2.0\n");
    drvui->frame_no = 0;
    drvui->nedges = 0;
    drvui->nsphere = drvui->npoly = drvui->nbond = drvui->nplane = drvui->ncone = 1;
    drvui->nbplane = drvui->nsurf = drvui->natprop = 1;
    if (drvui->modulated <= 0) {
	for (i = 0; i < 3; i++) {	/* initialize first element of subsystem array */
	    for (m = 0; m < 3; m++)
		drvui->subsys_fact[0][i][m] = 0.0f;
	    drvui->subsys_fact[0][i][i] = 1.0f;
	}
    }
// Initialize POV output limits

    if (!drvui->automation) {
	POV_Max[0] = POV_Max[1] = POV_Max[2] = -99999.0f;
	POV_Min[0] = POV_Min[1] = POV_Min[2] = 99999.0f;
    }
    natom = 0;
    NvertM = 1;			// make master list empty
    drvui->nmag = 0;
    drvui->nlabel = 1;

// start loop through frames

    while (drvui->frame_no < drvui->max_frame) {
	drvui->frame_no++;
	NvertM = 1;
	fprintf (drvui->fcns, "\n doing frame %d of %d\n", drvui->frame_no,
		 drvui->max_frame);
	fprintf (drvui->flout, "\n doing frame %d of %d\n", drvui->frame_no,
		 drvui->max_frame);
	clipflag = 0;
	packflag = 0;
	domolcomp = 0;
	drvui->mol_d = 0.0f;
	drvui->mag_matrix[0][0] = 1.0f;
	drvui->mag_matrix[1][0] = 0.0f;
	drvui->mag_matrix[2][0] = 0.0f;
	drvui->mag_matrix[2][1] = 0.0f;
	drvui->mag_matrix[1][2] = 0.0f;
	drvui->mag_matrix[0][1] = 0.0f;
	drvui->mag_matrix[0][2] = 0.0f;
	drvui->mag_matrix[1][1] = 1.0f;
	drvui->mag_matrix[2][2] = 1.0f;
	origin[0] = drvui->X_Origin;
	origin[1] = drvui->Y_Origin;
	origin[2] = drvui->Z_Origin;
	boxlim[0] = drvui->X_Boxlim;
	boxlim[1] = drvui->Y_Boxlim;
	boxlim[2] = drvui->Z_Boxlim;
	docell = Unit_Cell;
	drvui->do_ellipsoids = Color_Warning = 0;
	rad_cell = 0.02f;
	drvui->Sphere_Mult = drvui->SpMult;
	drvui->Bond_Mult = drvui->BndMult;
	slabmode = 0;
	drvui->n_ellips = 1;
	if (drvui->frame_no != 1) {
	    drvui->fpin = fopen (drvui->Cur_Temp, "r");
	}
	get_input (0);

// Calculate Grand rotation matrix

	if (keepmatrix == 0) {
	    Rotq = XYZ_Rot_to_Q (xrot, yrot, zrot);
	    crystal->calculate (rm);	// get rotation matrix
	    for (m = 0; m < 3; m++) {	//  and copy to G_Rot
		G_Rot[0][m] = rm[m];
		G_Rot[1][m] = rm[m + 4];
		G_Rot[2][m] = rm[m + 8];
	    }
	}
// generate all possible atoms to get the size of the drawing

	if (!s_vert) {
	    if ((s_vert = (float *) zalloc ((unsigned) (3 * (long) (drvui->verts_alloc
								    *
								    sizeof (float))))) ==
		NULL) {
		Error_Box ("Unable to get initial s_vert allocation");
		goto out;
	    }
	}
	if (drvui->vert_occ)
	    free (drvui->vert_occ);
	if (!
	    (drvui->vert_occ =
	     (float *) zalloc ((unsigned) (drvui->verts_alloc * sizeof (float))))) {
	    Error_Box ("Unable to get site occupancy allocation");
	    goto out;
	}
	get_atom_id ();		// rework atom id numbers to pointers
	if (slabmode)
	    generate_slab ();	// initialize cutout box if needed
	if (drvui->modulated > 0) {
	    float save_limits[6];

	    int j, atom;

	    fprintf (drvui->flout,
		     "\nTransformed (modulated) contents of 2x2x2 unit-cells corresponding to\n"
		     " the 'simulated' P1 positions of JANA 2000\n\n"
		     "  Atom       x'       y'       z'   site occ.   Sym. Op\n\n");
	    for (i = 0; i < 3; i++) {
		save_limits[i] = drvui->frames[drvui->frame_no].cryst_lim[i];
		drvui->frames[drvui->frame_no].cryst_lim[i] = -1.0f;	/* set limits from -1 to +1 in all directions */
		save_limits[i + 3] = drvui->frames[drvui->frame_no].cryst_lim[i + 3];
		drvui->frames[drvui->frame_no].cryst_lim[3 + i] = 1.0f;
	    }
	    build_box_contents ();	//  get all atoms in the 2x2x2 region
	    for (i = 0; i < 6; i++) {
		drvui->frames[drvui->frame_no].cryst_lim[i] = save_limits[i];	/* restore original limits */
	    }
	    for (atom = 0; atom < natom; atom++) {
		for (j = 1; j < NvertM; j++) {
		    int sym, k, m;

		    if (atom != drvui->atom_no[j] / 100)
			continue;
		    sym = drvui->atom_no[j] % 100;
		    k = 0;
		    for (m = 0; m < 3; m++) {
			if (xypos[3 * j + m] < -1.0f)
			    k = 1;
			if (xypos[3 * j + m] > 1.0f)
			    k = 1;
		    }
		    if (!k)
			fprintf (drvui->flout,
				 " %c%c%c%c%3d %8.5f %8.5f %8.5f %8.5f %7d\n",
				 drvui->atoms[atom].atom_l[0],
				 drvui->atoms[atom].atom_l[1],
				 drvui->atoms[atom].atom_l[2],
				 drvui->atoms[atom].atom_l[3],
				 drvui->atoms[atom].sv_atom_n, xypos[3 * j],
				 xypos[3 * j + 1], xypos[3 * j + 2], drvui->vert_occ[j],
				 sym + 1);
		}
	    }
	    NvertM = 1;		// make master list empty and reset limits
	    if (!drvui->automation) {
	        POV_Max[0] = POV_Max[1] = POV_Max[2] = -99999.0f;
	        POV_Min[0] = POV_Min[1] = POV_Min[2] = 99999.0f;
	    }
	}
	build_box_contents ();	// complete the master atom list
	fprintf (drvui->fcns, "\nMaster atom list contains %5d atoms.\n", NvertM - 1);
	fprintf (drvui->flout, "\nMaster atom list contains %5d atoms.\n", NvertM - 1);

// allocate space for the various arrays (at least 14 cell vertices plus first atom - thus extra 20)

	if (o_vert)
	    free (o_vert);
	if (!
	    (o_vert =
	     (float *) zalloc ((unsigned) (6 * (NvertM + 20) * sizeof (float))))) {
	    Error_Box ("Unable to get o_vert allocation");
	    goto out;
	}
	if (o_vert_nm)
	    free (o_vert_nm);
	if (!
	    (o_vert_nm =
	     (float *) zalloc ((unsigned) (6 * (NvertM + 20) * sizeof (float))))) {
	    Error_Box ("Unable to get o_vert_nm allocation");
	    goto out;
	}
	if (vert_sym_no)
	    free (vert_sym_no);
	if (!(vert_sym_no = (int *) zalloc ((unsigned) (2 * NvertM * sizeof (int))))) {
	    Error_Box ("Unable to get vert_sym_no allocation");
	    goto out;
	}
	if (vert_sym_nos)
	    free (vert_sym_nos);
	if (!(vert_sym_nos = (int *) zalloc ((unsigned) (2 * NvertM * sizeof (int))))) {
	    Error_Box ("Unable to get vert_sym_nos allocation");
	    goto out;
	}
	if (drvui->orig_atom_no)
	    free (drvui->orig_atom_no);
	if (!
	    (drvui->orig_atom_no =
	     (int *) zalloc ((unsigned) (2 * NvertM * sizeof (int))))) {
	    Error_Box ("Unable to get orig_atom_no allocation");
	    goto out;
	}
	if (!(poly_list = (int *) zalloc ((unsigned) (4 * (NvertM + 20)) * sizeof (int)))) {
	    Error_Box ("Unable to get poly_list allocation");
	    goto out;
	}
	if (!(vertex_list = (int *) zalloc ((unsigned) (4 * NvertM) * sizeof (int)))) {
	    Error_Box ("Unable to get vertex_list allocation");
	    goto out;
	}
	analyze_bonds ();	// print the bond distance table
	fprintf (drvui->fcns, "\nThe 'pov' and 'wrl' files contain:\n");
	if (docell != 0 || natom == 0) {	// if unit-cell corners drawn, add them for scaling
	    float vert[3];

	    int saved_Nvert;

	    int nocell;

	    for (nocell = 0; nocell < drvui->no_subsys; nocell++) {
		saved_Nvert = NvertM;
		for (i = 0; i <= 2; ++i)
		    vert[i] = 0.0f;
		sub_add_vert (vert, nocell);	// place 0,0,0 in list
		vert[2] = 1.0f;
		sub_add_vert (vert, nocell);	// place 0,0,1 in list
		vert[1] = 1.0f;
		sub_add_vert (vert, nocell);	// place 0,1,1 in list
		vert[2] = 0.0f;
		sub_add_vert (vert, nocell);	// place 0,1,0 in list
		vert[0] = 1.0f;
		sub_add_vert (vert, nocell);	// place 1,1,0 in list
		vert[1] = 0.0f;
		sub_add_vert (vert, nocell);	// place 1,0,0 in list
		vert[2] = 1.0f;
		sub_add_vert (vert, nocell);	// place 1,0,1 in list
		vert[1] = 1.0f;
		sub_add_vert (vert, nocell);	// place 1,1,1 in list
		NvertM = saved_Nvert;	// reset list count 
	    }
	}
	strcpy (col_bg_p, drvui->col_bg);
	strcpy (col_bg_v, drvui->col_bg);
	Transform_VRML_Color (col_bg_v);
	Transform_POV_Color (col_bg_p);
	temp = POV_Max[0] - POV_Min[0];
	if (POV_Max[1] - POV_Min[1] > temp)
	    temp = POV_Max[1] - POV_Min[1];
	if (POV_Max[2] - POV_Min[2] > temp)
	    temp = POV_Max[2] - POV_Min[2];
	if (drvui->voidflag != 0)
	    calculate_voids ();
	if (drvui->frame_no == 1) {	// only do for first frame
	    Scale = 15.0f * temp / Magnification;
	    Text_Size = 0.0005f * Scale;
	    if (doVrml) {
		if (Vrml2) {
		    if (X3D)
			fprintf (drvui->fpoutv,"#X3D V3.0 utf8\n\nPROFILE Immersive\n\n");
		    else
			fprintf (drvui->fpoutv, "#VRML V2.0 utf8\n\n");
		} else {
		    fprintf (drvui->fpoutv, "#VRML V1.0 ascii\n\n");
		}
		fprintf (drvui->fpoutv, "# Scene file created by DRAWxtl %s from %s\n",
			 VERSION, drvui->Cur_File);
		fprintf (drvui->fpoutv,
			 "# with command line options: -b %6.2f %6.2f %6.2f -o %5.2f %5.2f %5.2f\n",
			 boxlim[0], boxlim[1], boxlim[2], origin[0], origin[1],
			 origin[2]);
		fprintf (drvui->fpoutv,
			 "#  -p %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f -v %5.2f %5.2f %5.2f\n",
			 drvui->frames[drvui->frame_no].cryst_lim[0],
			 drvui->frames[drvui->frame_no].cryst_lim[3],
			 drvui->frames[drvui->frame_no].cryst_lim[1],
			 drvui->frames[drvui->frame_no].cryst_lim[4],
			 drvui->frames[drvui->frame_no].cryst_lim[2],
			 drvui->frames[drvui->frame_no].cryst_lim[5], xrot, yrot, zrot);

// Output VRML header

		if (!Vrml2) {
		    fprintf (drvui->fpoutv,
			     "DEF BackgroundColor Info { string \"%s\" }\n", col_bg_v);
		    fprintf (drvui->fpoutv, "Separator {\n");
		    fprintf (drvui->fpoutv, " Material {\n");
		    fprintf (drvui->fpoutv, "  ambientColor    0.1 0.1 0.1\n");
		    fprintf (drvui->fpoutv, "  diffuseColor    1 1 1\n");
		    fprintf (drvui->fpoutv, "  specularColor   0.9 0.9 0.9\n");
		    fprintf (drvui->fpoutv, "  shininess       0.8\n");
		    fprintf (drvui->fpoutv, " }\n");
		}

		if (!M_cameras && Vrml2) {	// Generate Orthographic Projection for VRML97
		    fprintf (drvui->fpoutv, "NavigationInfo {type \"EXAMINE\"}\n");
		    fprintf (drvui->fpoutv, "Background { skyColor %s }\n", col_bg_v);
		    fprintf (drvui->fpoutv,
			     "Viewpoint { fieldOfView 0.1 position 0. 0. %8.3f }\n",
			     Scale);
		    fprintf (drvui->fpoutv, "DirectionalLight {direction 0 0 -1}\n");
		    fprintf (drvui->fpoutv, "DirectionalLight {direction 0 0 1}\n");
		} else {
		    if (!Vrml2) {	// simulate orthographic projection for VRML1.0
			float fscale, fangle;

			if (!M_cameras) {
			    fscale = Scale * (10.0f / 7.0f);
			    fangle = 0.07853f;
			} else {
			    fscale = Scale / 7.0f;
			    fangle = 0.7853f;
			}
			fprintf (drvui->fpoutv, "PerspectiveCamera {\n");
			fprintf (drvui->fpoutv, "  position        0 0 %8.3f\n", fscale);
			fprintf (drvui->fpoutv, "  orientation     0 0 1 0\n");
			fprintf (drvui->fpoutv, "  focalDistance %8.3f\n", fscale);
			fprintf (drvui->fpoutv, "  heightAngle   %8.5f\n", fangle);
			fprintf (drvui->fpoutv, "}\n");
			fprintf (drvui->fpoutv,
				 "DirectionalLight{ direction 0 0 -1 intensity 1.0}\n");
			fprintf (drvui->fpoutv,
				 "DirectionalLight{ direction 0 0 1 intensity 1.0}\n");
		    } else {
			fprintf (drvui->fpoutv, "NavigationInfo {type \"EXAMINE\"}\n");
			fprintf (drvui->fpoutv, "Background { skyColor %s }\n", col_bg_v);
			fprintf (drvui->fpoutv, "Viewpoint { position 0. 0. %8.3f }\n",
				 Scale / 7.);
			fprintf (drvui->fpoutv, "DirectionalLight {direction 0 0 -1}\n");
			fprintf (drvui->fpoutv, "DirectionalLight {direction 0 0 1}\n");
		    }
		}

		if (Vrml2) {
		    fprintf (drvui->fpoutv, "Transform { rotation 0 0 1 %f\n",
			     zrot / RAD);
		} else {
		    fprintf (drvui->fpoutv, "Rotation {");
		    fprintf (drvui->fpoutv, " rotation 0 0 1 %f}\n", zrot / RAD);
		}
		if (Vrml2) {
		    fprintf (drvui->fpoutv, " children Transform {");
		    fprintf (drvui->fpoutv, " rotation 0 1 0 %f\n", yrot / RAD);
		} else {
		    fprintf (drvui->fpoutv, " Rotation {");
		    fprintf (drvui->fpoutv, " rotation 0 1 0 %f}\n", yrot / RAD);
		}
		if (Vrml2) {
		    fprintf (drvui->fpoutv, " children Transform {");
		    fprintf (drvui->fpoutv, " rotation 1 0 0 %f\n", xrot / RAD);
		} else {
		    fprintf (drvui->fpoutv, " Rotation {");
		    fprintf (drvui->fpoutv, " rotation 1 0 0 %f}\n", xrot / RAD);
		}

		if (Vrml2)
		    fprintf (drvui->fpoutv, "  children [\n");
	    }

	    Scale *= 0.9f;
	    float cpx = (POV_Max[0] + POV_Min[0]) / 2.0f;

	    float cpy = (POV_Max[1] + POV_Min[1]) / 2.0f;

// Output POV header
	    if (doPOV) {
		fprintf (drvui->fpoutp,
			 "/*\n\nPOVxtl-Persistence of Vision Ray Tracer\n");
		fprintf (drvui->fpoutp, "Crystal Structure Drawing Program %s\n\n",
			 VERSION);
		fprintf (drvui->fpoutp,
			 "Copyright 1996, 1997, 2000-2011 by Larry Finger, Martin Kroeker, and Brian Toby\n");
		fprintf (drvui->fpoutp, "   (all rights reserved).\n\n");
		fprintf (drvui->fpoutp, "Scene file created from %s \n", drvui->Cur_File);
		fprintf (drvui->fpoutp,
			 "with command line options: -b %6.2f %6.2f %6.2f -o %5.2f %5.2f %5.2f \n",
			 boxlim[0], boxlim[1], boxlim[2], origin[0], origin[1],
			 origin[2]);
		fprintf (drvui->fpoutp,
			 "-p %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f -v %5.2f %5.2f %5.2f -m %5.2f\n*/\n\n",
			 drvui->frames[drvui->frame_no].cryst_lim[0],
			 drvui->frames[drvui->frame_no].cryst_lim[3],
			 drvui->frames[drvui->frame_no].cryst_lim[1],
			 drvui->frames[drvui->frame_no].cryst_lim[4],
			 drvui->frames[drvui->frame_no].cryst_lim[2],
			 drvui->frames[drvui->frame_no].cryst_lim[5], xrot, yrot, zrot,
			 Magnification);

		fprintf (drvui->fpoutp, "#if (version>3.65)\n#version 3.6;\n#end\n\n");
		if (drvui->Stereo == 1)
		    fprintf (drvui->fpoutp, "#version unofficial stereopov 0.2;\n");
		fprintf (drvui->fpoutp, "#include \"colors.inc\"\n");
		fprintf (drvui->fpoutp, "#include \"chars.inc\"\n");

		if (drvui->Stereo == 2) {
		    fprintf (drvui->fpoutp, "#include \"transforms.inc\"\n");
		    fprintf (drvui->fpoutp, "#include \"math.inc\"\n\n");
		    fprintf (drvui->fpoutp, "/* Stereo camera macro adapted from Jaime Vives Piqueres, www.ignorancia.org */\n");
		    if (!M_cameras) {
			fprintf (drvui->fpoutp, "#macro meshcam_orthographic(res_x, res_y, c_up_length)\n");
			fprintf (drvui->fpoutp, "#local c_direction=1;\n");
			fprintf (drvui->fpoutp, "#local px_inc=c_up_length*(res_x/res_y)/res_x;\n");
			fprintf (drvui->fpoutp, "#local centering=<-res_x*.5*px_inc+px_inc*.5,res_y*.5*px_inc-px_inc-px_inc*2/3,0>;\n");
			fprintf (drvui->fpoutp, "#local tr_base_1=< 0, px_inc*2/3, 0>;\n");
			fprintf (drvui->fpoutp, "#local tr_base_2=< px_inc/2, -px_inc/3, 0>;\n");
			fprintf (drvui->fpoutp, "#local tr_base_3=<-px_inc/2, -px_inc/3, 0>;\n");
			fprintf (drvui->fpoutp, "mesh{\n#local row_count=0;\n#while (row_count<res_y)\n");
			fprintf (drvui->fpoutp, " #local col_count=0;\n #local d_y=row_count*px_inc;\n #while (col_count<res_x)\n");
			fprintf (drvui->fpoutp, "  #local d_x=col_count*px_inc;\n");
			fprintf (drvui->fpoutp, "  #local tr_vertex_1=tr_base_1+<d_x,-d_y,c_direction>+centering;\n");
			fprintf (drvui->fpoutp, "  #local tr_vertex_2=tr_base_2+<d_x,-d_y,c_direction>+centering;\n");
			fprintf (drvui->fpoutp, "  #local tr_vertex_3=tr_base_3+<d_x,-d_y,c_direction>+centering;\n");
		    } else {
			fprintf (drvui->fpoutp, "#macro meshcam(res_x, res_y, c_angle)\n");
			fprintf (drvui->fpoutp, "#local c_up=1*y;\n");
			fprintf (drvui->fpoutp, "#local c_right=(res_x/res_y)*x;\n");
			fprintf (drvui->fpoutp, "#local c_direction=0.5*vlength(c_right)/tan(radians(c_angle)/2);\n");
			fprintf (drvui->fpoutp, "#local sph_lens=sphere{0,c_direction inverse}\n");
			fprintf (drvui->fpoutp, "#local px_inc=(res_x/res_y)/res_x;\n");
			fprintf (drvui->fpoutp, "#local centering=<-res_x*.5*px_inc+px_inc*.5,res_y*.5*px_inc-px_inc-px_inc*2/3,0>;\n");
			fprintf (drvui->fpoutp, "#local tr_base_1=< 0, px_inc*2/3, 0>;\n");
			fprintf (drvui->fpoutp, "#local tr_base_2=< px_inc/2, -px_inc/3, 0>;\n");
			fprintf (drvui->fpoutp, "#local tr_base_3=<-px_inc/2, -px_inc/3, 0>;\n");
	        	fprintf (drvui->fpoutp, "mesh{\n#local row_count=0;\n#while (row_count<res_y)\n");
			fprintf (drvui->fpoutp, " #local col_count=0;\n #local d_y=row_count*px_inc;\n #while (col_count<res_x)\n");
			fprintf (drvui->fpoutp, "  #local d_x=col_count*px_inc;\n  #local Norm=<0,0,0>;\n");
			fprintf (drvui->fpoutp, "  #local tr_vertex_1=trace(sph_lens,0,tr_base_1+<d_x,-d_y,c_direction>+centering, Norm);\n");
			fprintf (drvui->fpoutp, "  #local tr_vertex_2=trace(sph_lens,0,tr_base_2+<d_x,-d_y,c_direction>+centering, Norm);\n");
			fprintf (drvui->fpoutp, "  #local tr_vertex_3=trace(sph_lens,0,tr_base_3+<d_x,-d_y,c_direction>+centering, Norm);\n");
		    }
		    fprintf (drvui->fpoutp, "  triangle{tr_vertex_1,tr_vertex_2,tr_vertex_3}\n");
		    fprintf (drvui->fpoutp, "  #local col_count=col_count+1;\n  #end\n #local row_count=row_count+1;\n");
		    fprintf (drvui->fpoutp, "#end\n}\n#end\n");
		    fprintf (drvui->fpoutp, "#macro meshcam_placement(c_lo, c_la)\n transform{\n");
		    fprintf (drvui->fpoutp, "  transform{Reorient_Trans(z,<c_la.x,0,c_la.z>-<c_lo.x,0,c_lo.z>)}\n");
		    fprintf (drvui->fpoutp, "  transform{Reorient_Trans(<c_la.x,0,c_la.z>-<c_lo.x,0,c_lo.z>,c_la-c_lo)}\n");
		    fprintf (drvui->fpoutp, "  translate c_lo\n }\n#end\n");
		}

		if (drvui->ambient > 0.) {
		    fprintf (drvui->fpoutp, "\ndefault {\n");
		    fprintf (drvui->fpoutp, "        finish {\n");
		    fprintf (drvui->fpoutp, "                ambient %5.2f\n",
			     drvui->ambient);
		    fprintf (drvui->fpoutp, "                diffuse %5.2f\n",
			     drvui->diffuse);
		    fprintf (drvui->fpoutp, "                specular %5.2f\n",
			     drvui->specular);
		    fprintf (drvui->fpoutp, "                roughness %5.2f\n",
			     drvui->roughness);
		    fprintf (drvui->fpoutp, "               }\n");
		    fprintf (drvui->fpoutp, "        }\n\n");
		}

		if (drvui->Stereo <2) {
		    if (!M_cameras) {
			fprintf (drvui->fpoutp, "camera{ orthographic\n");
		    } else {
			fprintf (drvui->fpoutp, "camera{\n");
		    }
		    fprintf (drvui->fpoutp, " location< %f ,%f, %f>\n", cpx, cpy,
			    1.6f * Scale);
		    if (!M_cameras) {
			fprintf (drvui->fpoutp, " up <0,%f,0>\n", 0.1 * Scale);
			fprintf (drvui->fpoutp, " right <%f,0,0>\n", -0.1 * Scale);
		    } else {
			fprintf (drvui->fpoutp, " up <0, 0.1, 0>\n");
			fprintf (drvui->fpoutp, " right <-0.1, 0, 0>\n");
		    }
		    fprintf (drvui->fpoutp, " look_at <%f,%f,0>\n", cpx, cpy);
		    if (drvui->Stereo == 1) {
			fprintf (drvui->fpoutp, " stereo_base %f\n", drvui->stereo_base);
			if (drvui->cross_eyed)
			    fprintf (drvui->fpoutp, " cross_eyed\n");
		    }
		    fprintf (drvui->fpoutp, "}\n");
		} else {
		    if (!M_cameras)
			fprintf (drvui->fpoutp, "#declare c_location=<%f,%f,%f>;\n",cpx,cpy,-1.6*Scale);
		    else 
			fprintf (drvui->fpoutp, "#declare c_location=<%f,%f,%f>;\n",cpx,cpy,-0.1*Scale);
		    fprintf (drvui->fpoutp, "#declare c_look_at=<%f,%f,0>;\n", cpx,cpy);
		    fprintf (drvui->fpoutp, "#declare c_stereo_base=%f;\n", drvui->stereo_base);
		    if (!M_cameras) {
			fprintf (drvui->fpoutp, "#declare c_up_length=%f;\n", 0.1*Scale);
			fprintf (drvui->fpoutp, "#declare camera_mesh=meshcam_orthographic(image_width*.5, image_height, c_up_length)\n");
		    } else {
			fprintf (drvui->fpoutp, "#declare c_angle=54;\n");
			fprintf (drvui->fpoutp, "#declare camera_mesh=meshcam(image_width*.5, image_height, c_angle)\n");
		    }
		    fprintf (drvui->fpoutp, "#declare c_location_left =c_location+c_stereo_base;\n");
		    fprintf (drvui->fpoutp, "#declare c_location_right=c_location-c_stereo_base;\n");
		    fprintf (drvui->fpoutp, "camera{\n");
		    fprintf (drvui->fpoutp, " mesh_camera{ 1 2\n");
		    if (drvui->cross_eyed) {
			fprintf (drvui->fpoutp, "  mesh{camera_mesh\n   meshcam_placement(c_location_right,c_look_at)\n  }\n");
			fprintf (drvui->fpoutp, "  mesh{camera_mesh\n   meshcam_placement(c_location_left,c_look_at)\n  }\n");
		    } else {
			fprintf (drvui->fpoutp, "  mesh{camera_mesh\n   meshcam_placement(c_location_left,c_look_at)\n  }\n");
			fprintf (drvui->fpoutp, "  mesh{camera_mesh\n   meshcam_placement(c_location_right,c_look_at)\n  }\n");
		    }
		    fprintf (drvui->fpoutp, " }\n location <0,0,1.01>\n}\n\n");
		} 

		fprintf (drvui->fpoutp, "background{color %s}\n", col_bg_p);
		if (drvui->Stereo == 2)
		    fprintf (drvui->fpoutp, "light_source {< 0, 0, %f>\n", -2.*Scale);
		else
		    fprintf (drvui->fpoutp, "light_source {< 0, 0, %f>\n", Scale);
		fprintf (drvui->fpoutp, "color red 4.0 green 4.0 blue 4.0\n}\n");
		fprintf (drvui->fpoutp, "object{\n union{\n");
		fprintf (drvui->fpoutp, "\n");
		fprintf (drvui->fpoutp, "#declare xrot = %5.1f;\n", xrot);
		fprintf (drvui->fpoutp, "#declare yrot = %5.1f;\n", yrot);
		fprintf (drvui->fpoutp, "#declare zrot = %5.1f;\n", zrot);
		fprintf (drvui->fpoutp, "#declare xmove = %5.1f;\n", drvui->Trans[0]);
		fprintf (drvui->fpoutp, "#declare ymove = %5.1f;\n", drvui->Trans[1]);
		fprintf (drvui->fpoutp, "\n");
	    }
	    if (doAsy) {
		fprintf (drvui->fpouta, "// Crystal Structure Drawing Program DRAWxtl %s\n\n",
			 VERSION);
		fprintf (drvui->fpouta,
			 "// Copyright 1996, 1997, 2000-2011 by Larry Finger, Martin Kroeker, and Brian Toby\n");
		fprintf (drvui->fpouta, "//   (all rights reserved).\n\n");
		fprintf (drvui->fpouta, "// Scene file created from %s \n", drvui->Cur_File);
		fprintf (drvui->fpouta,
			 "// with command line options: -b %6.2f %6.2f %6.2f -o %5.2f %5.2f %5.2f \n",
			 boxlim[0], boxlim[1], boxlim[2], origin[0], origin[1],
			 origin[2]);
		fprintf (drvui->fpouta,
			 "// -p %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f -v %5.2f %5.2f %5.2f -m %5.2f\n\n\n",
			 drvui->frames[drvui->frame_no].cryst_lim[0],
			 drvui->frames[drvui->frame_no].cryst_lim[3],
			 drvui->frames[drvui->frame_no].cryst_lim[1],
			 drvui->frames[drvui->frame_no].cryst_lim[4],
			 drvui->frames[drvui->frame_no].cryst_lim[2],
			 drvui->frames[drvui->frame_no].cryst_lim[5], xrot, yrot, zrot,
			 Magnification);
		fprintf (drvui->fpouta, "import three;\n");
		fprintf (drvui->fpouta, "unitsize(0.1cm);\n");
		fprintf (drvui->fpouta, "size(250,250);\n");
		fprintf (drvui->fpouta, "viewportmargin=(1cm,1cm);\n");
		if (!M_cameras) 
		    fprintf (drvui->fpouta, "currentprojection=orthographic (%5.2f,%5.2f,%5.2f,target=(%5.2f,%5.2f,0.),up=Y);\n", cpx,cpy,1.6f*Scale,cpx,cpy);
		else
		    fprintf (drvui->fpouta, "currentprojection=perspective (%5.2f,%5.2f,%5.2f,target=(%5.2f,%5.2f,0.),up=Y);\n", cpx,cpy,1.6f*Scale,cpx,cpy);
		fprintf (drvui->fpouta,"\npicture pic;\n");
		fprintf (drvui->fpouta,"transform3 t=rotate(%.2f,X)*rotate(%.2f,Y)*rotate(%.2f,Z);\n",
			xrot,yrot,zrot);
		fprintf (drvui->fpouta,"\n");

	    }
	    glNewList (drvui->crystalDL, GL_COMPILE);
	} else {
	    if (doPOV)
		fprintf (drvui->fpoutp, "\n/************ Frame %d **********/\n",
			 drvui->frame_no);
	    if (doVrml)
		fprintf (drvui->fpoutv, "\n############# Frame %d ##########\n",
			 drvui->frame_no);
	    if (doAsy)
		fprintf (drvui->fpouta, "\n//############# Frame %d ##########\n",
			 drvui->frame_no);

	}			// only do for drvui->frame_no == 1

	do_drawing ();		// generate all elements of the drawing for this frame
	free (vertex_list);
	vertex_list = NULL;
	free (poly_list);
	poly_list = NULL;
    }				// while drvui->frame_no < drvui->max_frame
    drvui->frame_no = drvui->max_frame;
    glEndList ();
    if (gl_size == 0.0f) {
	gl_size = max (POV_Max[1] - POV_Min[1], POV_Max[0] - POV_Min[0]);
	gl_size = max (gl_size, POV_Max[2] - POV_Min[2]);
    }

    if (doPOV) {
	if (drvui->noshadow == 1)
	    fprintf (drvui->fpoutp, " } no_shadow no_reflection\n");
	else
	    fprintf (drvui->fpoutp, " }\n");
	fprintf (drvui->fpoutp, "rotate <xrot, 0, 0>\n");
	fprintf (drvui->fpoutp, "rotate <0, yrot, 0>\n");
	fprintf (drvui->fpoutp, "rotate <0, 0, zrot>\n");
	fprintf (drvui->fpoutp, "translate <xmove, ymove, 0>\n");
	fprintf (drvui->fpoutp, "}\n");
	if (ShowMapLegend == 1 && drvui->frames[drvui->frame_no].slice > 1) {
	    char label[30];
	    float glr,glg,glb;
	    for (int j=0;j<256;j++){
		float rho = Map_Info.rhomn+0.0039f*(float)j*(Map_Info.rhomx-Map_Info.rhomn);
		colorramp (rho,&glr,&glg,&glb);
		fprintf (drvui->fpoutp,"cylinder { <%5.3f,%5.3f,%5.3f>,<%5.3f,%5.3f,%5.3f>,0.01\n",
			-Scale/20.+0.5,Scale/20.-2.5+0.0078f*(float)j,0.,
			-Scale/20.+0.8,Scale/20.-2.5+0.0078f*(float)j,0.);
		fprintf(drvui->fpoutp," texture{pigment{color red %5.3f green %5.3f blue %5.3f} finish{diffuse 0.25}}\n}\n",glr,glg,glb);
	    }
	    fprintf(drvui->fpoutp,"cylinder { <%5.3f,%5.3f,%5.3f>,<%5.3f,%5.3f,%5.3f>,0.01\n",
		    -Scale/20.+0.4,Scale/20.-2.5+0.0078f*256.*(0.f-Map_Info.rhomn)/(Map_Info.rhomx-Map_Info.rhomn),
		    0.0f, -Scale/20.+0.9,Scale/20.-2.5+0.0078f*256.*(0.f-Map_Info.rhomn)/(Map_Info.rhomx-Map_Info.rhomn),
		    0.0f);
	    fprintf(drvui->fpoutp," texture{pigment{color Black}}\n}\n");
	    for (int i=0;i<6;i++) {
		float fw=0.2f*(float)i*(Map_Info.rhomx-Map_Info.rhomn)+Map_Info.rhomn;
		sprintf(label,"% 5.3f",fw);

		float size = drvui->label_scale * 2.0f * Text_Size;
		fprintf (drvui->fpoutp, "  text { ttf \"crystal.ttf\",\"%s\" 0.15,0\n", label);
		fprintf (drvui->fpoutp, "    scale <%4.2f,%4.2f,%4.2f>\n", size, size, size);
		fprintf (drvui->fpoutp,
			"    translate <%8.5f,%8.5f,%8.5f> pigment{color Black}\n",
			-Scale/20.+1.,Scale/20.-2.5+0.4f*(float)i,0.0f);
		fprintf (drvui->fpoutp, "  }\n");
	    }
	}
    }
    if (doVrml) {
	if (Vrml2)
	    fprintf (drvui->fpoutv, "]\n");
	fprintf (drvui->fpoutv, " }\n");
	if (Vrml2) {
	    fprintf (drvui->fpoutv, " }\n");
	    fprintf (drvui->fpoutv, " }\n");
	}
    }
    if (doAsy) {
	if (ShowMapLegend == 1 && drvui->frames[drvui->frame_no].slice > 1) {
	    char label[30];
	    float glr,glg,glb;
	    for (int j=0;j<256;j++){
		float rho = Map_Info.rhomn+0.0039f*(float)j*(Map_Info.rhomx-Map_Info.rhomn);
		colorramp (rho,&glr,&glg,&glb);
		fprintf (drvui->fpouta," draw(pic, (%8.5f,%8.5f,%8.5f)--(%8.5f,%8.5f,%8.5f),\n",
			-Scale/20.+0.5,Scale/20.-2.5+0.0078f*(float)j,0.,
			-Scale/20.+0.8,Scale/20.-2.5+0.0078f*(float)j,0.);
		fprintf (drvui->fpouta,"rgb(%4.2f,%4.2f,%4.2f) );\n",glr,glg,glb);
	    }
	    fprintf (drvui->fpouta," draw(pic, (%8.5f,%8.5f,%8.5f)--(%8.5f,%8.5f,%8.5f),\n",
		    -Scale/20.+0.4,Scale/20.-2.5+0.0078f*256.*(0.f-Map_Info.rhomn)/(Map_Info.rhomx-Map_Info.rhomn),
		    0.0f, -Scale/20.+0.9,Scale/20.-2.5+0.0078f*256.*(0.f-Map_Info.rhomn)/(Map_Info.rhomx-Map_Info.rhomn),
		    0.0f);
	    fprintf (drvui->fpouta,"rgb(%4.2f,%4.2f,%4.2f) );\n",0.0f,0.0f,0.0f);
	    for (int i=0;i<6;i++) {
		float fw=0.2f*(float)i*(Map_Info.rhomx-Map_Info.rhomn)+Map_Info.rhomn;
		sprintf(label,"% 5.3f",fw);

		fprintf (drvui->fpouta, " label(pic, \"%s\",(%8.5f,%8.5f,%8.5f));\n",
			label, -Scale/20.+1.,Scale/20.-2.5+0.4f*(float)i,0.0f);
	    }
	}
	fprintf (drvui->fpouta,"add(t*pic);\n");
    }

    fprintf (drvui->flout, "End of File.\n");
    if (doVrml)
	fclose (drvui->fpoutv);
    if (doPOV)
	fclose (drvui->fpoutp);
    if (doAsy)
	fclose (drvui->fpouta);
    fclose (drvui->flout);
    fclose (drvui->fcns);
  out:
    in_progress = 0;
}				// end of Generate_Drawing

/* ************************************************************** */
/* ************************************************************** */

void
generate_slab (void)
{
// generate the slab parameters from the input data

    float csal, snal, csbe, snbe, csga, snga;

    float rot_mat[3][3];

    float bmat[3][3];

    float vert[3];

    int i, j;

    for (i = 3; i <= 5; ++i)
	if (drvui->slab_con[i] <= 0)
	    drvui->slab_con[i] = 90.0f;
    for (i = 0; i <= 2; ++i)
	for (j = 0; j <= 2; ++j)
	    bmat[j][i] = 0.0f;	// initialize matrix
    bmat[0][0] = 1.0f;
    csal = (float) cos (drvui->slab_con[3] / RAD);
    snal = (float) sin (drvui->slab_con[3] / RAD);
    csbe = (float) cos (drvui->slab_con[4] / RAD);
    snbe = (float) sin (drvui->slab_con[4] / RAD);
    csga = (float) cos (drvui->slab_con[5] / RAD);
    snga = (float) sin (drvui->slab_con[5] / RAD);
    bmat[0][1] = csga;
    bmat[1][1] = snga;
    bmat[0][2] = csbe;
    bmat[1][2] = (csal - csbe * csga) / snga;
    bmat[2][2] = (float) sqrt (1.0 - bmat[0][2] * bmat[0][2] - bmat[1][2] * bmat[1][2]);
    for (i = 0; i <= 2; ++i)
	for (j = 0; j <= 2; ++j)
	    bmat[j][i] *= drvui->slab_con[i];

    csal = (float) cos (drvui->slab_rot[0] / RAD);
    snal = (float) sin (drvui->slab_rot[0] / RAD);
    csbe = (float) cos (drvui->slab_rot[1] / RAD);
    snbe = (float) sin (drvui->slab_rot[1] / RAD);
    csga = (float) cos (drvui->slab_rot[2] / RAD);
    snga = (float) sin (drvui->slab_rot[2] / RAD);
    rot_mat[0][0] = csbe * csga;
    rot_mat[0][1] = snal * snbe * csga + csal * snga;
    rot_mat[0][2] = -csal * snbe * csga + snal * snga;
    rot_mat[1][0] = -csbe * snga;
    rot_mat[1][1] = -snal * snbe * snga + csal * csga;
    rot_mat[1][2] = csal * snbe * snga + snal * csga;
    rot_mat[2][0] = snbe;
    rot_mat[2][1] = -snal * csbe;
    rot_mat[2][2] = csal * csbe;

/* 0 0 0 */
    for (i = 0; i <= 2; ++i)
	vert[i] = (float) 0.0;
    slabx1 = 0.;
    slaby1 = 0.;
    slabz1 = 0.;

/* 0 0 1 */
    vert[2] = (float) 1.0;
    slabx2 = 0.;
    for (i = 0; i <= 2; ++i)
	slabx2 += (float) (bmat[0][i] * vert[i]);
    slaby2 = 0.;
    for (i = 0; i <= 2; ++i)
	slaby2 += (float) (bmat[1][i] * vert[i]);
    slabz2 = 0.;
    for (i = 0; i <= 2; ++i)
	slabz2 += (float) (bmat[2][i] * vert[i]);

/* 0 1 0 */
    vert[1] = (float) 1.0;
    vert[2] = (float) 0.0;
    slabx3 = 0.;
    for (i = 0; i <= 2; ++i)
	slabx3 += (float) (bmat[0][i] * vert[i]);
    slaby3 = 0.;
    for (i = 0; i <= 2; ++i)
	slaby3 += (float) (bmat[1][i] * vert[i]);
    slabz3 = 0.;
    for (i = 0; i <= 2; ++i)
	slabz3 += (float) (bmat[2][i] * vert[i]);

/* 1 0 0 */
    vert[0] = (float) 1.0;
    vert[1] = (float) 0.0;
    slabx4 = 0.;
    for (i = 0; i <= 2; ++i)
	slabx4 += (float) (bmat[0][i] * vert[i]);
    slaby4 = 0.;
    for (i = 0; i <= 2; ++i)
	slaby4 += (float) (bmat[1][i] * vert[i]);
    slabz4 = 0.;
    for (i = 0; i <= 2; ++i)
	slabz4 += (float) (bmat[2][i] * vert[i]);


    slabx1 -= drvui->slab_off[0];
    slabx2 -= drvui->slab_off[0];
    slabx3 -= drvui->slab_off[0];
    slabx4 -= drvui->slab_off[0];
    slaby1 -= drvui->slab_off[1];
    slaby2 -= drvui->slab_off[1];
    slaby3 -= drvui->slab_off[1];
    slaby4 -= drvui->slab_off[1];

    slabz1 -= drvui->slab_off[2];
    slabz2 -= drvui->slab_off[2];
    slabz3 -= drvui->slab_off[2];
    slabz4 -= drvui->slab_off[2];

    vert[0] = slabx1 * rot_mat[0][0] + slaby1 * rot_mat[1][0] + slabz1 * rot_mat[2][0];
    vert[1] = slabx1 * rot_mat[0][1] + slaby1 * rot_mat[1][1] + slabz1 * rot_mat[2][1];
    vert[2] = slabx1 * rot_mat[0][2] + slaby1 * rot_mat[1][2] + slabz1 * rot_mat[2][2];
    slabx1 = vert[0];
    slaby1 = vert[1];
    slabz1 = vert[2];
    vert[0] = slabx2 * rot_mat[0][0] + slaby2 * rot_mat[1][0] + slabz2 * rot_mat[2][0];
    vert[1] = slabx2 * rot_mat[0][1] + slaby2 * rot_mat[1][1] + slabz2 * rot_mat[2][1];
    vert[2] = slabx2 * rot_mat[0][2] + slaby2 * rot_mat[1][2] + slabz2 * rot_mat[2][2];
    slabx2 = vert[0];
    slaby2 = vert[1];
    slabz2 = vert[2];
    vert[0] = slabx3 * rot_mat[0][0] + slaby3 * rot_mat[1][0] + slabz3 * rot_mat[2][0];
    vert[1] = slabx3 * rot_mat[0][1] + slaby3 * rot_mat[1][1] + slabz3 * rot_mat[2][1];
    vert[2] = slabx3 * rot_mat[0][2] + slaby3 * rot_mat[1][2] + slabz3 * rot_mat[2][2];
    slabx3 = vert[0];
    slaby3 = vert[1];
    slabz3 = vert[2];
    vert[0] = slabx4 * rot_mat[0][0] + slaby4 * rot_mat[1][0] + slabz4 * rot_mat[2][0];
    vert[1] = slabx4 * rot_mat[0][1] + slaby4 * rot_mat[1][1] + slabz4 * rot_mat[2][1];
    vert[2] = slabx4 * rot_mat[0][2] + slaby4 * rot_mat[1][2] + slabz4 * rot_mat[2][2];
    slabx4 = vert[0];
    slaby4 = vert[1];
    slabz4 = vert[2];

    slabv[0] = slabx1;
    slabv[1] = slaby1;
    slabv[2] = slabz1;
    slabv[3] = slabx2;
    slabv[4] = slaby2;
    slabv[5] = slabz2;
    slabv[6] = slabx3;
    slabv[7] = slaby3;
    slabv[8] = slabz3;
    slabv[9] = slabx2 + (slabx3 - slabx1);
    slabv[10] = slaby2 + (slaby3 - slaby1);
    slabv[11] = slabz2 + (slabz3 - slabz1);
    slabv[12] = slabx4;
    slabv[13] = slaby4;
    slabv[14] = slabz4;
    slabv[15] = slabx3 + (slabx4 - slabx1);
    slabv[16] = slaby3 + (slaby4 - slaby1);
    slabv[17] = slabz3 + (slabz4 - slabz1);
    slabv[18] = slabx2 + (slabx4 - slabx1);
    slabv[19] = slaby2 + (slaby4 - slaby1);
    slabv[20] = slabz2 + (slabz4 - slabz1);
    slabv[21] = slabx2 + (slabx4 - slabx1) + slabx3 - slabx1;
    slabv[22] = slaby2 + (slaby4 - slaby1) + slaby3 - slaby1;
    slabv[23] = slabz2 + (slabz4 - slabz1) + slabz3 - slabz1;
}

/* ************************************************************** */
/* ************************************************************** */

void
show_slab (void)
{
// show the outlines of the slab when slabmode == 2
    glPushMatrix ();
    glDisable (GL_LIGHTING);
    glLineWidth (3);
    glColor3f (1., 0., 0.);
    glBegin (GL_LINES);

    glVertex3f (slabx1, slaby1, slabz1);
    glVertex3f (slabx2, slaby2, slabz2);

    glVertex3f (slabx3, slaby3, slabz3);
    glVertex3f (slabx1, slaby1, slabz1);

    glVertex3f (slabx3, slaby3, slabz3);
    glVertex3f (slabx2 + (slabx3 - slabx1), slaby2 + (slaby3 - slaby1),
		slabz2 + (slabz3 - slabz1));

    glVertex3f (slabx2, slaby2, slabz2);
    glVertex3f (slabx2 + (slabx3 - slabx1), slaby2 + (slaby3 - slaby1),
		slabz2 + (slabz3 - slabz1));

    glVertex3f (slabx1, slaby1, slabz1);
    glVertex3f (slabx4, slaby4, slabz4);


    glVertex3f (slabx3, slaby3, slabz3);
    glVertex3f (slabx3 + (slabx4 - slabx1), slaby3 + (slaby4 - slaby1),
		slabz3 + (slabz4 - slabz1));

    glVertex3f (slabx4, slaby4, slabz4);
    glVertex3f (slabx3 + (slabx4 - slabx1), slaby3 + (slaby4 - slaby1),
		slabz3 + (slabz4 - slabz1));

    glVertex3f (slabx2, slaby2, slabz2);
    glVertex3f (slabx2 + (slabx4 - slabx1), slaby2 + (slaby4 - slaby1),
		slabz2 + (slabz4 - slabz1));


    glVertex3f (slabx2 + (slabx4 - slabx1), slaby2 + (slaby4 - slaby1),
		slabz2 + (slabz4 - slabz1));
    glVertex3f (slabx4, slaby4, slabz4);

    glVertex3f (slabx2 + (slabx4 - slabx1), slaby2 + (slaby4 - slaby1),
		slabz2 + (slabz4 - slabz1));
    glVertex3f (slabx2 + (slabx4 - slabx1) + slabx3 - slabx1,
		slaby2 + (slaby4 - slaby1) + slaby3 - slaby1,
		slabz2 + (slabz4 - slabz1) + slabz3 - slabz1);


    glVertex3f (slabx3 + (slabx4 - slabx1), slaby3 + (slaby4 - slaby1),
		slabz3 + (slabz4 - slabz1));
    glVertex3f (slabx2 + (slabx4 - slabx1) + slabx3 - slabx1,
		slaby2 + (slaby4 - slaby1) + slaby3 - slaby1,
		slabz2 + (slabz4 - slabz1) + slabz3 - slabz1);

    glVertex3f (slabx2 + (slabx3 - slabx1), slaby2 + (slaby3 - slaby1),
		slabz2 + (slabz3 - slabz1));
    glVertex3f (slabx2 + (slabx4 - slabx1) + slabx3 - slabx1,
		slaby2 + (slaby4 - slaby1) + slaby3 - slaby1,
		slabz2 + (slabz4 - slabz1) + slabz3 - slabz1);
    glEnd ();
    glPopMatrix ();
    glEnable (GL_LIGHTING);
    glLineWidth (1);
}

/* ************************************************************** */
/* ************************************************************** */

void
label_cell (void)
{
// generate the labels corresponding to the unit-cell axes
    int i, j;

    int used[4];

    int old_type = 1;

    if (!Labels)
	return;
    if (drvui->lat_con[0] > 0. && drvui->lat_con[1] > 0.) {

// initialize the saved_x's  - first pass only
	if (!drvui->labels_inited) {
	    drvui->saved_x_label[0][0][0] = -0.5f / drvui->lat_con[0];
	    drvui->saved_x_label[0][0][1] = -0.5f / drvui->lat_con[1];
	    drvui->saved_x_label[0][0][2] = -0.5f / drvui->lat_con[2];
	    drvui->saved_x_label[0][3][2] = 1.0f + 0.5f / drvui->lat_con[2];
	    drvui->saved_x_label[1][3][2] = 1.4f / drvui->lat_con[2];
	    drvui->saved_x_label[0][2][1] = 1.0f + 0.5f / drvui->lat_con[1];
	    drvui->saved_x_label[1][2][1] = 1.4f / drvui->lat_con[1];
	    drvui->saved_x_label[0][1][0] = 1.0f + 0.5f / drvui->lat_con[0];
	    drvui->saved_x_label[1][1][0] = 1.4f / drvui->lat_con[0];
	    drvui->labels_inited = 1;
	    for (i = 0; i < 4; i++)
		drvui->triple[i] = 0;
	}
	if (fabs (offset[0]) + fabs (offset[1]) + fabs (offset[2]) < 0.01)
	    Locate_Triple ();
	for (i = 0; i < 4; i++)
	    used[i] = 0;

/* find any labels that are already in the list and save their current position */

	if (drvui->nlabel > 1) {
	    for (i = 1; i < drvui->nlabel; i++) {
		if (!strcmp (drvui->labels[i].label_label, "o")) {
		    drvui->triple[0] = 0;
		    for (j = 0; j < 3; j++)
			drvui->saved_x_label[0][0][j] = drvui->labels[i].label_x[j];
		    used[0] = 1;
		    old_type = 0;
		}
		if (!strcmp (drvui->labels[i].label_label, "triple_vect")) {
		    drvui->triple[0] = i;
		    for (j = 0; j < 3; j++)
			drvui->saved_x_label[1][0][j] = drvui->labels[i].label_x[j];
		    used[0] = 1;
		    old_type = 1;
		}
	    }
	    for (i = 1; i < drvui->nlabel; i++) {
		if (!strcmp (drvui->labels[i].label_label, "a")) {
		    drvui->triple[1] = i;
		    for (j = 0; j < 3; j++)
			drvui->saved_x_label[old_type][1][j] =
			    drvui->labels[i].label_x[j];
		    used[1] = 1;
		}
		if (!strcmp (drvui->labels[i].label_label, "b")) {
		    drvui->triple[2] = i;
		    for (j = 0; j < 3; j++)
			drvui->saved_x_label[old_type][2][j] =
			    drvui->labels[i].label_x[j];
		    used[2] = 1;
		}
		if (!strcmp (drvui->labels[i].label_label, "c")) {
		    drvui->triple[1] = i;
		    for (j = 0; j < 3; j++)
			drvui->saved_x_label[old_type][3][j] =
			    drvui->labels[i].label_x[j];
		    used[3] = 1;
		}
	    }
	}

/* following section adds labels to the list for those not already supplied */

	if (!used[0]) {
	    drvui->labels[drvui->nlabel].label_fn = 1;
	    for (j = 0; j < 3; j++) {
		drvui->labels[drvui->nlabel].label_x[j] =
		    drvui->saved_x_label[Display_axes][0][j];
	    }
	    if (!Display_axes) {
		drvui->triple[0] = 0;
		strcpy (drvui->labels[drvui->nlabel++].label_label, "o");
		check_dynamic_storage ();
	    } else {
		drvui->triple[0] = drvui->nlabel;
		strcpy (drvui->labels[drvui->nlabel++].label_label, "triple_vect");
		check_dynamic_storage ();
	    }
	}
	if (!used[1]) {
	    drvui->labels[drvui->nlabel].label_fn = 1;
	    for (j = 0; j < 3; j++) {
		drvui->labels[drvui->nlabel].label_x[j] =
		    drvui->saved_x_label[Display_axes][1][j];
	    }
	    drvui->triple[1] = drvui->nlabel;
	    strcpy (drvui->labels[drvui->nlabel++].label_label, "a");
	    check_dynamic_storage ();
	}
	if (!used[2]) {
	    drvui->labels[drvui->nlabel].label_fn = 1;
	    for (j = 0; j < 3; j++) {
		drvui->labels[drvui->nlabel].label_x[j] =
		    drvui->saved_x_label[Display_axes][2][j];
	    }
	    drvui->triple[2] = drvui->nlabel;
	    strcpy (drvui->labels[drvui->nlabel++].label_label, "b");
	    check_dynamic_storage ();
	}
	if (!used[3]) {
	    drvui->labels[drvui->nlabel].label_fn = 1;
	    for (j = 0; j < 3; j++) {
		drvui->labels[drvui->nlabel].label_x[j] =
		    drvui->saved_x_label[Display_axes][3][j];
	    }
	    drvui->triple[3] = drvui->nlabel;
	    strcpy (drvui->labels[drvui->nlabel++].label_label, "c");
	    check_dynamic_storage ();
	}

/* The following section will copy the appropriate saved values into the active label parameters */

	for (i = 1; i < drvui->nlabel; i++) {
	    if (Display_axes) {
		if (!strcmp (drvui->labels[i].label_label, "o")) {
		    strcpy (drvui->labels[i].label_label, "triple_vect");
		    drvui->triple[0] = i;
		    for (j = 0; j < 3; j++) {
			drvui->labels[i].label_x[j] = drvui->saved_x_label[1][0][j];
		    }
		}
	    } else {
		if (!strcmp (drvui->labels[i].label_label, "triple_vect")) {
		    strcpy (drvui->labels[i].label_label, "o");
		    drvui->triple[0] = 0;
		    for (j = 0; j < 3; j++) {
			drvui->labels[i].label_x[j] = drvui->saved_x_label[0][0][j];
		    }
		}
	    }
	    if (!strcmp (drvui->labels[i].label_label, "c")) {
		drvui->triple[3] = i;
		for (j = 0; j < 3; j++) {
		    drvui->labels[i].label_x[j] =
			drvui->saved_x_label[Display_axes][3][j];
		}
	    }
	    if (!strcmp (drvui->labels[i].label_label, "b")) {
		drvui->triple[2] = i;
		for (j = 0; j < 3; j++) {
		    drvui->labels[i].label_x[j] =
			drvui->saved_x_label[Display_axes][2][j];
		}
	    }
	    if (!strcmp (drvui->labels[i].label_label, "a")) {
		drvui->triple[1] = i;
		for (j = 0; j < 3; j++) {
		    drvui->labels[i].label_x[j] =
			drvui->saved_x_label[Display_axes][1][j];
		}
	    }
	}
    }
}
