// $Id: cursor.cxx 1003 2010-08-23 21:45:23Z martin $
//
// cursor.cxx - routine for DRAWxtl V5.5 - the GUI version 
//
// Coded using the FLTK 1.1.6 widget set
//
//     Larry W. Finger, Martin Kroeker and Brian Toby
//
// This module includes cursor support routines
//
// routines contained within this file:
//
// display_cursor_text - output the cursor position into the status line
// draw_cursor - put the cursor on the screen
// find_atom - find atom nearest the cursor position
// find_proj_atom - find atom that projects closest to the cursor
// move_cursor - move the cursor position

#include "drawxtl.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "draw_ext.h"
#include <FL/x.H>
#include "DRAWxtlViewUI.h"
#include "CrystalView.h"

#include "DRAWxtl_proto.h"

void
display_cursor_text (void)
{
    char Text[500];

    char text[100];

    char text2[100];

    char Text2[100];

    char Text3[100];

    if (!cur_show)
	return;

    strcpy (Text2, "");
    strcpy (Text3, "");
    if (!ReadFourMap)
	sprintf (Text, " Cursor: %.4f,%.4f,%.4f (step=%.2fA)",
		 cur_cen[0], cur_cen[1], cur_cen[2], drvui->cur_step);
    else
	sprintf (Text, " Cursor: %.4f,%.4f,%.4f (step=%.2fA). Rho=%.3f",
		 cur_cen[0], cur_cen[1], cur_cen[2], drvui->cur_step,
		 InterpolateMap (cur_cen[0], cur_cen[1], cur_cen[2]));
    if (cur_atom[0] > 0) {
	sprintf (Text2, " Atoms 1.%s", cur_name[0]);
    }
    if (cur_atom[1] > 0) {
	sprintf (text2, " 2.%s", cur_name[1]);
	strcat (Text2, text2);
	sprintf (Text3, "\n d12=%6.3fA", dist12);
	if (cur_atom[2] > 0) {
	    sprintf (text2, " 3.%s", cur_name[2]);
	    strcat (Text2, text2);
	    sprintf (text, ", d23=%6.3fA, a123=%5.2f°", dist23, ang123);
	    strcat (Text3, text);
	    if (cur_atom[3] > 0) {
		sprintf (text2, " 4.%s", cur_name[3]);
		strcat (Text2, text2);
		sprintf (text, ", d34=%6.3fA, a234=%5.2f°, tor1234=%5.2f°", dist34,
			 ang234, torsion_ang);
		strcat (Text3, text);
	    }
	}
    }
    strcat (Text, Text2);
    strcat (Text, Text3);
    drvui->Cursor_pos->value (Text);

    if (drvui->Cursor_posW) 
	update_cursor_window ();
}

void
draw_cursor (void)
{

    int i, j;

    float vert[3];

    if (!cur_show)
	return;

    nvert = 0;			/* initialize vertex list */
    for (i = 0; i <= 2; ++i) {
	for (j = 0; j < 3; j++)
	    vert[j] = cur_cen[j];
	vert[i] = cur_cen[i] - 0.5f / drvui->lat_con[i];
	add_vert_nc (vert);
	vert[i] = cur_cen[i] + 0.5f / drvui->lat_con[i];
	add_vert_nc (vert);
    }
    glPushMatrix ();
    glDisable (GL_LIGHTING);
    glColor3f (1.0, 0.0, 0.0);

    glBegin (GL_LINES);
    glVertex3f (s_vert[0], s_vert[1], s_vert[2]);
    glVertex3f (s_vert[3], s_vert[4], s_vert[5]);

    glVertex3f (s_vert[6], s_vert[7], s_vert[8]);
    glVertex3f (s_vert[9], s_vert[10], s_vert[11]);

    glVertex3f (s_vert[12], s_vert[13], s_vert[14]);
    glVertex3f (s_vert[15], s_vert[16], s_vert[17]);
    glEnd ();
    glEnable (GL_LIGHTING);

    glPopMatrix ();
}


int
find_atom (void)
{
    int i, j;

    double distance = 1000.0, d;

    int closest = 0;

    nvert = 0;
    for (i = 0; i < natom; i++)
	find_all_in_box (i);

    for (j = 0; j < nvert; j++) {
	d = (o_vert[3 * j] - cur_cen[0]) * (o_vert[3 * j] - cur_cen[0])
	    + (o_vert[3 * j + 1] - cur_cen[1]) * (o_vert[3 * j + 1] - cur_cen[1])
	    + (o_vert[3 * j + 2] - cur_cen[2]) * (o_vert[3 * j + 2] - cur_cen[2]);
	if (d < distance) {
	    distance = d;
	    closest = j;
	}
    }
    cur_cen[0] = o_vert[3 * closest];
    cur_cen[1] = o_vert[3 * closest + 1];
    cur_cen[2] = o_vert[3 * closest + 2];
    return closest;
}

int
find_proj_atom (int x, int y)
{
    int j, k;

    double distance = 1.0e15, d;

    int closest = 0;

    int clipped = 0;

    GLdouble winX, winY, winZ;

    int saved_nvert = nvert;

    for (j = 0; j < natom; ++j) {	// get all atoms in display box
	find_all_in_box (j);
    }
    for (j = saved_nvert; j < nvert; j++) {	// loop through all atoms
	clipped = 0;
	if (clipflag == 1) {
	    for (k = 0; k < 3; k++)
		if (o_vert[3 * j + k] < drvui->frames[drvui->frame_no].clip_lim[k] - 0.01f
		    || o_vert[3 * j + k] >
		    drvui->frames[drvui->frame_no].clip_lim[k + 3] + 0.01f)
		    clipped = 1;
	}
	if (clipped == 1)
	    continue;
	gluProject (s_vert[3 * j], s_vert[3 * j + 1], s_vert[3 * j + 2],
		    modelMatrix, projMatrix, viewport, &winX, &winY, &winZ);
	d = ((float) x - winX) * ((float) x - winX) +
	    ((float) y + winY - viewport[3]) * ((float) y + winY - viewport[3]);
	if (d < distance) {	// find closest atom (in projection
	    distance = d;
	    closest = j;
	}
    }
    nvert = saved_nvert;
    cur_cen[0] = o_vert[3 * closest];	// set cursor to atom coordinates
    cur_cen[1] = o_vert[3 * closest + 1];
    cur_cen[2] = o_vert[3 * closest + 2];

    return closest;
}

void
move_cursor (int axis, float inc_amt)
{
    if (!cur_show)
	return;

//    cur_name[0][0] = '\0';
    cur_cen[axis] += inc_amt * drvui->cur_step / drvui->lat_con[axis];
}

void
update_cursor_window (void)
{
    char string[100], string2[30];

    float centr2[3] = {0.0f, 0.0f, 0.0f}, centr3[3] = {0.0f, 0.0f, 0.0f}, centr4[3] = {0.0f,0.0f,0.0f};

    int lastatom, i;

    float dist1 = 0.0f, dist2 = 0.0f, dist3 = 0.0f;

    float ang1 = 0.0f, ang2 = 0.0f;

    drvui->Cursor_posW->clear();

    if (!cur_show) {
	sprintf (string, "Press 'c' to activate the cursor mode in the graphics area");
	drvui->Cursor_posW->add (string);
	return;
    }
  
    sprintf (string, "Cursor position:\t\t%.4f\t%.4f\t%.4f",
		 cur_cen[0], cur_cen[1], cur_cen[2]);
    drvui->Cursor_posW->add (string);
    sprintf (string, "Current step size\t%.2fA", drvui->cur_step);
    drvui->Cursor_posW->add (string);

    if (ReadFourMap) {
	sprintf (string, "Local rho value\t\t\t%.3f",
		 InterpolateMap (cur_cen[0], cur_cen[1], cur_cen[2]));
	drvui->Cursor_posW->add (string);
    }

    drvui->Cursor_posW->add ("\n");
    
    lastatom = -1;
    for (i = 0; i < 4; i++) 
	if (cur_atom[i] > 0) lastatom++;

    switch (lastatom) {
	case 3:
	    dist1=dist34;
	    dist2=dist23;
	    dist3=dist12;
	    ang1=ang234;
	    ang2=ang123;
	    break;
	case 2:
	    dist1=dist23;
	    dist2=dist12;
	    ang1=ang123;
            break;
	case 1:
	    dist1=dist12;
	    break;
	default:
	    break;
    }
    if (lastatom > -1) {
	sprintf (string, "Current atom:\t%4s\t%.4f\t%.4f\t%.4f", cur_name[lastatom],
		o_vert[3 * cur_atom[lastatom]], o_vert[3 * cur_atom[lastatom] + 1], 
		o_vert[3 * cur_atom[lastatom] + 2]);
	drvui->Cursor_posW->add (string);
    }
    strcpy (string,"\n");
    drvui->Cursor_posW->add (string);
    drvui->Cursor_posW->add ("@-");
    drvui->Cursor_posW->add (string);

    if (lastatom > 0) {
	sprintf (string, "Distance\t%4s-%4s\t\t\t%.3f",cur_name[lastatom],cur_name[lastatom-1],dist1);
	drvui->Cursor_posW->add (string);
    } 
    if (lastatom > 1) {
	sprintf (string, "Angle\t%4s-%4s\t-%4s\t\t%.3f",cur_name[lastatom],cur_name[lastatom-1],cur_name[lastatom-2],ang1);
	drvui->Cursor_posW->add (string);
    }
    if (cur_atom[3] > 0) {
	sprintf (string, "Dihedral angle\t%4s-%4s\t-%4s\t-%4s\t%.3f", cur_name[3],
		cur_name[2], cur_name[1], cur_name[0], torsion_ang);
	drvui->Cursor_posW->add (string);
	sprintf (string, "Distance\t%4s-%4s\t\t\t%.3f", cur_name[1], cur_name[2], dist23);
	drvui->Cursor_posW->add (string);
	sprintf (string, "Distance\t%4s-%4s\t\t\t%.3f", cur_name[1], cur_name[0], dist12);
	drvui->Cursor_posW->add (string);
	sprintf (string, "Angle\t%4s-%4s\t-%4s\t\t%.3f", cur_name[2], cur_name[1], cur_name[0],
		ang2);
	drvui->Cursor_posW->add (string);
    } else if (cur_atom[2] > 0) {
	sprintf (string, "Distance\t%4s-%4s\t\t\t%.3f", cur_name[1], cur_name[2], dist23);
	drvui->Cursor_posW->add (string);
    }

    strcpy (string, "\n");
    drvui->Cursor_posW->add (string);
    drvui->Cursor_posW->add ("@-");
    drvui->Cursor_posW->add (string);

    if (lastatom > 0) {
	sprintf (string, "Distance Matrix\t%4s\t%4s",
		cur_name[lastatom], cur_name[lastatom - 1]);
	if (lastatom > 1) { 
	    sprintf (string2, "\t%4s", cur_name[lastatom - 2]);
	    strcat (string, string2);
	}
        if (lastatom > 2) { 
	    sprintf (string2, "\t%4s", cur_name[lastatom - 3]);
	    strcat (string, string2);
	}
	drvui->Cursor_posW->add (string);
	sprintf (string, "%4s\t-\t%5.3f", cur_name[lastatom], dist1);
        if (lastatom > 1) { 
	    sprintf (string2, "\t%5.3f", dist (cur_atom[lastatom], cur_atom[lastatom - 2]));
	    strcat (string, string2);
	}
        if (lastatom > 2) {
	    sprintf(string2, "\t%5.3f", dist(cur_atom[lastatom],cur_atom[lastatom - 3]));
	    strcat (string, string2);
        }
	drvui->Cursor_posW->add (string);
	centr2[0] = (o_vert[3 * cur_atom[lastatom]] + o_vert[3 * cur_atom[lastatom - 1]]) / 2.0f; 
	centr2[1] = (o_vert[3 * cur_atom[lastatom] + 1] + o_vert[3 * cur_atom[lastatom - 1] + 1] ) / 2.0f;
	centr2[2] = (o_vert[3 * cur_atom[lastatom] + 2] + o_vert[3 * cur_atom[lastatom - 1] + 2] ) / 2.0f;
    } 

    if (lastatom > 1) {
	sprintf (string, "%4s\t%.3f\t-\t%.3f",
		cur_name[lastatom - 1],  dist1, dist2);
        if (lastatom > 2) {
	    sprintf (string2, "\t%5.3f", dist (cur_atom[lastatom - 1], cur_atom[lastatom - 3 ]));
	    strcat (string, string2);
        }
	drvui->Cursor_posW->add (string);
        
	sprintf (string, "%4s\t%.3f\t%.3f\t-",
		cur_name[lastatom - 2], dist (cur_atom[lastatom], cur_atom[lastatom - 2]),  dist (cur_atom[lastatom - 1], cur_atom[lastatom - 2 ]));
	if (lastatom > 2) {
	    sprintf (string2, "\t%.3f", dist3);
	    strcat (string, string2);
        }
	drvui->Cursor_posW->add (string);

	centr3[0] = (o_vert[3 * cur_atom[lastatom]] + o_vert[3 * cur_atom[lastatom - 1]] 
		    + o_vert[3 * cur_atom[lastatom - 2]]) / 3.0f;
	centr3[1] = (o_vert[3 * cur_atom[lastatom] + 1] + o_vert[3 * cur_atom[lastatom - 1] + 1] 
		    + o_vert[3 * cur_atom[lastatom - 2] + 1]) / 3.0f;
	centr3[2] = (o_vert[3 * cur_atom[lastatom] + 2] + o_vert[3 * cur_atom[lastatom - 1] + 2] 
		    + o_vert[3 * cur_atom[lastatom - 2] + 2]) / 3.0f;
    }
    if (lastatom > 2) {
	sprintf (string, "%4s\t%.3f\t%.3f\t%.3f\t-", cur_name[0],
		dist(cur_atom[lastatom], cur_atom[lastatom - 3]),
		dist(cur_atom[lastatom - 1], cur_atom[lastatom - 3]), dist3 );
	drvui->Cursor_posW->add (string);

	centr4[0] = (o_vert[3 * cur_atom[0]] + o_vert[3 * cur_atom[1]] + o_vert[3 * cur_atom[2]] 
		    + o_vert[3 * cur_atom[3]]) / 4.0f;
	centr4[1] = (o_vert[3 * cur_atom[0] + 1] + o_vert[3 * cur_atom[1] + 1]
		    + o_vert[3 * cur_atom[2] + 1] + o_vert[3 * cur_atom[3] + 1]) / 4.0f;
	centr4[2] = (o_vert[3 * cur_atom[0] + 2] + o_vert[3 * cur_atom[1] + 2]
		    + o_vert[3 * cur_atom[2] + 2] + o_vert[3 * cur_atom[3] + 2]) / 4.0f;
    }
    strcpy (string, "\n");
    drvui->Cursor_posW->add (string);
    drvui->Cursor_posW->add ("@-");
    drvui->Cursor_posW->add (string);
    if (cur_atom[1] > 0) {
	sprintf (string, "Center of top 2 atoms\t\t%.4f\t%.4f\t%.4f",
		centr2[0], centr2[1], centr2[2]);
	drvui->Cursor_posW->add (string);
    }
    if (cur_atom[2] > 0) {
	sprintf (string, "Center of top 3 atoms\t\t%.4f\t%.4f\t%.4f",
		centr3[0], centr3[1], centr3[2]);
	drvui->Cursor_posW->add (string);
    }
    if (cur_atom[3] > 0) {
	sprintf (string, "Center of all 4 atoms\t\t%.4f\t%.4f\t%.4f",
		centr4[0], centr4[1], centr4[2]);
	drvui->Cursor_posW->add (string);
    }
}
