// $Id: Edit2.cxx 1116 2011-02-26 14:24:12Z martin $
//
// Edit2.cxx - routine for DRAWxtl V5.5 - the GUI version 
//
// Coded using the FLTK 1.1.6 widget set
//
//     Larry W. Finger, Martin Kroeker and Brian Toby
//
// This module includes the majority of the edit screens for the GUI
//
// routines contained within this file:
//
//  Edit_Polyhedra_Close_cb - callback routine to dismiss Polyhedra Edit screen
//  Edit_Polyhedra_Save_cb - callback routine to update the polyhedral parameters
//  Edit_Spheres_cb - callback routine to load spheres edit screen
//  Edit_Spheres_Close_cb - callback routine to dismiss Sphere Edit screen
//  Edit_Spheres_Save_cb - callback routine to update the sphere parameters
//  Edit_STR_cb - callback routine to load Edit Str File screen
//  Edit_STR_Save_cb - callback routine to save "Edit Parameter" data
//  MapType_cb - callback routine to show or hide calculation button depending on map type  
//  Map_Info_cb - callback routine to load map information screen
//  Lone_Pair_Combo_cb - updates widgets when atom in combo box is changed
//  LonePair_Frame_Combo_cb - callback routine enterd when frame combo box changed
//  Modify_Bonds_cb - callback routine entered when text selected in main text window
//  Modify_Bonds_Distance_cb - callback routine when text selected in "Bond To" distance panel
//  Modify_LonePair_cb - callback routine when LonePair text buffer modified
//  Modify_Maps_cb - callback routine entered when text selected in main text window
//  Modify_Polyhedra_cb - callback routine when text selected in PolyhedraBuffer
//  Modify_Polyhedra_Distance_cb - callback routine when text selected in "Polyhedra To" edit panel
//  Modify_Spheres_cb - callback routine when text selected in "Sphere Parameters" edit panel
//  New_Arrow_Add_cb - callback routine entered when 'Add' button pushed on Arrow Edit Screen
//  New_Arrow_Input_cb - callback routine entered when any of the entries on the arrow add section are changed
//  New_Bond_Add_cb - callback routine entered when 'Add' button pushed on Bond Edit Screen
//  New_Bond_Input_cb - callback routine entered when any of the entries on the bond add line are changed
//  New_Lone_Pair_Add_cb - callback routine entered when lone pair cone added
//  New_Polyhedra_Add_cb - callback routine entered when 'Add' button pushed on Polyhedra Edit Screen
//  New_Polyhedra_Input_cb - callback routine entered when any of the entries on the polyhedra add line are changed
//  New_Sphere_Input_cb - callback routine entered when new sphere added
//  next_focus - routine to select next widget for focus when <TAB> is pressed
//  Polyhedra_Combo_cb - callback routine entered when atom in combo box is changed
//  Polyhedra_Frame_Combo_cb - callback entered when frame combo is changed
//  Sphere_Combo_cb - callback when Sphere combo box clicked
//  Sphere_Frame_Combo_cb - callback entered when frame combo is changed
//  View_Console_cb - callback routine to view console listing
//  View_Cursor_cb - callback routine to view cursor feedback window
//  View_File_cb - callback routine to view any file
//  View_Listing_cb - callback to view listing file
//  View_List_Close_cb - callback to close str file editor window
//  View_Listing_Close_cb - callback to close Listing windows
//  View_POV_cb - callback to view the POV file

#include "drawxtl.h"
#include "DRAWxtlViewUI.h"
#include "gl2ps.h"
#include "EditView.h"
#include "Ellipsoids.h"
#include "draw_ext.h"
#if WIN32
#include <direct.h>
#else
#include <unistd.h>
#define _chdir chdir
#endif
#include <ctype.h>

#include "DRAWxtl_proto.h"

static int zero2 = 0;

static int three = 3;

extern char Edit_title[128];

#ifdef WIN32
const char *flu_file_chooser (const char *message, const char *pattern,
			      const char *filename);
#endif

void
Edit_Polyhedra_Close_cb (void)
{
    Polyhedra->Polyhedra_Edit_Window->hide ();
    Restore_Working_Copy ();	// undo any changes
    Generate_Drawing (1);	// regenerate
}

void
Edit_Polyhedra_Save_cb (Fl_Button *, int *save)
{
// callback routine when 'save' or 'apply' button is pressed on the Edit Polyhhedra screen
    char atom1[5], atom2[5], color[40], widget[16382], type[3];

    char tmp[20];

    int i = 0;

    int np = 0;

    unsigned int j, k;

    int l;

    float d1, d2, transp;

    int Frame_No = 1;

    char *selection;

    if (drvui->max_frame > 1) {
	Frame_No = atoi (Polyhedra->Frame_No->value ());
	for (i = 1, j = 1; i < drvui->npoly; i++) {	// copy parameters for other frames to
	    if (drvui->polyhedra[i].poly_fn != Frame_No) {	//   start of list
		if ((int) j != i) {
		    drvui->polyhedra[j].poly_fn = drvui->polyhedra[i].poly_fn;
		    strncpy (drvui->polyhedra[j].poly_col, drvui->polyhedra[i].poly_col,
			     40);
		    strncpy (drvui->polyhedra[j].poly_l, drvui->polyhedra[i].poly_l, 4);
		    strncpy (drvui->polyhedra[j].poly_t, drvui->polyhedra[i].poly_t, 4);
		    drvui->polyhedra[j].poly_min = drvui->polyhedra[i].poly_min;
		    drvui->polyhedra[j].poly_size = drvui->polyhedra[i].poly_size;
		}
		j++;
	    }
	}
	drvui->npoly = j;
	for (i = 1, j = 1; i < drvui->nplane; i++) {	// copy parameters for other frames to
	    if (drvui->planes[i].plane_fn != Frame_No) {	//   start of list
		if ((int) j != i) {
		    drvui->planes[j].plane_fn = drvui->planes[i].plane_fn;
		    strncpy (drvui->planes[j].plane_col, drvui->planes[i].plane_col, 40);
		    strncpy (drvui->planes[j].plane_l, drvui->planes[i].plane_l, 4);
		    drvui->planes[j].plane_size = drvui->planes[i].plane_size;
		}
		j++;
	    }
	}
	drvui->nplane = j;
	i = drvui->npoly - 1;
    }
    selection = Polyhedra->PolyhedraBuffer->text ();
    strcpy (widget, selection);
    free (selection);
    drvui->npoly = 0;
    drvui->nplane = 0;
    if (strlen (widget) < 10) {
	strcpy (widget, "");
	Polyhedra->PolyhedraBuffer->text (widget);
    } else {
	while (strlen (widget) > 10) {
	    i++;
	    sscanf (widget, " %s", type);
	    strcpy (atom2, "");
	    transp = 0.;
	    if (!strncmp (type, "PS", 2) || !strncmp (type, "SH", 2)
		|| !strncmp (type, "PL", 2)) {
		sscanf (widget, "%s %s %f %f %s %s %f", type, atom1, &d1, &d2, color, tmp,
			&transp);
		strcpy (atom2, "");
		strcpy (drvui->polyhedra[i].poly_t, atom2);
	    } else {
		sscanf (widget, "%s %s %s %f %f %s %s %f", type, atom1, atom2, &d1, &d2,
			color, tmp, &transp);
		strcpy (drvui->polyhedra[i].poly_t, atom2);
	    }
	    trim_string (color, 40);
	    if (!strlen (color))
		strcpy (color, "Gray20");
	    if (isdigit (color[0]) || color[0] == '.') {
		char color1[10], color2[10], color3[10];

		if (!strncmp (type, "PS", 2) || !strncmp (type, "SH", 2)
		    || !strncmp (type, "PL", 2)) {
		    sscanf (widget, "%s %s %f %f %s %s %s %s %f", type, atom1, &d1, &d2,
			    color1, color2, color3, tmp, &transp);
		    strcpy (atom2, "");
		    strcpy (drvui->polyhedra[i].poly_t, atom2);
		} else {
		    sscanf (widget, "%s %s %s %f %f %s %s %s %s %f", type, atom1, atom2,
			    &d1, &d2, color1, color2, color3, tmp, &transp);
		    strcpy (drvui->polyhedra[i].poly_t, atom2);
		}
		strcpy (color, color1);
		strcat (color, " ");
		strcat (color, color2);
		strcat (color, " ");
		strcat (color, color3);
		trim_string (color, 40);
	    }
	    if (transp > 0.0f) {
		sprintf (tmp, " filter %.3f", transp);
		strcat (color, tmp);
	    }
	    if (!strncmp (type, "PL", 2)) {
		np++;
		drvui->planes[np].plane_size = d2;
		strncpy (drvui->planes[np].plane_col, color, 25);
		strcpy (drvui->planes[np].plane_l, atom1);
		drvui->planes[np].plane_fn = Frame_No;
		drvui->nplane = np + 1;
	    } else {
		strncpy (drvui->polyhedra[i - np].poly_l, atom1, 4);
		drvui->polyhedra[i - np].poly_min = d1;
		drvui->polyhedra[i - np].poly_size = d2;
		strcpy (drvui->polyhedra[i - np].poly_col, color);
		drvui->polyhedra[i - np].poly_fn = Frame_No;
		drvui->npoly = i + 1 - np;
		for (l = 0; l < drvui->nedges; l++) {
		    if (!strcmp (drvui->polyedges[l].name, atom1)) {
			strcpy (drvui->polyhedra[i-np].poly_col_edge, drvui->polyedges[l].color);
			drvui->polyhedra[i-np].poly_rad_edge = drvui->polyedges[l].radius;
			break;
		    }
		}
	    }
	    for (j = 0; j < strlen (widget); j++) {
		if (widget[j] == '\n')
		    break;
	    }
	    for (j++, k = 0; j < strlen (widget); j++)
		widget[k++] = widget[j];
	    widget[k] = 0;
	}
    }
    strcpy (drvui->col_edge, Polyhedra->Def_Edge_Color->value ());
    drvui->rad_edge = (float) atof (Polyhedra->Def_Edge_Radius->value ());
    if (drvui->rad_edge > 0.0)
	edges = 2;
    drvui->Str_File_Changed = 1;
    Update_Str (0);
    Generate_Drawing (0);
    if (*save != 3) {
	Save_Working_Copy ();
	Polyhedra->PolyInstr2->hide ();
    }
    if (*save == 1) {
	Polyhedra->Polyhedra_Edit_Window->hide ();
    }
    Fl::redraw ();
}

void
Edit_Slab_cb (void)
{
// routine to create the edit slab screen
    char string[100];

    static int one = 1;

    int y = 50;

    static const Fl_Menu_Item slabmodes[] = {
	{"ignore", 0, NULL, (void *) 0, 0, 0, 0, 14, 56},
	{"apply", 0, NULL, (void *) 0, 0, 0, 0, 14, 56},
	{"show", 0, NULL, (void *) 0, 0, 0, 0, 14, 56},
	{0, 0, 0, 0, 0, 0, 0, 0, 0}
    };

    if (!strlen (drvui->CurFile->value ())) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }
    if (!Slabs) {
	Slabs = new SlabParam;	// new instance of the slab parameters
	Slabs->SlabWindow = new Fl_Window (50, 50, 400, 430, "Edit Slab Parameters");
	Slabs->SlabWindow->callback ((Fl_Callback *) Edit_Slab_Close_cb);
	{
	    Fl_Text_Display *o = new Fl_Text_Display (0, y, 400, 0, "Slab Dimensions");

	    o->box (FL_NO_BOX);
	    o->labelfont (1);
	    o->textcolor (1);
	}
	y += 25;
	Slabs->Slab_A = new Fl_Input (20, y, 100, 25, "a");
	Slabs->Slab_A->align (FL_ALIGN_TOP);
	Slabs->Slab_A->labelfont (1);

	Slabs->Slab_B = new Fl_Input (150, y, 100, 25, "b");
	Slabs->Slab_B->align (FL_ALIGN_TOP);
	Slabs->Slab_B->labelfont (1);

	Slabs->Slab_C = new Fl_Input (280, y, 100, 25, "c");
	Slabs->Slab_C->align (FL_ALIGN_TOP);
	Slabs->Slab_C->labelfont (1);

	y += 40;
	Slabs->Slab_Alpha = new Fl_Input (20, y, 100, 25, "alpha");
	Slabs->Slab_Alpha->align (FL_ALIGN_TOP);
	Slabs->Slab_Alpha->labelfont (1);

	Slabs->Slab_Beta = new Fl_Input (150, y, 100, 25, "beta");
	Slabs->Slab_Beta->align (FL_ALIGN_TOP);
	Slabs->Slab_Beta->labelfont (1);

	Slabs->Slab_Gamma = new Fl_Input (280, y, 100, 25, "gamma");
	Slabs->Slab_Gamma->align (FL_ALIGN_TOP);
	Slabs->Slab_Gamma->labelfont (1);

	y += 60;
	{
	    Fl_Text_Display *o = new Fl_Text_Display (0, y, 400, 0, "Slab Orientation");

	    o->box (FL_NO_BOX);
	    o->labelfont (1);
	    o->textcolor (1);
	}

	y += 30;
	Slabs->Slab_Off_X = new Fl_Input (20, y, 100, 25, "X Offset");
	Slabs->Slab_Off_X->align (FL_ALIGN_TOP);
	Slabs->Slab_Off_X->labelfont (1);

	Slabs->Slab_Off_Y = new Fl_Input (150, y, 100, 25, "Y Offset");
	Slabs->Slab_Off_Y->align (FL_ALIGN_TOP);
	Slabs->Slab_Off_Y->labelfont (1);

	Slabs->Slab_Off_Z = new Fl_Input (280, y, 100, 25, "Z Offset");
	Slabs->Slab_Off_Z->align (FL_ALIGN_TOP);
	Slabs->Slab_Off_Z->labelfont (1);

	y += 60;
	Slabs->Slab_Rot_X = new Fl_Input (20, y, 100, 25, "X Rotation");
	Slabs->Slab_Rot_X->align (FL_ALIGN_TOP);
	Slabs->Slab_Rot_X->labelfont (1);

	Slabs->Slab_Rot_Y = new Fl_Input (150, y, 100, 25, "Y Rotation");
	Slabs->Slab_Rot_Y->align (FL_ALIGN_TOP);
	Slabs->Slab_Rot_Y->labelfont (1);

	Slabs->Slab_Rot_Z = new Fl_Input (280, y, 100, 25, "Z Rotation");
	Slabs->Slab_Rot_Z->align (FL_ALIGN_TOP);
	Slabs->Slab_Rot_Z->labelfont (1);

	y += 60;
	Slabs->Slab_Mode = new Fl_Choice (175, y, 70, 25, "Slab Mode");
	Slabs->Slab_Mode->align (FL_ALIGN_TOP);
	Slabs->Slab_Mode->labelfont (1);
	Slabs->Slab_Mode->menu (slabmodes);

	Slabs->Slab_Mode->value (slabmode);

#if !defined (WIN32) && !defined (__APPLE__)
	Slabs->SlabWindow->icon ((char *) drvui->icon);
#endif
	y += 60;
	Fl_Button *r = new Fl_Button (65, y, 70, 25, "Close");

	r->tooltip ("Close this window and discard all changes.");
	r->callback ((Fl_Callback *) Edit_Slab_Close_cb);
	Fl_Button *s = new Fl_Button (265, y, 70, 25, "Save");

	s->tooltip
	    ("Apply current contents of top box to drawing, then close this window.");
	s->callback ((Fl_Callback *) Edit_Slab_Save_cb, &one);
	Fl_Button *a = new Fl_Button (165, y, 70, 25, "Apply");

	a->tooltip
	    ("Apply current contents of top box to drawing, but leave this window open.");
	a->callback ((Fl_Callback *) Edit_Slab_Save_cb, &zero2);
	Slabs->SlabWindow->end ();
    }
    sprintf (string, "%6.3f", drvui->slab_con[0]);
    Slabs->Slab_A->value (string);
    sprintf (string, "%6.3f", drvui->slab_con[1]);
    Slabs->Slab_B->value (string);
    sprintf (string, "%6.3f", drvui->slab_con[2]);
    Slabs->Slab_C->value (string);
    sprintf (string, "%6.3f", drvui->slab_con[3]);
    Slabs->Slab_Alpha->value (string);
    sprintf (string, "%6.3f", drvui->slab_con[4]);
    Slabs->Slab_Beta->value (string);
    sprintf (string, "%6.3f", drvui->slab_con[5]);
    Slabs->Slab_Gamma->value (string);
    sprintf (string, "%6.3f", drvui->slab_off[0]);
    Slabs->Slab_Off_X->value (string);
    sprintf (string, "%6.3f", drvui->slab_off[1]);
    Slabs->Slab_Off_Y->value (string);
    sprintf (string, "%6.3f", drvui->slab_off[2]);
    Slabs->Slab_Off_Z->value (string);
    sprintf (string, "%6.3f", drvui->slab_rot[0]);
    Slabs->Slab_Rot_X->value (string);
    sprintf (string, "%6.3f", drvui->slab_rot[1]);
    Slabs->Slab_Rot_Y->value (string);
    sprintf (string, "%6.3f", drvui->slab_rot[2]);
    Slabs->Slab_Rot_Z->value (string);
    Slabs->SlabWindow->show ();
}

void
Edit_Slab_Close_cb (void)
{
    Slabs->SlabWindow->hide ();
}

void
Edit_Slab_Save_cb (Fl_Button *, int *save)
{
// callback routine when 'save' or 'apply' button is pressed on the Edit Slab screen

    Omit->nomits = 0;
    drvui->slab_con[0] = (float) atof (Slabs->Slab_A->value ());
    drvui->slab_con[1] = (float) atof (Slabs->Slab_B->value ());
    drvui->slab_con[2] = (float) atof (Slabs->Slab_C->value ());
    drvui->slab_con[3] = (float) atof (Slabs->Slab_Alpha->value ());
    drvui->slab_con[4] = (float) atof (Slabs->Slab_Beta->value ());
    drvui->slab_con[5] = (float) atof (Slabs->Slab_Gamma->value ());
    drvui->slab_off[0] = (float) atof (Slabs->Slab_Off_X->value ());
    drvui->slab_off[1] = (float) atof (Slabs->Slab_Off_Y->value ());
    drvui->slab_off[2] = (float) atof (Slabs->Slab_Off_Z->value ());
    drvui->slab_rot[0] = (float) atof (Slabs->Slab_Rot_X->value ());
    drvui->slab_rot[1] = (float) atof (Slabs->Slab_Rot_Y->value ());
    drvui->slab_rot[2] = (float) atof (Slabs->Slab_Rot_Z->value ());
    slabmode = Slabs->Slab_Mode->value ();
    drvui->Str_File_Changed = 1;
    if (*save) {
	Slabs->SlabWindow->hide ();
    }
    Update_Str (0);
    Generate_Drawing (0);
    Fl::redraw ();
}

void
Edit_Spheres_cb (void)
{
// callback routine to load sphere edit screen
    char string[100];

    static int one = 1;

    static int two = 2;

    int i;

    int y_size = 120;

    int y = 25;

    if (!strlen (drvui->CurFile->value ())) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }
/*
    if (!natom) {         // No atoms, no spheres...
        Error_Box("This structure file does not contain any atoms.");
        return;
    }
*/
    Save_Working_Copy ();	// save the working str file
    if (!Spheres) {
	Spheres = new SphereParam;	// new instance of the spheres parameters
	Spheres->Sphere_Edit_Window =
	    new Fl_Window (100, 100, 460, 350, "Edit Sphere Parameters");
	Spheres->Sphere_Edit_Window->callback ((Fl_Callback *) Edit_Spheres_Close_cb);
	Spheres->Sphere_Output_Buffer = NULL;
	Spheres->SphereBuffer = NULL;
	if (drvui->max_frame > 1) {
	    Flu_Combo_List *o = Spheres->Frame_No =
		new Flu_Combo_List (180, y, 75, 25, "Frame No.");
	    o->align (FL_ALIGN_TOP);
	    o->callback (Sphere_Frame_Combo_cb);
	    o->labelfont (1);
	    for (i = 1; i <= drvui->max_frame; i++) {
		sprintf (string, "%d", i);
		o->list.add (string);
	    }
	    o->pop_height (20 * drvui->max_frame);
	    o->value ("1");
	    y += 40;
	    y_size -= 40;
	}
	Spheres->Sphere_Edit = new Fl_Text_Editor (25, y, 410, y_size,
						   "Atom     Size     Color                         ");
	Spheres->SphereBuffer = new Fl_Text_Buffer;
	Spheres->SphereBuffer->add_modify_callback ((Fl_Text_Modify_Cb) Modify_Spheres_cb,
						    (void *) NULL);
	y += 5 + y_size;
	Spheres->SphereInstr = new Fl_Output (25, y, 410, 0, "Press 'Add' to replace "
					      "selected line - 'Remove' to delete it");
	Spheres->SphereInstr->hide ();
	Spheres->SphereInstr->align (FL_ALIGN_BOTTOM);
	Spheres->SphereInstr1 = new Fl_Output (25, y, 410, 0,
					       "Highlight or double-click text above to edit");
	Spheres->SphereInstr1->hide ();
	Spheres->SphereInstr1->align (FL_ALIGN_BOTTOM);
	Spheres->Sphere_Edit->textfont (FL_COURIER);
	Spheres->Sphere_Edit->textsize (12);
	Spheres->Sphere_Edit->buffer (Spheres->SphereBuffer);
	Spheres->Sphere_Edit->labelfont (FL_COURIER_BOLD);
	y += 40;
	Flu_Combo_List *o = Spheres->Sphere_Combo =
	    new Flu_Combo_List (55, y, 100, 25, "Atom");
	o->align (FL_ALIGN_TOP);
	o->callback (Sphere_Combo_cb);
	o->labelfont (1);
	Fl_Input *os = Spheres->New_Sphere_Size = new Fl_Input (165, y, 50, 25, "Size");

	os->align (FL_ALIGN_TOP);
	os->callback ((Fl_Callback *) New_Sphere_Input_cb);
	os->labelfont (1);
	Flu_Combo_List *ot = Spheres->New_Sphere_Color =
	    new Flu_Combo_List (225, y, 160, 25, "Color");
	Load_Color_Combo (ot);
	ot->align (FL_ALIGN_TOP);
	ot->callback ((Fl_Callback *) New_Sphere_Input_cb);
	ot->labelfont (1);
	y += 30;
	Fl_Button *om = Spheres->New_Sphere_Add = new Fl_Button (105, y, 70, 25, "Add");

	om->callback ((Fl_Callback *) New_Sphere_Add_cb, &one);
	om->tooltip ("When active, press to transfer data in boxes to window above");
	om->deactivate ();
	Fl_Button *mm = Spheres->New_Sphere_Remove =
	    new Fl_Button (195, y, 70, 25, "Remove");
	mm->callback ((Fl_Callback *) New_Sphere_Add_cb, &zero2);
	mm->tooltip ("When active, press to remove highlighted line.");
	mm->deactivate ();
	Fl_Button *pm = Spheres->New_Sphere_Convert =
	    new Fl_Button (285, y, 70, 25, "Convert");
	pm->callback ((Fl_Callback *) New_Sphere_Add_cb, &two);
	pm->tooltip
	    ("When active, press to convert atom in boxes from sphere to ellipsoid");
	pm->deactivate ();
	y += 40;
	Spheres->SphereInstr2 = new Fl_Output (25, y, 410, 0,
					       "Changes are temporary until \"Apply\" or \"Save\" is pressed.");
	Spheres->SphereInstr2->hide ();
	Spheres->SphereInstr2->align (FL_ALIGN_BOTTOM);
	y += 30;
	Fl_Button *r = new Fl_Button (105, y, 70, 25, "Close");

	r->tooltip ("Close this window and discard all changes.");
	r->callback ((Fl_Callback *) Edit_Spheres_Close_cb);
	Fl_Button *s = new Fl_Button (285, y, 70, 25, "Save");

	s->tooltip
	    ("Apply current contents of top box to drawing, then close this window.");
	s->callback ((Fl_Callback *) Edit_Spheres_Save_cb, &one);
	Fl_Button *a = new Fl_Button (195, y, 70, 25, "Apply");

	a->tooltip
	    ("Apply current contents of top box to drawing, but leave this window open.");
	a->callback ((Fl_Callback *) Edit_Spheres_Save_cb, &zero2);
	Spheres->Sphere_Edit_Window->end ();
#if !defined (WIN32) && !defined (__APPLE__)
	Spheres->Sphere_Edit_Window->icon ((char *) drvui->icon);
#endif
    }
    Sphere_Frame_Combo_cb (NULL, NULL);
    Spheres->Sphere_Edit_Window->show ();
}

void
Edit_Spheres_Close_cb (void)
{
    drvui->destroy |= SPHERE;
    Spheres->Sphere_Edit_Window->hide ();
    Restore_Working_Copy ();	// undo any changes
    Generate_Drawing (1);	// regenerate
}

void
Edit_Spheres_Save_cb (Fl_Button *, int *save)
{
// callback routine when 'save' or 'apply' button is pressed on the Edit Spheres screen
// save == 0 - apply the changes and update saved working copy
// save == 1 - save the changes and close the window
// save == 3 - apply the changes but do not update saved working copy
    char atom1[5], color[40], widget[16382];

    int j, k;

    int i = 0, n;

    float d1;

    int Frame_No = 1;

    if (drvui->max_frame > 1) {
	Frame_No = atoi (Spheres->Frame_No->value ());
	for (i = 1, j = 1; i < drvui->nsphere; i++) {	// copy parameters for other frames to
	    if (drvui->spheres[i].sphere_fn != Frame_No) {	//   start of list
		if ((int) j != i) {
		    drvui->spheres[j].sphere_fn = drvui->spheres[i].sphere_fn;
		    strcpy (drvui->spheres[j].sphere_col, drvui->spheres[i].sphere_col);
		    strcpy (drvui->spheres[j].sphere_l, drvui->spheres[i].sphere_l);
		    drvui->spheres[j].sphere_size = drvui->spheres[i].sphere_size;
		    drvui->spheres[j].sphere_n = drvui->spheres[i].sphere_n;
		}
		j++;
	    }
	}
	i = j - 1;
    }
    char *selection = Spheres->SphereBuffer->text ();

    strcpy (widget, selection);
    free (selection);
    if (strlen (widget) < 10) {
	drvui->nsphere = i + 1;	// sum over all frames or zero
	strcpy (widget, "");
	Spheres->SphereBuffer->text (widget);
    } else {
	while (strlen (widget) > 10) {
	    i++;
	    char filter1[10], filter2[10], numb[5];

	    memset (filter1, 0, 10);
	    (void) sscanf (widget, "%s %s %f %s %s %s", atom1, numb, &d1, color, filter1,
			   filter2);
	    strcpy (drvui->spheres[i].sphere_l, atom1);
	    drvui->spheres[i].sphere_size = d1;
	    if (!strcmp (filter1, "filter")) {
		strcat (color, " ");
		strcat (color, filter1);
		strcat (color, " ");
		strcat (color, filter2);
	    }
	    strcpy (drvui->spheres[i].sphere_col, color);
	    drvui->spheres[i].sphere_fn = Frame_No;
	    if (strstr (numb, "*"))
		n = -1;
	    else
		(void) sscanf (numb, "%d", &n);
	    drvui->spheres[i].sphere_n = n;
	    drvui->nsphere = i + 1;
	    for (j = 0; j < (int) strlen (widget); j++) {
		if (widget[j] == '\n')
		    break;
	    }
	    for (j++, k = 0; j < (int) strlen (widget); j++)
		widget[k++] = widget[j];
	    widget[k] = 0;
	}
    }
    drvui->Str_File_Changed = 1;
    if (*save != 3) {
	Save_Working_Copy ();
	Spheres->SphereInstr2->hide ();
    }
    if (*save == 1) {
	drvui->destroy |= SPHERE;
	Spheres->Sphere_Edit_Window->hide ();
    }
    Update_Str (0);
    Generate_Drawing (0);
    Fl::redraw ();
}

void
Edit_STR_cb (Fl_Menu_ *, void *arg)	// routine to edit the str file
{
    int r;
    const char *tempname;

    static int one = 1;

    static int mone = -1;

    if (textwindow) {
	textwindow->hide ();
	textwindow->~Fl_Window ();
	delete (textbuf);
	delete (textwindow);
    }
    if (!strlen (drvui->Cur_File)) {
	r =fl_choice ("A Structure File must be selected first.","Ok","Create new",NULL);
	if (!r) return;
	tempname=(fl_input("Enter filename","unnamed.str"));
	if (!tempname) return;
	strcpy(drvui->Cur_File,tempname);
    }
// callback routine to display Edit Str File screen
    Edit_Str_Type = arg;
    if (arg) {
	if (strlen (drvui->Cur_File) < 2)
	    return;
	strcpy (Edit_title, "Edit Original Version of STR File");
    } else {
	if (strlen (drvui->Cur_Temp) < 2)
	    return;
	strcpy (Edit_title, "Edit Working Copy of STR File");
	Update_Str (0);		// make str file current
    }
    Edit_changed = 0;
    Edit_loading = 1;
    textwindow = new Fl_Window (20, 20, 480, 550, Edit_title);
    textwindow->begin ();
    textwindow->callback ((Fl_Callback *) Edit_STR_Close_cb);
    Fl_Text_Editor *display = new Fl_Text_Editor (0, 0, 480, 510);

    display->textfont (FL_COURIER);
    display->textsize (12);
    textbuf = new Fl_Text_Buffer;
    display->buffer (textbuf);
    textbuf->add_modify_callback (Edit_Changed_cb, textwindow);
    textbuf->call_modify_callbacks ();

    Fl_Button *o = new Fl_Button (50, 515, 80, 30, "Close");

    o->tooltip ("Close this window and discard all changes.");
    o->callback ((Fl_Callback *) Edit_STR_Close_cb);
    Fl_Button *oo = new Fl_Button (150, 515, 80, 30, "Apply");

    oo->tooltip ("Apply current contents to drawing, but leave this window open.");
    oo->callback ((Fl_Callback *) Edit_STR_Save_cb, &mone);
    Fl_Button *po = new Fl_Button (250, 515, 80, 30, "Save");

    po->tooltip ("Apply current contents to drawing, then close this window.");
    po->callback ((Fl_Callback *) Edit_STR_Save_cb, &zero2);
    Fl_Button *op = new Fl_Button (350, 515, 80, 30, "Save As");

    op->callback ((Fl_Callback *) Edit_STR_Save_cb, &one);
    op->tooltip ("Rewrite current contents to new file, and close this window.");
    if (arg) {
	textbuf->loadfile (drvui->Cur_File);
    } else {
	textbuf->loadfile (drvui->Cur_Temp);
    }
    Edit_loading = 0;
#if !defined (WIN32) && !defined (__APPLE__)
    textwindow->icon ((char *) drvui->icon);
#endif
    textwindow->end ();
    textwindow->resizable (textwindow);
    textwindow->show ();
}

void
Edit_STR_Close_cb (void)	// callback to destruct editor window
{
// callback to close the STR edit screen
    if (Edit_changed) {		// buffer changed - query user
	int r = fl_choice ("The current file has been modified, but not saved.\n"
			   "Would you like to save it now?",
			   "Cancel", "Yes", "No");

	if (r == 1) {
	    drvui->Str_File_Changed = 1;
	    Edit_STR_Save_cb (NULL, &zero2);
	    return;
	}
	if (r == 0)
	    return;
    }
    drvui->destroy |= TEXT1;
    textwindow->hide ();
}

void
Edit_STR_Save_cb (Fl_Button *, int *action)	// callback to save edited str file
{
// callback routine to save edit parameter data
    FILE *file;

    char *selection;

    drvui->Str_File_Changed = 1;
    if (*action > 0) {		// true for 'Save As'
#if defined(WIN32)
	char drive[_MAX_DRIVE];

	char dir[_MAX_DIR];

	char fname[_MAX_FNAME];

	char ext[_MAX_EXT];

	const char *newfile =
	    flu_file_chooser ("Select New Name for Data/Experiment File",
			      "*.str", drvui->Cur_File);

	if (newfile) {
	    _splitpath (newfile, drive, dir, fname, ext);	//Windows code
	    strcpy (drvui->Cur_Dir, drive);	// Drive letter
	    strcat (drvui->Cur_Dir, dir);	//   and directory
	    strcpy (drvui->Cur_File, fname);	// copy file name
	    strcat (drvui->Cur_File, ext);	// and add extension
#else
	int i, k = 0;

	char *newfile =
	    fl_file_chooser ("Select New Name for Data/Experiment File", "*.str",
			     drvui->Cur_File);

	if (newfile) {
	    strcpy (drvui->Cur_File, newfile);
	    strcpy (drvui->Cur_Dir, newfile);
	    for (i = strlen (drvui->Cur_Dir); i > 0; --i) {	// Find final / in file name
		if (drvui->Cur_Dir[i - 1] == '/') {
		    drvui->Cur_Dir[i] = 0;
		    break;
		}
	    }
	    for (i = strlen (drvui->Cur_Dir); newfile[i] != 0; i++) {
		drvui->Cur_File[k++] = newfile[i];	// copy file name to Cur_File
	    }
	    drvui->Cur_File[k] = 0;
#endif
	    chdir (drvui->Cur_Dir);	// switch to new directory
	    if ((file = fopen (drvui->Cur_File, "r"))) {
		if (!fl_choice ("Selected file exists. Do you wish to overwrite it?",
				"No", "Yes", NULL)) {
		    return;
		}
	    }
	    WriteConfig ();
	}
    }				// end of *action > 0
    if (Edit_Str_Type || *action > 0) {
	file = fopen (drvui->Cur_File, "w");	// open new or existing file for writing
	if (!file) {
	    Error_Box ("Unable to open STR file for writing");
	    return;
	}
    } else {
	file = fopen (drvui->Cur_Temp, "w");
	if (!file) {
	    Error_Box ("Unable to open TMP file for writing");
	    return;
	}
    }

    selection = textbuf->text ();
    fprintf (file, "%s", selection);	// copy text from window to file
    free (selection);
    fclose (file);		// close the file
    Edit_changed = 0;
    if (Edit_Str_Type) {
	Process_Inp (2);
    } else {
	int l,m;
	drvui->fpin = fopen (drvui->Cur_Temp, "r");	// set up to process tmp file
	natom = 0;
	Omit->nomits = 0;
	drvui->frame_no = 1;
	drvui->modulated = 0;
	drvui->no_mod_vectors = 0;
	drvui->no_site_displace = 0;
	drvui->no_site_U_terms = 0;
	drvui->no_site_occ = 0;
	drvui->no_cell_vec = 0;

	drvui->sys=0;               /* reset crystal system/symmetry and */
	drvui->no_subsys = 1;       /* initialize the subsystem variables */
	for (l = 0; l < 3; l++) {
	    for (m = 0; m < 3; m++)
		drvui->subsys_fact[0][l][m] = 0.0f;
	    drvui->subsys_fact[0][l][l] = 1.0f;     /* set to identity */
	}
	drvui->subsys_ref_volume = 1.0f;
	drvui->subsys_vol[0] = 1.0;
	ShowMapLegend = 0;
	read_inp (2);
	make_bmat (drvui->sys, drvui->lat_con, drvui->b_mat, drvui->ginv, drvui->rec_lat_con);

	fclose (drvui->fpin);
    }
    Generate_Drawing (0);
    if (*action >= 0) {
	drvui->destroy |= TEXT1;
	textwindow->hide ();
    }
    Fl::redraw ();
}

void
Lone_Pair_Combo_cb (Fl_Widget *, void *)
{
// load other widgets when the atom in the combo box is changed
    if (strlen (LonePairs->Height->value ()) == 0)
	return;
    if (strlen (LonePairs->Radius1->value ()) == 0)
	return;
    if (strlen (LonePairs->Radius2->value ()) == 0)
	return;
    if (strlen (LonePairs->Number->value ()) == 0)
	return;
    if (strlen (LonePairs->LonePair_Color->value ()) == 0)
	return;
    LonePairs->LonePair_Add->activate ();
}

void
LonePair_Frame_Combo_cb (Fl_Widget *, void *)
{
// routine called when the frame number is changed on the Bonds Edit Screen
    int i, j;

    char atom1[5];

    char widget[16382];

    char atoms[100][5];

    char color[40];

    int no;

    float radius1, radius2;

    float d1;

    char string[100];

    int Frame_No = 1;

    if (drvui->max_frame > 1)
	Frame_No = atoi (LonePairs->Frame_No->value ());
    int nlist = Get_Unique_Atoms (atoms, Frame_No);	// get unique atom names in sorted list

    LonePairs->LonePair_Combo->list.clear ();	// clear out the old names
    for (j = 0; j < nlist; j++) {	// add atom names in alphabetic order
	LonePairs->LonePair_Combo->list.add (atoms[j]);
    }
    LonePairs->LonePair_Combo->value (atoms[0]);	// load first one in window
    Lone_Pair_Combo_cb (NULL, NULL);
    widget[0] = '\0';
    for (i = 1; i < drvui->ncone; i++) {
	if (drvui->cones[i].cone_fn != Frame_No)
	    continue;		//skip if not for this frame
	strncpy (atom1, drvui->cones[i].cone_l1, 4);
	atom1[4] = 0;
	for (j = 3; j >= 0; --j) {
	    if (atom1[j] == ' ')
		atom1[j] = 0;
	}
	no = drvui->cones[i].numlonepairs;
	radius1 = drvui->cones[i].cone_min;
	radius2 = drvui->cones[i].cone_max;
	d1 = drvui->cones[i].cone_height;
	strncpy (color, drvui->cones[i].col_cone, 40);
	trim_string (color, 40);
	if (!strlen (color))
	    strcpy (color, "Gray20");
	if (color[strlen (color) - 1] < ' ')
	    color[strlen (color) - 1] = 0;
	sprintf (string, "  %4s     %2d%9.3f%9.3f%9.3f     %s\n", atom1, no, d1,
		 radius1, radius2, color);
	strcat (widget, string);
    }
    LonePairs->LonePairBuffer->text (widget);
    LonePairs->LonePairInst1->show ();
    LonePairs->LonePairInst->hide ();
}

void
Main_Frame_Combo_cb (Fl_Widget *, void *)
{
// callback routine when Frame Number is changed on the main page
    int frame = atoi (drvui->Frame_No->value ());

    drvui->X_Min->value (drvui->frames[frame].cryst_lim[0]);
    drvui->X_Max->value (drvui->frames[frame].cryst_lim[3]);
    drvui->Y_Min->value (drvui->frames[frame].cryst_lim[1]);
    drvui->Y_Max->value (drvui->frames[frame].cryst_lim[4]);
    drvui->Z_Min->value (drvui->frames[frame].cryst_lim[2]);
    drvui->Z_Max->value (drvui->frames[frame].cryst_lim[5]);
    if (clipflag) {
	drvui->X_Min_clip->value (drvui->frames[frame].clip_lim[0]);
	drvui->Y_Min_clip->value (drvui->frames[frame].clip_lim[1]);
	drvui->Z_Min_clip->value (drvui->frames[frame].clip_lim[2]);
	drvui->X_Max_clip->value (drvui->frames[frame].clip_lim[3]);
	drvui->Y_Max_clip->value (drvui->frames[frame].clip_lim[4]);
	drvui->Z_Max_clip->value (drvui->frames[frame].clip_lim[5]);
    }
}

void
MapType_cb (void)
{
    if (!strcmp (Maps->MapType->value (), "CIF FoFc - fcf")
	|| !strcmp (Maps->MapType->value (), "JANA FoFc - m80")) {
	Maps->MapCalc->show ();
	Maps->MapCalcType->show ();
	Maps->Resolution->show ();
    } else {
	Maps->MapCalc->hide ();
	Maps->MapCalcType->hide ();
	Maps->Resolution->hide ();
    }
}

void
Map_Info_cb (void)
{
    char string[600];

    char buf_string[16384];

    static int six = 6, zero = 0;

// callback to display map header information
    if (helpwindow6) {
//        helpwindow6->show();
//        return;
	Fl::delete_widget (helpwindow6);
    }

    if (Map_Info.info_valid == 0)	// no map yet, so read it
	Edit_Maps_Save_cb ((Fl_Button *) zero, &zero);
    if (Map_Info.info_valid == 0)	// if still no data, return
	return;
    helpwindow6 = new Fl_Window (200, 180, 620, 360, "Map Info");
    helpwindow6->begin ();
    helpwindow6->callback ((Fl_Callback *) View_Help_Close_cb, &six);
    int y = 40;

    Fl_Text_Display *box1 =
	new Fl_Text_Display (0, y, 600, 0, "DRAWxtl V5.5 Map File Information");
    box1->box (FL_NO_BOX);
    box1->labelsize (24);
    box1->labelcolor ((Fl_Color) 1);
    y += 20;
    strcpy (string, "Map Title: ");
    strcat (string, Map_Info.title);
    Fl_Multiline_Output *box2 = new Fl_Multiline_Output (20, y, 580, 280, "");

    box2->box (FL_NO_BOX);
    box2->labelsize (14);
    box2->labelcolor ((Fl_Color) 186);
    box2->cursor_color (FL_BACKGROUND_COLOR);
    box2->color (FL_BACKGROUND_COLOR);
    strcpy (buf_string, string);
    strcat (buf_string, "\n\n");
    sprintf (string, "Map Cell Parameters: %7.3f %7.3f %7.3f %7.2f %7.2f %7.2f\n\n",
	     Map_Info.lat_con[0], Map_Info.lat_con[1], Map_Info.lat_con[2],
	     Map_Info.lat_con[3], Map_Info.lat_con[4], Map_Info.lat_con[5]);
    strcat (buf_string, string);
    sprintf (string, "Minimum and Maximum density: %8.1f %8.1f\n\n", Map_Info.rhomn,
	     Map_Info.rhomx);
    strcat (buf_string, string);
    sprintf (string, "Map Steps parallel to a, b, c: %5d %5d %5d\n\nSteps per A: %i\n\n",
	     Map_Info.map_int[0], Map_Info.map_int[1], Map_Info.map_int[2], Map_Info.res);
    strcat (buf_string, string);
    if (!drvui->modulated) {
	sprintf (string, "Map Calculated Region:            Axis     Min       Max\n"
		 "                                                    x: %8.3f %8.3f\n"
		 "                                                    y: %8.3f %8.3f\n"
		 "                                                    z: %8.3f %8.3f\n",
		 Map_Info.xlim[0], Map_Info.xlim[1],
		 Map_Info.ylim[0], Map_Info.ylim[1], Map_Info.zlim[0], Map_Info.zlim[1]);
    } else {
	sprintf (string,
		 "Map Calculated Region:            Axis     Min       Max     Axis    Min       Max\n"
		 "                                                    x: %8.3f %8.3f        x4: %8.3f %8.3f\n"
		 "                                                    y: %8.3f %8.3f        x5: %8.3f %8.3f\n"
		 "                                                    z: %8.3f %8.3f        x6: %8.3f %8.3f\n",
		 Map_Info.xlim[0], Map_Info.xlim[1], Map_Info.x4lim[0], Map_Info.x4lim[1],
		 Map_Info.ylim[0], Map_Info.ylim[1], Map_Info.x5lim[0], Map_Info.x5lim[1],
		 Map_Info.zlim[0], Map_Info.zlim[1], Map_Info.x6lim[0],
		 Map_Info.x6lim[1]);
    }
    strcat (buf_string, string);
    y += 250;
    Fl_Button *o = new Fl_Button (270, y, 80, 30, "Close");

    o->callback ((Fl_Callback *) View_Help_Close_cb, &six);
#if !defined (WIN32) && !defined (__APPLE__)
    helpwindow6->icon ((char *) drvui->icon);
#endif
    box2->value (buf_string);
    helpwindow6->end ();
    helpwindow6->show ();
}

void
Modify_Arrow_cb (Fl_Widget *, void *)
{
    char value[8], color[40];

    float length, diam;

    float pos[3], comp[3];

    int start, end;

    const char *selection;

    if (!arrows->ArrowBuffer->selected ()) {
	arrows->AddButton->deactivate ();
	arrows->RemoveButton->deactivate ();
	return;
    }
    memset (color, 0, 40);
    arrows->ArrowBuffer->selection_position (&start, &end);
    selection = arrows->ArrowBuffer->line_text (start);
    arrows->ArrowInstr1->hide ();
    arrows->ArrowInstr->show ();
    sscanf (selection, "%f %f %f %f %f %f %f %f %s", &pos[0], &pos[1], &pos[2],
	    &comp[0], &comp[1], &comp[2], &length, &diam, color);
    free ((char *) selection);
    if (strlen (color) && isalpha (color[0])) {
	sprintf (value, "%4.3f", pos[0]);
	arrows->Px->value (value);
	sprintf (value, "%4.3f", pos[1]);
	arrows->Py->value (value);
	sprintf (value, "%4.3f", pos[2]);
	arrows->Pz->value (value);
	sprintf (value, "%4.3f", comp[0]);
	arrows->Cx->value (value);
	sprintf (value, "%4.3f", comp[1]);
	arrows->Cy->value (value);
	sprintf (value, "%4.3f", comp[2]);
	arrows->Cz->value (value);
	sprintf (value, "%4.3f", length);
	arrows->Length->value (value);
	sprintf (value, "%4.3f", diam);
	arrows->Diameter->value (value);
	arrows->Color->value (color);
	arrows->AddButton->activate ();
	arrows->RemoveButton->activate ();
	New_Arrow_Input_cb (NULL, NULL);
    }
    return;
}

void
Modify_Bonds_cb (Fl_Widget *, void *)
{
    char from[10], color[40], type[10], to[10], value1[20];

    float mind, maxd, diam;

    int i, start, end;

    int itype;

    const char *selection;

    if (!Bonds->BondBuffer->selected ()) {
	Bonds->New_Bond_Add->deactivate ();
	Bonds->New_Bond_Remove->deactivate ();
	return;
    }
    memset (color, 0, 40);
    Bonds->BondBuffer->selection_position (&start, &end);
    selection = Bonds->BondBuffer->line_text (start);
    Bonds->BondInstr1->hide ();
    Bonds->BondInstr->show ();
    (void) sscanf (selection, "%s %s %s %f %f %f %s", type, from, to, &diam, &mind, &maxd,
		   color);
    itype = 0;
    if (!strncmp (type, "dash", 4)) {
	i = sscanf (selection, "%s %d %s %s %f %f %f %s", type, &itype, from, to, &diam, &mind, &maxd,
		   color);
	if ( i < 4 ) {
	    itype = 5;
	    (void) sscanf (selection, "%s %s %s %f %f %f %s", type, from, to, &diam, &mind, &maxd,
		   color);
	}
    }
    free ((char *) selection);
    Bonds->New_Bond_From->value (from);
    Bonds->New_Bond_To->value (to);
    sprintf (value1, "%2.3f", mind);
    Bonds->New_Bond_Min->value (value1);
    sprintf (value1, "%2.3f", maxd);
    Bonds->New_Bond_Max->value (value1);
    sprintf (value1, "%2.3f", diam);
    Bonds->New_Bond_Dia->value (value1);
    Bonds->New_Bond_Color->value (color);
    if (itype > 0) {
	Bonds->New_Bond_Style->set ();
	if (itype != 5) {
	    sprintf (value1, "%3d", itype);
	    Bonds->New_Bond_Dashes->value(value1);
	}
    }
    Bonds->New_Bond_Add->activate ();
    Bonds->New_Bond_Remove->activate ();
}

void
Modify_Bonds_Distance_cb (Fl_Widget *, void *)
{
    char atom[10], value[8];

    float dist;

    int start, end;

    const char *selection;

    if (!Bonds->Bond_Output_Buffer->selected ()) {
	Bonds->New_Bond_To->value ("");
	Bonds->New_Bond_Add->deactivate ();
	Bonds->New_Bond_Remove->deactivate ();
	return;
    }
    memset (atom, 0, 10);
    Bonds->Bond_Output_Buffer->selection_position (&start, &end);
    selection = Bonds->Bond_Output_Buffer->line_text (start);
    sscanf (selection, "%s %f", atom, &dist);
    free ((char *) selection);
    if (strlen (atom) && isalpha (atom[0])) {
	Bonds->New_Bond_To->value (atom);
	sprintf (value, "%2.3f", .1);
	Bonds->New_Bond_Dia->value (value);
	sprintf (value, "%2.3f", dist - .1);
	Bonds->New_Bond_Min->value (value);
	sprintf (value, "%2.3f", dist + .1);
	Bonds->New_Bond_Max->value (value);
	Bonds->New_Bond_Color->value ("Red");
	Bonds->New_Bond_Add->activate ();
	Bonds->New_Bond_Remove->activate ();
	New_Bond_Input_cb (NULL, NULL);
    } else {
	Bonds->New_Bond_To->value ("");
	Bonds->New_Bond_Add->deactivate ();
	Bonds->New_Bond_Remove->deactivate ();
    }
    return;
}

void
Modify_Ellipsoids_cb (Fl_Widget *, void *)
{
    char atom[10], color[40], string[10];

    int no, i;

    int start, end;

    const char *selection;

    if (!ellipsoids->ColorInputBuf->selected ()) {
	ellipsoids->Atom_Combo->value ("");
	ellipsoids->Color_Combo->value ("");
	ellipsoids->New_Ellipsoid_Add->deactivate ();
	ellipsoids->New_Ellipsoid_Convert->deactivate ();
	return;
    }
    memset (atom, 0, 10);
    memset (color, 0, 40);
    ellipsoids->ColorInputBuf->selection_position (&start, &end);
    selection = ellipsoids->ColorInputBuf->line_text (start);
    sscanf (selection, "%s %s %s", atom, string, color);
    free ((char *) selection);
    while (strlen (atom) < 4) {
	strcat (atom, " ");
    }
    if (strlen (atom) && isalpha (atom[0])) {
	if (strstr (string, "*"))
	    sprintf (string, "%4s *", atom);
	else {
	    sscanf (string, "%d", &no);
	    sprintf (string, "%4s%2d", atom, no);
	}
	ellipsoids->Atom_Combo->value (string);
	ellipsoids->Color_Combo->value (color);
	ellipsoids->New_Ellipsoid_Add->activate ();
	ellipsoids->New_Ellipsoid_Remove->activate ();
	set_tf_status ();	//  set the appropriate flags
	for (i = 0; i < natom; i++) {
	    if (drvui->atoms[i].atom_fn != drvui->frame_no)
		continue;
	    if (check_atom_name (drvui->atoms[i].atom_l, atom)
		&& drvui->atoms[i].TF_status == 1) {
		ellipsoids->New_Ellipsoid_Convert->activate ();
		break;
	    }
	}
	New_Ellipsoid_Input_cb (NULL, NULL);
    }
    return;
}

void
Modify_LonePair_cb (Fl_Widget *, void *)
{
    char atom[10], value[8], color[40];

    float d1, d2, height;

    int start, end, no;

    const char *selection;

    if (!LonePairs->LonePairBuffer->selected ()) {
	LonePairs->LonePair_Combo->value ("");
	return;
    }
    memset (atom, 0, 10);
    memset (color, 0, 40);

    LonePairs->LonePairBuffer->selection_position (&start, &end);
    selection = LonePairs->LonePairBuffer->line_text (start);
    LonePairs->LonePairInst1->hide ();
    LonePairs->LonePairInst->show ();
    sscanf (selection, " %s %d %f %f %f %s", atom, &no, &height, &d1, &d2, color);
    free ((char *) selection);
    if (strlen (atom) && isalpha (atom[0])) {
	LonePairs->LonePair_Combo->value (atom);
	sprintf (value, "%2.3f", d1);
	LonePairs->Radius1->value (value);
	sprintf (value, "%2.3f", d2);
	LonePairs->Radius2->value (value);
	sprintf (value, "%2.3f", height);
	LonePairs->Height->value (value);
	LonePairs->LonePair_Color->value (color);
	sprintf (value, "%d", no);
	LonePairs->Number->value (value);
	LonePairs->LonePair_Add->activate ();
	LonePairs->LonePair_Remove->activate ();
    } else {
	LonePairs->LonePair_Combo->value ("");
	LonePairs->Radius1->value ("");
	LonePairs->Radius2->value ("");
	LonePairs->Height->value ("");
	LonePairs->LonePair_Color->value ("");
	LonePairs->Number->value ("");
	LonePairs->LonePair_Add->deactivate ();
	LonePairs->LonePair_Remove->deactivate ();
    }
}

void
Modify_Maps_cb (Fl_Widget *, void *)
{
    char type[8], color[40], value[10];

    float level;

    float step, top;

    int start, end;

    const char *selection;

    if (!Maps->MapsBuffer->selected ()) {
	return;
    }
    memset (type, 0, 8);
    memset (color, 0, 40);
    Maps->MapsBuffer->selection_position (&start, &end);
    selection = Maps->MapsBuffer->line_text (start);
    Maps->MapsInstr1->hide ();
    Maps->MapsInstr->show ();
    if (!drvui->Fourier2d)
	sscanf (selection, " %f %s %s", &level, type, color);
    else
	sscanf (selection, " %f %f %f %s", &level, &step, &top, color);
    free ((char *) selection);
    if (strlen (color) && isalpha (color[0])) {
	sprintf (value, "%2.3f", level);
	Maps->Level->value (value);
	if (drvui->Fourier2d) {
	    sprintf (value, "%2.3f", step);
	    Maps->Step->value (value);
	    sprintf (value, "%2.3f", top);
	    Maps->Top->value (value);
	} else {
	    Maps->Type->value (type);
	}
	Maps->Color->value (color);
	Maps->Add_Button->activate ();
	Maps->Remove_Button->activate ();
    } else {
	Maps->Color->value ("");
	if (drvui->Fourier2d) {
	    Maps->Step->value ("");
	    Maps->Top->value ("");
	} else
	    Maps->Type->value ("");
	Maps->Level->value ("");
	Maps->Add_Button->activate ();
	Maps->Remove_Button->deactivate ();
    }
}

void
Modify_Polyhedra_cb (Fl_Widget *, void *)
{
    char atom[10], color[40], type[10], to[10], value1[20];

    float mind, maxd, transp;

    int start, end;

    int itype = 0;

    const char *selection;

    if (!Polyhedra->PolyhedraBuffer->selected ()) {
	Polyhedra->PolyInstr->hide ();
	Polyhedra->New_Polyhedra_Add->deactivate ();
	Polyhedra->New_Polyhedra_Remove->deactivate ();
	return;
    }
    memset (atom, 0, 10);
    memset (color, 0, 40);
    memset (to, 0, 10);
    transp = 0.0f;
    mind = 0.0f;
    maxd = 0.0f;
    Polyhedra->PolyhedraBuffer->selection_position (&start, &end);
    selection = Polyhedra->PolyhedraBuffer->line_text (start);
    Polyhedra->PolyInstr1->hide ();
    Polyhedra->PolyInstr->show ();
    sscanf (selection, "%s", type);
    if (!strncmp (type, "PS", 2)) {
	(void) sscanf (selection, "%s %s %f %f %s %*s %f", type, atom, &mind, &maxd,
		       color, &transp);
	if (isdigit (color[0]) || color[0] == '.') {
	    char color1[10], color2[10], color3[10];

	    (void) sscanf (selection, "%s %s %f %f %s %s %s %*s %f", type, atom, &mind,
			   &maxd, color1, color2, color3, &transp);
	    strcpy (color, color1);
	    strcat (color, " ");
	    strcat (color, color2);
	    strcat (color, " ");
	    strcat (color, color3);
	}
	itype = 0;
    } else if (!strncmp (type, "PV", 2)) {
	(void) sscanf (selection, "%s %s %s %f %f %s %*s %f", type, atom, to, &mind,
		       &maxd, color, &transp);
	if (isdigit (color[0]) || color[0] == '.') {
	    char color1[10], color2[10], color3[10];

	    (void) sscanf (selection, "%s %s %s %f %f %s %s %s %*s %f", type, atom, to,
			   &mind, &maxd, color1, color2, color3, &transp);
	    strcpy (color, color1);
	    strcat (color, " ");
	    strcat (color, color2);
	    strcat (color, " ");
	    strcat (color, color3);
	}
	itype = 1;
    } else if (!strncmp (type, "SH", 2)) {
	(void) sscanf (selection, "%s %s %f %f %s %*s %f", type, atom, &mind, &maxd,
		       color, &transp);
	if (isdigit (color[0]) || color[0] == '.') {
	    char color1[10], color2[10], color3[10];

	    (void) sscanf (selection, "%s %s %f %f %s %s %s %*s %f", type, atom, &mind,
			   &maxd, color1, color2, color3, &transp);
	    strcpy (color, color1);
	    strcat (color, " ");
	    strcat (color, color2);
	    strcat (color, " ");
	    strcat (color, color3);
	}
	itype = 2;
    } else if (!strncmp (type, "PL", 2)) {
	(void) sscanf (selection, "%s %s %f %f %s %*s %f", type, atom, &mind, &maxd,
		       color, &transp);
	if (isdigit (color[0]) || color[0] == '.') {
	    char color1[10], color2[10], color3[10];

	    (void) sscanf (selection, "%s %s %f %f %s %s %s %*s %f", type, atom, &mind,
			   &maxd, color1, color2, color3, &transp);
	    strcpy (color, color1);
	    strcat (color, " ");
	    strcat (color, color2);
	    strcat (color, " ");
	    strcat (color, color3);
	}
	itype = 3;
    }
    free ((char *) selection);
    Polyhedra->New_Polyhedra_From->value (atom);
    Polyhedra->Polyhedra_Combo->value (atom);
    Polyhedra_Combo_cb (NULL, NULL);
    Polyhedra->New_Polyhedra_To->value (to);
    sprintf (value1, "%2.3f", mind);
    Polyhedra->New_Polyhedra_Min->value (value1);
    sprintf (value1, "%2.3f", maxd);
    Polyhedra->New_Polyhedra_Max->value (value1);
    Polyhedra->New_Polyhedra_Color->value (color);
    if (transp > 0.) {
	sprintf (value1, "%.3f", transp);
	Polyhedra->New_Polyhedra_Transp->value (value1);
    }
    switch (itype) {
    case 0:
	Polyhedra->Polysz->set ();
	Polyhedra->Polyvert->clear ();
	Polyhedra->Polyshell->clear ();
	Polyhedra->Plane->clear ();
	break;
    case 1:
	Polyhedra->Polyvert->set ();
	Polyhedra->Polysz->clear ();
	Polyhedra->Polyshell->clear ();
	Polyhedra->Plane->clear ();
	break;
    case 2:
	Polyhedra->Polyvert->clear ();
	Polyhedra->Polyshell->set ();
	Polyhedra->Polysz->clear ();
	Polyhedra->Plane->clear ();
	break;
    case 3:
	Polyhedra->Plane->set ();
	Polyhedra->Polyvert->clear ();
	Polyhedra->Polysz->clear ();
	Polyhedra->Polyshell->clear ();
    }
    Polyhedra->New_Polyhedra_Add->activate ();
    Polyhedra->New_Polyhedra_Remove->activate ();
    return;
}

void
Modify_Polyhedra_Distance_cb (Fl_Widget *, void *)
{
// callback routine entered when a selection is made in the polyhedra bond-distance table
    char atom[10], value[8];

    float dist;

    int start, end;

    const char *selection;

    if (!Polyhedra->Polyhedra_Output_Buffer->selected ()) {
	Polyhedra->New_Polyhedra_To->value ("");
	Polyhedra->New_Polyhedra_Add->deactivate ();
	return;
    }
    memset (atom, 0, 10);
    Polyhedra->Polyhedra_Output_Buffer->selection_position (&start, &end);
    selection = Polyhedra->Polyhedra_Output_Buffer->line_text (start);
    sscanf (selection, "%s %f", atom, &dist);
    if (strlen (atom) && isalpha (atom[0])) {
	Polyhedra->New_Polyhedra_To->value (atom);
	sprintf (value, "%2.3f", dist - .2);
	Polyhedra->New_Polyhedra_Min->value (value);
	sprintf (value, "%2.3f", dist + .2);
	Polyhedra->New_Polyhedra_Max->value (value);
	Polyhedra->New_Polyhedra_Color->value ("Red");
	Polyhedra->New_Polyhedra_Add->activate ();
	New_Polyhedra_Input_cb (NULL, NULL);
    } else {
	Polyhedra->New_Polyhedra_To->value ("");
	Polyhedra->New_Polyhedra_Add->deactivate ();
	Polyhedra->New_Polyhedra_Remove->deactivate ();
    }
    free ((char *) selection);
    return;
}

void
Modify_Spheres_cb (Fl_Widget *, void *)
{
    char atom[10], value[8], color[40], string[40];

    float dist;

    int start, end, i, n;

    const char *selection;

    if (!Spheres->SphereBuffer->selected ()) {
	Spheres->Sphere_Combo->value ("");
	Spheres->New_Sphere_Size->value ("");
	Spheres->New_Sphere_Add->deactivate ();
	Spheres->New_Sphere_Convert->deactivate ();
	return;
    }
    memset (atom, 0, 10);
    memset (color, 0, 40);
    Spheres->SphereBuffer->selection_position (&start, &end);
    selection = Spheres->SphereBuffer->line_text (start);
    Spheres->SphereInstr1->hide ();
    Spheres->SphereInstr->show ();
    if (strlen (selection) == 0)
	return;			// user clicked on empty line by mistake
    sscanf (selection, "%s %s %f %s", atom, string, &dist, color);
    while (strlen (atom) < 4)
	strcat (atom, " ");
    free ((char *) selection);
    if (strstr (string, "*"))
	n = -1;
    else
	(void) sscanf (string, "%d", &n);
    if (strlen (atom) && isalpha (atom[0])) {
	strcpy (string, atom);
	if (n < 0)
	    strcpy (value, " *");
	else
	    sprintf (value, "%2d", n);
	strcat (string, value);
	Spheres->Sphere_Combo->value (string);
	sprintf (value, "%2.3f", dist);
	Spheres->New_Sphere_Size->value (value);
	Spheres->New_Sphere_Color->value (color);
	Spheres->New_Sphere_Add->activate ();
	Spheres->New_Sphere_Remove->activate ();
	for (i = 1; i < drvui->n_ellips; i++) {
	    if (check_atom_name (drvui->ellips[i].ellips_l, atom)
		&& (drvui->ellips[i].ellips_n == n || n == -1)
		&& drvui->ellips[i].ell_type == 1)
		Spheres->New_Sphere_Convert->activate ();
	}
	New_Sphere_Input_cb (NULL, NULL);
    }
    return;
}

void
New_Arrow_Add_cb (class Fl_Widget *, int *action)
{
// Add new arrow to display list on Arrow Edit Screen
    char string[100], color[40];

    int start, end;

    char *selection = NULL;

    arrows->ArrowInstr2->show ();
    if (arrows->ArrowBuffer->selected ()) {
	arrows->ArrowBuffer->selection_position (&start, &end);
	selection = arrows->ArrowBuffer->line_text (start);
	if (!*action)
	    arrows->ArrowBuffer->remove (arrows->ArrowBuffer->line_start (start),
					 arrows->ArrowBuffer->line_end (end) + 1);
    }

    if (*action) {
	strcpy (color, arrows->Color->value ());
	trim_string (color, 40);
	if (!strlen (color))
	    strcpy (color, "Gray20");
	sprintf (string, "%6.3f %6.3f %6.3f %10.3f %6.3f %6.3f %10.3f %6.3f   %s\n",
		 atof (arrows->Px->value ()), atof (arrows->Py->value ()),
		 atof (arrows->Pz->value ()), atof (arrows->Cx->value ()),
		 atof (arrows->Cy->value ()), atof (arrows->Cz->value ()),
		 atof (arrows->Length->value ()), atof (arrows->Diameter->value ()),
		 color);
	if (selection)
	    arrows->ArrowBuffer->replace (arrows->ArrowBuffer->line_start (start),
					  arrows->ArrowBuffer->line_end (end) + 1,
					  string);
	else
	    arrows->ArrowBuffer->append (string);
    }
    if (selection)
	free ((char *) selection);
    arrows->Px->value ("");
    arrows->Py->value ("");
    arrows->Pz->value ("");
    arrows->Cx->value ("");
    arrows->Cy->value ("");
    arrows->Cz->value ("");
    arrows->Color->value ("");
    arrows->Diameter->value ("");
    arrows->Length->value ("");
    arrows->AddButton->deactivate ();
    arrows->RemoveButton->deactivate ();
    arrows->ArrowInstr->hide ();
    arrows->ArrowInstr1->show ();
    Edit_Arrow_Save_cb (NULL, &three);
}

void
New_Arrow_Input_cb (class Fl_Widget *, void *)
{
// callback routine to make 'Add' new bond button active whenever all fields are non-blank
    if (!strlen (arrows->Px->value ()))
	return;
    if (!strlen (arrows->Py->value ()))
	return;
    if (!strlen (arrows->Pz->value ()))
	return;
    if (!strlen (arrows->Cx->value ()))
	return;
    if (!strlen (arrows->Cy->value ()))
	return;
    if (!strlen (arrows->Cz->value ()))
	return;
    if (!strlen (arrows->Length->value ()))
	return;
    if (!strlen (arrows->Diameter->value ()))
	return;
    if (!strlen (arrows->Color->value ()))
	return;
    arrows->AddButton->activate ();
}

void
New_Bond_Add_cb (class Fl_Widget *, int *action)
{
// Add new bond to display list on Bond Edit Screen
    char string[100], type[5];

    float dia, dmin, dmax;

    int start, end, numd;

    char *selection = NULL;

    Bonds->BondInstr2->show ();
    if (Bonds->BondBuffer->selected ()) {
	Bonds->BondBuffer->selection_position (&start, &end);
	selection = Bonds->BondBuffer->line_text (start);
	if (*action)
	    Bonds->BondBuffer->remove (Bonds->BondBuffer->line_start (start),
				       Bonds->BondBuffer->line_end (end) + 1);
    }
    if (!*action) {
	dia = (float) atof (Bonds->New_Bond_Dia->value ());
	dmin = (float) atof (Bonds->New_Bond_Min->value ());
	dmax = (float) atof (Bonds->New_Bond_Max->value ());
	if (Bonds->New_Bond_Style->value ()) {
	    strcpy (type, "dash");
	    numd = atoi (Bonds->New_Bond_Dashes->value ());
	    if (numd > 0) 
		sprintf (string, "%4s %3d %4s   %4s  %9.3f %9.3f %9.3f     %s\n",
			type, numd, Bonds->New_Bond_From->value (), Bonds->New_Bond_To->value (), dia,
			dmin, dmax, Bonds->New_Bond_Color->value ());
	    else
		sprintf (string, "%4s     %4s   %4s  %9.3f %9.3f %9.3f     %s\n",
			type, Bonds->New_Bond_From->value (), Bonds->New_Bond_To->value (), dia,
			dmin, dmax, Bonds->New_Bond_Color->value ());
	} else {
	    strcpy (type, "bond");
	    sprintf (string, "%4s     %4s   %4s  %9.3f %9.3f %9.3f     %s\n",
		 type, Bonds->New_Bond_From->value (), Bonds->New_Bond_To->value (), dia,
		 dmin, dmax, Bonds->New_Bond_Color->value ());
	}
	if (selection)
	    Bonds->BondBuffer->replace (Bonds->BondBuffer->line_start (start),
					Bonds->BondBuffer->line_end (end) + 1, string);
	else
	    Bonds->BondBuffer->append (string);
    }
    if (selection)
	free ((char *) selection);
    Bonds->New_Bond_To->value ("");
    Bonds->New_Bond_Dia->value ("");
    Bonds->New_Bond_Min->value ("");
    Bonds->New_Bond_Max->value ("");
    Bonds->New_Bond_Color->value ("");
    Bonds->New_Bond_Style->value (0);
    Bonds->New_Bond_Dashes->value ("");
    Bonds->New_Bond_Add->deactivate ();
    Bonds->New_Bond_Remove->deactivate ();
    Bonds->BondInstr1->show ();
    Bonds->BondInstr->hide ();
    Edit_Bond_Save_cb (NULL, &three);
}

void
New_Ellipsoid_Add_cb (class Fl_Widget *, int *action)
{
    int start, end;

    char string[40], string2[40], atom[5], color[40], num[4];

    int no, i = 0;

    char *selection = NULL;

    ellipsoids->Instr2->show ();
    strcpy (string, ellipsoids->Atom_Combo->value ());
    strcpy (color, ellipsoids->Color_Combo->value ());
    if (ellipsoids->ColorInputBuf->selected ()) {
	ellipsoids->ColorInputBuf->selection_position (&start, &end);
	selection = ellipsoids->ColorInputBuf->line_text (start);
	if (*action != 1)
	    ellipsoids->ColorInputBuf->remove (ellipsoids->ColorInputBuf->
					       line_start (start),
					       ellipsoids->ColorInputBuf->line_end (end) +
					       1);
    }
    ellipsoids->New_Ellipsoid_Add->deactivate ();
    ellipsoids->New_Ellipsoid_Convert->deactivate ();
    ellipsoids->New_Ellipsoid_Remove->deactivate ();
    ellipsoids->Instr1->show ();
    ellipsoids->Instr->hide ();
    if (*action == 1) {		// add
	sscanf (string, "%s %s", atom, num);
	if (strstr (num, "*"))
	    sprintf (string, "%4s   *   %s\n", atom, color);
	else {
	    sscanf (num, "%d", &no);
	    sprintf (string, "%4s %3d   %s\n", atom, no, color);
	}
	if (selection)
	    ellipsoids->ColorInputBuf->replace (ellipsoids->ColorInputBuf->
						line_start (start),
						ellipsoids->ColorInputBuf->
						line_end (end) + 1, string);
	else
	    ellipsoids->ColorInputBuf->append (string);
    } else if (*action == 2) {	// convert selection to sphere
	while (1) {		// delete selected item from combo list
	    if (ellipsoids->Atom_Combo->list.text (++i) == NULL)
		break;
	    strcpy (string2, ellipsoids->Atom_Combo->list.text (i));
	    if (!strcmp (string, string2)) {
		ellipsoids->Atom_Combo->list.remove (i);
		break;
	    }
	}
	if (!selection)
	    return;		// happens when the atom is not in the ellipsoid list,
	// but the user tried to force its conversion anyway
	sscanf (selection, "%s %s %s", atom, num, color);
	if (strstr (num, "*"))
	    no = -1;
	else
	    sscanf (num, "%d", &no);
	while (strlen (atom) < 4)
	    strcat (atom, " ");
	for (i = 1; i < drvui->n_ellips; i++) {
	    if (check_atom_name (atom, drvui->ellips[i].ellips_l) &&
		(no == drvui->ellips[i].ellips_n || no == -1)) {
		drvui->ellips[i].ell_type = 1;	// no longer an ellipse
		drvui->spheres[drvui->nsphere].sphere_fn = 1;	// it goes into frame 1
		strcpy (drvui->spheres[drvui->nsphere].sphere_l, atom);	// same name as ellipse
		if (no == -1)
		    drvui->spheres[drvui->nsphere].sphere_n = no;
		else
		    drvui->spheres[drvui->nsphere].sphere_n = drvui->ellips[i].ellips_n;
		drvui->spheres[drvui->nsphere].sphere_size =
		    (drvui->ellips[i].ellips_RMS[0]
		     + drvui->ellips[i].ellips_RMS[1]
		     + drvui->ellips[i].ellips_RMS[2]) * drvui->Ellipsoid_Scale / 3.0f;	//size
		strcpy (drvui->spheres[drvui->nsphere].sphere_col, color);
		i = drvui->n_ellips;
	    }
	}
	if (!drvui->spheres) {
	    drvui->spheres =
		(struct sphere_struct *) zalloc (20 * sizeof (struct sphere_struct));
	    if (!drvui->spheres) {
		Error_Box ("Unable to allocate storage for spheres.");
		return;
	    }
	}
	if (++drvui->nsphere > drvui->nsphere_alloc) {
	    drvui->spheres =
		(struct sphere_struct *) realloc (drvui->spheres,
						  20 * sizeof (struct sphere_struct));
	    if (!drvui->spheres) {
		Error_Box ("Unable to allocate storage for spheres.");
		exit (0);
	    }
	    drvui->nsphere_alloc += 20;
	}
	drvui->Str_File_Changed = 1;
	Update_Str (0);
	Generate_Drawing (0);
    }
    if (selection)
	free ((char *) selection);
    ellipsoids->Atom_Combo->value ("");
    ellipsoids->Color_Combo->value ("");
    Edit_Ellipsoid_Save_cb (NULL, &three);
}

void
New_Ellipsoid_Input_cb (class Fl_Widget *, void *)
{
// callback routine to make 'Add' new ellipsoid button active whenever all required fields are non-blank
    if (!strlen (ellipsoids->Atom_Combo->value ()))
	return;
    if (!strlen (ellipsoids->Color_Combo->value ()))
	return;
    ellipsoids->New_Ellipsoid_Add->activate ();
    ellipsoids->New_Ellipsoid_Convert->activate ();
}

void
New_Bond_Input_cb (class Fl_Widget *, void *)
{
// callback routine to make 'Add' new bond button active whenever all fields are non-blank
    if (!strlen (Bonds->New_Bond_To->value ()))
	return;
    if (!atof (Bonds->New_Bond_Dia->value ()))
	return;
    if (!atof (Bonds->New_Bond_Min->value ()))
	return;
    if (!atof (Bonds->New_Bond_Max->value ()))
	return;
    if (!strlen (Bonds->New_Bond_Color->value ()))
	return;
    Bonds->New_Bond_Add->activate ();
}

void
New_Map_Add_cb (class Fl_Widget *, int *action)
{
// Add new map contour to display list on Map Edit Screen
    char string[100];

    char value[20];

    float level;

    int start, end;

    char *selection = NULL;

    Maps->MapsInstr2->show ();
    if (Maps->MapsBuffer->selected ()) {
	Maps->MapsBuffer->selection_position (&start, &end);
	selection = Maps->MapsBuffer->line_text (start);
	Maps->MapsBuffer->remove (Maps->MapsBuffer->line_start (start),
				  Maps->MapsBuffer->line_end (end) + 1);
    }
    if (*action) {
	strcpy (value, Maps->Level->value ());
	level = (float) atof (value);
	if (drvui->Fourier2d) {
	    float step, top;

	    strcpy (value, Maps->Step->value ());
	    step = (float) atof (value);
	    strcpy (value, Maps->Top->value ());
	    top = (float) atof (value);
	    sprintf (string, "%7.3f %7.3f %7.3f %s\n", level, step, top,
		     Maps->Color->value ());
	} else {
	    sprintf (string, "%8.3f    %6s     %s\n", level, Maps->Type->value (),
		     Maps->Color->value ());
	}
	Maps->MapsBuffer->append (string);
    }
    if (selection)
	free ((char *) selection);
    Maps->Level->value ("");
    if (drvui->Fourier2d) {
	Maps->Step->value ("");
	Maps->Top->value ("");
    }
    Maps->Add_Button->deactivate ();
    Maps->Remove_Button->deactivate ();
    Maps->MapsInstr1->show ();
    Maps->MapsInstr->hide ();
    Edit_Maps_Save_cb (NULL, &three);
}

void
New_Map_Input_cb (Fl_Widget *, void *)
{
// callback routine to make 'Add' new map level button active whenever all fields are non-blank
    if (!strlen (Maps->Color->value ()))
	return;
    if (!strlen (Maps->Level->value ()))
	return;
    if (drvui->Fourier2d) {
	if (!strlen (Maps->Step->value ()))
	    return;
	if (!strlen (Maps->Top->value ()))
	    return;
    } else if (!strlen (Maps->Type->value ()))
	return;
    Maps->Add_Button->activate ();
}

void
New_Polyhedra_Add_cb (class Fl_Widget *, int *action)
{
// Add new polyhedron to display list on Polyhedra Edit Screen
    char string[100];

    float dmin, dmax, transp;

    char temp[10];

    char atom1[5], atom2[5], type[3] = "PS";

    int start, end, i;

    char *selection = NULL;

    Polyhedra->PolyInstr2->show ();
    if (Polyhedra->PolyhedraBuffer->selected ()) {
	Polyhedra->PolyhedraBuffer->selection_position (&start, &end);
	selection = Polyhedra->PolyhedraBuffer->line_text (start);
	if (!*action)
	    Polyhedra->PolyhedraBuffer->remove (Polyhedra->PolyhedraBuffer->
						line_start (start),
						Polyhedra->PolyhedraBuffer->
						line_end (end) + 1);
    }
    if (*action) {
	strcpy (temp, Polyhedra->New_Polyhedra_Min->value ());
	dmin = (float) atof (temp);
	dmax = (float) atof (Polyhedra->New_Polyhedra_Max->value ());
	transp = (float) atof (Polyhedra->New_Polyhedra_Transp->value ());

	strcpy (atom1, Polyhedra->New_Polyhedra_From->value ());
	strcpy (atom2, Polyhedra->New_Polyhedra_To->value ());
	if (Polyhedra->Polysz->value ()) {
	    strcpy (type, "PS");
	    strcpy (atom2, "");
	    dmin = 0.005f;
	}
	if (Polyhedra->Polyvert->value ()) {
	    strcpy (type, "PV");
	    dmin = 0.005f;
	}
	if (Polyhedra->Polyshell->value ()) {
	    strcpy (type, "SH");
	    strcpy (atom2, "");
	}
	if (Polyhedra->Plane->value ()) {
	    strcpy (type, "PL");
	    strcpy (atom2, "");
	    dmin = 0.005f;
	}
	if (transp > 0.0f)
	    sprintf (string, "  %2s      %4s    %4s %9.3f %9.3f     %s filter %.3f\n",
		     type, atom1, atom2, dmin, dmax,
		     Polyhedra->New_Polyhedra_Color->value (), transp);
	else
	    sprintf (string, "  %2s      %4s    %4s %9.3f %9.3f     %s\n", type, atom1,
		     atom2, dmin, dmax, Polyhedra->New_Polyhedra_Color->value ());
	if (selection)
	    Polyhedra->PolyhedraBuffer->replace (Polyhedra->PolyhedraBuffer->
						 line_start (start),
						 Polyhedra->PolyhedraBuffer->
						 line_end (end) + 1, string);
	else
	    Polyhedra->PolyhedraBuffer->append (string);
	if (Polyhedra->Edge_Color->value()) {
	    int found = 0;
	    for ( i = 0; i<drvui->nedges; i++) {
		if (!strcmp (drvui->polyedges[i].name, atom1)) {
		    strcpy (drvui->polyedges[i].color, Polyhedra->Edge_Color->value());
		    drvui->polyedges[i].radius = (float) atof (Polyhedra->Edge_Radius->value());
		    found = 1;
		    break;
		}
	    }
	    if (found == 0) {
		trim_string (atom1, 4);
		strcpy (drvui->polyedges[drvui->nedges].name, atom1);
		strcpy (drvui->polyedges[drvui->nedges].color, Polyhedra->Edge_Color->value());
		drvui->polyedges[drvui->nedges++].radius = (float) atof (Polyhedra->Edge_Radius->value());
		check_dynamic_storage ();
	    }
	}
    }
    if (selection)
	free ((char *) selection);
    Polyhedra->New_Polyhedra_To->value ("");
    Polyhedra->New_Polyhedra_Min->value ("");
    Polyhedra->New_Polyhedra_Max->value ("");
    Polyhedra->New_Polyhedra_Color->value ("");
    Polyhedra->New_Polyhedra_Transp->value ("");
    Polyhedra->New_Polyhedra_Add->deactivate ();
    Polyhedra->New_Polyhedra_Remove->deactivate ();
    Polyhedra->PolyInstr1->show ();
    Edit_Polyhedra_Save_cb (NULL, &three);
}

void
New_Lone_Pair_Add_cb (class Fl_Widget *, int *action)
{
// Add/remove lone pair in display list on Lone-Pair Edit Screen
    char string[100];

    float height, radius1, radius2;

    int start, end;

    int no;

    char *selection = NULL;

    LonePairs->LonePairInst2->show ();
    if (LonePairs->LonePairBuffer->selected ()) {
	LonePairs->LonePairBuffer->selection_position (&start, &end);
	selection = LonePairs->LonePairBuffer->line_text (start);
	if (*action != 1)
	    LonePairs->LonePairBuffer->remove (LonePairs->LonePairBuffer->
					       line_start (start),
					       LonePairs->LonePairBuffer->line_end (end) +
					       1);
    }
    if (*action == 1) {
	height = (float) atof (LonePairs->Height->value ());
	radius1 = (float) atof (LonePairs->Radius1->value ());
	radius2 = (float) atof (LonePairs->Radius2->value ());
	no = atoi (LonePairs->Number->value ());

	sprintf (string, "  %4s     %2d%9.3f%9.3f%9.3f     %s\n",
		 LonePairs->LonePair_Combo->value (), no, height, radius1, radius2,
		 LonePairs->LonePair_Color->value ());
	if (selection)
	    LonePairs->LonePairBuffer->replace (LonePairs->LonePairBuffer->
						line_start (start),
						LonePairs->LonePairBuffer->
						line_end (end) + 1, string);
	else
	    LonePairs->LonePairBuffer->append (string);
    }
    if (selection)
	free ((char *) selection);
    LonePairs->Number->value ("");
    LonePairs->Height->value ("");
    LonePairs->Radius1->value ("");
    LonePairs->Radius2->value ("");
    LonePairs->LonePair_Color->value ("");
    LonePairs->LonePair_Add->deactivate ();
    LonePairs->LonePair_Remove->deactivate ();
    LonePairs->LonePairInst1->show ();
    LonePairs->LonePairInst->hide ();
    Edit_Lone_Pair_Save_cb (NULL, &three);
}

void
New_Polyhedra_Input_cb (class Fl_Widget *, void *)
{
// callback routine to make 'Add' new polyhedra button active whenever all required fields are non-blank
    if (Polyhedra->Polyvert->value () && !strlen (Polyhedra->New_Polyhedra_To->value ()))
	return;
    if (Polyhedra->Polyshell->value ()
	&& atof (Polyhedra->New_Polyhedra_Min->value ()) == 0.0)
	return;
    if (atof (Polyhedra->New_Polyhedra_Max->value ()) == 0.0)
	return;
    if (!strlen (Polyhedra->New_Polyhedra_Color->value ()))
	return;
    Polyhedra->New_Polyhedra_Add->activate ();
}

void
New_Sphere_Add_cb (class Fl_Widget *, int *action)
{
// * action == 0 - delete a sphere
//          == 1 - add new sphere 
//          == 2 - convert existing sphere to ellipsoid (only if Uij's available)

    char string[100];

    int start, end, i, n;

    float size;

    char atom1[10], num[4], color[40];

    char *selection = NULL;

    size = (float) atof (Spheres->New_Sphere_Size->value ());

    Spheres->SphereInstr2->show ();
    strcpy (atom1, Spheres->Sphere_Combo->value ());
    if (Spheres->SphereBuffer->selected ()) {
	Spheres->SphereBuffer->selection_position (&start, &end);
	selection = Spheres->SphereBuffer->line_text (start);
	if (*action != 1)
	    Spheres->SphereBuffer->remove (Spheres->SphereBuffer->line_start (start),
					   Spheres->SphereBuffer->line_end (end) + 1);
    }
    if (*action == 1) {
	sprintf (string, "%6s %9.3f     %s\n", atom1, size,
		 Spheres->New_Sphere_Color->value ());
	if (selection)
	    Spheres->SphereBuffer->replace (Spheres->SphereBuffer->line_start (start),
					    Spheres->SphereBuffer->line_end (end) + 1,
					    string);
	else
	    Spheres->SphereBuffer->append (string);
    } else if (*action == 2) {
	if (selection && strlen (selection) > 0) {
	    sscanf (selection, "%s %s %f %s", atom1, num, &size, color);
	    if (strstr (num, "*"))
		n = -1;
	    else
		(void) sscanf (num, "%d", &n);
	} else {
	    strcpy (string, atom1);
	    sscanf (string, "%s %d", atom1, &n);
	    if (strstr (string, "*"))
		n = -1;
	    strcpy (color, Spheres->New_Sphere_Color->value ());
	}
	while (strlen (atom1) < 4)
	    strcat (atom1, " ");
	for (i = 1; i < drvui->n_ellips; i++) {
	    if (check_atom_name (atom1, drvui->ellips[i].ellips_l)
		&& drvui->ellips[i].ell_type < 1000 && (drvui->ellips[i].ellips_n == n
							|| n == -1)) {
		drvui->ellips[i].ell_type = 1001;
		strcpy (drvui->ellips[i].ellips_col, color);
		drvui->do_ellipsoids = 1;
	    }
	}
    }
    free ((char *) selection);
    Spheres->New_Sphere_Size->value ("");
    Spheres->New_Sphere_Color->value ("");
    Spheres->Sphere_Combo->value ("");
    Spheres->New_Sphere_Add->deactivate ();
    Spheres->New_Sphere_Remove->deactivate ();
    Spheres->New_Sphere_Convert->deactivate ();
    Spheres->SphereInstr->hide ();
    Spheres->SphereInstr1->show ();
    Edit_Spheres_Save_cb (NULL, &three);
}

void
New_Sphere_Input_cb (class Fl_Widget *, void *)
{
    int i, n;

    char atom[5], number[5];

// callback routine to make 'Add' new sphere button active whenever all required fields are non-blank
    if (atof (Spheres->New_Sphere_Size->value ()) == 0.0)
	return;
    if (!strlen (Spheres->New_Sphere_Color->value ()))
	return;
    Spheres->New_Sphere_Add->activate ();
    (void) sscanf (Spheres->Sphere_Combo->value (), "%s %s", atom, number);
    if (strstr (number, "*")) {
	n = -1;
    } else {
	(void) sscanf (number, "%d", &n);
    }
    while (strlen (atom) < 4)
	strcat (atom, " ");
    for (i = 1; i < drvui->n_ellips; i++) {
	if (check_atom_name (atom, drvui->ellips[i].ellips_l)
	    && (drvui->ellips[i].ellips_n == n || n == -1)
	    && drvui->ellips[i].ell_type == 1) {
	    Spheres->New_Sphere_Convert->activate ();
	}
    }
}

void
next_focus (void)
{
// routine to change focus to next widget when <TAB> is pushed
    Fl_Widget *o = Fl::focus ();	// get pointer to widget with focus

    if (o == drvui->Origin_X) {	// from Origin_X go to Origin_Y
	drvui->Origin_Y->take_focus ();
	drvui->Origin_Y->position (0);	// highlight text
	drvui->Origin_Y->mark (strlen (drvui->Origin_Y->value ()));
	return;
    }
    if (o == drvui->Origin_Y) {
	drvui->Origin_Z->take_focus ();
	drvui->Origin_Z->position (0);
	drvui->Origin_Z->mark (strlen (drvui->Origin_Z->value ()));
	return;
    }
    if (o == drvui->Origin_Z) {
	drvui->X_Rot->take_focus ();
	drvui->X_Rot->position (0);
	drvui->X_Rot->mark (strlen (drvui->X_Rot->value ()));
	return;
    }
    if (o == drvui->X_Rot) {
	drvui->Y_Rot->take_focus ();
	drvui->Y_Rot->position (0);
	drvui->Y_Rot->mark (strlen (drvui->Y_Rot->value ()));
	return;
    }
    if (o == drvui->Y_Rot) {
	drvui->Z_Rot->take_focus ();
	drvui->Z_Rot->position (0);
	drvui->Z_Rot->mark (strlen (drvui->Z_Rot->value ()));
	return;
    }
    if (o == drvui->Z_Rot) {
	drvui->X_Min->take_focus ();
	return;
    }
// set up the tab loop on the parameter edit screen
    if (edtprm) {
	if (o == edtprm->Magnification) {
	    edtprm->List->take_focus ();
	    edtprm->List->position (0);
	    edtprm->List->mark (strlen (edtprm->List->value ()));
	}
	if (o == edtprm->List) {
	    edtprm->DepthCue->take_focus ();
	    edtprm->DepthCue->position (0);
	    edtprm->DepthCue->mark (strlen (edtprm->DepthCue->value ()));
	}
	if (o == edtprm->DepthCue) {
	    edtprm->Mol_Comp_Dist->take_focus ();
	    edtprm->Mol_Comp_Dist->position (0);
	    edtprm->Mol_Comp_Dist->mark (strlen (edtprm->Mol_Comp_Dist->value ()));
	}
	if (o == edtprm->Mol_Comp_Dist) {
	    edtprm->Cell_Edge_Width->take_focus ();
	    edtprm->Cell_Edge_Width->position (0);
	    edtprm->Cell_Edge_Width->mark (strlen (edtprm->Cell_Edge_Width->value ()));
	}
	if (o == edtprm->Cell_Edge_Width) {
	    edtprm->Poly_Limit->take_focus ();
	    edtprm->Poly_Limit->position (0);
	    edtprm->Poly_Limit->mark (strlen (edtprm->Poly_Limit->value ()));
	}
	if (o == edtprm->Poly_Limit) {
	    edtprm->Phong_Refl->take_focus ();
	    edtprm->Phong_Refl->position (0);
	    edtprm->Phong_Refl->mark (strlen (edtprm->Phong_Refl->value ()));
	}
	if (o == edtprm->Phong_Refl) {
	    edtprm->Phong_Size->take_focus ();
	    edtprm->Phong_Size->position (0);
	    edtprm->Phong_Size->mark (strlen (edtprm->Phong_Size->value ()));
	}
	if (o == edtprm->Phong_Size) {
	    edtprm->Ambient_Finish->take_focus ();
	    edtprm->Ambient_Finish->position (0);
	    edtprm->Ambient_Finish->mark (strlen (edtprm->Ambient_Finish->value ()));
	}
	if (o == edtprm->Ambient_Finish) {
	    edtprm->Diffuse_Finish->take_focus ();
	    edtprm->Diffuse_Finish->position (0);
	    edtprm->Diffuse_Finish->mark (strlen (edtprm->Diffuse_Finish->value ()));
	}
	if (o == edtprm->Diffuse_Finish) {
	    edtprm->Specular_Finish->take_focus ();
	    edtprm->Specular_Finish->position (0);
	    edtprm->Specular_Finish->mark (strlen (edtprm->Specular_Finish->value ()));
	}
	if (o == edtprm->Specular_Finish) {
	    edtprm->Finish_Roughness->take_focus ();
	    edtprm->Finish_Roughness->position (0);
	    edtprm->Finish_Roughness->mark (strlen (edtprm->Finish_Roughness->value ()));
	}
	if (o == edtprm->Finish_Roughness) {
	    edtprm->Magnification->take_focus ();
	    edtprm->Magnification->position (0);
	    edtprm->Magnification->mark (strlen (edtprm->Magnification->value ()));
	}
    }
// set up the tab loop on the ellipsoid edit screen
    if (ellipsoids) {
	if (o == ellipsoids->Probability) {
	    ellipsoids->Axis_Width->take_focus ();
	    ellipsoids->Axis_Width->position (0);
	    ellipsoids->Axis_Width->mark (strlen (ellipsoids->Axis_Width->value ()));
	}
	if (o == ellipsoids->Axis_Width) {
	    ellipsoids->Axis_Color->take_focus ();
//            ellipsoids->Axis_Color->position(0);
//            ellipsoids->Axis_Color->mark(strlen(ellipsoids->Axis_Color->value()));
	}
	if (o == ellipsoids->Axis_Color) {
	    ellipsoids->Cutout_Color->take_focus ();
//            ellipsoids->Cutout_Color->position(0);
//            ellipsoids->Cutout_Color->mark(strlen(ellipsoids->Cutout_Color->value()));
	}
	if (o == ellipsoids->Cutout_Color) {
	    ellipsoids->Probability->take_focus ();
	    ellipsoids->Probability->position (0);
	    ellipsoids->Probability->mark (strlen (ellipsoids->Probability->value ()));
	}
    }
}

void
Polyhedra_Combo_cb (Fl_Widget *, void *)
{
// update bond distance table when the atom in the combo box is changed
    const char *atom;
    char value[10];
    int i;

    if (!drvui->table)
	drvui->table = (char *) zalloc (20480 * sizeof (char));
    atom = Polyhedra->Polyhedra_Combo->value ();	// get name from combo box
    Polyhedra->New_Polyhedra_From->value (atom);	// put it in the 'From' location
    Load_Bond_Data (atom, drvui->table);
    Polyhedra->Polyhedra_Output_Buffer->text (drvui->table);
    for (i = 0; i<drvui->nedges; i++) {
	if (!strcmp (atom, drvui->polyedges[i].name)) {
		Polyhedra->Edge_Color->value( (const char*) drvui->polyedges[i].color);
		sprintf (value, "%.3f", drvui->polyedges[i].radius);
		Polyhedra->Edge_Radius->value(value);
		return;
	}
    }
    Polyhedra->Edge_Color->value("");
    Polyhedra->Edge_Radius->value("");
}

void
Polyhedra_Frame_Combo_cb (Fl_Widget *, void *)
{
// routine called when the frame number is changed on the Polyhedra Edit Screen
    int i, j;

    char atom1[5];

    char atom2[5];

    char widget[16382];

    char atoms[100][5];

    char color[40];

    char string[100];

    int Frame_No = 1;

    if (drvui->max_frame > 1)
	Frame_No = atoi (Polyhedra->Frame_No->value ());
    int nlist = Get_Unique_Atoms (atoms, Frame_No);	// get unique atom names in sorted list

    Polyhedra->Polyhedra_Combo->list.clear ();	// clear out the old names
    for (j = 0; j < nlist; j++) {	// add atom names in alphabetic order
	Polyhedra->Polyhedra_Combo->list.add (atoms[j]);
    }
    Polyhedra->Polyhedra_Combo->value (atoms[0]);	// load first one in window
    Polyhedra_Combo_cb (NULL, NULL);
    widget[0] = '\0';
    for (i = 1; i < drvui->npoly; i++) {
	if (drvui->polyhedra[i].poly_fn != Frame_No)
	    continue;		//skip if not for this frame
	float d1, d2;

	char type[3] = "PS";

	memset (atom1, 0, 5);
	memset (atom2, 0, 5);
	strncpy (atom1, drvui->polyhedra[i].poly_l, 4);
	strncpy (atom2, drvui->polyhedra[i].poly_t, 4);
	if (strlen (atom2))
	    strcpy (type, "PV");
	d1 = drvui->polyhedra[i].poly_min;
	if (d1 > 0.005)
	    strcpy (type, "SH");
	d2 = drvui->polyhedra[i].poly_size;
	strcpy (color, drvui->polyhedra[i].poly_col);
	if (color[strlen (color) - 1] < ' ')
	    color[strlen (color) - 1] = 0;
	if (strlen (color) > 25)
	    color[25] = 0;
	trim_string (atom1, 5);
	if (strlen (atom2))
	    trim_string (atom2, 5);
	sprintf (string, "  %2s      %4s    %4s %9.3f %9.3f     %s\n", type, atom1, atom2,
		 d1, d2, color);
	strcat (widget, string);
    }
    for (i = 1; i < drvui->nplane; i++) {
	float d1, d2;

	char type[3] = "PL";

	strncpy (atom1, drvui->planes[i].plane_l, 4);
	strcpy (atom2, "");
	atom1[4] = 0;
	for (j = 3; j >= 0; --j) {
	    if (atom1[j] == ' ')
		atom1[j] = 0;
	}
	d1 = 0.;
	d2 = drvui->planes[i].plane_size;
	strcpy (color, drvui->planes[i].plane_col);
	if (color[strlen (color) - 1] < ' ')
	    color[strlen (color) - 1] = 0;
	if (strlen (color) > 25)
	    color[25] = 0;
	sprintf (string, "  %s      %4s    %4s %9.3f %9.3f     %s\n", type, atom1, atom2,
		 d1, d2, color);
	strcat (widget, string);
    }
    Polyhedra->PolyhedraBuffer->text (widget);
    Polyhedra->PolyInstr1->show ();
    Polyhedra->PolyInstr->hide ();
}

void
Sphere_Combo_cb (Fl_Widget *, void *)
{
// load size and color when the Sphere Combo is changed

    Spheres->New_Sphere_Size->value ("0.4");
    Spheres->New_Sphere_Color->value ("Red");
    Spheres->New_Sphere_Add->activate ();
/*
    for (i=1; i<drvui->n_ellips; i++) {
        if (check_atom_name(atom1,drvui->ellips_l[i]) && drvui->ell_type[i] < 1000 && drvui->ell_type[i] >= 0) {
            Spheres->New_Sphere_Convert->activate();
            i = drvui->n_ellips;
        }
    }
*/
}

struct combo_atoms_struct
{
    char combo_atoms[7];
};

void
Sphere_Frame_Combo_cb (Fl_Widget *, void *)
{
// routine called when the frame number is changed on the Spheres Edit Screen
    int i, j, no_ad, n, n_combo = 0;

    char atom[5];

    char tstring[5];

    struct combo_atoms_struct *combo;

    char widget[16382];

    char color[40];

    char string[100];

    int Frame_No = 1;

    combo =
	(combo_atoms_struct *) zalloc (drvui->atom_alloc *
				       sizeof (struct combo_atoms_struct));
    if (!combo) {
	Error_Box ("Unable to allocate space for combo atoms.\n");
	exit (0);
    }
    if (drvui->max_frame > 1)
	Frame_No = atoi (Spheres->Frame_No->value ());
    Spheres->Sphere_Combo->list.add ("");	// add a dummy at start
    for (j = 0; j < natom; j++) {
	memset (atom, 0, 5);
	if (drvui->atoms[j].atom_fn != Frame_No)
	    continue;
	strcpy (atom, drvui->atoms[j].atom_l);
	strcpy (string, drvui->atoms[j].atom_l);
	n = drvui->atoms[j].sv_atom_n;;
	no_ad = 0;
	for (i = 1; i < drvui->n_ellips; i++) {
	    if (check_atom_name (atom, drvui->ellips[i].ellips_l)
		&& drvui->ellips[i].ell_type > 1000 && drvui->ellips[i].ellips_n == n) {
		no_ad = 1;
		i = drvui->n_ellips;
	    }
	}
	if (!no_ad) {
	    sprintf (tstring, "%2d", n);
	    strcat (string, tstring);
	    Spheres->Sphere_Combo->list.add (string);
	    n = 1;
	    for (i = 0; i < n_combo; i++) {
		if (check_atom_name (atom, combo[i].combo_atoms))
		    n = 0;
	    }
	    if (n) {
		strcpy (string, atom);
		strcat (string, " *");
		strcpy (combo[n_combo++].combo_atoms, string);
	    }
	}
    }
    for (i = 0; i < n_combo; i++)
	Spheres->Sphere_Combo->list.insert (i + 2, combo[i].combo_atoms);
    Spheres->Sphere_Combo->value ("");	// load blank in window
    Sphere_Combo_cb (NULL, NULL);
    widget[0] = '\0';
    for (i = 1; i < drvui->nsphere; i++) {	// process sphere info
	if (drvui->spheres[i].sphere_fn != Frame_No)
	    continue;		//skip if not for this frame
	float d1;

	d1 = drvui->spheres[i].sphere_size;
	strcpy (atom, drvui->spheres[i].sphere_l);
	while (strlen (atom) < 4)
	    strcat (atom, " ");
	strcpy (color, drvui->spheres[i].sphere_col);
	if (drvui->spheres[i].sphere_n > 0) {
	    sprintf (string, "%4s%2d %9.3f     %s\n", atom, drvui->spheres[i].sphere_n,
		     d1, color);
	} else {
	    sprintf (string, "%4s * %9.3f     %s\n", atom, d1, color);
	}
	strcat (widget, string);
    }
    Spheres->SphereBuffer->text (widget);
    if (strlen (widget) > 5)
	Spheres->SphereInstr1->show ();
    free (combo);
}

void
View_Console_cb (void)		// callback to list the 'cns' file
{				// see View_Listing_cb for detailed comments
    static int one = 1;

// callback routine to view the console listing
    if (!strlen (drvui->Cur_File)) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }
    if (!listwindow1) {
	listwindow1 = new Fl_Window (50, 50, 680, 450, "Console Output");
	listwindow1->begin ();
	listwindow1->callback ((Fl_Callback *) View_Listing_Close_cb, &one);
	Fl_Text_Editor *display = new Fl_Text_Editor (0, 0, 680, 410);

	textbuf1 = new Fl_Text_Buffer;
	display->buffer (textbuf1);
	display->textfont (FL_COURIER);
	Fl_Button *o = new Fl_Button (300, 415, 80, 30, "Close");

	o->tooltip ("Close this window and discard all changes.");
	o->callback ((Fl_Callback *) View_Listing_Close_cb, &one);
	listwindow1->end ();
#if !defined (WIN32) && !defined (__APPLE__)
	listwindow1->icon ((char *) drvui->icon);
#endif
    }
    textbuf1->loadfile (drvui->Cur_Console);
    listwindow1->show ();
}

void
View_Cursor_cb (void)		// callback to show a dedicated window for
{				// Graphics Cursor feedback
    static int five = 5;

    if (!strlen (drvui->Cur_File)) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }

    if (!listwindow5) {
	listwindow5 = new Fl_Window (50, 50, 550, 550, "Geometry Information");
	listwindow5->begin ();
	listwindow5->callback ((Fl_Callback *) View_Listing_Close_cb, &five);

	drvui->Cursor_posW = new Fl_Browser (15, 15, 520, 450);
	static int bwidths[] = { 160,100,80,80,80,0 };
        static const Fl_Menu_Item resetmodes[] = {
	    {"never", 0, NULL, (void *) 0, 0, 0, 0, 14, 56},
	    {"after 2 atoms", 0, NULL, (void *) 0, 0, 0, 0, 14, 56},
	    {"after 3 atoms", 0, NULL, (void *) 0, 0, 0, 0, 14, 56},
	    {"after 4 atoms", 0, NULL, (void *) 0, 0, 0, 0, 14, 56},
	    {0, 0, 0, 0, 0, 0, 0, 0, 0}
	};
	drvui->Cursor_posW->column_widths (bwidths);
	drvui->Cursor_posW->column_char ('\t');
	drvui->Cursor_posW->type (FL_MULTI_BROWSER);
	drvui->Cursor_posW->box (FL_FRAME_BOX);
	drvui->Cursor_posW->color (FL_BACKGROUND_COLOR);

	update_cursor_window();

	listwindow5->resizable(drvui->Cursor_posW);

	Fl_Choice *p = drvui->Cursor_reset =
	    new Fl_Choice (220, 480, 125, 25, "Reset list ");
	p->align (FL_ALIGN_LEFT);
	p->callback (Cursor_Reset_Combo_cb);
	p->labelfont (1);
	p->menu (resetmodes);
	p->value (0);

	Fl_Button *o = new Fl_Button (240, 515, 80, 30, "Close");

	o->tooltip ("Close this window");
	o->callback ((Fl_Callback *) View_Listing_Close_cb, &five);
	listwindow5->end ();
#if !defined (WIN32) && !defined (__APPLE__)
	listwindow5->icon ((char *) drvui->icon);
#endif
    }
    listwindow5->show ();
}

void
View_File_cb (void)		// callback to view any file
{				// see View_Listing_cb for detailed comments
    static int three = 3;

    char *newfile;

    static char label[100];

// callback routine to view the any file
    if (!listwindow3) {
	listwindow3 = new Fl_Window (50, 50, 680, 450, "");
	listwindow3->begin ();
	listwindow3->callback ((Fl_Callback *) View_Listing_Close_cb, &three);
	Fl_Text_Editor *display = new Fl_Text_Editor (0, 0, 680, 410);

	textbuf3 = new Fl_Text_Buffer;
	display->buffer (textbuf3);
	display->textfont (FL_COURIER);
	Fl_Button *o = new Fl_Button (300, 415, 80, 30, "Close");

	o->tooltip ("Close this window.");
	o->callback ((Fl_Callback *) View_Listing_Close_cb, &three);
	listwindow3->end ();
#if !defined (WIN32) && !defined (__APPLE__)
	listwindow3->icon ((char *) drvui->icon);
#endif
    }
    newfile = fl_file_chooser ("Select File to View", "*.*", NULL, 1);
    if (newfile) {
	textbuf3->loadfile (newfile);
	sprintf (label, "File View: %s", newfile);
	listwindow3->label (label);
	listwindow3->show ();
    }
}

void
View_Listing_cb (void)		// callback to list the '.out' file
{
    static int two = 2;

    if (!strlen (drvui->Cur_File)) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }
// callback to view listing file
    if (!listwindow2) {
	listwindow2 = new Fl_Window (50, 50, 760, 450, "Listing View");
	listwindow2->begin ();
	listwindow2->callback ((Fl_Callback *) View_Listing_Close_cb, &two);
	Fl_Text_Display *display = new Fl_Text_Display (0, 0, 760, 410);

	display->textfont (FL_COURIER);
	textbuf2 = new Fl_Text_Buffer;
	display->buffer (textbuf2);
	Fl_Button *o = new Fl_Button (340, 415, 80, 30, "Close");

	o->tooltip ("Close this window and discard all changes.");
	o->callback ((Fl_Callback *) View_Listing_Close_cb, &two);
	listwindow2->end ();
#if !defined (WIN32) && !defined (__APPLE__)
	listwindow2->icon ((char *) drvui->icon);
#endif
    }
    textbuf2->loadfile (drvui->Cur_Listing);
    listwindow2->show ();
}

void
View_Listing_Close_cb (Fl_Window *, int *arg)	// callback to hide listing windows
{
    switch (*arg) {
    case 1:
	listwindow1->hide ();
	break;
    case 2:
	listwindow2->hide ();
	break;
    case 3:
	listwindow3->hide ();
	break;
    case 4:
	listwindow4->hide ();
	break;
    case 5:
	listwindow5->hide ();
	break;
    }
}

void
View_POV_cb (void)
{
// Callback routine to view the POV file
    char cmd[512], incpath[256], errfile[1024];

    int i;

    if (!strlen (drvui->CurFile->value ())) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }
    Update_Str (0);		// Update str file
    Generate_Drawing (0);	// Update picture
    Fl::flush ();		// Clear the menu bleed through

    incpath[0] = '\0';
    if (strlen (drvui->POV_Include) > 10)
	strcpy (incpath, drvui->POV_Include);	// make copy of colors.inc full path name

#ifdef WIN32
    strcpy (cmd, "\"\"");	// build the command string
    strcat (cmd, drvui->POV_Path);	//   for Windows
    strcat (cmd, "\"");
    strcat (cmd, " ");
    strcat (cmd, drvui->POV_Options);
    strcat (cmd, " ");
    strcat (cmd, " +I\"");
#else
    strcpy (cmd, drvui->POV_Path);	// build command string
    strcat (cmd, " ");		//    for Linux
    strcat (cmd, drvui->POV_Options);
    strcat (cmd, " ");
    if (strlen (incpath) > 0) {
	strcat (cmd, " +HI");
	strcat (cmd, incpath);
	strcat (cmd, " ");
    }
#endif
    strcat (cmd, drvui->Cur_Root);
    strcat (cmd, ".pov");

    strcpy (errfile, drvui->Cur_Root);
    strcat (errfile, ".err");
#ifdef WIN32
    strcat (cmd, "\"\"");
#else
    strcat (cmd, " +GF");
    strcat (cmd, errfile);
#endif
    FILE *ef = fopen (errfile, "w");

    fprintf (ef, "An error occurred running the POVRAY program.\n");
    fprintf (ef, "The most probable cause is that the path is not correct.\n");
    fprintf (ef, "Check settings in File/Configure/POV menu\n");
    fclose (ef);

    Fl::flush ();
#ifdef WIN32
    drvui->mainWindow->hide ();
#endif
    i = system (cmd);		// call the POV program
#ifdef WIN32
    drvui->mainWindow->show ();
#endif
    if (i != 0) {
	if (!listwindow4) {
	    static int four = 4;

	    listwindow4 = new Fl_Window (50, 50, 760, 450, "POV error messages");
	    listwindow4->begin ();
	    listwindow4->callback ((Fl_Callback *) View_Listing_Close_cb, &four);
	    Fl_Text_Display *display = new Fl_Text_Display (0, 0, 760, 410);

	    display->textfont (FL_COURIER);
	    textbuf4 = new Fl_Text_Buffer;
	    display->buffer (textbuf4);
	    Fl_Button *o = new Fl_Button (340, 415, 80, 30, "Close");

	    o->tooltip ("Close this window.");
	    o->callback ((Fl_Callback *) View_Listing_Close_cb, &four);
	    listwindow4->end ();
#if !defined (WIN32) && !defined (__APPLE__)
	    listwindow4->icon ((char *) drvui->icon);
#endif
	}
	textbuf4->loadfile (errfile);
	listwindow4->show ();
//        Error_Box("An error occurred running the POVRAY program.\nCheck settings in File/Configure/POV menu\nand contents of POV output file PROJECT.err.");
    } else
	unlink (errfile);
}

void
Write_Map_cb (void)
{
    int i, j, k, ijk;

    FILE *mapout;

    char *newfile =
	fl_file_chooser ("Select Output 'grd'", "*.*", Maps->Filename->value (), 1);
    if (newfile) {
	if (!(mapout = fopen (newfile, "w"))) {
	    Error_Box ("Unable to open selected output file.");
	    return;
	}
	if (!FourierPt) {
	    Error_Box ("Fourier Map Not Available.");
	    fclose (mapout);
	    return;
	}
	fprintf (mapout, "%s\n", Map_Info.title);
	fprintf (mapout, "%f %f %f %f %f %f\n", Map_Info.lat_con[0],
		 Map_Info.lat_con[1], Map_Info.lat_con[2], Map_Info.lat_con[3],
		 Map_Info.lat_con[4], Map_Info.lat_con[5]);
	fprintf (mapout, "%d %d %d\n", Map_Info.map_int[0], Map_Info.map_int[1],
		 Map_Info.map_int[2]);

	ijk = 0;
/* Write the Rho values */
	for (i = 0; i < Map_Info.map_int[0]; i++) {
	    for (j = 0; j < Map_Info.map_int[1]; j++) {
		for (k = 0; k < Map_Info.map_int[2]; k++) {
		    fprintf (mapout, " %f\n", FourierPt[ijk++]);
		}
	    }
	}
	fclose (mapout);
	Maps->Filename->value (newfile);
	FourierMapType = 1;	// new type
	Maps->MapType->value ("GSAS - grd");
    }
}

void
Dump_View_cb (void)
{
    FILE *fp;

    char file[256];

    int state = GL2PS_OVERFLOW, buffsize = 0;

    int format = GL2PS_EPS;

    int sort;

    int options =
	GL2PS_DRAW_BACKGROUND | GL2PS_NO_BLENDING | GL2PS_OCCLUSION_CULL | GL2PS_BEST_ROOT
	| GL2PS_SILENT | GL2PS_TIGHT_BOUNDING_BOX;
    int nbcol = 16;

    GLint viewport[4];

    float cpx, cpy, cpz;

    float m[16];
    GLfloat mat_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat light_ambient[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    GLfloat light_specular[] = { 0.1f, 0.1f, 0.1f, 1.0f };
    GLfloat light_diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat mat_shininess[] = { 0.8f };
    GLfloat light_position[] = { 0.0f, 1.0f, 1.0f, 0.0f };
    GLfloat light_direction[] = { 0.0f, -1.0f, 0.0f };

    if (!strlen (drvui->CurFile->value ())) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }
//   if (ReadFourMap)
//     sort = GL2PS_SIMPLE_SORT;
//   else
    sort = GL2PS_BSP_SORT;


    Update_Str (0);		// Update str file
    Generate_Drawing (0);	// Update picture

    glGetIntegerv (GL_VIEWPORT, viewport);
    strcpy (file, drvui->Cur_Root);
    strcat (file, ".ps");

    fp = fopen (file, "wb");

    while (state == GL2PS_OVERFLOW) {
	buffsize += 1024 * 1024;
	gl2psBeginPage (file, "DRAWxtl 5.5", viewport, format, sort, options,
			GL_RGBA, 0, NULL, nbcol, nbcol, nbcol, buffsize, fp, file);

	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	if (M_cameras == 0) {
	    glOrtho (-gl_size, gl_size, -gl_size, gl_size, -10000., 10000.);
	} else {
	    gluPerspective (17., 1., 0.01, 1000.);
	}
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity ();
	cpx = (POV_Max[0] + POV_Min[0]) / 2.0f;
	cpy = (POV_Max[1] + POV_Min[1]) / 2.0f;
	cpz = (POV_Max[2] + POV_Min[2]) / 2.0f;
//    glShadeModel (GL_SMOOTH);
	glMaterialfv (GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv (GL_FRONT, GL_SHININESS, mat_shininess);
	glLightfv (GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv (GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv (GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightModeli (GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	glEnable (GL_LIGHTING);
	glEnable (GL_LIGHT0);
	glEnable (GL_COLOR_MATERIAL);
	glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable (GL_DEPTH_TEST);
	glDepthFunc (GL_LESS);
	glClearColor (drvui->glback[0], drvui->glback[1], drvui->glback[2], 0.0f);
	glClear ((GLbitfield) (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT));
//    glLoadIdentity();
//    glEnable(GL_BLEND);
//    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	gluLookAt (cpx, cpy, Scale * .50,	// camera position
		   cpx, cpy, -1.0,	// camera lookat point
		   0.0f, 1.0f, 0.0f);	// camera "up" vector
	glLightfv (GL_LIGHT0, GL_POSITION, light_position);
	glLightfv (GL_LIGHT0, GL_SPOT_DIRECTION, light_direction);
	glColor3f (0.0f, 0.0f, 0.0f);
	glPushMatrix ();
//    glTranslatef (gl_pos_x, gl_pos_y, gl_pos_z);
	glTranslatef (drvui->Trans[0], drvui->Trans[1], drvui->Trans[2]);
	quaternion_to_rotmatrix (&Rotq, m);

	glMultMatrixf (m);
	glTranslatef (-cpx, -cpy, -cpz);
	draw_cursor ();
	glCallList (drvui->crystalDL);
	for (drvui->frame_no = 1; drvui->frame_no <= drvui->max_frame; drvui->frame_no++)
	    generate_gl_texts ();
	glPopMatrix ();
	drvui->frame_no = drvui->max_frame;


	state = gl2psEndPage ();
    }

    if (state > 2)
	Error_Box ("Generation of Postscript file failed.");

    fclose (fp);

}

void
Cursor_Reset_Combo_cb (Fl_Widget *, void *)
{
    drvui->cur_reset = drvui->Cursor_reset->value ();
    if (drvui->cur_reset == 0) drvui->cur_reset --;
}

void
Edit_Surfaces_cb (void)
{
// Callback routine to show the Surface Parameter screen and load the widgets on that page
    char string[128];

    static int zero = 0;

    static int one = 1;

    int y;

    int i;

    int wh;

    char widget[2048];

    char type[6];

    if (!strlen (drvui->CurFile->value ())) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }
    Save_Working_Copy ();
    if (!Surf) {
	char title[30];

	    strcpy (title, "Edit Surface parameters");
	Surf = new SurfParam;	// new instance of the surface parameters
	
	wh = 535;
	if (drvui->max_frame > 1) 
	    wh += 60;
	Surf->Surfaces_Edit_Window = new Fl_Window (50, 50, 500, wh, title);

	Surf->Surfaces_Edit_Window->callback ((Fl_Callback *) Edit_Surfaces_Close_cb);
	y = 20;
	if (drvui->max_frame > 1) {
	    Surf->Frame_No =
		new Flu_Combo_List (220, y, 75, 25, "Frame No.");
	    Surf->Frame_No->align (FL_ALIGN_TOP);
	    Surf->Frame_No->callback (Surface_Frame_Combo_cb);
	    Surf->Frame_No->labelfont (1);
	    for (i = 1; i <= drvui->max_frame; i++) {
		sprintf (string, "%d", i);
		Surf->Frame_No->list.add (string);
	    }
	    Surf->Frame_No->pop_height (20 * drvui->max_frame);
	    Surf->Frame_No->value ("1");
            y += 50;
	}
	Fl_Box *ob = new Fl_Box (10,y-20,480,y+200,"Bader AIM theory surfaces");
	ob->box(FL_ENGRAVED_BOX);
	ob->align(FL_ALIGN_TOP|FL_ALIGN_INSIDE);
	y += 15;
 
	strcpy (widget, "");
	Surf->AimSurfBuffer = new Fl_Text_Buffer;
	Surf->AimSurfBuffer->add_modify_callback ((Fl_Text_Modify_Cb) Modify_AimSurfaces_cb,
					       (void *) NULL);
	Fl_Text_Editor *o =
		new Fl_Text_Editor (50, y, 420, 80, "Atom       surfFile         Type       Color     ");
	o->labelfont (FL_COURIER_BOLD);
	o->textfont (FL_COURIER);
	o->buffer (Surf->AimSurfBuffer);

	Surf->AimSurfInstr = new Fl_Output (20, y + 100, 450, 0, "Press 'Add' to replace "
					 "selected line - 'Remove' to delete it");
	Surf->AimSurfInstr->hide ();
	Surf->AimSurfInstr->align (FL_ALIGN_BOTTOM);
	Surf->AimSurfInstr1 = new Fl_Output (20, y + 85, 450, 0, "Highlight text above "
					  "or double click to edit line");
	Surf->AimSurfInstr1->align (FL_ALIGN_BOTTOM);
	y += 110;
	Flu_Combo_List *oc = Surf->AimSurf_Combo =
	    new Flu_Combo_List (45, y, 75, 25, "Atom");
	oc->align (FL_ALIGN_BOTTOM);
//	oc->callback (AimSurf_Combo_cb);

	Surf->AimFile = new Fl_Input (120, y, 90, 25, "Filename");
	Surf->AimFile->align (FL_ALIGN_BOTTOM);
	Surf->AimFile->callback ((Fl_Callback *) New_AimSurf_Input_cb);
	Surf->AimSurfType = new Flu_Combo_List (220, y, 70, 25, "Type");
	Surf->AimSurfType->align (FL_ALIGN_BOTTOM);
	Surf->AimSurfType->callback ((Fl_Callback *) New_AimSurf_Input_cb);
	Surf->AimSurfType->list.add ("mesh");
	Surf->AimSurfType->list.add ("solid");
	Surf->AimSurfType->list.add ("dots");
	Surf->AimSurfType->pop_height (40);

	Surf->AimSurfColor = new Flu_Combo_List (300, y, 160, 25, "Color");
	Load_Color_Combo (Surf->AimSurfColor);
	Surf->AimSurfColor->align (FL_ALIGN_TOP);
	Surf->AimSurfColor->callback ((Fl_Callback *) New_AimSurf_Input_cb);
	Surf->AimSurfColor->align (FL_ALIGN_BOTTOM);
	y += 45;

	Fl_Button *om;
	Fl_Button *mm;

	om = Surf->Add_Button = new Fl_Button (180, y, 70, 25, "Add");
	om->callback ((Fl_Callback *) New_AimSurf_Add_cb, &one);
	om->tooltip ("When active, press to transfer data in boxes to window above");
	om->deactivate ();
	mm = Surf->Remove_Button = new Fl_Button (270, y, 70, 25, "Remove");
	mm->callback ((Fl_Callback *) New_AimSurf_Add_cb, &zero);
	mm->tooltip ("When active, press to remove highlighted line.");
	mm->deactivate ();

	Fl_Box *obb = new Fl_Box (10, y + 35, 480, 300, "Contact Surfaces and Void Volumes");
	obb->box(FL_ENGRAVED_BOX);
	obb->align(FL_ALIGN_TOP|FL_ALIGN_INSIDE);

	Surf->SurfType = new Flu_Combo_List (40, y + 60, 180, 25, "Method");
	Surf->SurfType->align (FL_ALIGN_BOTTOM);
//	Surf->SurfType->callback ((Fl_Callback *) New_Surf_Input_cb);
	Surf->SurfType->list.add ("No surface");
	Surf->SurfType->list.add ("Grid search");
	Surf->SurfType->list.add ("MSMS (run now)");
	Surf->SurfType->list.add ("MSMS (precalculated)");
	Surf->SurfType->list.add ("Surface area test");
	Surf->SurfType->pop_height (40);
	Surf->SurfColor = new Flu_Combo_List (300, y + 60, 160, 25, "Color");
	Load_Color_Combo (Surf->SurfColor);
	Surf->SurfColor->list.add("byatom");
	Surf->SurfColor->align (FL_ALIGN_TOP);
	Surf->SurfColor->align (FL_ALIGN_BOTTOM);
	Surf->Probe = new Fl_Input (120, y + 100, 30, 25, "Probe radius:");
	Surf->Probe->align (FL_ALIGN_LEFT);
	Surf->GridX = new Fl_Input (220, y + 100, 40, 25, "Gridsize:");
	Surf->GridX->align (FL_ALIGN_LEFT);
	Surf->GridY = new Fl_Input (270, y + 100, 40, 25, "x");
	Surf->GridY->align (FL_ALIGN_LEFT);
	Surf->GridZ = new Fl_Input (320, y + 100, 40, 25, "x");
	Surf->GridZ->align (FL_ALIGN_LEFT);

	y += 150;

	Surf->SurfBuffer = new Fl_Text_Buffer;
	Surf->SurfBuffer->add_modify_callback ((Fl_Text_Modify_Cb) Modify_Surfaces_cb,
					       (void *) NULL);
	Fl_Text_Editor *os =
		new Fl_Text_Editor (140, y, 200, 60, "Atom       Radius");
	os->labelfont (FL_COURIER_BOLD);
	os->textfont (FL_COURIER);
	os->buffer (Surf->SurfBuffer);

	Surf->SurfInstr = new Fl_Output (20, y + 80, 450, 0, "Press 'Add' to replace "
					 "selected line - 'Remove' to delete it");
	Surf->SurfInstr->hide ();
	Surf->SurfInstr->align (FL_ALIGN_BOTTOM);
	Surf->SurfInstr1 = new Fl_Output (20, y + 60, 450, 0, "Highlight text above "
					  "or double click to edit line");
	Surf->SurfInstr1->align (FL_ALIGN_BOTTOM);
	y += 80;
	Flu_Combo_List *osc = Surf->Surf_Combo =
	    new Flu_Combo_List (145, y, 75, 25, "Atom");
	osc->align (FL_ALIGN_BOTTOM);
//	osc->callback (Surf_Combo_cb);

	Surf->Radius = new Fl_Input (220, y, 90, 25, "Radius");
	Surf->Radius->align (FL_ALIGN_BOTTOM);
	Surf->Radius->callback ((Fl_Callback *) New_Surf_Input_cb);

	y += 40;
	Fl_Button *osm;

	Fl_Button *msm;

	osm = Surf->Add_Button2 = new Fl_Button (180, y, 70, 20, "Add");
	osm->callback ((Fl_Callback *) New_Radius_Add_cb, &one);
	osm->tooltip ("When active, press to transfer data in boxes to window above");
	osm->deactivate ();
	msm = Surf->Remove_Button2 = new Fl_Button (270, y, 70, 20, "Remove");
	msm->callback ((Fl_Callback *) New_Radius_Add_cb, &zero);
	msm->tooltip ("When active, press to remove highlighted line.");
	msm->deactivate ();

	y += 10;

	Surf->AimSurfInstr2 = new Fl_Output (25, y, 450, 0,
					  "Changes are temporary until \"Apply\" or \"Save\" is pressed.");
	Surf->AimSurfInstr2->hide ();
	Surf->AimSurfInstr2->align (FL_ALIGN_BOTTOM);
	y += 20;

	Fl_Button *r = new Fl_Button (125, y, 70, 25, "Close");

	r->tooltip ("Close this window and discard all changes.");
	r->callback ((Fl_Callback *) Edit_Surfaces_Close_cb);
	Fl_Button *s = new Fl_Button (325, y, 70, 25, "Save");

	s->tooltip
	    ("Apply current contents of top box to drawing, then close this window.");
	s->callback ((Fl_Callback *) Edit_Surfaces_Save_cb, &one);
	Fl_Button *a = new Fl_Button (225, y, 70, 25, "Apply");

	a->tooltip
	    ("Apply current contents of top box to drawing, but leave this window open.");
	a->callback ((Fl_Callback *) Edit_Surfaces_Save_cb, &zero);
#if !defined (WIN32) && !defined (__APPLE__)
	Surf->Surfaces_Edit_Window->icon ((char *) drvui->icon);
#endif
    }
    strcpy (widget, "");
    if (drvui->nsurf > 1) 
        for (i = 1; i < drvui->nsurf; i++) {	// fill in the widgets

	    switch (drvui->surftype[i]) {
		case 0:
		default:
		strcpy (type, "mesh ");
	        break;
		case 1:
		strcpy (type, "solid");
		break;
		case 2:
		strcpy (type, "dots ");
		break;
	    }
	    sprintf (string, "%4s%2d    %s      %s     %s\n",
		     drvui->surfatom[i],drvui->surfnum[i], drvui->surffile[i], type,
		     drvui->surfcolor[i]);
	strcat (widget, string);
    }
    Surf->AimSurfBuffer->text (widget);

    Surf->AimSurfInstr1->show ();
    Surface_Frame_Combo_cb(NULL,NULL);

    if (drvui->voidflag==0) {
	Surf->SurfType->value("No surface");
    } else {
	char tmpstr[20];
	if (drvui->voidflag ==1) Surf->SurfType->value("Grid search");
	if (drvui->voidflag ==2) Surf->SurfType->value("MSMS (run now)");
	if (drvui->voidflag ==-2) Surf->SurfType->value("MSMS (precalculated)");
	if (drvui->voidflag ==3) Surf->SurfType->value("Surface area test");
	Surf->SurfColor->value(drvui->voidcolor);
	if (drvui->probesize > 0.) {
		sprintf(tmpstr,"%.2f",drvui->probesize);
		Surf->Probe->value(tmpstr);
	}
	sprintf (tmpstr,"%d",drvui->voidgrid[0]); 
	Surf->GridX->value(tmpstr);
	sprintf (tmpstr,"%d",drvui->voidgrid[1]); 
	Surf->GridY->value(tmpstr);
	sprintf (tmpstr,"%d",drvui->voidgrid[2]); 
	Surf->GridZ->value(tmpstr);

    strcpy (widget, "");
    if (drvui->natprop > 1) 
        for (i = 1; i < drvui->natprop; i++) {	// fill in the widgets
	    if (drvui->atprops[i].atprop_n != -1) 
		sprintf (string, "%4s%2d    %f\n",
			 drvui->atprops[i].atprop_l,drvui->atprops[i].atprop_n,drvui->atprops[i].radius);
	    else
		sprintf (string, "%4s *    %f\n",
			 drvui->atprops[i].atprop_l,drvui->atprops[i].radius);
	    strcat (widget, string);
    }
    Surf->SurfBuffer->text (widget);

    }
    Surf->Surfaces_Edit_Window->end ();
    Surf->Surfaces_Edit_Window->show ();
}

void
Edit_Surfaces_Close_cb (void)
{
    Surf->Surfaces_Edit_Window->~Fl_Window ();	// this window needs to be deleted
    delete (Surf->Surfaces_Edit_Window);	// not just killed (2d/3d nature might change)
    delete (Surf->AimSurfBuffer);
    delete (Surf);
    Surf = NULL;
    Restore_Working_Copy ();	// undo any changes
    Generate_Drawing (1);	// regenerate
}

void
Edit_Surfaces_Save_cb (Fl_Button *, int *save)
{
// callback routine when 'save' or 'apply' button is pressed on the Edit Maps screen
    char type[5], number[5], color[40], widget[16382];
    char filename[128];
    char atom[5];
    char t_name[5];
    int num;
    float radius;

    unsigned int j, k;

    int kk;

    int i = 0;

    int Frame_No = 1;

    char *selection = Surf->AimSurfBuffer->text ();

    if (drvui->max_frame > 1)
	Frame_No = atoi (Surf->Frame_No->value ());

    strcpy (widget, selection);
    free (selection);

    drvui->nsurf = 1;

    if (strlen (widget) < 10) {
	strcpy (widget, "");
	Surf->AimSurfBuffer->text (widget);
    } else {
	while (strlen (widget) > 10) {
	    i++;
		sscanf (widget, "%s %d %s %s %s", atom, &num, filename, type, color);
	    trim_string (color, 40);
	    if (!strlen (color))
		strcpy (color, "Gray20");
		if (strncmp (type, "mesh", 4) == 0)
		    drvui->surftype[i] = 0;
		else if (strncmp (type, "solid", 5) == 0)
		    drvui->surftype[i] = 1;
		else
		    drvui->surftype[i] = 2;

	    strcpy (drvui->surffile[i], filename);
	    strcpy (drvui->surfcolor[i], color);
            strcpy (drvui->surfatom[i], atom);
            drvui->surfnum[i] = num;
 
	    drvui->nsurf++;
	    for (j = 0; j < strlen (widget); j++) {
		if (widget[j] == '\n')
		    break;
	    }
	    for (j++, k = 0; j < strlen (widget); j++)
		widget[k++] = widget[j];
	    widget[k] = 0;
	}
    }
    
    if (!strncmp( Surf->SurfType->value(), "No",2) ) {
	drvui->voidflag = 0;
    } else if (!strncmp ( Surf->SurfType->value(), "Grid", 4) ) {
	drvui->voidflag = 1;
    } else if (!strncmp ( Surf->SurfType->value(), "MSMS", 4) ) {
	drvui->voidflag = 2;
        if (!strstr(Surf->SurfType->value(), "now") ) 
	    drvui->voidflag = -2;
    } else {
	drvui->voidflag = 3;
    }

    if ( drvui->voidflag != 0 ) {
        drvui->probesize = (float)atof(Surf->Probe->value());
        drvui->voidgrid[0] = atoi(Surf->GridX->value());
        drvui->voidgrid[1] = atoi(Surf->GridY->value());
        drvui->voidgrid[2] = atoi(Surf->GridZ->value());
        strcpy(drvui->voidcolor,Surf->SurfColor->value());

	selection = Surf->SurfBuffer->text ();

	strcpy (widget, selection);
	free (selection);

	if (strlen (widget) < 10) {
//	    drvui->voidflag = 0;
	    strcpy (widget, "");
	    Surf->SurfBuffer->text (widget);
	} else {
	    drvui->natprop=1;
	    while (strlen (widget) > 10) {
		i++;
		sscanf (widget, "%s %s %f", atom, number, &radius);
		if (number[0]== '*') {
		    num = -1;
		} else {
		    (void) sscanf (number, "%d", &num);
		}
		drvui->atprops[drvui->natprop].atprop_n=num;
		strcpy(drvui->atprops[drvui->natprop].atprop_l,atom);
		drvui->atprops[drvui->natprop].radius=radius;
		drvui->natprop++;
		for (kk = 0; kk < natom; kk++) {       /* find this atom in list */
		    if (check_atom_name (t_name, drvui->atoms[kk].atom_l)) {
			if ((num == -1) || (num == drvui->atoms[kk].atom_n)) {
			    drvui->atoms[kk].radius = radius;
			}
		    }
		}
		for (j = 0; j < strlen (widget); j++) {
		    if (widget[j] == '\n')
			break;
		}
		for (j++, k = 0; j < strlen (widget); j++)
		    widget[k++] = widget[j];
		widget[k] = 0;
	    }
	}
    }
    drvui->Str_File_Changed = 1;
    if (*save != 3) {
	Save_Working_Copy ();
	Surf->AimSurfInstr2->hide ();
    }
    if (*save == 1) {
	Fl::delete_widget (Surf->Surfaces_Edit_Window);
	delete (Surf);
	Surf = NULL;
    }
    Update_Str (0);
    Generate_Drawing (0);
    Fl::redraw ();
}

void
Modify_AimSurfaces_cb (Fl_Widget *, void *)
{
    char type[8], color[40], value[10],tstring[10];
    char atom[5],filename[128];
    int num;

    int start, end;

    const char *selection;

    if (!Surf->AimSurfBuffer->selected ()) {
	return;
    }
    memset (type, 0, 8);
    memset (color, 0, 40);
    Surf->AimSurfBuffer->selection_position (&start, &end);
    selection = Surf->AimSurfBuffer->line_text (start);
    Surf->AimSurfInstr1->hide ();
    Surf->AimSurfInstr->show ();
	sscanf (selection, " %s %2d %s %s %s", atom, &num, filename, type, color);
    free ((char *) selection);
    if (strlen (color) && isalpha (color[0])) {
	strcpy(tstring,atom);
	if (strlen(tstring)<4) strcat(tstring," ");
	if (strlen(tstring)<4) strcat(tstring," ");
	if (strlen(tstring)<4) strcat(tstring," ");
	sprintf(value,"%2d",num);
	strcat(tstring,value);
	Surf->AimSurf_Combo->value (tstring);
	Surf->AimFile->value (filename);
	Surf->AimSurfType->value (type);
	Surf->AimSurfColor->value (color);
	Surf->Add_Button->activate ();
	Surf->Remove_Button->activate ();
    } else {
	Surf->AimSurf_Combo->value ("");
	Surf->AimFile->value ("");
	Surf->AimSurfColor->value ("");
	Surf->AimSurfType->value ("");
	Surf->Add_Button->activate ();
	Surf->Remove_Button->deactivate ();
    }
}

void
Modify_Surfaces_cb (Fl_Widget *, void *)
{
    char type[8], value[10],tstring[10];
    char atom[5];
    double rad;

    int start, end;

    const char *selection;

    if (!Surf->SurfBuffer->selected ()) {
	return;
    }
    memset (atom, 0, 5);
    memset (type, 0, 8);
    Surf->SurfBuffer->selection_position (&start, &end);
    selection = Surf->SurfBuffer->line_text (start);
    Surf->SurfInstr1->hide ();
    Surf->SurfInstr->show ();
	sscanf (selection, "%s %s %lf", atom, type, &rad);
    free ((char *) selection);
    
	strcpy(tstring,atom);
    if (strlen(tstring) && rad >0.) {
	if (strlen(tstring)<4) strcat(tstring," ");
	if (strlen(tstring)<4) strcat(tstring," ");
	if (strlen(tstring)<4) strcat(tstring," ");
	strcat(tstring,type);
	Surf->Surf_Combo->value (tstring);
	sprintf(value,"%.2f",rad);
	Surf->Radius->value(value);
	Surf->Add_Button2->activate ();
	Surf->Remove_Button2->activate ();
    } else {
	Surf->Surf_Combo->value ("");
	Surf->Radius->value ("");
	Surf->Add_Button2->activate ();
	Surf->Remove_Button2->deactivate ();
    }
}

void
New_AimSurf_Input_cb (Fl_Widget *, void *)
{
// callback routine to make 'Add' new map level button active whenever all fields are non-blank
    if (!strlen (Surf->AimSurfColor->value ()))
	return;
    if (!strlen (Surf->AimSurfType->value ()))
	return;
    Surf->Add_Button->activate ();
}

void
New_Surf_Input_cb (Fl_Widget *, void *)
{
// callback routine to make 'Add' new map level button active whenever all fields are non-blank
    if (!strlen (Surf->Surf_Combo->value ()))
	return;
    if (!strlen (Surf->Radius->value ()))
	return;
    Surf->Add_Button2->activate ();
}

void
New_AimSurf_Add_cb (class Fl_Widget *, int *action)
{
    char string[100];

    int start, end;

    char *selection = NULL;

    Surf->AimSurfInstr2->show ();
    if (Surf->AimSurfBuffer->selected ()) {
	Surf->AimSurfBuffer->selection_position (&start, &end);
	selection = Surf->AimSurfBuffer->line_text (start);
	Surf->AimSurfBuffer->remove (Surf->AimSurfBuffer->line_start (start),
				  Surf->AimSurfBuffer->line_end (end) + 1);
    }
    if (*action) {
	    sprintf (string, "%s    %s      %s      %s\n", Surf->AimSurf_Combo->value(),Surf->AimFile->value(), 
		     Surf->AimSurfType->value (), Surf->AimSurfColor->value ());
	Surf->AimSurfBuffer->append (string);
    }
    if (selection)
	free ((char *) selection);
    Surf->Add_Button->deactivate ();
    Surf->Remove_Button->deactivate ();
    Surf->AimSurfInstr1->show ();
    Surf->AimSurfInstr->hide ();
    Edit_Surfaces_Save_cb (NULL, &three);
}

void
New_Radius_Add_cb (class Fl_Widget *, int *action)
{
    char string[100];

    int start, end;

    char *selection = NULL;

    Surf->AimSurfInstr2->show ();
    if (Surf->SurfBuffer->selected ()) {
	Surf->SurfBuffer->selection_position (&start, &end);
	selection = Surf->SurfBuffer->line_text (start);
	Surf->SurfBuffer->remove (Surf->SurfBuffer->line_start (start),
				  Surf->SurfBuffer->line_end (end) + 1);
    }
    if (*action) {
	    sprintf (string, "%s    %s\n", Surf->Surf_Combo->value(),Surf->Radius->value());
	Surf->SurfBuffer->append (string);
    }
    if (selection)
	free ((char *) selection);
    Surf->Add_Button2->deactivate ();
    Surf->Remove_Button2->deactivate ();
    Surf->SurfInstr1->show ();
    Surf->SurfInstr->hide ();
    Edit_Surfaces_Save_cb (NULL, &three);
}

//void AimSurf_Combo_cb(Fl_Widget*, void*)
//{}
//void Surf_Combo_cb(Fl_Widget*, void*)
//{}

void
Surface_Frame_Combo_cb (Fl_Widget *, void *)
{
    int i, j, n;

    char atom[5];

    char tstring[5];

    char string[100];

    int Frame_No = 1;

    if (drvui->max_frame > 1)
	Frame_No = atoi (Surf->Frame_No->value ());

    Surf->AimSurf_Combo->list.clear ();	// clear out the old names
    Surf->AimSurf_Combo->list.add ("");	// add a dummy at start

    Surf->Surf_Combo->list.clear ();	// clear out the old names
    Surf->Surf_Combo->list.add ("");	// add a dummy at start

    for (j = 0; j < natom; j++) {
	memset (atom, 0, 5);
	if (drvui->atoms[j].atom_fn != Frame_No)
	    continue;
	strcpy (string, drvui->atoms[j].atom_l);
	n = drvui->atoms[j].sv_atom_n;;
	    sprintf (tstring, "%2d", n);
	    strcat (string, tstring);
	    Surf->AimSurf_Combo->list.add (string);
	    Surf->Surf_Combo->list.add (string);
    }
    for (j = 0; j < natom; j++) {
	memset (atom, 0, 5);
	if (drvui->atoms[j].atom_fn != Frame_No)
	    continue;
        i=0;
	strcpy (string, drvui->atoms[j].atom_l);
        for (n=0;n<j;n++) {
		if (!strcmp(string,drvui->atoms[n].atom_l)) {
		i=1;
            	continue;
		}
	}
	if (i==1) continue;
	    sprintf (tstring, "*");
	    if (strlen(string)<4) strcat(string," ");
	    if (strlen(string)<4) strcat(string," ");
	    if (strlen(string)<4) strcat(string," ");
	    strcat (string, tstring);
	    Surf->Surf_Combo->list.add (string);
    }
    Surf->AimSurf_Combo->value ("");	// load blank in window
    Surf->Surf_Combo->value ("");	// load blank in window
}

