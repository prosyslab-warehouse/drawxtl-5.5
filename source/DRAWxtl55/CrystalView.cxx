// $Id: CrystalView.cxx 1103 2011-01-15 16:14:13Z martin $
//
// CrystalView.cxx - main routine for DRAWxtl V5.5 - the GUI version 
//
// Coded using the FLTK 1.1.6 widget set
//
//     Larry W. Finger, Martin Kroeker and Brian Toby
//
// This module includes many of the support routines for the GUI
//
// routines contained within this file:
//
//  main - the entry point
//  AngleInRange - places angle in range -180 to 180
//  Clear_Last_Omit_cb - callback routine to remove last entry in 'omit' item list
//  CrystalView::CrystalView constructor
//  CrystalView::draw - draw routine for CrystalView class
//  Destroy_Open_Windows - call the destructor for all open windows
//  Error_Box - display error box
//  Exit_cb - the callback routine from the "Exit" menu button pressed or when main window is closed
//  ImportDataFile_cb - Callback from Import External Format menu item
//  Include_Cutouts_cb - callback routine when 'Include Cutouts' check box is changed
//  Load_Bond_Data - get bond data from out file into array
//  LoadConfig - load configuration file
//  Max_Min_cb - callback routine to readout Max- and Min- sliders
//  Offset_cb - callback routine when Origin values are changed
//  pick_box - routine to pick the corners of the slab box overlay
//  process_hits - routine to process the omit list
//  Process_Inp - read the DRAWxtl input file
//  Restore_Working_Copy - restores working copy of str file from saved version
//  Rotation_cb - callback routine to readout rotation widgets
//  Save_Current_cb - Callback from Save Current menu button
//  Save_Working_Copy - Copies the working version of the str file into a "save" file
//  SelectDataFile_cb - Callback from Select Data File menu item
//  show_slab_ovl - overlay draw routine for the temporary slab outline 
//  start_picking - routine to pick objects to be omitted from drawing
//  update_box - routine to process the GL hits from pick_box
//  Update_Str - routine to update the 'str' file from widget contents
//  WriteConfig - write updated configuration file
//  XYZ_Rot_to_Q - X, Y, Z rotations to quaternion
//  moveto_atom - routine to pick an atom and move the crosshair to it (UNUSED)
//  pick_label - routine to pick a labeltext and return its number

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <FL/x.H>
#include <FL/Fl_Tooltip.H>
#include "drawxtl.h"
#include "DRAWxtlViewUI.h"
#include "EditView.h"
#include "Ellipsoids.h"
#include "draw_gbl.h"

extern int Block_CIF;

extern int volatile w_width;

extern volatile int w_height;

extern int w_st_x;

extern int w_st_y;

#ifdef WIN32
#include <io.h>
#include <direct.h>
#define AMBLEV 0.2f
#define snprintf _snprintf
#else
#include <unistd.h>
#define _chdir chdir
#if !defined(__APPLE__)
#include <X11/xpm.h>
# endif
#define AMBLEV 0.2f
#endif

// function prototypes

#include "DRAWxtl_proto.h"

#ifdef WIN32
const char *flu_file_chooser (const char *, const char *, const char *);
char Configure_file[] = { ".drawxtlrc" };

char foldersymbol = '\\';
#else
char Configure_file[] = { "~/.drawxtlrc" };

char foldersymbol = '/';
#endif

#ifdef FREEGLUT24
struct freeglut_compat
{				// ugly hack to access internal data structures of freeglut 2.4
    int n0;
    int n1;
    GLboolean n2;
    int n3;
    int n4;
    GLboolean n5;
    unsigned int n6;
    GLboolean Initialised;
};

extern freeglut_compat fgState;
#endif

int
main (int argc, char **argv)
{
    int i = 0;

#if !defined (WIN32) && !defined (__APPLE__)
#include "DRAWxtl.xpm"
#endif

    LoadConfig (0);		/* quick read to get saved size only */
    drvui = new DRAWxtlViewUI;
    drvui->destroy = 0;
    drvui->max_frame = 0;
    drvui->origin1_flag = 0;
    drvui->Stereo = 0;
    drvui->cross_eyed = 1;
    drvui->stereo_base = 0.05f;
    drvui->atom_no = NULL;
    drvui->atom_so = NULL;
    drvui->orig_atom_no = NULL;
    drvui->table = NULL;
    drvui->msgbuffer = (char*) zalloc(sizeof(char));
    drvui->vert_occ = NULL;
    drvui->label_scale = 1.0f;
    drvui->triple[0] = 0;
    drvui->labels_inited = 0;
    drvui->auto_ellipse = 0;
    drvui->nmag_alloc = 0;
    drvui->nsurf = 0;
    drvui->crystalDL = 0;
    drvui->voidmap = NULL;
    drvui->voidflag = 0;
    drvui->voiddata1 = NULL;
    drvui->voiddata2 = NULL;
    drvui->Cursor_posW = NULL;
    drvui->cur_reset = -1;
    init_dynamic_storage ();
    memset (drvui->saved_x_label, 0, 24 * sizeof (float));
    memset (drvui->slab_con, 0, sizeof (drvui->slab_con));
    offset[0] = offset[1] = offset[2] = 0.0f;
    drvui->glback[0] = drvui->glback[1] = drvui->glback[2] = 1.0f;
    vzero (drvui->Trans);
    Fl_Tooltip::font (FL_COURIER_BOLD);
    Fl_Tooltip::size (12);
    Fl_Tooltip::color (23);
    strcpy (drvui->Cur_File, "");
    strcpy (drvui->Cur_Root, "");
    drvui->Str_File_Changed = 0;
    strcpy (drvui->ProgramPath, argv[0]);
    Omit = new OmitParam;
    drvui->b_mat[0][0] = drvui->b_mat[1][1] = drvui->b_mat[2][2] = 1.;
    drvui->b_mat[0][1] = drvui->b_mat[1][0] = 0.;
    drvui->b_mat[0][2] = drvui->b_mat[2][0] = 0.;
    drvui->b_mat[1][2] = drvui->b_mat[2][1] = 0.;

#ifdef WIN32
    drvui->mainWindow->icon ((char *) LoadIcon (fl_display, MAKEINTRESOURCE (IDI_ICON1)));
#else
    fl_open_display ();
#if !defined(__APPLE__)
    XpmCreatePixmapFromData (fl_display, DefaultRootWindow (fl_display), drawxtl_xpm,
			     &drvui->icon, &shapemask, NULL);
    drvui->mainWindow->icon ((char *) drvui->icon);
#endif
#endif
    Fl::args (argc, argv, i);

#ifdef FREEGLUT24
    fgState.Initialised = 1;	// convince freeglut 2.4 that we did call glutInit()  
#endif

    LoadConfig (1);
    if (i < argc) {
	if (strcmp (argv[i], "-h") && strcmp (argv[i], "-?")) {
	    strcpy (drvui->Cur_File, argv[i]);	// input file specified on start line
	    if (!strchr (drvui->Cur_File, foldersymbol)) {
		getcwd (drvui->Cur_Dir, 1023);
	    } else {
		strcpy (drvui->Cur_Dir, drvui->Cur_File);
		char *end = strrchr (drvui->Cur_Dir, foldersymbol);

		end++;
		*end = '\0';

// path may be relative to the cwd, go there to obtain the absolute path
		chdir (drvui->Cur_Dir);
		getcwd (drvui->Cur_Dir, 1023);

// remove path component from filename
		char *start = strrchr (drvui->Cur_File, foldersymbol);

		start++;
		memmove (drvui->Cur_File, start, strlen (start) + 1);
	    }

	    if (strstr (drvui->Cur_File, ".cif")
		&& strlen (strstr (drvui->Cur_File, ".cif")) == 4) {	// filename ends in .cif
		char tmp_file[256], string[256], newfile[256];

		FILE *newstr, *inp;

		static int one = 1;

		strcpy (newfile, drvui->Cur_File);
		strcpy (tmp_file, drvui->Cur_File);	// copy original cif filename
		strcat (tmp_file, ".str");	// and add extension

		newstr = fopen (tmp_file, "r");
		if (!(inp = fopen (newfile, "r"))) {
		    sprintf (string, "The file you selected ('%s') cannot be read\n"
			     "Do you wish to continue?", newfile);
		    if (fl_choice ("%s", "No", "Yes", NULL, string)) {
			inp = fopen (tmp_file, "w");
			fclose (inp);
			Edit_STR_cb (NULL, &one);
		    }
		} else {
		    strcpy (drvui->Cur_File, tmp_file);
		    WriteConfig ();	// update configuration file
		    drvui->CurFile->value (drvui->Cur_File);	// update main screen widgets
		    drvui->CurDir->value (drvui->Cur_Dir);
		    newstr = fopen (drvui->Cur_File, "w");
		    fprintf (newstr, "titl imported\n");
		    Block_CIF = 0;
		    fprintf (newstr, "import cif %s\n", newfile);
		    drvui->auto_ellipse = 1;
		    fprintf (newstr, "box 0.02 Black\n");
		    fprintf (newstr, "background White\n");
		    fprintf (newstr, "view 0. 0. 0.\n");
		    fprintf (newstr, "pack  -0.05 1.05 -0.05 1.05 -0.05 1.05\n");
		    fprintf (newstr, "end\n");
		    fclose (newstr);
		    fclose (inp);
		}
	    }
	    strcpy (drvui->Cur_Root, drvui->Cur_File);
	    trim_string (drvui->Cur_Root, 256);
	    int j;

	    for (j = strlen (drvui->Cur_Root); j > 0; --j) {
		if (drvui->Cur_Root[j] == '.') {
		    drvui->Cur_Root[j] = 0;
		    break;
		}
	    }
	}
    }
    WriteConfig ();
    drvui->mainWindow->show (argc, argv);
    return Fl::run ();
}

CrystalView::~CrystalView ()
{
}

double
AngleInRange (double angle)
{
// routine to return an angle in the range -180 to 180
    if (angle > 180.0)
	angle -= 180.0;
    if (angle < -180.0)
	angle += 180.0;
    if (angle < 0.0) {
	angle = int (100.0 * angle - 0.5) / 100.0;
    } else {
	angle = int (100.0 * angle + 0.5) / 100.0;
    }
    return angle;
}

void
Clean_Up_Files (void)
{
// delete the "frm", "save", "tmp" and "cns" files
    int r;

    char string[512], tmp[50];

    if (!strlen (drvui->Cur_Root))
	return;

    for (r = 1; r <= drvui->max_frame; r++) {	// delete the frm* files
	strcpy (string, drvui->Cur_Root);
	sprintf (tmp, ".frm%d", r);
	strcat (string, tmp);
	unlink (string);
    }
    strcpy (string, drvui->Cur_Root);
    strcat (string, ".save");
    unlink (string);		// delete the "save" file
    unlink (drvui->Cur_Temp);	// delete the tmp file
    unlink (drvui->Cur_Console);	// delete the cns file
}

void
Clear_Last_Omit_cb (Fl_Button *, void *)
{
    if (Omit->nomits > 0) {
	drvui->Str_File_Changed = 1;
	--Omit->nomits;
	if (!Omit->nomits) {
	    edtprm->ClearLastOmit->deactivate ();
	    edtprm->ClearOmit->deactivate ();
	}
    }
    Update_Str (0);		// update the 'str' file
    Generate_Drawing (1);	// regenerate the drawing
    Fl::redraw ();		// update the screen
}

CrystalView::CrystalView (int x, int y, int w, int h, const char *l)
    :
Tb_Window (x, y, w, h, l)
{
// CrystalView constructor - subclass of Tb_Window - use default constructor
}

void
Destroy_Open_Windows (void)
{
// call the destructor for all open windows and delete the database
    if (arrows) {
	//arrows->ArrowWindow->~Fl_Window ();
	delete (arrows->ArrowWindow);
	delete (arrows->ArrowBuffer);
	delete (arrows);
	arrows = NULL;
    }
    if (Bonds) {
//        delete(Bonds->Bond_Output_Buffer);
//        delete(Bonds->BondBuffer);
	//Bonds->Bond_Edit_Window->~Fl_Window ();
	delete (Bonds);
	Bonds = NULL;
    }
    if (Configure) {
	//Configure->ConfigWindow->~Fl_Window ();
	delete (Configure->ConfigWindow);
	delete (Configure);
	Configure = NULL;
    }
    if (MiscConfigure) {
	//MiscConfigure->MiscConfigWindow->~Fl_Window ();
	delete (MiscConfigure->MiscConfigWindow);
	delete (MiscConfigure);
	MiscConfigure = NULL;
    }
    if (MSMSConfigure) {
	//MSMSConfigure->MSMSConfigWindow->~Fl_Window ();
	delete (MSMSConfigure->MSMSConfigWindow);
	delete (MSMSConfigure);
	MSMSConfigure = NULL;
    }
    if (edtprm) {
	//edtprm->editWindow->~Fl_Window ();
	delete (edtprm->editWindow);
	delete (edtprm);
	edtprm = NULL;
    }
    if (ellipsoids) {
	//ellipsoids->Ellips_Window->~Fl_Window ();
	delete (ellipsoids->Ellips_Window);
	delete (ellipsoids->ColorInputBuf);
	delete (ellipsoids);
	ellipsoids = NULL;
    }
    if (LonePairs) {
	//LonePairs->LonePair_Edit_Window->~Fl_Window ();
	delete (LonePairs->LonePair_Edit_Window);
	delete (LonePairs->LonePairBuffer);
	delete (LonePairs);
	LonePairs = NULL;
    }
    if (Maps) {
	//Maps->Maps_Edit_Window->~Fl_Window ();
	delete (Maps->Maps_Edit_Window);
	delete (Maps->MapsBuffer);
	delete (Maps);
	Maps = NULL;
    }
    if (Modparms) {
// we cannot simply call the destructor here as this menu contains spinners
	Fl::delete_widget (Modparms->Mods_Edit_Window);
	delete (Modparms);
	Modparms = NULL;
    }
    if (MiscConfigure) {
	//MiscConfigure->MiscConfigWindow->~Fl_Window ();
	delete (MiscConfigure->MiscConfigWindow);
	delete (MiscConfigure);
	MiscConfigure = NULL;
    }
    if (Polyhedra) {
	//Polyhedra->Polyhedra_Edit_Window->~Fl_Window ();
	delete (Polyhedra->Polyhedra_Edit_Window);
	delete (Polyhedra->PolyhedraBuffer);
	delete (Polyhedra->Polyhedra_Output_Buffer);
	delete (Polyhedra);
	Polyhedra = NULL;
    }
    if (Slabs) {
	//Slabs->SlabWindow->~Fl_Window ();
	delete (Slabs->SlabWindow);
	delete (Slabs);
	Slabs = NULL;
    }
    if (Spheres) {
	//Spheres->Sphere_Edit_Window->~Fl_Window ();
	delete (Spheres->Sphere_Edit_Window);
	if (Spheres->SphereBuffer) {
	    delete (Spheres->SphereBuffer);
	    Spheres->SphereBuffer = NULL;
	}
	delete (Spheres);
	Spheres = NULL;
    }
    if (listwindow1) {
	//listwindow1->~Fl_Window ();
	if (textbuf1) {
	    delete (textbuf1);
	    textbuf1 = NULL;
	}
	delete (listwindow1);
	listwindow1 = NULL;
    }
    if (listwindow2) {
	//listwindow2->~Fl_Window ();
	if (textbuf2) {
	    delete (textbuf2);
	    textbuf2 = NULL;
	}
	delete (listwindow2);
	listwindow2 = NULL;
    }
    if (listwindow3) {
	//listwindow3->~Fl_Window ();
	if (textbuf3) {
	    delete (textbuf3);
	    textbuf3 = NULL;
	}
	delete (listwindow3);
	listwindow3 = NULL;
    }
    if (helpwindow) {
	//helpwindow->~Fl_Window ();
	if (helpbuf) {
	    delete (helpbuf);
	    helpbuf = NULL;
	}
	delete (helpwindow);
	helpwindow = NULL;
    }
    if (helpwindow1) {
	//helpwindow1->~Fl_Window ();
	if (helpbuf1) {
	    delete (helpbuf1);
	    helpbuf1 = NULL;
	}
	delete (helpwindow1);
	helpwindow1 = NULL;
    }
    if (helpwindow2) {
	//helpwindow2->~Fl_Window ();
	delete (helpwindow2);
	helpwindow2 = NULL;
    }
    if (helpwindow3) {
	//helpwindow3->~Fl_Window ();
	delete (helpwindow3);
	helpwindow3 = NULL;
    }
    if (helpwindow4) {
	//helpwindow4->~Fl_Window ();
	delete (helpwindow4);
	helpwindow4 = NULL;
    }
    if (helpwindow5) {
	//helpwindow5->~Fl_Window ();
	delete (helpwindow5);
	helpwindow5 = NULL;
    }
    if (helpwindow6) {
	Fl::delete_widget (helpwindow6);
//        helpwindow6->~Fl_Window();
//        delete(helpwindow6);
	helpwindow6 = NULL;
    }
    if (textwindow) {
	//textwindow->~Fl_Window ();
	if (textbuf) {
	    delete (textbuf);
	    textbuf = NULL;
	}
	delete (textwindow);
	textwindow = NULL;
    }
    drvui->destroy = 0;
}

void
CrystalView::draw (void)
{
// routine to draw crystal object to GL Window
    static int first_time = 1;

    float cpx, cpy, cpz;

    float m[16];
    GLfloat mat_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat light_ambient[] = { AMBLEV, AMBLEV, AMBLEV, 1.0f };
    GLfloat light_specular[] = { 0.1f, 0.1f, 0.1f, 1.0f };
    GLfloat light_diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat mat_shininess[] = { 0.8f };
    GLfloat light_position[] = { 0.0f, 1.0f, 1.0f, 0.0f };
    GLfloat light_direction[] = { 0.0f, -1.0f, 0.0f };
    float ratio = (float) w () / (float) h ();

    if (!valid ()) {
	glLoadIdentity ();
	glViewport (0, 0, w (), h ());
	glEnable (GL_DEPTH_TEST);
    }
    if (first_time) {
	first_time = 0;
	drvui->crystalDL = glGenLists (1);
	drvui->frame_no = 1;
	if (strncmp (drvui->LoadOnStartup, "yes", 3) == 0 || strlen (drvui->Cur_Root)) {
	    Process_Inp (2);
	    Generate_Drawing (0);
	    Rotq = XYZ_Rot_to_Q (Xrot, Yrot, Zrot);
	    if (drvui->max_frame > 1) {
		Add_Frame_Main ();
		Main_Frame_Combo_cb (NULL, NULL);
	    } else {
		drvui->Frame_No->hide ();
	    }
	}
    }
    glMatrixMode (GL_PROJECTION);

    glLoadIdentity ();
    if (M_cameras == 0) {
	if (ratio <= 1.)
	    glOrtho (-gl_size, gl_size, -gl_size / ratio, gl_size / ratio, -10000.,
		     10000.);
	else
	    glOrtho (-gl_size * ratio, gl_size * ratio, -gl_size, gl_size, -10000.,
		     10000.);
    } else {
	gluPerspective (17, ratio, 0.01, 1000.);	// view angle, aspect,near/far clip
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
    glClearColor (drvui->glback[0], drvui->glback[1], drvui->glback[2], 0.0f);
    glClear ((GLbitfield) (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT));
    gluLookAt (cpx, cpy, Scale * .50,	// camera position
	       cpx, cpy, -1.0,	// camera lookat point
	       0.0f, 1.0f, 0.0f);	// camera "up" vector
    glColor3f (0.0f, 0.0f, 0.0f);
    glLightfv (GL_LIGHT0, GL_POSITION, light_position);
    glLightfv (GL_LIGHT0, GL_SPOT_DIRECTION, light_direction);
    glPushMatrix ();
    //   glTranslatef (gl_pos_x, gl_pos_y, gl_pos_z);
    Tb_Window::calculate (m);	// convert quaternion to matrix

    Yrot = -asin (m[2]) * RAD;	// convert to equivalent rotations
    double cosy = cos (Yrot / RAD);

    if (fabs (cosy) > 1.0e-3) {
	Xrot = atan2 (m[6] / cosy, m[10] / cosy) * RAD;
	Zrot = atan2 (m[1] / cosy, m[0] / cosy) * RAD;
    } else {
	Zrot = 0.0;
	Xrot = atan2 (-m[9], m[5]) * RAD;
    }
    char string[20];

    sprintf (string, "%6.2f", Zrot);	// load widgets on main screen
    drvui->Z_Rot->value (string);
    sprintf (string, "%6.2f", Yrot);
    drvui->Y_Rot->value (string);
    sprintf (string, "%6.2f", Xrot);
    drvui->X_Rot->value (string);

//    glRotatef ((float)Xrot, (float) 1., (float) 0., (float) 0.);
//    glRotatef ((float)Yrot, (float) 0., (float) 1., (float) 0.);
//    glRotatef ((float)Zrot, (float) 0., (float) 0., (float) 1.);
    glMultMatrixf (m);		// this line equivalent to 3 above
//    glTranslatef (-cpx,-cpy,-cpz); //this can mess up initial positioning w.r.t POV
    draw_cursor ();
    glCallList (drvui->crystalDL);
    int saved_frame_no = drvui->frame_no;

    for (drvui->frame_no = 1; drvui->frame_no <= drvui->max_frame; drvui->frame_no++)
	generate_gl_texts ();
    if (ShowMapLegend) 
	MapLegend();
    drvui->frame_no = saved_frame_no;
    glGetIntegerv (GL_VIEWPORT, viewport);	// save matrices for later
    glGetDoublev (GL_MODELVIEW_MATRIX, modelMatrix);
    glGetDoublev (GL_PROJECTION_MATRIX, projMatrix);
    glPopMatrix ();
    if (drvui->destroy) {	// close windows that are no longer needed
	if (drvui->destroy & LIST1) {
	    listwindow1->~Fl_Window ();
	    listwindow1 = NULL;
	    drvui->destroy &= !LIST1;
	}
	if (drvui->destroy & LIST2) {
	    listwindow2->~Fl_Window ();
	    listwindow2 = NULL;
	    if (textbuf2) {
		delete (textbuf2);
		textbuf2 = NULL;
	    }
	    drvui->destroy &= !LIST2;
	}
	if (drvui->destroy & LIST3) {
	    listwindow3->~Fl_Window ();
	    listwindow3 = NULL;
	    if (textbuf3) {
		delete (textbuf3);
		textbuf3 = NULL;
	    }
	    drvui->destroy &= !LIST3;
	}
	if (drvui->destroy & SPHERE) {
	    Spheres->Sphere_Edit_Window->~Fl_Window ();
	    Spheres->Sphere_Edit_Window = NULL;
	    if (Spheres->SphereBuffer) {
		delete (Spheres->SphereBuffer);
		Spheres->SphereBuffer = NULL;
	    }
	    if (Spheres->Sphere_Output_Buffer) {
		delete (Spheres->Sphere_Output_Buffer);
		Spheres->Sphere_Output_Buffer = NULL;
	    }
	    delete (Spheres);
	    Spheres = NULL;
	    drvui->destroy &= !SPHERE;
	}
	if (drvui->destroy & ELLIPSOID) {
	    ellipsoids->Ellips_Window->~Fl_Window ();
	    ellipsoids->Ellips_Window = NULL;
	    delete (ellipsoids->ColorInputBuf);
	    delete (ellipsoids);
	    ellipsoids = NULL;
	    drvui->destroy &= !ELLIPSOID;
	}
	if (drvui->destroy & LONEPAIR) {
	    LonePairs->LonePair_Edit_Window->~Fl_Window ();
	    delete (LonePairs->LonePairBuffer);
	    delete (LonePairs);
	    LonePairs = NULL;
	    drvui->destroy &= !LONEPAIR;
	}
	if (drvui->destroy & TEXT1) {
	    textwindow->~Fl_Window ();
	    textwindow = NULL;
	    delete (textbuf);
	    textbuf = NULL;
	    drvui->destroy &= !TEXT1;
	}
    }
    if (drvui->origin1_flag)
	drvui->Origin1_Msg->show ();
    else
	drvui->Origin1_Msg->hide ();
}

void
CrystalView::draw_overlay (void)
{
    static float cpx, cpy, cpz;

    float m[16];


    if (slabmode < 2)
	return;

    float ratio = 1.0f * (float) w () / (float) h ();

    glPushMatrix ();
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    if (M_cameras == 0) {
	if (ratio <= 1.)
	    glOrtho (-gl_size, gl_size, -gl_size / ratio, gl_size / ratio, -10000.,
		     10000.);
	else
	    glOrtho (-gl_size * ratio, gl_size * ratio, -gl_size, gl_size, -10000.,
		     10000.);
    } else {
	gluPerspective (17, ratio, 0.01, 1000.);
    }

    glMatrixMode (GL_MODELVIEW);
    Tb_Window::calculate (m);

    glMultMatrixf (m);
    cpx = (POV_Max[0] + POV_Min[0]) / 2.0f;
    cpy = (POV_Max[1] + POV_Min[1]) / 2.0f;
    cpz = (POV_Max[2] + POV_Min[2]) / 2.0f;
    glTranslatef (-cpx, -cpy, -cpz);

    show_slab_ovl ();
    glPopMatrix ();
}

void
Error_Box (const char *message)
{
    int yy = 100;

// routine for Error message windows

    if (errorbox) {
	errorbox->~Fl_Window ();
//    Fl::delete_widget(errorbox);
	errorbox = NULL;
    }

    drvui->msgbuffer= (char*) realloc(drvui->msgbuffer,(strlen(drvui->msgbuffer)+strlen(message)+2)*sizeof(char));
    strcat(drvui->msgbuffer,message);
    if (strstr (message, "\n"))
	yy = 140;
    errorbox = new Fl_Window (200, 180, 500, yy, "Error");
    errorbox->begin ();
    errorbox->callback ((Fl_Callback *) Error_Box_cb);

    int y = 20;

    Fl_Multiline_Output *box = new Fl_Multiline_Output (0, 10, 500, yy - 50, "");

    box->box (FL_NO_BOX);
    box->labelsize (14);
    box->labelcolor ((Fl_Color) 1);
    box->color (FL_GRAY);
    box->cursor_color (FL_BACKGROUND_COLOR);
    box->value (drvui->msgbuffer);
    y += 15 + yy - 80;
    Fl_Button *o = new Fl_Button (210, y, 80, 30, "Close");

    o->callback ((Fl_Callback *) Error_Box_cb);
#if !defined (WIN32) && !defined (__APPLE__)
    errorbox->icon ((char *) drvui->icon);
#endif
    errorbox->end ();
    errorbox->show ();
    strcat(drvui->msgbuffer,"\n");
//    errorbox->set_modal();
}

void
Error_Box_cb (Fl_Widget * w, void *)
{
// callback to hide error message window
    w->window ()->hide ();
    free (drvui->msgbuffer);
    drvui->msgbuffer = (char*) zalloc (sizeof (char));
}

void
Exit_cb (void)
{
// callback routine when exit menuitem is clicked or when main window is closed
    int r;

    if (strlen (drvui->Cur_File) && drvui->Str_File_Changed) {
	r = fl_choice ("The current str file has not been saved.\n"
		       "What action should be taken?", "Cancel", "Save", "Discard");
	if (!r)
	    return;
	if (r == 1)
	    Update_Str (1);
    }
    WriteConfig ();
    Clean_Up_Files ();		// delete all the temporary files
    Destroy_Open_Windows ();
    if (FourierPt)
	free (FourierPt);
    free (crystal);
    delete (Omit);
    if (s_vert)
	free (s_vert);
    if (o_vert)
	free (o_vert);
    if (o_vert_nm)
	free (o_vert_nm);
    if (drvui->vert_occ)
	free (drvui->vert_occ);
    if (xypos)
	free (xypos);
    if (xypos_nm)
	free (xypos_nm);
    if (vert_sym_no)
	free (vert_sym_no);
    if (vert_sym_nos)
	free (vert_sym_nos);
    if (vertex_list)
	free (vertex_list);
    if (poly_list)
	free (poly_list);
    if (drvui->atom_no)
	free (drvui->atom_no);
    if (drvui->atom_so)
	free (drvui->atom_so);
    if (drvui->orig_atom_no)
	free (drvui->orig_atom_no);
    if (drvui->modulate_x)
	free (drvui->modulate_x);
    if (drvui->modulate_3x)
	free (drvui->modulate_3x);
    if (drvui->modulate_3t)
	free (drvui->modulate_3t);
    if (drvui->atoms)
	free (drvui->atoms);
    if (drvui->ellips)
	free (drvui->ellips);
    if (drvui->fourier)
	free (drvui->fourier);
    if (drvui->labels)
	free (drvui->labels);
    if (drvui->bplanes)
	free (drvui->bplanes);
    if (drvui->planes)
	free (drvui->planes);
    if (drvui->polyhedra)
	free (drvui->polyhedra);
    if (drvui->polyedges)
	free (drvui->polyedges);
    if (drvui->spheres)
	free (drvui->spheres);
    if (drvui->atprops)
	free (drvui->atprops);
    if (drvui->modulate_gbl)
	free (drvui->modulate_gbl);
    if (drvui->frames)
	free (drvui->frames);
    if (drvui->nsurf > 1) {
	for (r = 1; r < drvui->nsurf; r++) {
	    free (drvui->surfx[r]);
	    free (drvui->surfy[r]);
	    free (drvui->surfz[r]);
	}
    }
    if (drvui->voidmap) {
	for (int j = 0; j < drvui->voidgrid[0]; j++) {
	    for (int k = 0; k < drvui->voidgrid[1]; k++)
		free (drvui->voidmap[j][k]);
	    free (drvui->voidmap[j]);
	}
	free (drvui->voidmap);
    }
    if (!strncmp (drvui->LoadOnStartup, "maybe", 5)) {
	strcpy (drvui->LoadOnStartup, "yes");
	WriteConfig ();
    }
    Fl::delete_widget (drvui->mainWindow);
    free (drvui->table);
    free (drvui->arrows);
    free (drvui->cones);
    free (drvui->bonds);
    free (drvui->msgbuffer);
    delete (drvui);
    exit (0);
}

int
Get_Unique_Atoms (char atoms[100][5], int Frame_No)
{
// routine to return the unique atom names in current frame
    int j;

    int mlist = 0;

    char atom1[5];

    FILE *aout;

    char string[256];

    strcpy (string, drvui->Cur_Root);
    strcat (string, ".frm");
    sprintf (atom1, "%d", Frame_No);
    strcat (string, atom1);
    if (!(aout = fopen (string, "r"))) {
	Error_Box ("Unable to open frm file.");
	return 0;
    }
    strcpy (string, "");
    while (!feof (aout)) {
	fgets (string, 100, aout);
	memset (atom1, 0, 5);
	if (strlen (string) < 4)
	    break;
	(void) sscanf (string, "%s", atom1);
	atom1[4] = 0;
	for (j = 3; j >= 0; --j) {
	    if (atom1[j] == ' ')
		atom1[j] = 0;
	}
	strncpy (atoms[mlist++], atom1, 5);	// add new one to end of list

	for (j = 0; j < mlist - 1; j++) {	// if not unique, truncate list
	    if (!strncmp (atoms[j], atom1, 4))
		--mlist;
	}
    }
    for (j = 0; j < mlist - 1; j++) {
	int k;

	for (k = j + 1; k < mlist; k++) {	// sort atom names
	    if (strncmp (atoms[k], atoms[j], 4) < 0) {
		strncpy (atom1, atoms[j], 4);
		strncpy (atoms[j], atoms[k], 4);
		strncpy (atoms[k], atom1, 4);
	    }
	}
    }
    fclose (aout);
    return mlist;
}

void
ImportDataFile_cb (Fl_Widget *, void *arg)
{
// Callback routine to import a structure file
    int r;

    int i;

    FILE *inp;

    FILE *newstr;

    static int one = 1;

    const char *phase;

    char string[2048];

    char tmp_dir[1024];

    char tmp_file[1024];

    char newfile2[1024];

    char Input_string[256];
    char string0[] = { "Select CIF Data File" };
    char string1[] = { "Select FDAT Data File" };
    char string2[] = { "Select GSAS Data File" };
    char string3[] = { "Select SCHAKAL Data File" };
    char string4[] = { "Select SHELX Data File" };
    char string5[] = { "Select WIEN2k Data File" };
    char string6[] = { "Select DISCUS Data File" };
    char string7[] = { "Select FULLPROF Data File" };
    char string8[] = { "Select EXCITING Data File" };

    switch ((long) arg) {
    case 0:
	strcpy (Input_string, string0);
	break;
    case 1:
	strcpy (Input_string, string1);
	break;
    case 2:
	strcpy (Input_string, string2);
	break;
    case 3:
	strcpy (Input_string, string3);
	break;
    case 4:
	strcpy (Input_string, string4);
	break;
    case 5:
	strcpy (Input_string, string5);
	break;
    case 6:
	strcpy (Input_string, string6);
	break;
    case 7:
	strcpy (Input_string, string7);
	break;
    case 8:
	strcpy (Input_string, string8);
	break;
    }
    if (strlen (drvui->Cur_File) && drvui->Str_File_Changed) {
	r = fl_choice ("The current str file has not been saved.\n"
		       "What action should be taken?", "Cancel", "Save", "Discard");
	if (!r)
	    return;
	if (r == 1)
	    Update_Str (1);
    }
#if defined(WIN32)
    char drive[_MAX_DRIVE];

    char dir[_MAX_DIR];

    char fname[_MAX_FNAME];

    char ext[_MAX_EXT];

    const char *newfile = flu_file_chooser (Input_string, "*", drvui->Cur_File);

    if (newfile) {
	_splitpath (newfile, drive, dir, fname, ext);	//Windows code
	strcpy (tmp_dir, drive);	// Drive letter
	strcat (tmp_dir, dir);	//   and directory
	strcpy (tmp_file, fname);	// copy file name
	strcat (tmp_file, ".str");	// and add extension
#else
    int k = 0;

    const char *newfile = fl_file_chooser (Input_string, "*", drvui->Cur_File);

    if (newfile) {
	strcpy (tmp_file, newfile);
	strcpy (tmp_dir, newfile);
	for (i = strlen (tmp_dir); i > 0; --i) {	// Find final / in file name
	    if (tmp_dir[i - 1] == '/') {
		tmp_dir[i] = 0;
		break;
	    }
	}
	for (i = strlen (tmp_dir); newfile[i] != 0; i++) {
	    tmp_file[k++] = newfile[i];	// copy file name to tmp_file
	}
	tmp_file[k] = 0;
	strcat (tmp_file, ".str");
#endif
	chdir (tmp_dir);	// switch to that directory
	newstr = fopen (tmp_file, "r");
	if (newstr) {
	    snprintf (string, 2048,
		      "The implied output file ('%s') will be overwritten.\n"
		      "Do you wish to continue?", tmp_file);
	    fclose (newstr);
	    if (!fl_choice ("%s", "No", "Yes", NULL, string)) {
		chdir (drvui->Cur_Dir);	// restore the original directory
		return;
	    }
	}
	if (!(inp = fopen (newfile, "r"))) {
	    snprintf (string, 2048, "The file you selected ('%s') cannot be read\n"
		      "Do you wish to continue?", newfile);
	    if (fl_choice ("%s", "No", "Yes", NULL, string)) {
		inp = fopen (tmp_file, "w");
		fclose (inp);
		Edit_STR_cb (NULL, &one);
	    }
	} else {
	    drvui->Str_File_Changed = 0;
	    Clean_Up_Files ();	// delete the temporary files,
	    Destroy_Open_Windows ();	// input dialogs and maps 
	    if (FourierPt) {	// associated with the previous
		free (FourierPt);	// structure
		FourierPt = NULL;
	    }
	    if (drvui->nsurf > 1) {
		for (i = 1; i < drvui->nsurf; i++) {
		    free (drvui->surfx[i]);
		    free (drvui->surfy[i]);
		    free (drvui->surfz[i]);
		    drvui->surfx[i] = NULL;
		    drvui->surfy[i] = NULL;
		    drvui->surfz[i] = NULL;
		}
	    }
	    strcpy (drvui->Cur_File, tmp_file);
	    strcpy (drvui->Cur_Dir, tmp_dir);
	    WriteConfig ();	// update configuration file
	    drvui->CurFile->value (drvui->Cur_File);	// update main screen widgets
	    drvui->CurDir->value (drvui->Cur_Dir);
	    newstr = fopen (drvui->Cur_File, "w");
	    if (!newstr) {
		snprintf (string, 256,
			  "The implied output file ('%s') could not be created.",
			  drvui->Cur_File);
		Error_Box (string);
		char *newfile =
		    fl_file_chooser ("Please choose a new location in a writable directory", "*.str",
				     drvui->Cur_File);
		if (!newfile || !(newstr=fopen(newfile,"w"))) {
		    fclose (inp);
		    return;
		} else {
		    strcpy(drvui->Cur_File,newfile);
		    strcpy(drvui->Cur_Dir,newfile);
		    for (i = strlen (drvui->Cur_Dir); i > 0; --i) {	// Find final / in file name
	    		if (drvui->Cur_Dir[i - 1] == '/') {
			    drvui->Cur_Dir[i] = 0;
			    break;
			}
		    }
		    WriteConfig ();	// update configuration file
		    drvui->CurFile->value (drvui->Cur_File);	// update main screen widgets
		    drvui->CurDir->value (drvui->Cur_Dir);
		}
	    }
	    fprintf (newstr, "titl imported\n");

	    if (strlen (newfile) > 80) {
// remove path component from filename if the full name plus import instruction 
// would exceed our maximum input line length of 100 bytes
		const char *start = strrchr (newfile, '/');

		if (!start) {
		    start = strrchr (newfile, '\\');	// look for the other Windows folder marker
		    if (!start) {
			Error_Box
			    ("Unable to find folder/directory markers in file string\n"
			     "Is your input file name really longer than 80 characters?\n"
			     "If so, please rename it.");
			fclose (inp);
			fclose (newstr);
			return;
		    }
		}
		start++;
		memmove (newfile2, start, strlen (start) + 1);
		if (strlen (newfile2) > 80) {
		    Error_Box
			("Your file name is longer than 80 characters. Please rename it.");
		    fclose (inp);
		    fclose (newstr);
		    return;
		}
	    } else {
		strcpy (newfile2, newfile);
	    }

	    switch ((long) arg) {
	    case 0:		// CIF file
		Block_CIF = 0;
//		fprintf (newstr, "import cif %s\n", newfile2);
		fprintf (newstr, "inline cif\n");
		while (!feof (inp)) {
		    fgets (string, 133, inp);
		    fprintf (newstr, "%s", string);
		}
		fprintf (newstr, "# End of data for %s\n", newfile2);
		drvui->auto_ellipse =
		    fl_choice ("Should ellipsoids be displayed automatically?", "No",
			       "Yes", NULL);
		break;
	    case 1:		// FDAT file
		fprintf (newstr, "import fdat %s\n", newfile2);
		break;
	    case 2:		// GSAS file
		phase = fl_input ("Enter GSAS phase number:", "1");
		if (!phase)
		    fprintf (newstr, "import gsas %s 1\n", newfile2);
		else
		    fprintf (newstr, "import gsas %s %s\n", newfile2, phase);
		break;
	    case 3:		// SCHAKAL file
		fprintf (newstr, "import schakal %s\n", newfile2);
		break;
	    case 4:		// process SHELX file
		drvui->auto_ellipse =
		    fl_choice ("Should ellipsoids be displayed automatically?", "No",
			       "Yes", NULL);
		fprintf (newstr, "inline shelx\n");
		while (!feof (inp)) {
		    fgets (string, 133, inp);
		    fprintf (newstr, "%s", string);
		    if (!strncmp (string, "HKLF", 4) || !strncmp (string, "hklf", 4)
			|| !strncmp (string, "END", 3) || !strncmp (string, "end", 3))
			break;
		}
		if (drvui->autolabel == 1)
		    drvui->autolabel = 2;
		break;
	    case 5:		// WIEN file
		fprintf (newstr, "import wien2k %s\n", newfile2);
		break;
	    case 6:		// DISCUS file
		fprintf (newstr, "import discus %s\n", newfile2);
		break;
	    case 7:		// FULLPROF file
		phase = fl_input ("Enter FULLPROF phase number:", "1");
		if (!phase)
		    fprintf (newstr, "import pcr %s 1\n", newfile2);
		else
		    fprintf (newstr, "import pcr %s %s\n", newfile2, phase);
		break;
	    case 8:		// Exciting file
		fprintf (newstr, "import exciting %s\n", newfile2);
		break;
	    }
	    fprintf (newstr, "box 0.02 Black\n");
	    fprintf (newstr, "background White\n");
	    fprintf (newstr, "view 0. 0. 0.\n");
	    fprintf (newstr, "pack  -0.05 1.05 -0.05 1.05 -0.05 1.05\n");
	    fprintf (newstr, "end\n");
	    fclose (newstr);
	    fclose (inp);
	}
	Process_Inp (2);
	Rotq = XYZ_Rot_to_Q (Xrot, Yrot, Zrot);
	Add_Frame_Main ();
	Generate_Drawing (0);	// generate this structure
	Fl::redraw ();		//  and draw it
    }
}

void
Include_Cutouts_cb (void)
{
// callback routine when Use_Cutouts checkbox is changed
    char string[128];

    drvui->Str_File_Changed = 1;
    if (ellipsoids->Use_Cutouts->value ()) {
	ellipsoids->Cutout_Color->activate ();
	drvui->El_Cutout = 1;
	strcpy (string, ellipsoids->Cutout_Color->value ());
	if (!strlen (string))
	    strcpy (string, (char *) "Gray20");
	ellipsoids->Cutout_Color->value (string);
    } else {
	drvui->El_Cutout = 0;
	ellipsoids->Cutout_Color->deactivate ();
    }
}

void
Load_Bond_Data (const char *atom, char table[20480])
{
// get sorted bond data for 'atom' from out file into character array table
    FILE *flin;

    int i;

    float f1, f2, f3, d;

    char filename[256];

    char string[100];

    char string1[20], string2[40], string3[20];

    char atom1[5];

    int nfound = 0;

    char atoms[1000][5];

    float d_list[1000];

    memset (string1, 0, 20);
    memset (string2, 0, 20);
    memset (string3, 0, 20);
    strcpy (filename, drvui->Cur_Root);	// generate name of '.out' file
    strcat (filename, ".out");
    if (!(flin = fopen (filename, "r"))) {
	strcpy (table, "Generated Bond Table\n  not available!");
	return;
    }

    table[0] = 0;
    strcpy (string, "    ");
    while (strncmp (string, " Bond Distance", 14)) {	// find bond-distance table
	if (!fgets (string, 100, flin)) {
	    strcpy (table, "Generated Bond Table\n  not available!");
	    fclose (flin);
	    return;
	}
    }
    while (strncmp (string, "End of File.", 12) && strncmp (string, " doing frame", 12)) {	// bond-distances over when end of file is reached
	int j;

	if (!fgets (string, 100, flin)) {
	    strcpy (table, "Generated Bond Table\n  not available!");
	    fclose (flin);
	    return;
	}
	sscanf (string, "%s %s %s", string1, string2, string3);
	if ((!strncmp (string1, "Distances", 9))
	    && (!strncmp (string3, atom, strlen (string3)))) {
	    if (!fgets (string, 100, flin)) {
		if (!strlen (table)) {
		    strcpy (table, "Generated Bond Table\n  not available!");
		}
		fclose (flin);
		return;
	    }
	    if (!fgets (string, 100, flin)) {
		if (!strlen (table)) {
		    strcpy (table, "Generated Bond Table\n  not available!");
		}
		fclose (flin);
		return;
	    }
	    while (strlen (string) > 10 && nfound < 1000) {
		sscanf (string, "%s %d %s %f %f %f %f", string1, &j, string2, &f1, &f2,
			&f3, &d);
		strncpy (atoms[nfound], string1, 4);
		atoms[nfound][4] = 0;
		d_list[nfound++] = d;
		if (!fgets (string, 100, flin)) {
		    if (!strlen (table)) {
			strcpy (table, "Generated Bond Table\n  not available!");
		    }
		    fclose (flin);
		    return;
		}
	    }
	}
    }
    for (i = 0; i < nfound - 1; i++) {
	int k;

	float tmp;

	for (k = i + 1; k < nfound; k++) {	// sort distances
	    if (d_list[k] < d_list[i]) {
		tmp = d_list[i];
		d_list[i] = d_list[k];
		d_list[k] = tmp;
		strncpy (atom1, atoms[i], 4);
		strncpy (atoms[i], atoms[k], 4);
		strncpy (atoms[k], atom1, 4);
	    }
	}
    }
    for (i = 0; i < nfound; i++) {
	sprintf (string, "%4s  %8.3f\n", atoms[i], d_list[i]);
	strcat (table, string);
    }

    if (!strlen (table)) {
	strcpy (table, "Generated Bond Table\n  not available!");
    }
    fclose (flin);
}

void
LoadConfig (bool full_check)
{
// load configuration file. If full_check is true, load all data, else just saved screen size
    FILE *fname;

    char filename[256];

    char temp[256];

    char stereo[256];

    char autolabel[256];

#ifdef WIN32
    char profile[100];

    OSVERSIONINFOEX osvi;

    memset (&osvi, 0, sizeof (OSVERSIONINFOEX));
    osvi.dwOSVersionInfoSize = sizeof (OSVERSIONINFOEX);
    GetVersionEx ((OSVERSIONINFO *) & osvi);
    if (osvi.dwMajorVersion == 4)
	strcpy (profile, "c:");
    else
	strcpy (profile, getenv ("USERPROFILE"));
    strcpy (filename, profile);
    strcat (filename, "\\");
    strcat (filename, Configure_file);
#else
    fl_filename_expand (filename, Configure_file);
#endif
    fname = fopen (filename, "r");
    if (!full_check) {
	if (fname) {
	    int i;

	    for (i = 0; i < 15; i++)
		fgets (temp, 256, fname);
	    if (strlen (temp) > 10) {
		int xpos, ypos, width, height;

		if (sscanf (temp, "%d %d %d %d", &xpos, &ypos, &width, &height) == 4) {
		    w_st_x = xpos;
		    w_st_y = ypos;
		    w_width = width;
		    w_height = height;
		}
	    }
	    fclose (fname);
	}
	return;
    }
    if (!fname) {
	strcpy (drvui->Cur_Dir,
		"Unable to open Configuration - Please use 'Configure' menu items.");
	strcpy (drvui->Cur_File, "");
	strcpy (drvui->DefaultFinish, "finish 0.70 0.30 0.08 0.01");
#ifdef WIN32			// Allow for different windows defaults
	strcpy (drvui->POV_Options, "+W600 +H600 +D +FC +A0.3 +P");
//      strcpy(drvui->POV_Include,"c:\\colors.inc");   // Is there a good default for windows?
//      strcpy(drvui->POV_Path,"povray.exe");
#else
	strcpy (drvui->POV_Options, "+W600 +H600 +D +FC +A0.3 +P");
#ifdef __APPLE__
	strcpy (drvui->POV_Include, "/sw/share/povray-3.5/include/colors.inc");
	strcpy (drvui->POV_Path, "/sw/bin/povray");
#else
	strcpy (drvui->POV_Include, "/usr/share/povray-3.7/include/colors.inc");
	strcpy (drvui->POV_Path, "/usr/bin/povray");
#endif
#endif
    } else {
	fgets (drvui->DRAWxtl_Path, 256, fname);
	trim_string (drvui->DRAWxtl_Path, 256);
	fgets (drvui->POV_Path, 256, fname);
	trim_string (drvui->POV_Path, 256);
	fgets (drvui->VRML_Path, 256, fname);
	trim_string (drvui->VRML_Path, 256);
	fgets (drvui->EditName, 256, fname);
	trim_string (drvui->EditName, 256);
	fgets (drvui->FileViewName, 256, fname);
	trim_string (drvui->FileViewName, 256);
	fgets (drvui->POV_Options, 256, fname);
	trim_string (drvui->POV_Options, 256);
	fgets (drvui->Cur_File, 256, fname);
	trim_string (drvui->Cur_File, 256);
	fgets (drvui->Cur_Dir, 256, fname);
	trim_string (drvui->Cur_Dir, 256);
	fgets (drvui->POV_Include, 256, fname);
	trim_string (drvui->POV_Include, 256);
	fgets (drvui->LoadOnStartup, 10, fname);
	trim_string (drvui->LoadOnStartup, 10);
	fgets (drvui->DefaultFinish, 256, fname);
	trim_string (drvui->DefaultFinish, 256);
	fgets (temp, 256, fname);	// skip over executable path
	fgets (autolabel, 256, fname);
	fgets (stereo, 256, fname);
	trim_string (stereo, 256);
	fgets (temp, 256, fname);	// skip over window dimensions
	temp[0] = '\0';
	fgets (temp, 256, fname);
	trim_string (temp, 256);
	strcpy (drvui->MSMS_Path, "");
	fgets (drvui->MSMS_Path, 256, fname);
	trim_string (drvui->MSMS_Path, 256);
	strcpy (drvui->Mencoder_Path, "");
	fgets (drvui->Mencoder_Path, 256, fname);
	trim_string (drvui->Mencoder_Path, 256);
	strcpy (drvui->FFmpeg_Path, "");
	fgets (drvui->FFmpeg_Path, 256, fname);
	trim_string (drvui->FFmpeg_Path, 256);
	fclose (fname);
    }
    if (!strncmp (drvui->LoadOnStartup, "maybe", 5)) {
	Error_Box ("Previous attempt to load a structure on startup failed.\n"
		   "LoadOnStartup cleared - Use Configure menu to reinstate.");
	strcpy (drvui->LoadOnStartup, "");
    } else if (!strlen (drvui->LoadOnStartup)) {
	strcpy (drvui->LoadOnStartup, "maybe");
    }

    sscanf (drvui->DefaultFinish, "%*s %f %f %f %f",
	    &drvui->ambient, &drvui->diffuse, &drvui->specular, &drvui->roughness);

    if (!strncmp (autolabel, "autolabel", 9))
	drvui->autolabel = 1;
    else
	drvui->autolabel = 0;
    if (strlen (stereo) > 1) {
	sscanf (stereo, "%d %d %f", &drvui->Stereo, &drvui->cross_eyed,
		&drvui->stereo_base);
    }
    if (strlen (temp) > 10) {
	sscanf (temp, "povray %d vrml %d asy %d", &doPOV, &doVrml, &doAsy);
    } else {
	doPOV = doVrml = doAsy = 1;
    }
    WriteConfig ();
}

void
Max_Min_cb (void)
{
// callback routine to readout Max - Min spinners
    static int in_progress;

    if (in_progress) {
#ifndef WIN32
	fprintf (stderr, "calculation already in progress\n");
#endif
	return;
    }
    in_progress = 1;
    int frame = atoi (drvui->Frame_No->value ());

    if (!frame)
	frame = 1;
    Omit->nomits = 0;
    xmin = drvui->X_Min->value ();
    drvui->frames[frame].cryst_lim[0] = (float) xmin;
    xmax = drvui->X_Max->value ();
    drvui->frames[frame].cryst_lim[3] = (float) xmax;
    ymin = drvui->Y_Min->value ();
    drvui->frames[frame].cryst_lim[1] = (float) ymin;
    ymax = drvui->Y_Max->value ();
    drvui->frames[frame].cryst_lim[4] = (float) ymax;
    zmin = drvui->Z_Min->value ();
    drvui->frames[frame].cryst_lim[2] = (float) zmin;
    zmax = drvui->Z_Max->value ();
    drvui->frames[frame].cryst_lim[5] = (float) zmax;
    if (drvui->Use_Clipping->value () != 0) {	// test check button
	clipflag = 1;		// check button set
	// save spinner widget values
	drvui->frames[frame].clip_lim[0] =
	    (float) max (drvui->X_Min_clip->value (), drvui->X_Min->value ());
	// MUST be within the min-max range
	drvui->frames[frame].clip_lim[0] =
	    (float) min (drvui->frames[frame].clip_lim[0], drvui->X_Max->value ());
	drvui->frames[frame].clip_lim[1] =
	    (float) max (drvui->Y_Min_clip->value (), drvui->Y_Min->value ());
	drvui->frames[frame].clip_lim[1] =
	    (float) min (drvui->frames[frame].clip_lim[1], drvui->Y_Max->value ());
	drvui->frames[frame].clip_lim[2] =
	    (float) max (drvui->Z_Min_clip->value (), drvui->Z_Min->value ());
	drvui->frames[frame].clip_lim[2] =
	    (float) min (drvui->frames[frame].clip_lim[2], drvui->Z_Max->value ());
	drvui->frames[frame].clip_lim[3] =
	    (float) min (drvui->X_Max_clip->value (), drvui->X_Max->value ());
	drvui->frames[frame].clip_lim[3] =
	    (float) max (drvui->frames[frame].clip_lim[3], drvui->X_Min->value ());
	drvui->frames[frame].clip_lim[4] =
	    (float) min (drvui->Y_Max_clip->value (), drvui->Y_Max->value ());
	drvui->frames[frame].clip_lim[4] =
	    (float) max (drvui->frames[frame].clip_lim[4], drvui->Y_Min->value ());
	drvui->frames[frame].clip_lim[5] =
	    (float) min (drvui->Z_Max_clip->value (), drvui->Z_Max->value ());
	drvui->frames[frame].clip_lim[5] =
	    (float) max (drvui->frames[frame].clip_lim[5], drvui->Z_Min->value ());
	drvui->X_Min_clip->value (drvui->frames[frame].clip_lim[0]);	// reload widgets with tested value
	drvui->Y_Min_clip->value (drvui->frames[frame].clip_lim[1]);
	drvui->Z_Min_clip->value (drvui->frames[frame].clip_lim[2]);
	drvui->X_Max_clip->value (drvui->frames[frame].clip_lim[3]);
	drvui->Y_Max_clip->value (drvui->frames[frame].clip_lim[4]);
	drvui->Z_Max_clip->value (drvui->frames[frame].clip_lim[5]);

	drvui->X_Min_clip->activate ();	// make widgets active
	drvui->X_Max_clip->activate ();
	drvui->Y_Min_clip->activate ();
	drvui->Y_Max_clip->activate ();
	drvui->Z_Min_clip->activate ();
	drvui->Z_Max_clip->activate ();
    } else {
	clipflag = 0;		// check button off - clipping disabled
	drvui->X_Min_clip->deactivate ();	// deactivate widgets
	drvui->X_Max_clip->deactivate ();
	drvui->Y_Min_clip->deactivate ();
	drvui->Y_Max_clip->deactivate ();
	drvui->Z_Min_clip->deactivate ();
	drvui->Z_Max_clip->deactivate ();
#if 0
	drvui->X_Min_clip->value (drvui->X_Min->value ());	//load widgets with min or max value
	drvui->X_Max_clip->value (drvui->X_Max->value ());
	drvui->Y_Min_clip->value (drvui->Y_Min->value ());
	drvui->Y_Max_clip->value (drvui->Y_Max->value ());
	drvui->Z_Min_clip->value (drvui->Z_Min->value ());
	drvui->Z_Max_clip->value (drvui->Z_Max->value ());
#endif
    }
    drvui->Str_File_Changed = 1;
    Update_Str (0);
    Generate_Drawing (0);
    Fl::redraw ();
    in_progress = 0;
}

void
Offset_cb (void)
{
// callback routine when Origin values are changed
    const char *string;

    Omit->nomits = 0;
    drvui->Str_File_Changed = 1;
    string = drvui->Origin_X->value ();
    origin[0] = (float) atof (string);
    string = drvui->Origin_Y->value ();
    origin[1] = (float) atof (string);
    string = drvui->Origin_Z->value ();
    origin[2] = (float) atof (string);
    Update_Str (0);
    Generate_Drawing (0);
    Fl::redraw ();
}

int
pick_box (int x, int y, int w, int h)
{
    GLint viewport[4];

    int hits;

    float cpx, cpy, cpz;

    float m[16];

    extern void show_slab (void);

    float ratio = 1.0f * (float) w / (float) h;

    memset (selectBuf, 0, BUFSIZE);
    glSelectBuffer (BUFSIZE, selectBuf);
    glRenderMode (GL_SELECT);

    glMatrixMode (GL_PROJECTION);
    glPushMatrix ();
    glLoadIdentity ();
    glGetIntegerv (GL_VIEWPORT, viewport);

    gluPickMatrix (x, viewport[3] - y, 1, 1, viewport);

    if (M_cameras == 0) {
	if (ratio <= 1.)
	    glOrtho (-gl_size, gl_size, -gl_size / ratio, gl_size / ratio, -10000.,
		     10000.);
	else
	    glOrtho (-gl_size * ratio, gl_size * ratio, -gl_size, gl_size, -10000.,
		     10000.);
    } else {
	gluPerspective (17, ratio, 0.01, 1000.);
    }

    glMatrixMode (GL_MODELVIEW);

    glTranslatef (drvui->Trans[0], drvui->Trans[1], drvui->Trans[2]);

    quaternion_to_rotmatrix (&Rotq, m);

    glMultMatrixf (m);

    cpx = (POV_Max[0] + POV_Min[0]) / 2.0f;
    cpy = (POV_Max[1] + POV_Min[1]) / 2.0f;
    cpz = (POV_Max[2] + POV_Min[2]) / 2.0f;

    glTranslatef (-cpx, -cpy, -cpz);

    glInitNames ();
    glPushName (~0);

    show_slab_ovl ();

    glMatrixMode (GL_PROJECTION);
    glMatrixMode (GL_MODELVIEW);
    glFlush ();
    glPopMatrix ();
    glFlush ();
    hits = glRenderMode (GL_RENDER);
    if (!hits)
	return (0);
    else
	hits = update_box (hits, selectBuf);

    return (hits);
}

void
process_hits (int hits, GLuint buffer[])
{
    int i;

    GLuint names, *ptr, minZ, *ptrNames = 0, numberOfNames;

    numberOfNames = 0;

    ptr = (GLuint *) buffer;
    minZ = 0xffffffff;
    for (i = 0; i < hits; i++) {
	names = *ptr;
	ptr++;
	if (*ptr < minZ && *(ptr + 2) != (GLuint) - 1) {
	    numberOfNames = names;
	    minZ = *ptr;
	    ptrNames = ptr + 2;
	}
	ptr += names + 2;
    }
    if (numberOfNames == 0)
	return;

    ptr = ptrNames;

    Omit->omit1[Omit->nomits] = *ptr;
    ptr++;
    Omit->omit2[Omit->nomits++] = *ptr;

    if (Omit->nomits == 1 && edtprm) {
	edtprm->ClearLastOmit->activate ();
	edtprm->ClearOmit->activate ();
    }
    Generate_Drawing (1);	/* regenerate drawing */
}

void
Process_Inp (int i)
{
// routine to process the DRAWxtl input file
    char string[60];

    Init_DRAWxtl ();		// initialize variables
    drvui->CurDir->value (drvui->Cur_Dir);
    if (strlen (drvui->Cur_File) < 2) {
	drvui->X_Min->value (-0.1);
	drvui->X_Max->value (1.1);
	drvui->Y_Min->value (-0.1);
	drvui->Y_Max->value (1.1);
	drvui->Z_Min->value (-0.1);
	drvui->Z_Max->value (1.1);
	xmin = -.1;
	xmax = 1.1;
	ymin = -.1;
	ymax = 1.1;
	zmin = -.1;
	zmax = 1.1;
	drvui->X_Origin = 0.5;
	drvui->Y_Origin = 0.5;
	drvui->Z_Origin = 0.5;
	drvui->X_Boxlim = 10.;
	drvui->Y_Boxlim = 10.;
	drvui->Z_Boxlim = 10.;
	drvui->Origin_X->value ("0.5");
	drvui->Origin_Y->value ("0.5");
	drvui->Origin_Z->value ("0.5");
	return;
    }
    memset (string, 0, 60);
    strcpy (string, drvui->Cur_File);
    if (string[strlen (string) - 1] == '_')
	string[strlen (string) - 1] = 0;
    drvui->CurFile->value (string);

    _chdir (drvui->Cur_Dir);

    drvui->fpin = fopen (drvui->Cur_File, "r");	// open 'str' file
    if (!drvui->fpin) {
	Error_Box ("Cannot open structure file !");
	strcpy (drvui->Cur_Temp, "");
	drvui->X_Min->value (-0.1);
	drvui->X_Max->value (1.1);
	drvui->Y_Min->value (-0.1);
	drvui->Y_Max->value (1.1);
	drvui->Z_Min->value (-0.1);
	drvui->Z_Max->value (1.1);
	xmin = -.1;
	xmax = 1.1;
	ymin = -.1;
	ymax = 1.1;
	zmin = -.1;
	zmax = 1.1;
	drvui->X_Origin = 0.5;
	drvui->Y_Origin = 0.5;
	drvui->Z_Origin = 0.5;
	drvui->X_Boxlim = 10.;
	drvui->Y_Boxlim = 10.;
	drvui->Z_Boxlim = 10.;
	drvui->Origin_X->value ("0.5");
	drvui->Origin_Y->value ("0.5");
	drvui->Origin_Z->value ("0.5");

	return;
    } else {
	drvui->frame_no = 1;
	while (1) {
	    read_inp (i);	// read the 'str' file
	    make_bmat (drvui->sys, drvui->lat_con, drvui->b_mat, drvui->ginv, drvui->rec_lat_con);	/* create the lattice metric */
	    rewind (drvui->fpin);
	    drvui->frame_no++;
	    if (drvui->frame_no > drvui->max_frame)
		break;
	}
	drvui->frame_no = drvui->max_frame;
	fclose (drvui->fpin);	// close str file
	strcpy (drvui->Cur_Listing, drvui->Cur_File);
	trim_string (drvui->Cur_Listing, 256);
	int j;

	for (j = strlen (drvui->Cur_Listing); j > 0; --j) {
	    if (drvui->Cur_Listing[j] == '.') {
		drvui->Cur_Listing[j] = 0;
		break;
	    }
	}
	strcpy (drvui->Cur_Root, drvui->Cur_Listing);
	strcpy (drvui->Cur_Console, drvui->Cur_Listing);
	strcpy (drvui->Cur_Temp, drvui->Cur_Listing);
	strcat (drvui->Cur_Listing, ".out");
	strcat (drvui->Cur_Console, ".cns");
	strcat (drvui->Cur_Temp, ".tmp");
	Fl_Text_Buffer *atextbuf = new Fl_Text_Buffer;

	atextbuf->loadfile (drvui->Cur_File);	// get text of str file
	if ( atextbuf->savefile (drvui->Cur_Temp)) {	// output it to temp file
	    Error_Box("Unable to create temporary str file - disk full or directory not writable");

            char *newfile =
		fl_file_chooser ("Please save this file in a writable directory", "*.str",
				drvui->Cur_File);
	    if (newfile) {
		atextbuf->savefile(newfile);

#if defined(WIN32)
		char drive[_MAX_DRIVE];

		char dir[_MAX_DIR];

		char fname[_MAX_FNAME];

		char ext[_MAX_EXT];

		_splitpath (newfile, drive, dir, fname, ext);       //Windows code
		strcpy (drvui->Cur_Dir, drive);     // Drive letter
		strcat (drvui->Cur_Dir, dir);       //   and directory
		strcpy (drvui->Cur_File, fname);    // copy file name
		strcat (drvui->Cur_File, ext);      // and add extension
#else
		strcpy(drvui->Cur_File,newfile);
		strcpy (drvui->Cur_Dir,drvui->Cur_Listing);
		for (j = strlen (drvui->Cur_Dir); j > 0; --j) {
		    if (drvui->Cur_Dir[j] == '/') {
			drvui->Cur_Dir[j] = 0;
			break;
		    }
		}
#endif

		strcpy (drvui->Cur_Listing, drvui->Cur_File);
		trim_string (drvui->Cur_Listing, 256);
		int j;

		for (j = strlen (drvui->Cur_Listing); j > 0; --j) {
		    if (drvui->Cur_Listing[j] == '.') {
			drvui->Cur_Listing[j] = 0;
			break;
		    }
		}
		strcpy (drvui->Cur_Root, drvui->Cur_Listing);
		strcpy (drvui->Cur_Console, drvui->Cur_Listing);
		strcpy (drvui->Cur_Temp, drvui->Cur_Listing);
		strcat (drvui->Cur_Listing, ".out");
		strcat (drvui->Cur_Console, ".cns");
		strcat (drvui->Cur_Temp, ".tmp");
		atextbuf->savefile(drvui->Cur_Temp);
		drvui->CurFile->value (drvui->Cur_File);
		drvui->CurDir->value (drvui->Cur_Dir);
	    }
	}
	delete (atextbuf);	// free the buffer
    }
    int frame = atoi (drvui->Frame_No->value ());

    if (frame == 0)
	frame = 1;
/* set limits of search regions */
    if (packflag == 0 && boxflag == 0) {	/* neither given */
	drvui->frames[frame].cryst_lim[0] = origin[0] - 0.5f;
	drvui->frames[frame].cryst_lim[1] = origin[1] - 0.5f;
	drvui->frames[frame].cryst_lim[2] = origin[2] - 0.5f;	/* set to generate one cell */
	drvui->frames[frame].cryst_lim[3] = origin[0] + 0.5f;
	drvui->frames[frame].cryst_lim[4] = origin[1] + 0.5f;
	drvui->frames[frame].cryst_lim[5] = origin[2] + 0.5f;
	packflag = 1;
    }
    xmin = drvui->frames[frame].cryst_lim[0];	// transfer DRAWxtl variables to widgets
    drvui->X_Min->value (xmin);
    xmax = drvui->frames[frame].cryst_lim[3];
    drvui->X_Max->value (xmax);
    ymin = drvui->frames[frame].cryst_lim[1];
    drvui->Y_Min->value (ymin);
    ymax = drvui->frames[frame].cryst_lim[4];
    drvui->Y_Max->value (ymax);
    zmin = drvui->frames[frame].cryst_lim[2];
    drvui->Z_Min->value (zmin);
    zmax = drvui->frames[frame].cryst_lim[5];
    drvui->Z_Max->value (zmax);
    if (clipflag == 1) {	// clip active - load clip widgets
	drvui->Use_Clipping->value (1);
	drvui->X_Min_clip->value (drvui->frames[frame].clip_lim[0]);
	drvui->Y_Min_clip->value (drvui->frames[frame].clip_lim[1]);
	drvui->Z_Min_clip->value (drvui->frames[frame].clip_lim[2]);
	drvui->X_Max_clip->value (drvui->frames[frame].clip_lim[3]);
	drvui->Y_Max_clip->value (drvui->frames[frame].clip_lim[4]);
	drvui->Z_Max_clip->value (drvui->frames[frame].clip_lim[5]);

	drvui->X_Min_clip->activate ();	// activate them
	drvui->X_Max_clip->activate ();
	drvui->Y_Min_clip->activate ();
	drvui->Y_Max_clip->activate ();
	drvui->Z_Min_clip->activate ();
	drvui->Z_Max_clip->activate ();
    } else {			// no clipping - load widgets with min and max values
	drvui->Use_Clipping->value (0);	//  and deactivate them
	drvui->X_Min_clip->deactivate ();
	drvui->X_Max_clip->deactivate ();
	drvui->Y_Min_clip->deactivate ();
	drvui->Y_Max_clip->deactivate ();
	drvui->Z_Min_clip->deactivate ();
	drvui->Z_Max_clip->deactivate ();
	drvui->X_Min_clip->value (drvui->X_Min->value ());
	drvui->X_Max_clip->value (drvui->X_Max->value ());
	drvui->Y_Min_clip->value (drvui->Y_Min->value ());
	drvui->Y_Max_clip->value (drvui->Y_Max->value ());
	drvui->Z_Min_clip->value (drvui->Z_Min->value ());
	drvui->Z_Max_clip->value (drvui->Z_Max->value ());
    }
    Xrot = xrot;
    sprintf (string, "%6.2f", Xrot);
    drvui->X_Rot->value (string);
    sprintf (string, "%6.2f", Yrot);
    Yrot = yrot;
    drvui->Y_Rot->value (string);
    Zrot = zrot;
    sprintf (string, "%6.2f", Zrot);
    drvui->Z_Rot->value (string);
    sprintf (string, "%6.2f", origin[0]);
    drvui->Origin_X->value (string);
    sprintf (string, "%6.2f", origin[1]);
    drvui->Origin_Y->value (string);
    sprintf (string, "%6.2f", origin[2]);
    drvui->Origin_Z->value (string);
    if (!Vrml2)
	drvui->Generate_VRML1->set ();	// VRML1 checkbox is set
    else
	drvui->Generate_VRML1->clear ();
    if (!doVrml)
	drvui->Generate_VRML1->deactivate ();
    if (!doPOV)
	DRAWxtlViewUI::drawxtl_menu[39].deactivate ();
    if (!M_cameras)
	drvui->Orthographic_View->set ();	// orthographic view checkbox set
    else
	drvui->Orthographic_View->clear ();
    if (Display_axes)
	drvui->Show_Vector_Triple->set ();	// Show vector triple checkbox set
    else
	drvui->Show_Vector_Triple->clear ();

    Update_Str (0);
}

void
Progress_Window (int type, const char *label, float max)
{

// Routine to handle a progress bar. The first argument defines the type of call
//      type == -1, open a window. The second argument gives the label, 3rd the maximum.
//      type == -2, close the window. 2nd and 3rd arguments ignored.
//      type == 0, 3rd argument contains progress argument. 2nd ignored.

    static Fl_Window *prog_window;

    static Fl_Progress *prog_bar;

    switch (type) {

    case -1:
	prog_window = new Fl_Window (400, 300, 250, 50, "Please Wait");
	prog_window->begin ();
	prog_bar = new Fl_Progress (20, 10, 200, 30, label);
	prog_bar->minimum (0.0f);
	prog_bar->maximum (max);
#if !defined (WIN32) && !defined (__APPLE__)
	prog_window->icon ((char *) drvui->icon);
#endif
	prog_window->end ();
	prog_window->show ();
	break;
    case -2:
	prog_window->~Fl_Window ();
	delete(prog_window);
	break;
    default:
	prog_bar->value (max);
    }
    Fl::check ();
}

void
Restore_Working_Copy (void)
{
// restore working str file from saved version
    Fl_Text_Buffer *text = new Fl_Text_Buffer;

    char savefile[256];

    strcpy (savefile, drvui->Cur_Root);
    strcat (savefile, ".save");
    text->loadfile (savefile);
    text->savefile (drvui->Cur_Temp);
    delete (text);
}

void
Rotation_cb (void)
{
// callback routine to readout Rotation widgets
    Xrot = atof (drvui->X_Rot->value ());
    Xrot = AngleInRange (Xrot);
    xrot = (float) Xrot;
    Yrot = atof (drvui->Y_Rot->value ());
    Yrot = AngleInRange (Yrot);
    yrot = (float) Yrot;
    Zrot = atof (drvui->Z_Rot->value ());
    Zrot = AngleInRange (Zrot);
    zrot = (float) Zrot;
    Update_Str (0);
    Rotq = XYZ_Rot_to_Q (Xrot, Yrot, Zrot);
    Fl::redraw ();
}

void
Save_Current_cb (void)
{
// callback routine to save the current .str and .pov under a new name
    char oldpovfile[1024], newpovfile[1024];

    char old_dir[1024], old_file[1024], new_file[1024];

    int old_state;

#if defined(WIN32)
    char drive[_MAX_DRIVE];

    char dir[_MAX_DIR];

    char fname[_MAX_FNAME];

    char ext[_MAX_EXT];

    const char *newfile;
#else
    int i;

    char *newfile;
#endif

    if (!strlen (drvui->CurDir->value ())) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }
    strcpy (old_dir, drvui->Cur_Dir);
    strcpy (old_file, drvui->Cur_File);
    strcpy (oldpovfile, drvui->Cur_Root);
    strcat (oldpovfile, ".pov");

#if defined(WIN32)
    newfile = flu_file_chooser ("Select New Name for Current .STR/.POV Files",
				"*", drvui->Cur_Dir);
    if (newfile) {
	_splitpath (newfile, drive, dir, fname, ext);	//Windows code
	strcpy (drvui->Cur_Dir, drive);	// Drive letter
	strcat (drvui->Cur_Dir, dir);	//   and directory
	strcpy (drvui->Cur_File, fname);	// copy file name
    } else
	return;
#else
    newfile =
	fl_file_chooser ("Select New Name for Current .STR/.POV Files", "*",
			 drvui->Cur_Dir);
    if (newfile) {
	strcpy (drvui->Cur_File, newfile);
	strcpy (drvui->Cur_Dir, newfile);
	for (i = strlen (drvui->Cur_Dir); i > 0; --i) {	// Find final / in file name
	    if (drvui->Cur_Dir[i - 1] == '/') {
		drvui->Cur_Dir[i] = 0;
		break;
	    }
	}
	for (i = strlen (drvui->Cur_File); i > 0; --i) {	// find final . in file name
	    if (drvui->Cur_File[i - 1] == '.') {
		drvui->Cur_File[i - 1] = 0;
		break;
	    }
	}
    } else
	return;
#endif

    strcpy (newpovfile, drvui->Cur_File);
    strcat (newpovfile, ".pov");
    rename (oldpovfile, newpovfile);

    strcpy (new_file, drvui->Cur_File);
    strcat (new_file, ".str");
    strcpy (drvui->Cur_File, new_file);
    old_state = drvui->Str_File_Changed;
    Update_Str (1);
    drvui->Str_File_Changed = old_state;
    strcpy (drvui->Cur_File, old_file);
    strcpy (drvui->Cur_Dir, old_dir);
}

void
Save_Working_Copy (void)
{
// copy working version of str file into a "save" file
    Fl_Text_Buffer *text = new Fl_Text_Buffer;

    char savefile[256];

    strcpy (savefile, drvui->Cur_Root);
    strcat (savefile, ".save");
    text->loadfile (drvui->Cur_Temp);
    text->savefile (savefile);
    delete (text);
}

void
SelectDataFile_cb (void)
{
// Callback routine to select the data file
    int r, i;

    FILE *inp;

    static int one = 1;

    if (strlen (drvui->Cur_File) && drvui->Str_File_Changed) {
	r = fl_choice ("The current str file has not been saved.\n"
		       "What action should be taken?", "Cancel", "Save", "Discard");
	if (!r)
	    return;
	if (r == 1)
	    Update_Str (1);
    }

    drvui->Str_File_Changed = 0;

#if defined(WIN32)
    char drive[_MAX_DRIVE];

    char dir[_MAX_DIR];

    char fname[_MAX_FNAME];

    char ext[_MAX_EXT];

    const char *newfile = flu_file_chooser ("Select Data/Experiment File",
					    "*.str", drvui->Cur_File);

    if (newfile) {
	_splitpath (newfile, drive, dir, fname, ext);	//Windows code
	strcpy (drvui->Cur_Dir, drive);	// Drive letter
	strcat (drvui->Cur_Dir, dir);	//   and directory
	strcpy (drvui->Cur_File, fname);	// copy file name
	strcat (drvui->Cur_File, ext);	// and add extension
#else
    int k = 0;

    char *newfile =
	fl_file_chooser ("Select Data/Experiment File", "*.str", drvui->Cur_File);
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

	Clean_Up_Files ();	// get rid of the temporary files
	Destroy_Open_Windows ();

	if (FourierPt) {
	    free (FourierPt);
	    FourierPt = NULL;
	}
	for (i = 0; i < 4; i++) {
	    cur_atom[i] = 0;
	    strcpy (cur_name[i], "");
	}
	cur_show = 0;
	ReadFourMap = 0;
	strcpy (FourierFileName, "");
	drvui->Cursor_pos->value ("");
	drvui->labels_inited = 0;	// make sure the labels get reinited

	chdir (drvui->Cur_Dir);	// switch to that directory
	WriteConfig ();		// update configuration file
	drvui->CurFile->value (drvui->Cur_File);	// update main screen widgets
	drvui->CurDir->value (drvui->Cur_Dir);
	if (!(inp = fopen (drvui->Cur_File, "r"))) {
	    char string[256];

	    sprintf (string, "The file you selected ('%s') cannot be read\n"
		     "Do you wish to continue?", drvui->Cur_File);
	    if (fl_choice ("%s", "No", "Yes", NULL, string)) {
		inp = fopen (drvui->Cur_File, "w");
		fclose (inp);
		Edit_STR_cb (NULL, &one);
	    }
	} else {
	    fclose (inp);
	}
	vzero (drvui->Trans);
	drvui->origin1_flag = 0;
	gl_size = 0.0f;
	Process_Inp (2);
	Rotq = XYZ_Rot_to_Q (Xrot, Yrot, Zrot);
	Add_Frame_Main ();
	Fl::flush ();
	Generate_Drawing (0);	// generate this structure
	Fl::redraw ();		//  and draw it
    }
}

void
show_slab_ovl (void)
{
    glPushMatrix ();

    glDisable (GL_LIGHTING);
    glLineWidth (3);
    glColor3f (0., 1., 0.);
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
    glPushMatrix ();
    glLoadName (1);
    glTranslatef (slabx1, slaby1, slabz1);
    glutSolidCube (0.004 * Scale);
    glPopMatrix ();
    glPushMatrix ();
    glLoadName (2);
    glTranslatef (slabx2, slaby2, slabz2);
    glutSolidCube (0.004 * Scale);
    glPopMatrix ();
    glPushMatrix ();
    glLoadName (3);
    glTranslatef (slabx3, slaby3, slabz3);
    glutSolidCube (0.004 * Scale);
    glPopMatrix ();
    glPushMatrix ();
    glLoadName (4);
    glTranslatef (slabx4, slaby4, slabz4);
    glutSolidCube (0.004 * Scale);
    glPopMatrix ();

    glEnable (GL_LIGHTING);
    glPopMatrix ();
    glLineWidth (1);
}

void
start_picking (int x, int y, int w, int h)
{
    GLint viewport[4];

    int hits;

    float cpx, cpy, cpz;

    float m[16];

    int n, i;

    float ratio = 1.0f * (float) w / (float) h;

    memset (selectBuf, 0, BUFSIZE);
    glGetIntegerv (GL_VIEWPORT, viewport);

    glPushMatrix ();
    glSelectBuffer (BUFSIZE, selectBuf);
    glRenderMode (GL_SELECT);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();

    gluPickMatrix (x, viewport[3] - y, viewport[2] / 100, viewport[3] / 100, viewport);

    if (M_cameras == 0) {
	if (w <= h)
	    glOrtho (-gl_size, gl_size, -gl_size / ratio, gl_size / ratio, -10000.,
		     10000.);
	else
	    glOrtho (-gl_size * ratio, gl_size * ratio, -gl_size, gl_size, -10000.,
		     10000.);
    } else {
	gluPerspective (17, ratio, 0.01, 1000.);	// view angle, aspect,near/far clip
    }
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity ();

    glInitNames ();
    glPushName (~0);

    cpx = (POV_Max[0] + POV_Min[0]) / 2.0f;
    cpy = (POV_Max[1] + POV_Min[1]) / 2.0f;
    cpz = (POV_Max[2] + POV_Min[2]) / 2.0f;

    gluLookAt (cpx, cpy, Scale * .50,   // camera position
               cpx, cpy, -1.0,          // camera lookat point
               0.0f, 1.0f, 0.0f);       // camera "up" vector

    glTranslatef (drvui->Trans[0], drvui->Trans[1], drvui->Trans[2]);

    quaternion_to_rotmatrix (&Rotq, m);

    glMultMatrixf (m);

    glCallList (drvui->crystalDL);

    glPopMatrix ();

    hits = glRenderMode (GL_RENDER);
    glFlush ();

    if (hits > 0) {
	process_hits (hits, selectBuf);
    } else {
	hits = pick_label (x, y, Fl::w (), Fl::h ());
	if (hits > 0 && hits <= drvui->nlabel) {
	    if (hits != drvui->triple[0] && hits != drvui->triple[1]
		&& hits != drvui->triple[2] && hits != drvui->triple[3]) {
		for (n = hits; n < drvui->nlabel; n++) {	// remove label, but not a, b, c, or origin
		    drvui->labels[n].label_fn = drvui->labels[n + 1].label_fn;
		    strcpy (drvui->labels[n].label_label,
			    drvui->labels[n + 1].label_label);
		    for (i = 0; i < 3; i++)
			drvui->labels[n].label_x[i] = drvui->labels[n + 1].label_x[i];
		    for (i = 0; i < 4; i++)
			if (drvui->triple[i] == n + 1)
			    drvui->triple[i] = n;
		}
		drvui->nlabel--;
	    }
	    Update_Str (0);
	}
    }

    Fl::redraw ();

}

int
update_box (int hits, GLuint buffer[])
{
    int i, j;

    GLuint names, *ptr, minZ, *ptrNames = 0, numberOfNames;

    numberOfNames = 0;
    j = 0;

    ptr = (GLuint *) buffer;
    minZ = 0xffffffff;
    for (i = 0; i < hits; i++) {
	names = *ptr;
	ptr++;
	if (*ptr < minZ && *(ptr + 2) != (GLuint) - 1) {
	    numberOfNames = names;
	    minZ = *ptr;
	    ptrNames = ptr + 2;
	}
	ptr += names + 2;
    }

    if (numberOfNames == 0)
	return 0;

    ptr = ptrNames;
    j = *ptr;

    return (j);
}

void
Update_Objects (int Frame_No, FILE * out)
{
// Update the objects that belong to this frame
    int i;

    unsigned int j;

    char string[256];

// add bond lines
    if (drvui->nbond > 1) {
	for (i = 1; i < drvui->nbond; i++) {
	    char atom1[5], atom2[5];

	    if (drvui->bonds[i].bond_fn == Frame_No) {
		strncpy (string, drvui->bonds[i].col_bond, 39);
		string[39] = 0;
		for (j = 0; j < strlen (string); j++)
		    if (string[j] < ' ')
			string[j] = 0;
		strncpy (atom1, drvui->bonds[i].bond_l1, 4);
		strncpy (atom2, drvui->bonds[i].bond_l2, 4);
		atom1[4] = 0;
		atom2[4] = 0;
		if (drvui->bonds[i].bond_style == 0)
		    fprintf (out, "bond %s %s %6.3f %6.3f %6.3f %s\n", atom1,
			     atom2, drvui->bonds[i].bond_size, drvui->bonds[i].bond_min,
			     drvui->bonds[i].bond_max, string);
		else {
		    if (drvui->bonds[i].bond_style != 5)
			fprintf (out, "dash %d %s %s %6.3f %6.3f %6.3f %s\n",
				 drvui->bonds[i].bond_style, atom1,
				 atom2, drvui->bonds[i].bond_size,
				 drvui->bonds[i].bond_min, drvui->bonds[i].bond_max,
				 string);
		    else
			fprintf (out, "dash %s %s %6.3f %6.3f %6.3f %s\n", atom1,
				 atom2, drvui->bonds[i].bond_size,
				 drvui->bonds[i].bond_min, drvui->bonds[i].bond_max,
				 string);
		}
	    }
	}
    }
// add sphere lines
    if (drvui->nsphere > 1) {
	for (i = 1; i < drvui->nsphere; i++) {
	    if (drvui->spheres[i].sphere_fn == Frame_No) {
		if (drvui->spheres[i].sphere_n == -1) {
		    fprintf (out, "sphere %s %6.3f %s\n", drvui->spheres[i].sphere_l,
			     drvui->spheres[i].sphere_size, drvui->spheres[i].sphere_col);
		} else {
		    fprintf (out, "sphere %s %d %6.3f %s\n", drvui->spheres[i].sphere_l,
			     drvui->spheres[i].sphere_n, drvui->spheres[i].sphere_size,
			     drvui->spheres[i].sphere_col);
		}
	    }
	}
    }
// add magnetic arrow commands
    if (drvui->nmag > 0) {
	for (i = 0; i < drvui->nmag; i++) {
	    if (drvui->arrows[i].arrow_fn == Frame_No) {
		fprintf (out,
			 "arrow %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %s\n",
			 drvui->arrows[i].mag_xp[0], drvui->arrows[i].mag_xp[1],
			 drvui->arrows[i].mag_xp[2], drvui->arrows[i].mag_xc[0],
			 drvui->arrows[i].mag_xc[1], drvui->arrows[i].mag_xc[2],
			 drvui->arrows[i].arrow_length, drvui->arrows[i].arrow_diam,
			 drvui->arrows[i].col_arrow);
	    }
	}
// add magnetic transformation command
	fprintf (out, "mag_trans %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
		 drvui->mag_matrix[0][0], drvui->mag_matrix[0][1],
		 drvui->mag_matrix[0][2], drvui->mag_matrix[1][0],
		 drvui->mag_matrix[1][1], drvui->mag_matrix[1][2],
		 drvui->mag_matrix[2][0], drvui->mag_matrix[2][1],
		 drvui->mag_matrix[2][2]);
    }
// add various types of polyhedra commands
    if (drvui->npoly > 1) {
	for (i = 1; i < drvui->npoly; i++) {
	    if (drvui->polyhedra[i].poly_fn == Frame_No) {
		char atom1[5], atom2[5];

		strncpy (string, drvui->polyhedra[i].poly_col, 39);
		string[39] = 0;
		for (j = 0; j < strlen (string); j++)
		    if (string[j] < ' ')
			string[j] = 0;
		strncpy (atom2, drvui->polyhedra[i].poly_t, 4);
		strncpy (atom1, drvui->polyhedra[i].poly_l, 4);
		atom1[4] = 0;
		atom2[4] = 0;
		for (j = 0; j < strlen (atom1); j++)
		    if (atom1[j] <= ' ')
			atom1[j] = 0;
		for (j = 0; j < strlen (atom2); j++)
		    if (atom2[j] <= ' ')
			atom2[j] = 0;
		if (strlen (atom2)) {	// polyvert command if 'to' atom
		    fprintf (out, "polyvert %s %s %6.3f %s\n", atom1,
			     atom2, drvui->polyhedra[i].poly_size, string);
		} else if (drvui->polyhedra[i].poly_min > 0.005) {	// shell command if minimum distance specified
		    fprintf (out, "shell %s %6.3f %6.3f %s\n", atom1,
			     drvui->polyhedra[i].poly_min, drvui->polyhedra[i].poly_size,
			     string);
		} else {	// polysz command
		    fprintf (out, "polysz %s %6.3f %s\n", atom1,
			     drvui->polyhedra[i].poly_size, string);
		}
	    }
	}
// add polyedge commands 
	for (i = 1; i < drvui->npoly; i++) {
	    if (drvui->polyhedra[i].poly_fn == Frame_No) {
		if (drvui->polyhedra[i].poly_rad_edge > 0.0f) {
		    char atom1[5];

		    strncpy (atom1, drvui->polyhedra[i].poly_l, 4);
		    atom1[4] = 0;
		    for (j = 0; j < strlen (atom1); j++)
			if (atom1[j] <= ' ')
			    atom1[j] = 0;
		    fprintf (out, "polyedge %s %6.3f %s\n", atom1,
			     drvui->polyhedra[i].poly_rad_edge,
			     drvui->polyhedra[i].poly_col_edge);
		}
	    }
	}
    }
// add plane commands
    if (drvui->nplane > 1) {
	for (i = 1; i < drvui->nplane; i++) {
	    if (drvui->planes[i].plane_fn == Frame_No) {
		char atom1[5];

		strncpy (string, drvui->planes[i].plane_col, 39);
		string[39] = 0;
		for (j = 0; j < strlen (string); j++)
		    if (string[j] < ' ')
			string[j] = 0;
		strncpy (atom1, drvui->planes[i].plane_l, 4);
		atom1[4] = 0;
		for (j = 0; j < strlen (atom1); j++)
		    if (atom1[j] <= ' ')
			atom1[j] = 0;
		fprintf (out, "plane %s %6.3f %s\n", atom1, drvui->planes[i].plane_size,
			 string);
	    }
	}
    }
// add lonepair command(s)
    if (drvui->ncone > 1) {
	for (i = 1; i < drvui->ncone; i++) {
	    char atom[5];

	    if (drvui->cones[i].cone_fn == Frame_No) {
		strncpy (string, drvui->cones[i].col_cone, 39);
		string[39] = 0;
		for (j = 0; j < strlen (string); j++)
		    if (string[j] < ' ')
			string[j] = 0;
		strncpy (atom, drvui->cones[i].cone_l1, 4);
		atom[4] = 0;
		for (j = 0; j < strlen (atom); j++)
		    if (atom[j] <= ' ')
			atom[j] = 0;
		fprintf (out, "lonepair %s %d %6.3f %6.3f %6.3f %s\n", atom,
			 drvui->cones[i].numlonepairs, drvui->cones[i].cone_height,
			 drvui->cones[i].cone_min, drvui->cones[i].cone_max, string);
	    }
	}
    }
// add labeltext command(s)
    if (drvui->nlabel > 1) {
	for (i = 1; i < drvui->nlabel; i++) {
	    if (drvui->labels[i].label_fn == Frame_No) {
		if (!strcmp (drvui->labels[i].label_label, "triple_vect")) {
		    fprintf (out, "vectors %.3f %.3f %.3f\n", offset[0], offset[1],
			     offset[2]);
		} else {
		    fprintf (out, "labeltext %6.3f %6.3f %6.3f %s\n",
			     drvui->labels[i].label_x[0], drvui->labels[i].label_x[1],
			     drvui->labels[i].label_x[2], drvui->labels[i].label_label);
		}
	    }
	}
    }
// add pack commands for this frame
    if (packflag) {
	fprintf (out, "pack %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
		 drvui->frames[Frame_No].cryst_lim[0],
		 drvui->frames[Frame_No].cryst_lim[3],
		 drvui->frames[Frame_No].cryst_lim[1],
		 drvui->frames[Frame_No].cryst_lim[4],
		 drvui->frames[Frame_No].cryst_lim[2],
		 drvui->frames[Frame_No].cryst_lim[5]);
    }
// add clip commands for this frame
    if (clipflag) {
	fprintf (out, "clip %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
		 drvui->frames[Frame_No].clip_lim[0], drvui->frames[Frame_No].clip_lim[3],
		 drvui->frames[Frame_No].clip_lim[1], drvui->frames[Frame_No].clip_lim[4],
		 drvui->frames[Frame_No].clip_lim[2],
		 drvui->frames[Frame_No].clip_lim[5]);
    }
    if (Map_Info.info_valid) {
	if (drvui->modulated != 0)
	    fprintf (out, "mapregion %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
		     drvui->frames[Frame_No].map_lim[0], drvui->frames[Frame_No].map_lim[3],
		     drvui->frames[Frame_No].map_lim[1], drvui->frames[Frame_No].map_lim[4],
		     drvui->frames[Frame_No].map_lim[2], drvui->frames[Frame_No].map_lim[5],
		     drvui->frames[Frame_No].map_lim[6], drvui->frames[Frame_No].map_lim[7],
		     drvui->frames[Frame_No].map_lim[8]);
	else
	    fprintf (out, "mapregion %.4f %.4f %.4f %.4f %.4f %.4f\n", 
		     drvui->frames[Frame_No].map_lim[0], drvui->frames[Frame_No].map_lim[3],
		     drvui->frames[Frame_No].map_lim[1], drvui->frames[Frame_No].map_lim[4],
		     drvui->frames[Frame_No].map_lim[2], drvui->frames[Frame_No].map_lim[5]);
    	if (drvui->frames[Frame_No].slice >0) 
	    fprintf (out, "mapslice %.3f %.3f %.3f  %.3f %.3f %.3f  %d\n",
		     drvui->frames[Frame_No].mapslice[0], drvui->frames[Frame_No].mapslice[1], 
		     drvui->frames[Frame_No].mapslice[2], drvui->frames[Frame_No].mapnorm[0], 
		     drvui->frames[Frame_No].mapnorm[1], drvui->frames[Frame_No].mapnorm[2],
		     drvui->frames[Frame_No].slice);
    }

// add atomic property value lines
    if (drvui->natprop > 1) {
	for (i = 1; i < drvui->natprop; i++) {
	    if (drvui->atprops[i].atprop_fn == Frame_No) {
		if (drvui->atprops[i].atprop_n == -1) {
		    fprintf (out, "values %s * %6.3f\n", drvui->atprops[i].atprop_l,
			     drvui->atprops[i].radius);
		} else {
		    fprintf (out, "values %s %d %6.3f\n", drvui->atprops[i].atprop_l,
			     drvui->atprops[i].atprop_n, drvui->atprops[i].radius);
		}
	    }
	}
    }
}

void
Update_Str (int overwrite)
{
// routine to update the tmp file with the new values from the widgets
// if 'overwrite' is true, delete the str file and move (rename) the tmp file to str
    FILE *inp;

    FILE *out;

    char string[256];

    char temp_out[256];

    static const char *surftypes[3] = { "mesh", "solid", "dots"};

    int i;

    int Upd_magnification = 1;

    int Upd_depthcue = 1;

    int Upd_molcomp = 1;

    int Upd_ellipsoids = 1;

    int Upd_cutout = 1;

    int Upd_list = 1;

    int Upd_vrml = 1;

    int Upd_orthographic = 1;

    int Upd_axislines = 1;

    int Upd_box = 1;

    int Upd_edges = 1;

    int Upd_polytolerance = 1;

    int Upd_background = 1;

    int Frame_No;

    int Upd_finish = 1;

    int Upd_xyzoff = 1;

    int copy_end = 0;

    Frame_No = 1;
    inp = fopen (drvui->Cur_Temp, "r");	// open the existing file
    if (!inp)
	return;

    strcpy (temp_out, drvui->Cur_Temp);
    strcat (temp_out, "a");
    if (!(out = fopen (temp_out, "w"))) {	// open a temp file
	fclose (inp);
	Error_Box ("Unable to open temporary file.");
	return;
    }
    while (1) {
	fgets (string, 256, inp);
	trim_string (string, 256);
	if (feof (inp) != 0)
	    break;

	/* see if we need to preserve the terminal END of an inlined dataset */
	if (strstr (string, "inline") || strstr (string, "INLINE")) {
	    if (strstr (string, "shelx")
		|| strstr (string, "shakal")
		|| strstr (string, "schakal"))
		copy_end = 1;
	    if (strstr (string, "SHELX")
		|| strstr (string, "SHAKAL")
		|| strstr (string, "SCHAKAL"))
		copy_end = 1;
	}
	if (copy_end == 1) {	// for shelx data, hklf is as good as end
	    if (strstr (string, "HKLF") || strstr (string, "hklf"))
		copy_end = 0;
	}
	if (!strncmp (string, "end", 3) || !strncmp (string, "END", 3)) {
	    if (copy_end == 0)
		strcpy (string, "");
	    else
		copy_end = 0;
	}
	if (strncmp (string, "frame", 5) == 0) {
	    Update_Objects (Frame_No, out);	// update stuff for this frame
	    Frame_No++;
	}
//
// update lines containing information contained in widgets
//      magnification, depthcue, molcomp, ellipsoids, cutout, list, pack 
//      view, vrml, orthographic, vectors, axislines, box, edges
//      phong, origin, polytolerance, background, clip, bond, special
//      polysz, plane, polyvert, shell, sphere, lonepair, slab, noshadow
//      arrow, mag_trans, ellipcolor, dash, phong, mapcontour, mapcontour2d,
//      labeltext, mapread, mapregion, xyzoff, labelscale, import cif,
//      lookat -> view

	if (strncmp (string, "magnification", 13) == 0) {
	    if (fabs (Magnification - 1.0) > 0.001) {
		sprintf (string, "magnification %6.2f", Magnification);
		Upd_magnification = 0;
	    } else {
		strcpy (string, "");
	    }
	}
	if (strncmp (string, "depthcue", 8) == 0) {
	    if (DepthCue > 0.001) {
		sprintf (string, "depthcue %6.3f", DepthCue);
		Upd_depthcue = 0;
	    } else {
		strcpy (string, "");
	    }
	}
	if (strncmp (string, "molcomp", 7) == 0) {
	    if (drvui->mol_d > 0.0f) {
		sprintf (string, "molcomp %6.2f", drvui->mol_d);
		Upd_molcomp = 0;
	    } else
		strcpy (string, "");
	}
	if (strncmp (string, "ellipsoids", 10) == 0) {
	    if (drvui->do_ellipsoids)
		sprintf (string, "ellipsoids %6.2f", drvui->Ellipsoid_Prob);
	    else
		strcpy (string, "");
	    Upd_ellipsoids = 0;
	}
	if (strncmp (string, "cutout", 6) == 0) {
	    if (strlen (drvui->Cutout_color) != 0) {
		sprintf (string, "cutout %s", drvui->Cutout_color);
		Upd_cutout = 0;
	    } else
		strcpy (string, "");
	}
	if (strncmp (string, "list", 4) == 0) {
	    if (fabs (printdist - 3.5) > 0.001) {
		sprintf (string, "list %6.2f", printdist);
		Upd_list = 0;
	    } else
		strcpy (string, "");
	}
	if (strncmp (string, "pack", 4) == 0) {
	    strcpy (string, "");
	}
	if (strncmp (string, "slab", 4) == 0) {
	    strcpy (string, "");
	}
	if (strncmp (string, "view", 4) == 0) {
	    strcpy (string, "");
	}
	if (strncmp (string, "vrml", 4) == 0) {
	    Upd_vrml = 0;
	    if (Vrml2 == 0)
		sprintf (string, "vrml1");
	    else
		sprintf (string, "vrml97");
	}
	if (strncmp (string, "orthographic", 12) == 0) {
	    if (M_cameras == 0) {
		sprintf (string, "orthographic");
		Upd_orthographic = 0;
	    } else
		strcpy (string, "");
	}
	if (strncmp (string, "vectors", 7) == 0) {
	    strcpy (string, "");
	}
	if (strncmp (string, "polyedge", 8) == 0) {
	    strcpy (string, "");
	}
	if (strncmp (string, "axislines", 9) == 0) {
	    sprintf (string, "axislines %6.2f %s", drvui->Ellipaxis_width,
		     drvui->Ellipaxis_color);
	    Upd_axislines = 0;
	}
	if (strncmp (string, "box", 3) == 0) {
	    sprintf (string, "box %6.3f %s", rad_cell, drvui->col_cell);
	    Upd_box = 0;
	}
	if (strncmp (string, "edges", 5) == 0) {
	    if (edges) {
		sprintf (string, "edges %6.2f %s", drvui->rad_edge, drvui->col_edge);
		Upd_edges = 0;
	    } else
		strcpy (string, "");
	}
	if (strncmp (string, "origin", 6) == 0)
	    strcpy (string, "");
	if (strncmp (string, "polytolerance", 13) == 0) {
	    if (drvui->polylimit > 0.105) {
		sprintf (string, "polytolerance %6.2f", drvui->polylimit);
		Upd_polytolerance = 0;
	    } else
		strcpy (string, "");
	}
	if (strncmp (string, "background", 10) == 0) {
	    sprintf (string, "background %s", drvui->col_bg);
	    Upd_background = 0;
	}
	if (strncmp (string, "finish", 6) == 0) {
	    sprintf (string, "finish %6.2f %6.2f %6.2f %6.2f", drvui->ambient,
		     drvui->diffuse, drvui->specular, drvui->roughness);
	    Upd_finish = 0;
	}
	if (strncmp (string, "clip", 4) == 0)
	    strcpy (string, "");
	if (strncmp (string, "bond", 4) == 0)
	    strcpy (string, "");
	if (strncmp (string, "dash", 4) == 0)
	    strcpy (string, "");
	if (strncmp (string, "special", 7) == 0)
	    strcpy (string, "");
	if (strncmp (string, "polysz", 6) == 0)
	    strcpy (string, "");
	if (strncmp (string, "plane", 5) == 0)
	    strcpy (string, "");
	if (strncmp (string, "polyvert", 8) == 0)
	    strcpy (string, "");
	if (strncmp (string, "shell", 5) == 0)
	    strcpy (string, "");
	if (strncmp (string, "sphere", 6) == 0)
	    strcpy (string, "");
	if (strncmp (string, "lonepair", 8) == 0)
	    strcpy (string, "");
	if (strncmp (string, "noshadow", 8) == 0)
	    strcpy (string, "");
	if (strncmp (string, "arrow", 5) == 0)
	    strcpy (string, "");
	if (strncmp (string, "mag_trans", 9) == 0)
	    strcpy (string, "");
	if (strncmp (string, "ellipcolor", 10) == 0)
	    strcpy (string, "");
	if (strncmp (string, "phong", 5) == 0)
	    strcpy (string, "");
	if (strncmp (string, "mapcontour", 10) == 0)
	    strcpy (string, "");
	if (strncmp (string, "labeltext", 9) == 0)
	    strcpy (string, "");
	if (strncmp (string, "mapread", 7) == 0)
	    strcpy (string, "");
	if (strncmp (string, "mapregion", 9) == 0)
	    strcpy (string, "");
	if (strncmp (string, "mapslice", 8) == 0)
	    strcpy (string, "");
	if (strncmp (string, "xyzoff", 6) == 0) {
	    Upd_xyzoff = 0;
	    drvui->origin1_flag = 0;
	}
	if (strncmp (string, "labelscale", 10) == 0) {
	    strcpy (string, "");
	}
	if (strncmp (string, "lookat", 6) == 0) {
	    strcpy (string, "");
	}
	if (strncmp (string, "nolabels", 8) == 0) {
	    strcpy (string, "");
	}
	if (strncmp (string, "average", 7) == 0) {
	    strcpy (string, "");
	}
	if (strncmp (string, "phaseshift", 10) == 0) {
	    strcpy (string, "");
	}
	if (strncmp (string, "occupancy", 9) == 0) {
	    strcpy (string, "");
	}
	if (strncmp (string, "aimsurf", 7) == 0) {
	    strcpy (string, "");
	}
	if (strncmp (string, "voids", 5) == 0) {
	    strcpy (string, "");
	}
	if (strncmp (string, "values", 6) == 0) {
	    strcpy (string, "");
	}
	if (strncmp (string, "import", 6) == 0 && strstr (string, "cif")) {
	    char string2[5];

	    i = strlen (string) - 1;
	    if (isdigit (string[i])) {	// if number at end 
		for (int j = i; j > i - 6; j--) {
		    if (isspace (string[j])) {	// strip off number
			string[j] = '\0';
			break;
		    }
		    if (isalpha (string[j]))
			break;	// unless it is attached to filename
		}
	    }
	    sprintf (string2, " %d", Block_CIF);
	    strcat (string, string2);	// add block number for CIF
	}
	if (strncmp (string, "rem xyzoff", 10) == 0)
	    strcpy (string, "");
	if (strncmp (string, "rem - following lines indicate", 30) == 0)
	    strcpy (string, "");
	if ((i = strlen (string)) > 0)
	    if (string[i - 1] <= 13)
		string[i - 1] = 0;
	if ((i = strlen (string)) > 0) {
	    if (string[i - 1] <= 13)
		string[i - 1] = 0;
	    fprintf (out, "%s\n", string);
	}
    }
    fclose (inp);
    Update_Objects (Frame_No, out);	// update last frame stuff
//
// Add parameters that are or may be generated by widgets but NOT in input file
//
    if (drvui->noshadow)
	fprintf (out, "noshadow\n");
    if (Labels == 0)
	fprintf (out, "nolabels\n");
    if (fabs (Magnification - 1.0) > 0.001 && Upd_magnification)
	fprintf (out, "magnification %6.2f\n", Magnification);
    if (DepthCue > 0.001 && Upd_depthcue)
	fprintf (out, "depthcue %6.3f\n", DepthCue);
    if (drvui->mol_d > 0.0f && Upd_molcomp)
	fprintf (out, "molcomp %6.2f\n", drvui->mol_d);
    if (drvui->do_ellipsoids && Upd_ellipsoids)
	fprintf (out, "ellipsoids %6.2f\n", drvui->Ellipsoid_Prob);
    if (strlen (drvui->Cutout_color) != 0 && Upd_cutout)
	fprintf (out, "cutout %s\n", drvui->Cutout_color);
    if (fabs (printdist - 3.5) > 0.001 && Upd_list)
	fprintf (out, "list %6.2f\n", printdist);
    if (fabs (drvui->label_scale - 1.0) > 0.001)
	fprintf (out, "labelscale %.3f\n", drvui->label_scale);
    if (drvui->slab_con[0] > 0.)
	fprintf (out, "slab %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f "
		 "%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %d\n",
		 drvui->slab_con[0], drvui->slab_con[1], drvui->slab_con[2],
		 drvui->slab_con[3], drvui->slab_con[4], drvui->slab_con[5],
		 drvui->slab_off[0], drvui->slab_off[1], drvui->slab_off[2],
		 drvui->slab_rot[0], drvui->slab_rot[1], drvui->slab_rot[2], slabmode);
    fprintf (out, "view %6.1f %6.1f %6.1f\n", Xrot, Yrot, Zrot);
    if (Vrml2 == 0 && Upd_vrml == 1)
	fprintf (out, "vrml1\n");
    if (M_cameras == 0 && Upd_orthographic)
	fprintf (out, "orthographic\n");
    if ((drvui->Phong_Size != 1.0f) || (drvui->Phong_Value != 0.2f))
	fprintf (out, "phong %6.2f %6.2f\n", drvui->Phong_Value, drvui->Phong_Size);
    if (Upd_axislines)
	fprintf (out, "axislines %6.2f %s\n", drvui->Ellipaxis_width,
		 drvui->Ellipaxis_color);
    if (Upd_box)
	fprintf (out, "box %6.3f %s\n", rad_cell, drvui->col_cell);
    if (edges && Upd_edges)
	fprintf (out, "edges %6.2f %s\n", drvui->rad_edge, drvui->col_edge);
    if ((fabs (origin[0] - 0.5) > 0.0001) || (fabs (origin[1] - 0.5) > 0.0001) ||
	(fabs (origin[2] - 0.5) > 0.0001))
	fprintf (out, "origin %6.2f %6.2f %6.2f\n", origin[0], origin[1], origin[2]);
    if (drvui->polylimit > 0.105 && Upd_polytolerance)
	fprintf (out, "polytolerance %6.2f\n", drvui->polylimit);
    if (Upd_background)
	fprintf (out, "background %s\n", drvui->col_bg);
    if (Upd_finish)
	fprintf (out, "finish %6.2f %6.2f %6.2f %6.2f\n", drvui->ambient,
		 drvui->diffuse, drvui->specular, drvui->roughness);
    if (drvui->do_ellipsoids) {
	for (i = 1; i < drvui->n_ellips; i++) {
	    if (drvui->ellips[i].ell_type > 1000) {
		if (drvui->ellips[i].save_el_number != -1)
		    fprintf (out, "ellipcolor %s %d %s\n", drvui->ellips[i].ellips_l,
			     drvui->ellips[i].ellips_n, drvui->ellips[i].ellips_col);
		else {
		    int j, haveit;

		    haveit = 0;
		    for (j = 1; j < i; j++)
			if (check_atom_name
			    (drvui->ellips[i].ellips_l, drvui->ellips[j].ellips_l))
			    haveit = 1;
		    if (haveit == 0)
			fprintf (out, "ellipcolor %s * %s\n", drvui->ellips[i].ellips_l,
				 drvui->ellips[i].ellips_col);
		}
	    }
	}
    }
// add mapcontour and mapcontour2d commands
    if (drvui->numOfFourierContours) {
	for (i = 1; i <= drvui->numOfFourierContours; i++) {
	    char type[6];

	    if (drvui->Fourier2d) {
		fprintf (out, "mapcontour2d %.3f %.3f %.3f %s",
			 drvui->fourier[i].FourierContourLevel,
			 drvui->fourier[i].FourierContourStep,
			 drvui->fourier[i].FourierContourTop,
			 drvui->fourier[i].FourierContourColor);
		if (strlen(drvui->fourier[i].FourierBackColor))
		    fprintf(out," %s\n",drvui->fourier[i].FourierBackColor);
		else
		    fprintf(out, "\n");
	    } else {
		if (drvui->fourier[i].FourierContourSolid)
		    strcpy (type, "solid");
		else
		    strcpy (type, "mesh");
		fprintf (out, "mapcontour  %.3f %s %s\n",
			 drvui->fourier[i].FourierContourLevel, type,
			 drvui->fourier[i].FourierContourColor);
	    }
	}
    }
// add mapread command if needed
    switch (FourierMapType) {
    case 1:
	strcpy (string, "mapread grd ");
	break;
    case 2:
	strcpy (string, "mapread stf ");
	break;
    case 3:
	strcpy (string, "mapread w2k ");
	break;
    case 4:
	strcpy (string, "mapread vsp ");
	break;
    case 5:
	strcpy (string, "mapread flp ");
	break;
    case 6:
	strcpy (string, "mapread fcf ");
	break;
    case 7:
	strcpy (string, "mapread dn6 ");
	break;
    case 8:
	strcpy (string, "mapread m80 ");
	break;
    case 9:
	strcpy (string, "mapread exc ");
	break;
    case 10:
	strcpy (string, "mapread m81 ");
	break;
    case 11:
	strcpy (string, "mapread xsf ");
	break;
    default:
	strcpy (string, "");
	break;
    }
    if (strlen (string) > 0) {
	char res[20];
	strcat (string, FourierFileName);
	if (FourierMapType == 6 || FourierMapType == 8) {	// add calc type
	    if (Map_Info.map_type == 1)
		strcat (string, " Fc");
	    else if (Map_Info.map_type == 2)
		strcat (string, " Fo-Fc");
	    else if (Map_Info.map_type == 3)
		strcat (string, " 2Fo-Fc");
	    else if (Map_Info.map_type == 4)
		strcat (string, " Fo2");
	    else
		strcat (string, " Fo");
	}
	sprintf(res, " %i", Map_Info.res);
	strcat(string, res);
	fprintf (out, "%s\n", string);
    }
// add options and parameters related to modulation
    if (drvui->modulated != 0) {
	if (drvui->modulated < 0)
	    fprintf (out, "average\n");
	if (drvui->phaseshift[0] + drvui->phaseshift[1] + drvui->phaseshift[2] > 0.0)
	    fprintf (out, "phaseshift %.4f %.4f %.4f\n",
		     drvui->phaseshift[0], drvui->phaseshift[1], drvui->phaseshift[2]);
	for (i = 0; i < natom; i++)
	    if (drvui->atoms[i].min_occ > 0.)
		fprintf (out, "occupancy %s %d %.2f %.2f\n", drvui->atoms[i].atom_l,
			 drvui->atoms[i].sv_atom_n, drvui->atoms[i].occupancy,
			 drvui->atoms[i].min_occ);
    }
// add surface-related options and parameters
    if (drvui->nsurf >1 ) {
	for (i = 1; i < drvui->nsurf; i++) 
	    fprintf(out,"aimsurf %s %d %s %s %s\n", drvui->surfatom[i], drvui->surfnum[i],
		    drvui->surffile[i], surftypes[drvui->surftype[i]], drvui->surfcolor[i]);
    }
    if (drvui->voidflag != 0) {
	fprintf (out,"voids %d %.4f %d %d %d %s\n", drvui->voidflag, drvui->probesize,
		 drvui->voidgrid[0], drvui->voidgrid[1], drvui->voidgrid[2], drvui->voidcolor);
    }
// add special command (if needed)
    if (Omit->nomits > 0) {
	for (i = 0; i < Omit->nomits; i++) {
	    fprintf (out, "special %d %d\n", Omit->omit1[i], Omit->omit2[i]);
	}
    }
    if (Upd_xyzoff && drvui->origin1_flag)
	fprintf (out, "rem xyzoff %.3f %.3f %.3f\n", -drvui->origin_offset[0],
		 -drvui->origin_offset[1], -drvui->origin_offset[2]);
    if (fabs (drvui->Old_Xrot - Xrot) > 0.1)
	drvui->Str_File_Changed = 1;
    if (fabs (drvui->Old_Yrot - Yrot) > 0.1)
	drvui->Str_File_Changed = 1;
    if (fabs (drvui->Old_Zrot - Zrot) > 0.1)
	drvui->Str_File_Changed = 1;
    fprintf (out, "end\n");

    if (fclose (out) != 0) {
	Error_Box ("Unable to update the working file!");
	return;
    }

    if (!overwrite) {
	unlink (drvui->Cur_Temp);	// delete the input file
	if (rename (temp_out, drvui->Cur_Temp)) {	// rename temporary to working
	    Error_Box ("Unable to update working file!");
	    return;
	}
    } else {
	unlink (drvui->Cur_File);	// delete the str file
	if (rename (temp_out, drvui->Cur_File)) {	// rename tmp to str
	    Error_Box ("Unable to create revised str file!");
	    return;
	}
	drvui->Str_File_Changed = 0;
    }
}

void
WriteConfig (void)
{
// write updated configuration file
    FILE *fname;

    char filename[256];

#ifdef WIN32
    char profile[100];

    OSVERSIONINFOEX osvi;

    memset (&osvi, 0, sizeof (OSVERSIONINFOEX));
    osvi.dwOSVersionInfoSize = sizeof (OSVERSIONINFOEX);
    GetVersionEx ((OSVERSIONINFO *) & osvi);
    if (osvi.dwMajorVersion == 4)
	strcpy (profile, "c:");
    else
	strcpy (profile, getenv ("USERPROFILE"));
    strcpy (filename, profile);
    strcat (filename, "\\");
    strcat (filename, Configure_file);
#else
    fl_filename_expand (filename, Configure_file);
#endif
    fname = fopen (filename, "w");
    if (!fname) {
	return;
    }
    fprintf (fname, "%s\n%s\n", drvui->DRAWxtl_Path, drvui->POV_Path);
    fprintf (fname, "%s\n%s\n", drvui->VRML_Path, drvui->EditName);
    fprintf (fname, "%s\n%s\n", drvui->FileViewName, drvui->POV_Options);
    fprintf (fname, "%s\n%s\n", drvui->Cur_File, drvui->Cur_Dir);
    fprintf (fname, "%s\n", drvui->POV_Include);
    fprintf (fname, "%s\n", drvui->LoadOnStartup);
    fprintf (fname, "%s\n", drvui->DefaultFinish);
    fprintf (fname, "%s\n", drvui->ProgramPath);
    if (drvui->autolabel == 0)
	fprintf (fname, "noautolabel\n");
    else
	fprintf (fname, "autolabel\n");
    fprintf (fname, "%d %d %.3f\n", drvui->Stereo, drvui->cross_eyed, drvui->stereo_base);
    fprintf (fname, "%d %d %d %d\n", drvui->mainWindow->x (), drvui->mainWindow->y (),
	     drvui->mainWindow->w (), drvui->mainWindow->h ());
    fprintf (fname, "povray %d vrml %d asy %d\n", doPOV, doVrml, doAsy);
    fprintf (fname, "%s\n", drvui->MSMS_Path);
    fprintf (fname, "%s\n", drvui->Mencoder_Path);
    fprintf (fname, "%s\n", drvui->FFmpeg_Path);
    fclose (fname);
}

QUAT
XYZ_Rot_to_Q (double Xrot, double Yrot, double Zrot)
{
// routine to convert the X,Y,Z rotations into a quaternion
    float *axis;

    QUAT q1, q2, q3, qtemp, Rotq;

    axis = (float *) zalloc (3 * sizeof (float));
    axis[0] = 1.0f;
    axis_to_quaternion (axis, (float) (Xrot / RAD), &q1);
    axis[0] = axis[2] = 0.0;
    axis[1] = 1.0f;
    axis_to_quaternion (axis, (float) (Yrot / RAD), &q2);
    axis[2] = 1.0f;
    axis[1] = axis[0] = 0.0;
    axis_to_quaternion (axis, (float) (Zrot / RAD), &q3);
    qmult (&q3, &q2, &qtemp);
    qmult (&qtemp, &q1, &Rotq);
    free (axis);
    return Rotq;
}

#if 0
void
moveto_atom (int x, int y, int w, int h)
{
// picks the sphere at the current cursor position, moves the crosshair to it
// and updates distance and angle from previous atom(s)

    int i, j;

    GLint viewport[4];

    GLuint names, *ptr, minZ;

    GLuint *ptrNames = 0, numberOfNames = 0;

    int hits;

    int sphere, num;

    float cpx, cpy, cpz;

    float m[16];

    float dot;

    char atnum[5];

    int ellips;

    char cur_name_t[10];

    float ratio = 1.0f * (float) w / (float) h;

    memset (selectBuf, 0, BUFSIZE);
    glGetIntegerv (GL_VIEWPORT, viewport);
    glPushMatrix ();

    glSelectBuffer (BUFSIZE, selectBuf);
    glRenderMode (GL_SELECT);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();

    gluPickMatrix (x, viewport[3] - y, viewport[2] / 100, viewport[3] / 100, viewport);

    if (M_cameras == 0) {
	if (w <= h)
	    glOrtho (-gl_size, gl_size, -gl_size / ratio, gl_size / ratio, -10000.,
		     10000.);
	else
	    glOrtho (-gl_size * ratio, gl_size * ratio, -gl_size, gl_size, -10000.,
		     10000.);
    } else {
	gluPerspective (17, ratio, 0.01, 1000.);	// view angle, aspect,near/far clip
    }
    glMatrixMode (GL_MODELVIEW);

    glInitNames ();
    glPushName (~0);

    glTranslatef (drvui->Trans[0], drvui->Trans[1], drvui->Trans[2]);

    quaternion_to_rotmatrix (&Rotq, m);

    glMultMatrixf (m);
    cpx = (POV_Max[0] + POV_Min[0]) / 2.0f;
    cpy = (POV_Max[1] + POV_Min[1]) / 2.0f;
    cpz = (POV_Max[2] + POV_Min[2]) / 2.0f;


    glCallList (drvui->crystalDL);

    glPopMatrix ();

    hits = glRenderMode (GL_RENDER);
    glFlush ();

    if (hits == 0)
	return;

    ptr = (GLuint *) selectBuf;
    minZ = 0xffffffff;
    for (i = 0; i < hits; i++) {
	names = *ptr;
	ptr++;
	if (*ptr < minZ && *(ptr + 2) != (GLuint) - 1) {
	    numberOfNames = names;
	    minZ = *ptr;
	    ptrNames = ptr + 2;
	}
	ptr += names + 2;
    }
    if (numberOfNames == 0)
	return;

    ptr = ptrNames;

    sphere = *ptr;
    ptr++;
    num = *ptr;

    ellips = 0;

    if (sphere > 90000 && sphere / 100000 < drvui->n_ellips) {
	ellips = 1;
	sphere /= 100000;
    } else if (sphere >= drvui->nsphere)
	return;			// not a sphere


    nvert = 0;			// clear the vertex list
    for (j = 0; j < natom; ++j) {	// loop through atoms 
	if (drvui->atoms[j].atom_fn != drvui->frame_no)
	    continue;
	if ((drvui->atoms[j].atom_n & 255) == sphere
	    || ((drvui->atoms[j].atom_n >> 24) & 255) == sphere)
	    find_all_in_box (j);
    }

    if (nvert == 0)
	return;

    cur_cen[0] = o_vert[3 * num];
    cur_cen[1] = o_vert[3 * num + 1];
    cur_cen[2] = o_vert[3 * num + 2];
    if (cur_atom[3] > 0) {
	cur_atom[0] = cur_atom[1];
	cur_atom[1] = cur_atom[2];
	cur_atom[2] = cur_atom[3];
	strcpy (cur_name[0], cur_name[1]);
	strcpy (cur_name[1], cur_name[2]);
	strcpy (cur_name[2], cur_name[3]);
    }
    i = 0;
    if (ellips == 0) {
	for (j = 0; j < 4; j++) {
	    if (drvui->spheres[sphere].sphere_l[j] != ' ')
		cur_name_t[i++] = drvui->spheres[sphere].sphere_l[j];
	}
    } else {
	for (j = 0; j < 4; j++) {
	    if (drvui->ellips[sphere].ellips_l[j] != ' ')
		cur_name_t[i++] = drvui->ellips[sphere].ellips_l[j];
	}
    }
    cur_name_t[i] = '\0';
    sprintf (atnum, "%d", drvui->atoms[drvui->orig_atom_no[num]].sv_atom_n);
    strcat (cur_name_t, atnum);

    if (cur_atom[0] <= 0) {
	cur_atom[0] = i;
	strcpy (cur_name[0], cur_name_t);
    } else if (cur_atom[1] <= 0) {
	cur_atom[1] = i;
	dist12 = dist (cur_atom[0], cur_atom[1]);
	strcpy (cur_name[1], cur_name_t);
    } else if (cur_atom[2] <= 0) {
	cur_atom[2] = i;
	strcpy (cur_name[2], cur_name_t);
	dist12 = dist (cur_atom[0], cur_atom[1]);
	dist23 = dist (cur_atom[1], cur_atom[2]);
	if (dist12 == 0.0f || dist23 == 0.0f || cur_atom[2] == cur_atom[0]) {
	    ang123 = 0.0f;
	} else {
	    dot =
		dot0_3d (s_vert[3 * cur_atom[1]], s_vert[3 * cur_atom[1] + 1],
			 s_vert[3 * cur_atom[1] + 2], s_vert[3 * cur_atom[0]],
			 s_vert[3 * cur_atom[0] + 1], s_vert[3 * cur_atom[0] + 2],
			 s_vert[3 * cur_atom[2]], s_vert[3 * cur_atom[2] + 1],
			 s_vert[3 * cur_atom[2] + 2]);
	    float temp = dot / (dist12 * dist23);

	    if (temp > 1.0f)
		temp = 1.0f;
	    if (temp < -1.0f)
		temp = -1.0f;
	    ang123 = (float) (acos (temp) * RAD);
	}
    } else {
	float v1[3], v2[3], v0[3], p1[3], p2[3];

	int j;

	cur_atom[3] = i;
	strcpy (cur_name[3], cur_name_t);
	dist23 = dist (cur_atom[1], cur_atom[2]);
	dist34 = dist (cur_atom[2], i);
	for (j = 0; j < 3; j++) {
	    v0[j] = s_vert[3 * cur_atom[0] + j] - s_vert[3 * cur_atom[1] + j];
	    v1[j] = s_vert[3 * cur_atom[1] + j] - s_vert[3 * cur_atom[2] + j];
	    v2[j] = s_vert[3 * cur_atom[2] + j] - s_vert[3 * cur_atom[3] + j];
	}
	if (dist23 == 0.0f || dist34 == 0.0f)
	    ang234 = 0.0f;
	else {
	    dot = vdot (v1, v2);
	    float temp = -dot / (dist23 * dist34);

	    if (temp > 1.0f)
		temp = 1.0f;
	    if (temp < -1.0f)
		temp = -1.0f;
	    ang234 = (float) (acos (temp) * RAD);
	}
	vcross (v0, v1, p1);
	vcross (v2, v1, p2);
	if (vlength (p1) < 0.1f || vlength (p2) < 0.1f) {
	    torsion_ang = 0.0f;
	    return;
	}
	torsion_ang = -vdot (p1, p2) / (vlength (p1) * vlength (p2));
	if (torsion_ang < -1.0f)
	    torsion_ang = -1.0f;
	if (torsion_ang > 1.0f)
	    torsion_ang = 1.0f;
	torsion_ang = (float) (acos (torsion_ang) * RAD);
	if (vdot (p1, v2) > 0.0f)	// set the sign
	    torsion_ang *= -1.0f;
    }

    Fl::redraw ();
}
#endif

int
pick_label (int x, int y, int w, int h)
{
// picks the label at the current cursor position

    int i;

    GLint viewport[4];

    GLuint names, *ptr, minZ;

    GLuint *ptrNames = 0, numberOfNames = 0;

    int hits;

    int label;

    float m[16];

    int offset;

    float ratio = 1.0f * (float) w / (float) h;

    memset (selectBuf, 0, BUFSIZE);
    glGetIntegerv (GL_VIEWPORT, viewport);
    glPushMatrix ();

    glSelectBuffer (BUFSIZE, selectBuf);
    glRenderMode (GL_SELECT);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();

    offset = (int) (100.0f / gl_size);
    gluPickMatrix (x + offset, viewport[3] - y + offset, 25, 25, viewport);

    if (M_cameras == 0) {
	if (w <= h)
	    glOrtho (-gl_size, gl_size, -gl_size / ratio, gl_size / ratio, -10000.,
		     10000.);
	else
	    glOrtho (-gl_size * ratio, gl_size * ratio, -gl_size, gl_size, -10000.,
		     10000.);
    } else {
	gluPerspective (17, ratio, 0.01, 1000.);	// view angle, aspect,near/far clip
    }
    glMatrixMode (GL_MODELVIEW);

    glInitNames ();
    glPushName (~0);

    glTranslatef (drvui->Trans[0], drvui->Trans[1], drvui->Trans[2]);

    quaternion_to_rotmatrix (&Rotq, m);

    glMultMatrixf (m);

    for (drvui->frame_no = 1; drvui->frame_no <= drvui->max_frame; drvui->frame_no++)
	generate_gl_texts ();
    drvui->frame_no = drvui->max_frame;

    glPopMatrix ();

    hits = glRenderMode (GL_RENDER);
    glFlush ();

    if (hits == 0)
	return (0);

    ptr = (GLuint *) selectBuf;
    minZ = 0xffffffff;
    for (i = 0; i < hits; i++) {
	names = *ptr;
	ptr++;
	if (*ptr < minZ && *(ptr + 2) != (GLuint) - 1) {
	    numberOfNames = names;
	    minZ = *ptr;
	    ptrNames = ptr + 2;
	}
	ptr += names + 2;
    }
    if (numberOfNames == 0)
	return (0);

    ptr = ptrNames;

    label = *ptr;

    //fprintf(stderr,"picked label no %d\n",label);

    if (label >= drvui->nlabel)
	return (0);		// not a label ?


    //cur_cen[0]= drvui->label_x[label][0];
    //cur_cen[1]= drvui->label_x[label][1];
    //cur_cen[2]= drvui->label_x[label][2];

    //Fl::redraw();

    return (label);
}
