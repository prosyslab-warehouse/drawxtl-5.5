// $Id: Edit1.cxx 1125 2011-03-09 00:23:18Z larry $
//
// Edit1.cxx - routine for DRAWxtl V5.5 - the GUI version 
//
// Coded using the FLTK 1.1.6 widget set
//
//     Larry W. Finger, Martin Kroeker and Brian Toby
//
// This module includes the majority of the edit screens for the GUI
//
// routines contained within this file:
//
//  Add_Frame_Main - add frame widget to main screen
//  Automation_Abort_cb - callback routine to abort automation
//  Automation_Edit_cb - Main window for automation of modulated structures
//  Automation_Go_cb - callback for starting automation
//  Bond_Combo_cb - updates bond distance table when atom in combo box is changed
//  Bond_Frame_Combo_cb - callback routine entered when frame combo box changed
//  Check_Box_cb - processes main page check boxes (not clipflag)
//  Configure_cb - callback routine to load configuration screen
//  Configure_Close_cb - callback routine to cancel configuration file
//  Configure_Save_cb - callback routine to save configuration file
//  Configure_MSMS_cb - callback routine to create MSMS path entry window
//  Configure_MSMS_loc_cb - browser window for MSMS path
//  Configure_Close_MSMS_cb - callback to close MSMS configuration window without saving
//  Configure_Save_MSMS_cb - callback to close the MSMS configuration window
//  ConfigurePOV_cb - callback to select path of POV executable
//  ConfigurePOVOptions_cb - callback to select options for POV
//  Edit_Arrow_cb - callback routine to create Arrow Edit screen
//  Edit_Arrow_Close_cb - callback routine to dismiss Arrow Edit screen
//  Edit_Arrow_Save_cb - callback to save (or apply) Arrow Edit changes
//  Edit_Bond_cb - callback routine to load the bond editing screen
//  Edit_Bond_Close_cb - callback routine entered when the 'Cancel' Button is pressed
//  Edit_Bond_Save_cb - callback routine entered when the 'Save' Button is pressed
//  Edit_Ellipsoid_cb - callback routine to load ellipsoid edit screen
//  Edit_Ellipsoid_Close_cb - callback routine from cancel button on ellipsoids screen
//  Edit_Ellipsoid_Save_cb - callback routine from save button on ellipsoid screen
//  Edit_LonePair_cb - callback routine to load lone-pair edit screen
//  Edit_Lone_Pair_Close_cb - callback routine called when the 'Cancel' Button is pressed
//  Edit_Lone_Pair_Save_cb - callback routine called when the 'Save' Button is pressed
//  Edit_Parmeters_Close_cb - callback from cancel button on Edit Parameter Screen
//  Browse_Map_File_cb - callback to select a fourier map file
//  Edit_Maps_cb - callback routine to load map parameter edit screen
//  Edit_Maps_Close_cb - callback from 'Cancel' button
//  Edit_Maps_Save_cb - callback from 'Save' or 'Apply' button
//  Edit_Parmeters_cb - callback routine for edit parameters screen
//  Edit_Parmeters_Save_cb - callback from save or apply button on Edit Parameter Screen
//  Edit_Polyhedra_cb - callback routine to load polyhedral edit screen
//  Load_Color_Combo - load color combo list widget
//  Edit_Modparms_Close_cb - callback from 'Cancel' button on modulation screen
//  Edit_Modparms_Change_cb - callback for immediate feedback from the modulation t spinners
//  Edit_Modparms_Save_cb - callback from 'Save' or 'Apply' button on modulation screen
//  Edit_Modparms_cb - callback routine to load modulation parameters screen
//  Edit_Slice_cb - callback routine to load map slice screen
//  Edit_Slice_Close_cb - callback routine to close map slice screen
//  Edit_Slice_Save_cb - callback routine to save map slice screen parameters
//  Modify_Occ_cb - callback entered when a line in the occupancy threshold list is selected
//  Occ_Combo_cb - callback associated with the atom list combo on the modulation screen
//  New_Occ_Input_cb - callback triggered by changes to the occupancy entry fields
//  New_Occ_Add_cb - callback entered when Add button for occupancy thresholds is pressed 
//  Slice_Frame_Combo_cb - callback routine entered when the frame number combo on Edit Slice screen is changed

#include "drawxtl.h"
#include "DRAWxtlViewUI.h"
#include "EditView.h"
#include "Ellipsoids.h"
#include "draw_ext.h"
#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#define _chdir chdir
#endif
#include <ctype.h>

#include "DRAWxtl_proto.h"

static int zero = 0;

char Edit_title[128];

char Colors_Combo[] =
    { "\nRed\nGreen\nBlue\nYellow\nCyan\n"
"Magenta\nClear\nWhite\nBlack\nGray05\n"
"Gray10\nGray15\nGray20\nGray25\nGray30\n"
"Gray35\nGray40\nGray45\nGray50\nGray55\n"
"Gray60\nGray65\nGray70\nGray75\nGray80\n"
"Gray85\nGray90\nGray95\nDimGray\nDimGrey\n"
"Gray\nGrey\nLightGray\nLightGrey\nVLightGray\n"
"VLightGrey\nAquamarine\nBlueViolet\nBrown\nCadetBlue\n"
"Coral\nCornflowerBlue\nDarkGreen\nDarkOliveGreen\nDarkOrchid\n"
"DarkSlateBlue\nDarkSlateGray\nDarkSlateGrey\nDarkTurquoise\nFirebrick\n"
"ForestGreen\nGold\nGoldenrod\nGreenYellow\nIndianRed\n"
"Khaki\nLightBlue\nLightSteelBlue\nLimeGreen\nMaroon\n"
"MediumAquamarine\nMediumBlue\nMediumForestGreen\nMediumGoldenrod\nMediumOrchid\n"
"MediumSeaGreen\nMediumSlateBlue\nMediumSpringGreen\nMediumTurquoise\nMediumVioletRed\n"
"MidnightBlue\nNavy\nNavyBlue\nOrange\nOrangeRed\n"
"Orchid\nPaleGreen\nPink\nPlum\nSalmon\n"
"SeaGreen\nSienna\nSkyBlue\nSlateBlue\nSpringGreen\n"
"SteelBlue\nTan\nThistle\nTurquoise\nViolet\n"
"VioletRed\nWheat\nYellowGreen\nSummerSky\nRichBlue\n"
"Brass\nCopper\nBronze\nBronze2\nSilver\n"
"BrightGold\nOldGold\nFeldspar\nQuartz\nMica\n"
"NeonPink\nDarkPurple\nNeonBlue\nCoolCopper\nMandarinOrange\n"
"LightWood\nMediumWood\nDarkWood\nSpicyPink\nSemiSweetChoc\n"
"BakersChoc\nFlesh\nNewTan\nNewMidnightBlue\nVeryDarkBrown\n"
"DarkBrown\nDarkTan\nGreenCopper\nDkGreenCopper\nDustyRose\n"
"HuntersGreen\nScarlet\nMed_Purple\nLight_Purple\nVery_Light_Purple\n"
};

#ifdef WIN32
const char *flu_file_chooser (const char *message, const char *pattern,
			      const char *filename);
#endif

void
Add_Frame_Main (void)
{
// make frame number combo box on main screen hidden/shown as desired
    int i;

    char string[20];

    if (drvui->max_frame > 1) {
	drvui->Frame_No->list.clear ();
	drvui->Frame_No->show ();
	for (i = 1; i <= drvui->max_frame; i++) {
	    sprintf (string, "%d", i);
	    drvui->Frame_No->list.add (string);
	}
	drvui->Frame_No->pop_height (20 * drvui->max_frame);
	drvui->Frame_No->value ("1");
    } else {
	drvui->Frame_No->hide ();
    }
}

void
Arrow_Frame_Combo_cb (Fl_Widget *, void *)
{
// update arrow information table when the frame numver in the combo box is changed
    int i;

    int Frame_No = 1;

    if (!drvui->table)
	drvui->table = (char *) zalloc (20480 * sizeof (char));
    if (drvui->max_frame > 1)
	Frame_No = atoi (arrows->Frame_No->value ());
    strcpy (drvui->table, "");
    for (i = 0; i < drvui->nmag; i++) {	// fill in the widgets
	if (drvui->arrows[i].arrow_fn == Frame_No) {
	    char string[128];

	    sprintf (string, "%6.3f %6.3f %6.3f %10.3f %6.3f %6.3f %10.3f %6.3f   %s\n",
		     drvui->arrows[i].mag_xp[0], drvui->arrows[i].mag_xp[1],
		     drvui->arrows[i].mag_xp[2], drvui->arrows[i].mag_xc[0],
		     drvui->arrows[i].mag_xc[1], drvui->arrows[i].mag_xc[2],
		     drvui->arrows[i].arrow_length, drvui->arrows[i].arrow_diam,
		     drvui->arrows[i].col_arrow);
	    strcat (drvui->table, string);
	}
    }
    arrows->ArrowBuffer->text (drvui->table);
}

void Automation_Abort_cb(void)
{
// callback routine to abort a running automation

    int keep;

    drvui->automation = 0;
    keep = Automate->keeptemps->value();
    if (!keep) {
	float start, end, step;
	int nstep, no_steps;
	char POV_filename[200];
	char string[200],string2[20];

	strcpy(POV_filename, Automate->POV_Filename->value ());
	start = (float) atof (Automate->t_start->value());
	end = (float) atof (Automate->t_end->value());
	step = (float) atof (Automate->t_step->value());
	no_steps = (int)(((end - start) / step) + 0.5f) + 1;
	for (nstep = 1; nstep < no_steps; nstep++) {
	    strcpy(string, POV_filename);
	    sprintf(string2, "%.3d.pov", nstep);
	    strcat(string, string2);
            if (!unlink(string)) break;
	    strcpy(string, POV_filename);
	    if (Automate->Mencoder->value() )
		sprintf(string2,"%dd.tga",nstep);
	    else
		sprintf(string2,"%dd.png",nstep);
	    strcat(string, string2);
	    unlink(string);
	}
    }

//    Automate->Automation_Edit_Window->hide();
//    delete(Automate);
//    Automate = NULL;
}

void Automation_Edit_cb(void)
{
// callback routine to open automation window
    int y;
    char string[200];

    if (!drvui->modulated) {
	Error_Box("Automation requires a modulated structure");
	return;
    }
    if (Automate)
	delete(Automate);
    Automate = new AutomationParam;	// get new instance of automation parms
    Automate->Automation_Edit_Window =
	new Fl_Window (50, 50, 410, 350, "Modulated Structure Automation");
    Automate->Automation_Edit_Window->callback ((Fl_Callback *) Automation_Abort_cb);
#if !defined (WIN32) && !defined (__APPLE__)
    Automate->Automation_Edit_Window->icon ((char *) drvui->icon);
#endif
    y = 30;
    Automate->t_start = new Fl_Input(70, y, 70, 25, "t Start");
    Automate->t_start->align (FL_ALIGN_TOP);
    sprintf(string, "%.3f", drvui->phaseshift[0]);
    Automate->t_start->value(string);
    Automate->t_end = new Fl_Input(170, y, 70, 25, "t End");
    Automate->t_end->align (FL_ALIGN_TOP);
    Automate->t_end->value("1.000");
    Automate->t_step = new Fl_Input(270, y, 70, 25, "t Step");
    Automate->t_step->align (FL_ALIGN_TOP);
    Automate->t_step->value("0.001");
    y += 50;
    Automate->POV_Filename = new Fl_Input(30, y, 350, 25, "POV Basic File Name");
    Automate->POV_Filename->align (FL_ALIGN_TOP);
    strcpy(string, drvui->Cur_Dir);
#ifdef WIN32
    strcat(string, "\\");
#else
    strcat(string,"/");
#endif
    strcat(string, drvui->Cur_Root);
    Automate->POV_Filename->value (string);
    Automate->POV_Filename->tooltip("This name will have \"xxx.pov\" appended.");
    y += 30;

	Fl_Group *qq = new Fl_Group (70, y, 260, 90);

	qq->labelfont (1);
	qq->box (FL_THIN_UP_BOX);
	Automate->NoMovie = new Fl_Radio_Button (80, y + 10, 12, 12, "Prepare image files only");
	Automate->NoMovie->type (102);
	Automate->NoMovie->selection_color ((Fl_Color) 1);
	Automate->NoMovie->align (FL_ALIGN_RIGHT);
	Automate->NoMovie->set ();
	Automate->Mencoder = new Fl_Radio_Button (80, y + 30, 12, 12, "Create AVI using mencocder");
	Automate->Mencoder->type (102);
	Automate->Mencoder->selection_color ((Fl_Color) 1);
	Automate->Mencoder->align (FL_ALIGN_RIGHT);
	Automate->Ffmpeg = new Fl_Radio_Button (80, y + 50, 12, 12, "Create MPG using ffmpeg");
	Automate->Ffmpeg->type (102);
	Automate->Ffmpeg->selection_color ((Fl_Color) 1);
	Automate->Ffmpeg->align (FL_ALIGN_RIGHT);
	Automate->FfmpegG = new Fl_Radio_Button (80, y + 70, 12, 12, "Create animated GIF using ffmpeg");
	Automate->FfmpegG->type (102);
	Automate->FfmpegG->selection_color ((Fl_Color) 1);
	Automate->FfmpegG->align (FL_ALIGN_RIGHT);
	qq->end ();
    y += 100;
    Automate->keeptemps = new Fl_Check_Button(70, y, 70, 25, "Keep intermediate POV and image files");
    y += 50;
    Automate->width = new Fl_Input(70, y, 70, 25, "Image width");
    Automate->width->align (FL_ALIGN_TOP);
    Automate->height = new Fl_Input(170, y, 70, 25, "Image height");
    Automate->height->align (FL_ALIGN_TOP);
    Automate->fps = new Fl_Input(270, y, 70, 25, "Frames/sec");
    Automate->fps->align (FL_ALIGN_TOP);
    y += 50;
    Automate->Go = new Fl_Button(270, y, 70, 25, "Go");
    Automate->Go->callback((Fl_Callback *) Automation_Go_cb);
    Automate->Abort = new Fl_Button(170, y, 70, 25, "Abort");
    Automate->Abort->callback((Fl_Callback *) Automation_Abort_cb);
    Automate->Close = new Fl_Button(70, y, 70, 25, "Close");
    Automate->Close->callback((Fl_Callback *) Automation_Close_cb);
    Automate->Automation_Edit_Window->show();
}

void Automation_Go_cb(void)
{
// callback routine to start automation loop
    char POV_filename[200];
    char incpath[256];
    char cmdbase[512],cmd[512];
    char string[200];
    char string2[200];
    char format[20];
    int nstep, no_steps;
    int width,height,fps;
    int keep;
    int result;
    int dopov;
    int doasy;
    int dowrl;
    int num;
    FILE *filetmp;
    float start, end, step;

    doasy = doAsy;
    dopov = doPOV;
    dowrl = doVrml;
    doAsy = 0;
    doPOV = 1;
    doVrml = 0;
    strcpy(POV_filename, Automate->POV_Filename->value ());
    strcpy(string, POV_filename);
    strcat(string, ".pov");
    filetmp = fopen(string, "w");
    if (!filetmp) {
	Error_Box("Unable to write POV files in specified directory");
	return;
    }
    fclose(filetmp);
    unlink(string);	/* remove the temporary file */
    Progress_Window(-1, "Automation Progress", 100.0f);
    drvui->automation = 1;
    start = (float) atof (Automate->t_start->value());
    end = (float) atof (Automate->t_end->value());
    step = (float) atof (Automate->t_step->value());
    width = atoi (Automate->width->value());
    height = atoi (Automate->height->value());
    fps = atoi (Automate->fps->value());
    if (fps == 0) fps = 30;
    keep = Automate->keeptemps->value();
    no_steps = (int)(((end - start) / step) + 0.5f) + 1;
    num = (int)log10((double)no_steps)+1;
    incpath[0] = '\0';
    if (strlen (drvui->POV_Include) > 10)
        strcpy (incpath, drvui->POV_Include);   // make copy of colors.inc full path name
#ifdef WIN32
    strcpy (cmdbase, "\"\"");       // build the command string
    strcat (cmdbase, drvui->POV_Path);      //   for Windows
    strcat (cmdbase, "\"");
    strcat (cmdbase, " ");
    strcat (cmdbase, drvui->POV_Options);
    strcat (cmdbase, " ");
    strcat (cmdbase, " +I\"");
#else
    strcpy (cmdbase, drvui->POV_Path);      // build command string
    strcat (cmdbase, " ");          //    for Linux
    strcat (cmdbase, drvui->POV_Options); 
    sprintf (string2," -V -GA -GD -GF -GR -GS -GW -D -P +W%d +H%d",width,height);
    strcat (cmdbase,string2);
    if (Automate->Mencoder->value() )
	strcat(cmdbase," +FT");
    else
	strcat(cmdbase," +FN");
    strcat (cmdbase, " ");
    if (strlen (drvui->POV_Include) > 10) {
        strcat (cmdbase, " +HI");
        strcat (cmdbase, drvui->POV_Include);
        strcat (cmdbase, " ");
    }
#endif

    for (nstep = 1; nstep < no_steps; nstep++) {
	float progress;
	drvui->phaseshift[0] = start + (float)nstep * step;
	strcpy(string, POV_filename);
        sprintf(format,"%%.%dd.pov",num);
	sprintf(string2, (const char*)format, nstep);
	strcat(string, string2);
	if (Automate->Mencoder->value() )
	    sprintf(format," +O%%s%%.%dd.tga",num);
	else
	    sprintf(format," +O%%s%%.%dd.png",num);
	sprintf(string2,format,POV_filename,nstep);
	drvui->automate_name = string;
	progress = 100.0f * (float)nstep / (float)no_steps;
	Progress_Window(0, NULL, progress);
	Fl::redraw();
	Update_Str(0);
	Generate_Drawing(0);
	Fl::redraw();
	strcpy(cmd,cmdbase);
	strcat(cmd,string);
	strcat(cmd,string2);
	result=system(cmd);
        if (result != 0) drvui->automation=0;
	if (!keep) unlink(string);
	if (!drvui->automation)
	    break;
    }
    Progress_Window(-2, NULL, 0.0f);
    if (drvui->automation && !Automate->NoMovie->value()) { // if movie selected and run not already aborted
#ifdef WIN32
    strcpy (cmdbase, "\"\"");       // build the command string
    if (Automate->Mencoder->value() ) {
	sprintf(cmd,"%s mf://%s*.tga -mf w=%d:h=%d:fps=%d:type=tga -ovc lavc -lavcopts vcodec=mpeg4 -o %s.mpg", drvui->Mencoder_Path, POV_filename,width,height,fps,POV_filename);
//	sprintf(cmd,"mencoder mf://%s*.tga -mf w=%d:h=%d:fps=%d:type=tga -force-avi-aspect 1.0 -ovc copy -oac copy -o %s.mpg",POV_filename,width,height,fps,POV_filename);
    } else if (Automate->Ffmpeg->value() ) {
	sprintf(cmd,"%s -r %d -f image2 -i %s%%%dd.png -s %dx%d %s.mpg", drvui->FFmpeg_Path, fps, POV_filename, num, width,height,POV_filename);
    } else {
	sprintf(cmd,"%s -f image2 -pix_fmt rgb24 -r 1 -i %s%%%dd.png -s %dx%d -f gif %s.gif", drvui->FFmpeg_Path, POV_filename, num, width,height,POV_filename);
    }
    strcat (cmdbase,cmd);
    strcat (cmdbase, "\"\"");      
    result = system(cmdbase);
#else
    if (Automate->Mencoder->value() ) {
	sprintf(cmd,"%s mf://%s*.tga -mf w=%d:h=%d:fps=%d:type=tga -ovc lavc -lavcopts vcodec=mpeg4 -o %s.mpg", drvui->Mencoder_Path, POV_filename,width,height,fps,POV_filename);
//	sprintf(cmd,"mencoder mf://%s*.tga -mf w=%d:h=%d:fps=%d:type=tga -force-avi-aspect 1.0 -ovc copy -oac copy -o %s.mpg",POV_filename,width,height,fps,POV_filename);
    } else if (Automate->Ffmpeg->value() ) {
	sprintf(cmd,"%s -r %d -f image2 -i %s%%%dd.png -s %dx%d %s.mpg", drvui->FFmpeg_Path, fps, POV_filename, num, width,height,POV_filename);
    } else {
	sprintf(cmd,"%s -f image2 -pix_fmt rgb24 -r 1 -i %s%%%dd.png -s %dx%d -f gif %s.gif", drvui->FFmpeg_Path, POV_filename, num, width,height,POV_filename);
    }
    result = system(cmd);
#endif
    if (result != 0) Error_Box("Video encoding failed");
    }
    if (!keep) {
	if (Automate->Mencoder->value() )
	    sprintf(format,"%%.%dd.tga",num);
	else
	    sprintf(format,"%%.%dd.png",num);
	for (nstep = 1; nstep < no_steps; nstep++) {
	    strcpy(string, POV_filename);
	    sprintf(string2, format, nstep);
	    strcat(string, string2);
            unlink(string);
	}
    }
    drvui->automation = 0;
    doAsy = doasy;
    doPOV = dopov;
    doVrml = dowrl;
    Automate->Automation_Edit_Window->hide();
}

void Automation_Close_cb(void)
{
    Automate->Automation_Edit_Window->hide();
}

void
Bond_Combo_cb (Fl_Widget *, void *)
{
// update bond distance table when the atom in the combo box is changed
    const char *atom;

    if (!drvui->table)
	drvui->table = (char *) zalloc (20480 * sizeof (char));
    atom = Bonds->Bond_Combo->value ();	// get name from combo box
    Bonds->New_Bond_From->value (atom);	// put it in the 'From' location
    Load_Bond_Data (atom, drvui->table);
    Bonds->Bond_Output_Buffer->text (drvui->table);
}

void
Bond_Frame_Combo_cb (Fl_Widget *, void *)
{
// routine called when the frame number is changed on the Bonds Edit Screen
    int i, j;

    char type[5];

    char atom1[5];

    char atom2[5];

    char widget[16382];

    char atoms[100][5];

    char color[40];

    float rad, d1, d2;

    char string[100];

    int Frame_No = 1;

    if (drvui->max_frame > 1)
	Frame_No = atoi (Bonds->Frame_No->value ());
    int nlist = Get_Unique_Atoms (atoms, Frame_No);	// get unique atom names in sorted list

    Bonds->Bond_Combo->list.clear ();	// clear out the old names
    for (j = 0; j < nlist; j++) {	// add atom names in alphabetic order
	Bonds->Bond_Combo->list.add (atoms[j]);
    }
    Bonds->Bond_Combo->value (atoms[0]);	// load first one in window
    Bond_Combo_cb (NULL, NULL);
    widget[0] = '\0';
    for (i = 1; i < drvui->nbond; i++) {
	if (drvui->bonds[i].bond_fn != Frame_No)
	    continue;		//skip if not for this frame
	strncpy (atom1, drvui->bonds[i].bond_l1, 4);
	strncpy (atom2, drvui->bonds[i].bond_l2, 4);
	atom1[4] = 0;
	atom2[4] = 0;
	for (j = 3; j >= 0; --j) {
	    if (atom1[j] == ' ')
		atom1[j] = 0;
	}
	for (j = 3; j >= 0; --j) {
	    if (atom2[j] == ' ')
		atom2[j] = 0;
	}
	rad = drvui->bonds[i].bond_size;
	d1 = drvui->bonds[i].bond_min;
	d2 = drvui->bonds[i].bond_max;
	strncpy (color, drvui->bonds[i].col_bond, 25);
	if (color[strlen (color) - 1] < ' ')
	    color[strlen (color) - 1] = 0;
	if (strlen (color) > 25)
	    color[25] = 0;
	if (drvui->bonds[i].bond_style == 0) {
	    strcpy (type, "bond");
	    sprintf (string, "%4s     %4s   %4s  %9.3f %9.3f %9.3f     %s\n", type, atom1, atom2,
		 rad, d1, d2, color);
	} else {
	    strcpy (type, "dash");
	    sprintf (string, "%4s %3d %4s   %4s  %9.3f %9.3f %9.3f     %s\n", type, 
		drvui->bonds[i].bond_style, atom1, atom2, rad, d1, d2, color);
	}
	strcat (widget, string);
	Bonds->BondInstr1->show ();
    }
    Bonds->BondBuffer->text (widget);
}

void
Browse_Map_File_cb (void)
{
// callback routine to select a map file
#if defined(WIN32)
    const char *newfile =
	flu_file_chooser ("Select Map File", "*.*", Maps->Filename->value ());
#else
    char *newfile =
	fl_file_chooser ("Select Map File", "*.*", Maps->Filename->value (), 1);
#endif
    if (newfile) {
	Maps->Filename->value (newfile);
	ReadFourMap = 0;
	Map_Info.info_valid = 0;
    }
}

void
Check_Box_cb (void)
{
// callback routine entered when a main-screen checkbox is changed
    if (drvui->Generate_VRML1->value () != 0) {
	if (Vrml2)
	    drvui->Str_File_Changed = 1;
	Vrml2 = 0;
    } else {
	if (!Vrml2)
	    drvui->Str_File_Changed = 1;
	Vrml2 = 1;
    }
    if (drvui->Orthographic_View->value () != 0) {
	if (M_cameras)
	    drvui->Str_File_Changed = 1;
	M_cameras = 0;
	gl_size = max (POV_Max[1] - POV_Min[1], POV_Max[0] - POV_Min[0]);
    } else {
	if (!M_cameras)
	    drvui->Str_File_Changed = 1;
	M_cameras = 1;
    }
    if (drvui->Show_Vector_Triple->value () != 0) {
	if (!Display_axes)
	    drvui->Str_File_Changed = 1;
	Display_axes = 1;
	label_cell ();
    } else {
	if (Display_axes)
	    drvui->Str_File_Changed = 1;
	Display_axes = 0;
	label_cell ();
    }
    Update_Str (0);		// update the 'str' file
    Generate_Drawing (0);	// regenerate the drawing
    Fl::redraw ();		// update the screen
}

void
Configure_cb (void)
{
// callback routine to display configuration screen
    int y;

    char string[30];

    if (Configure) {
	Configure->ConfigWindow->show ();
	return;
    }
    Configure = new ConfigParm;	// get new instance of config parms
    Configure->ConfigWindow =
	new Fl_Window (50, 50, 640, 600, "POV Configuration Window");
    Configure->ConfigWindow->callback ((Fl_Callback *) Configure_Close_cb);
    y = 30;
    Configure->POVOptions = new Fl_Input (170, y, 300, 40, "POV Options");
    Configure->POVOptions->textcolor (1);
    Configure->POVOptions->value (drvui->POV_Options);
    Configure->POVOptions->
	tooltip
	("+Wnum and +Hnum for size, +D for preview, +P for pause at end, +A for antialiased lines, +F for filetype (FC=TGA,FN=PNG), +L for POV library path");
    y += 60;
    Configure->POVPath = new Fl_Input (170, y, 300, 40, "POV program");
    Configure->POVPath->textcolor (1);
    Configure->POVPath->value (drvui->POV_Path);
    Configure->POVPath->
	tooltip ("Enter the name and full path of the pov-ray executable here");
    Fl_Button *a = new Fl_Button (475, y + 5, 70, 25, "Browse");

    a->callback ((Fl_Callback *) ConfigurePOV_cb);
    y += 60;
    Configure->POVIncludePath = new Fl_Input (170, y, 300, 40, "POV Color File");
    Configure->POVIncludePath->textcolor (1);
    Configure->POVIncludePath->value (drvui->POV_Include);
    Fl_Button *b = new Fl_Button (475, y + 5, 70, 25, "Browse");

    b->callback ((Fl_Callback *) ConfigurePOVOptions_cb);
    y += 60;
    Configure->POVDefaultFinish = new Fl_Input (170, y, 300, 40, "POV Object Finish");
    Configure->POVDefaultFinish->textcolor (1);
    Configure->POVDefaultFinish->value (drvui->DefaultFinish);
    y += 60;
    Configure->Stereo_Base = new Fl_Input (170, y, 300, 40, "Stereo Base Separation");
    sprintf (string, "%.3f", drvui->stereo_base);
    Configure->Stereo_Base->value (string);
    Configure->Stereo_Base->textcolor (1);
    y += 60;
    Configure->Stereo =
	new Fl_Check_Button (150, y, 325, 25,
			     "Generate stereo POV pair - StereoPOV needed");
    if (drvui->Stereo == 1)
	Configure->Stereo->set ();
    y += 30;
    Configure->StereoMesh =
	new Fl_Check_Button (150, y, 325, 25,
			     "Generate stereo POV pair - POV 3.7 or newer needed");
    if (drvui->Stereo == 2)
	Configure->StereoMesh->set ();
    y += 30;
    Configure->CrossEyed =
	new Fl_Check_Button (150, y, 175, 25, "Use \"cross-eyed\" view");
    if (drvui->cross_eyed)
	Configure->CrossEyed->set ();

    y += 40;

    Configure->MencoderPath = new Fl_Input (170, y, 300, 40, "Mencoder program");
    Configure->MencoderPath->textcolor (1);
    Configure->MencoderPath->value (drvui->Mencoder_Path);
    Configure->MencoderPath->
	tooltip ("Enter the name and full path of the mencoder program for video encoding");
    Fl_Button *c = new Fl_Button (475, y + 5, 70, 25, "Browse");

    c->callback ((Fl_Callback *) ConfigureMencoder_cb);

    y += 60;
    Configure->FFmpegPath = new Fl_Input (170, y, 300, 40, "FFmpeg program");
    Configure->FFmpegPath->textcolor (1);
    Configure->FFmpegPath->value (drvui->FFmpeg_Path);
    Configure->FFmpegPath->
	tooltip ("Enter the name and full path of the ffmpeg program for video encoding");
    Fl_Button *d = new Fl_Button (475, y + 5, 70, 25, "Browse");

    d->callback ((Fl_Callback *) ConfigureFFmpeg_cb);

    y += 60;
    Fl_Button *o = new Fl_Button (330, y, 90, 30, "Save");

    o->callback ((Fl_Callback *) Configure_Save_cb);
    o->tooltip ("Apply current contents of top box to drawing, then close this window.");
    Fl_Button *p = new Fl_Button (220, y, 90, 30, "Cancel");

    p->tooltip ("Close this window and discard all changes.");
    p->callback ((Fl_Callback *) Configure_Close_cb);
#if !defined (WIN32) && !defined (__APPLE__)
    Configure->ConfigWindow->icon ((char *) drvui->icon);
#endif
    Configure->ConfigWindow->end ();
    Configure->ConfigWindow->show ();
}

void
Configure_Close_cb (void)
{
    Configure->ConfigWindow->hide ();
}

void
Configure_Save_cb (void)
{
// callback routine to save configuration file

    char string[100];

    strcpy (drvui->POV_Options, Configure->POVOptions->value ());
    strcpy (drvui->POV_Path, Configure->POVPath->value ());
    strcpy (drvui->POV_Include, Configure->POVIncludePath->value ());
    strcpy (drvui->DefaultFinish, Configure->POVDefaultFinish->value ());
    if (Configure->Stereo->value ())
	drvui->Stereo = 1;
    else if (Configure->StereoMesh->value ())
	drvui->Stereo = 2;
    else
	drvui->Stereo = 0;
    if (Configure->CrossEyed->value ())
	drvui->cross_eyed = 1;
    else
	drvui->cross_eyed = 0;
    strcpy (string, Configure->Stereo_Base->value ());
    sscanf (string, "%f", &drvui->stereo_base);
    strcpy (drvui->Mencoder_Path, Configure->MencoderPath->value ());
    strcpy (drvui->FFmpeg_Path, Configure->FFmpegPath->value ());
    WriteConfig ();
    Configure->ConfigWindow->hide ();
}

void
Configure_Misc_cb (void)
{
// callback routine to display configuration screen for miscellaneous parameters
    int y;

    if (MiscConfigure) {
	MiscConfigure->MiscConfigWindow->show ();
	return;
    }
    MiscConfigure = new ConfigMiscParm;	// get new instance of misc. config parms
    MiscConfigure->MiscConfigWindow =
	new Fl_Window (50, 50, 400, 250, "Misc. Configuration Window");
    MiscConfigure->MiscConfigWindow->callback ((Fl_Callback *) Configure_Misc_Close_cb);
    y = 30;
    MiscConfigure->LoadLast =
	new Fl_Check_Button (40, y, 275, 25, "Load most recently used file on startup");
    if (strncmp (drvui->LoadOnStartup, "yes", 3) == 0)
	MiscConfigure->LoadLast->set ();
    y += 30;
    MiscConfigure->AutoLabel =
	new Fl_Check_Button (40, y, 275, 25, "Label atoms when importing from SHELX");
    if (drvui->autolabel)
	MiscConfigure->AutoLabel->set ();
    y += 30;
    MiscConfigure->doVrml =
	new Fl_Check_Button (40, y, 275, 25, "Generate VRML file while rendering");
    if (doVrml)
	MiscConfigure->doVrml->set ();
    y += 30;
    MiscConfigure->doPOV =
	new Fl_Check_Button (40, y, 275, 25, "Generate POV input file while rendering");
    if (doPOV)
	MiscConfigure->doPOV->set ();
    y += 30;
    MiscConfigure->doAsy =
	new Fl_Check_Button (40, y, 275, 25, "Generate Asymptote input file while rendering");
    if (doAsy)
	MiscConfigure->doAsy->set ();
    y += 60;
    Fl_Button *n = new Fl_Button (100, y, 90, 30, "Cancel");

    n->tooltip ("Close this window and discard all changes.");
    n->callback ((Fl_Callback *) Configure_Misc_Close_cb);
    Fl_Button *o = new Fl_Button (210, y, 90, 30, "Save");

    o->callback ((Fl_Callback *) Configure_Misc_Save_cb);
    o->tooltip ("Apply current contents of top box to drawing, then close this window.");
#if !defined (WIN32) && !defined (__APPLE__)
    MiscConfigure->MiscConfigWindow->icon ((char *) drvui->icon);
#endif
    MiscConfigure->MiscConfigWindow->end ();
    MiscConfigure->MiscConfigWindow->show ();
}

void
Configure_Misc_Close_cb (void)
{
    MiscConfigure->MiscConfigWindow->hide ();
}

void
Configure_Misc_Save_cb (void)
{
    if (MiscConfigure->LoadLast->value ())
	strcpy (drvui->LoadOnStartup, "yes");
    else
	strcpy (drvui->LoadOnStartup, "no");
    if (MiscConfigure->AutoLabel->value ())
	drvui->autolabel = 1;
    else
	drvui->autolabel = 0;
    if (MiscConfigure->doVrml->value ()) {
	doVrml = 1;
	drvui->Generate_VRML1->activate ();
    } else {
	doVrml = 0;
	drvui->Generate_VRML1->deactivate ();
    }
    if (MiscConfigure->doPOV->value ()) {
	doPOV = 1;
	DRAWxtlViewUI::drawxtl_menu[39].activate ();
    } else {
	doPOV = 0;
	DRAWxtlViewUI::drawxtl_menu[39].deactivate ();
    }
    if (MiscConfigure->doAsy->value ()) {
	doAsy = 1;
    } else {
	doAsy = 0;
    }
    WriteConfig ();
    MiscConfigure->MiscConfigWindow->hide ();
    Generate_Drawing (0);
}

void
Configure_MSMS_cb (void)
{
    if (MSMSConfigure) {
	MSMSConfigure->MSMSConfigWindow->show ();
	return;
    }
    MSMSConfigure = new ConfigMSMSParm;	// get new instance
    MSMSConfigure->MSMSConfigWindow =
	new Fl_Window (50, 50, 640, 140, "MSMS Configuration Window");
    MSMSConfigure->MSMSConfigWindow->callback ((Fl_Callback *) Configure_Close_MSMS_cb);
    MSMSConfigure->MSMSPath = new Fl_Input (170, 30, 300, 40, "MSMS program");
    MSMSConfigure->MSMSPath->textcolor (1);
    MSMSConfigure->MSMSPath->value (drvui->MSMS_Path);
    MSMSConfigure->MSMSPath->
	tooltip ("Enter the name and full path of the msms executable here");
    Fl_Button *a = new Fl_Button (475, 35, 70, 25, "Browse");

    a->callback ((Fl_Callback *) Configure_MSMS_loc_cb);
    Fl_Button *o = new Fl_Button (330, 80, 90, 30, "Save");

    o->callback ((Fl_Callback *) Configure_Save_MSMS_cb);
    o->tooltip ("Apply current contents of top box to drawing, then close this window.");
    Fl_Button *p = new Fl_Button (220, 80, 90, 30, "Cancel");

    p->tooltip ("Close this window and discard all changes.");
    p->callback ((Fl_Callback *) Configure_Close_MSMS_cb);
#if !defined (WIN32) && !defined (__APPLE__)
    MSMSConfigure->MSMSConfigWindow->icon ((char *) drvui->icon);
#endif
    MSMSConfigure->MSMSConfigWindow->end ();
    MSMSConfigure->MSMSConfigWindow->show ();
}

void
Configure_MSMS_loc_cb (void)
{
// routine to configure MSMS exec path
#ifdef WIN32
    const char *newfile =
	flu_file_chooser ("Select MSMS exe file", "*.exe", drvui->MSMS_Path);
#else
    char *newfile = fl_file_chooser ("Select MSMS program", "*", drvui->MSMS_Path);
#endif
    if (newfile) {
	strcpy (drvui->MSMS_Path, newfile);
	MSMSConfigure->MSMSPath->value (drvui->MSMS_Path);
    }
}

void
ConfigurePOV_cb (void)
{
// routine to configure POV exec path
#ifdef WIN32
    const char *newfile =
	flu_file_chooser ("Select POVRAY exe file", "*.exe", drvui->POV_Path);
#else
    char *newfile = fl_file_chooser ("Select POVRAY program", "*", drvui->POV_Path);
#endif
    if (newfile) {
	strcpy (drvui->POV_Path, newfile);
	Configure->POVPath->value (drvui->POV_Path);
    }
}

void
ConfigurePOVOptions_cb (void)
{
// routine to configure POV include path
#ifdef WIN32
    const char *newfile =
	flu_file_chooser ("Select Path for POV colors.inc", "*.inc", drvui->POV_Include);
#else
    char *newfile =
	fl_file_chooser ("Select Path for POV 'colors.inc'", "*.inc", drvui->POV_Include);
#endif
    if (newfile) {
	strcpy (drvui->POV_Include, newfile);
	Configure->POVIncludePath->value (drvui->POV_Include);
    }
}

void
ConfigureMencoder_cb (void)
{
// routine to configure mencoder exec path
#ifdef WIN32
    const char *newfile =
	flu_file_chooser ("Select mencoder exe file", "*.exe", drvui->Mencoder_Path);
#else
    char *newfile = fl_file_chooser ("Select mencoder program", "*", drvui->Mencoder_Path);
#endif
    if (newfile) {
	strcpy (drvui->Mencoder_Path, newfile);
	Configure->MencoderPath->value (drvui->Mencoder_Path);
    }
}

void
ConfigureFFmpeg_cb (void)
{
// routine to configure FFmpeg exec path
#ifdef WIN32
    const char *newfile =
	flu_file_chooser ("Select FFmpeg exe file", "*.exe", drvui->FFmpeg_Path);
#else
    char *newfile = fl_file_chooser ("Select FFmpeg program", "*", drvui->FFmpeg_Path);
#endif
    if (newfile) {
	strcpy (drvui->FFmpeg_Path, newfile);
	Configure->FFmpegPath->value (drvui->FFmpeg_Path);
    }
}

void
Configure_Close_MSMS_cb (void)
{
    MSMSConfigure->MSMSConfigWindow->hide ();
}

void
Configure_Save_MSMS_cb (void)
{
// callback routine to save configuration file
    strcpy (drvui->MSMS_Path, MSMSConfigure->MSMSPath->value ());
    WriteConfig ();
    MSMSConfigure->MSMSConfigWindow->hide ();
}

void
Edit_Arrow_cb (void)
{
    char widget[16382];
    int y = 25;
    static int one = 1;
    int yy;

    if (!strlen (drvui->CurFile->value ())) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }
    Save_Working_Copy ();
    if (!arrows) {
	arrows = new ArrowParam;

	{
	    Fl_Window *o = arrows->ArrowWindow =
		new Fl_Window (0, 0, 490, 500, "Edit Arrow Parameters");
	    o->callback ((Fl_Callback *) Edit_Arrow_Close_cb);
#if !defined (WIN32) && !defined (__APPLE__)
	    o->icon ((char *) drvui->icon);
#endif
	    {
		Fl_Text_Display *o = new Fl_Text_Display (250, y, 0, 0,
							  "  Position                   Components            Length   Diam.  Color");

		o->labelfont (1);
	    }
	    {
		Fl_Text_Editor *o = new Fl_Text_Editor (20, y, 450, 80);

		o->labelfont (1);
		arrows->ArrowBuffer = new Fl_Text_Buffer;
		strcpy (widget, "");
		o->buffer (arrows->ArrowBuffer);
	    }
	    arrows->ArrowBuffer->add_modify_callback ((Fl_Text_Modify_Cb) Modify_Arrow_cb,
						      (void *) NULL);
	    arrows->ArrowInstr =
		new Fl_Output (20, y + 90, 450, 0,
			       "Press 'Add' to replace "
			       "selected line - 'Remove' to delete it");
	    arrows->ArrowInstr->hide ();
	    arrows->ArrowInstr->align (FL_ALIGN_BOTTOM);
	    arrows->ArrowInstr1 =
		new Fl_Output (20, y + 90, 450, 0,
			       "Highlight text above " "or double click to edit line");
	    arrows->ArrowInstr1->hide ();
	    arrows->ArrowInstr1->align (FL_ALIGN_BOTTOM);
	    y += 130;
	    {
		Fl_Text_Display *o = new Fl_Text_Display (244, y, 0, 0, "Position");

		o->labelfont (1);
	    }
	    y += 10;
	    {
		Fl_Input *o = arrows->Px = new Fl_Input (110, y, 70, 25, "x");

		o->labelfont (1);
		o->align (FL_ALIGN_TOP);
		o->callback ((Fl_Callback *) New_Arrow_Input_cb);
	    }
	    {
		Fl_Input *o = arrows->Py = new Fl_Input (210, y, 70, 25, "y");

		o->labelfont (1);
		o->align (FL_ALIGN_TOP);
		o->callback ((Fl_Callback *) New_Arrow_Input_cb);
	    }
	    {
		Fl_Input *o = arrows->Pz = new Fl_Input (310, y, 70, 25, "z");

		o->labelfont (1);
		o->align (FL_ALIGN_TOP);
		o->callback ((Fl_Callback *) New_Arrow_Input_cb);
	    }
	    if (drvui->max_frame > 1) {
		int i;

		char string[128];

		Flu_Combo_List *o = arrows->Frame_No =
		    new Flu_Combo_List (400, y, 75, 25, "Frame No.");
		o->align (FL_ALIGN_TOP);
		o->callback (Arrow_Frame_Combo_cb);
		o->labelfont (1);
		for (i = 1; i <= drvui->max_frame; i++) {
		    sprintf (string, "%d", i);
		    o->list.add (string);
		}
		o->pop_height (20 * drvui->max_frame);
		o->value ("1");
	    }
	    y += 50;
	    {
		Fl_Text_Display *o = new Fl_Text_Display (250, y, 0, 0, "Components");

		o->labelfont (1);
	    }
	    y += 10;
	    {
		Fl_Input *o = arrows->Cx = new Fl_Input (110, y, 70, 25, "x");

		o->labelfont (1);
		o->align (FL_ALIGN_TOP);
		o->callback ((Fl_Callback *) New_Arrow_Input_cb);
	    }
	    {
		Fl_Input *o = arrows->Cy = new Fl_Input (210, y, 70, 25, "y");

		o->labelfont (1);
		o->align (FL_ALIGN_TOP);
		o->callback ((Fl_Callback *) New_Arrow_Input_cb);
	    }
	    {
		Fl_Input *o = arrows->Cz = new Fl_Input (310, y, 70, 25, "z");

		o->labelfont (1);
		o->align (FL_ALIGN_TOP);
		o->callback ((Fl_Callback *) New_Arrow_Input_cb);
	    }
	    y += 45;
	    {
		Fl_Text_Display *o =
		    new Fl_Text_Display (75, y, 344, 0, "Arrow Size and Color");
		o->labelfont (1);
	    }
	    y += 15;
	    {
		Fl_Input *o = arrows->Length = new Fl_Input (90, y, 70, 25, "Length");

		o->labelfont (1);
		o->align (FL_ALIGN_TOP);
		o->callback ((Fl_Callback *) New_Arrow_Input_cb);
	    }
	    {
		Fl_Input *o = arrows->Diameter =
		    new Fl_Input (170, y, 70, 25, "Diameter");
		o->labelfont (1);
		o->align (FL_ALIGN_TOP);
		o->callback ((Fl_Callback *) New_Arrow_Input_cb);
	    }
	    {
		Flu_Combo_List *o = arrows->Color =
		    new Flu_Combo_List (250, y, 150, 25, "Color");
		o->labelfont (1);
		Load_Color_Combo (o);
		o->align (FL_ALIGN_TOP);
		o->callback ((Fl_Callback *) New_Arrow_Input_cb);
	    }
	    y += 30;
	    {
		Fl_Button *o = arrows->AddButton = new Fl_Button (165, y, 70, 25, "Add");

		o->down_box (FL_DOWN_BOX);
		o->labelfont (1);
		o->tooltip
		    ("When active, press to transfer data in boxes to window above");
		o->deactivate ();
		o->callback ((Fl_Callback *) New_Arrow_Add_cb, &one);
	    }
	    {
		Fl_Button *o = arrows->RemoveButton =
		    new Fl_Button (250, y, 70, 25, "Remove");
		o->down_box (FL_DOWN_BOX);
		o->labelfont (1);
		o->tooltip ("When active, press to remove highlighted line.");
		o->deactivate ();
		o->callback ((Fl_Callback *) New_Arrow_Add_cb, &zero);
	    }
	    y += 30;
	    arrows->ArrowInstr2 = new Fl_Output (25, y, 440, 0,
						 "Changes are temporary until \"Apply\" or \"Save\" is pressed.");
	    arrows->ArrowInstr2->hide ();
	    arrows->ArrowInstr2->align (FL_ALIGN_BOTTOM);
	    y += 20;
	    {
		Fl_Text_Display *o = new Fl_Text_Display (248, y, 0, 0,
							  "Matrix Between Magnetic and Nuclear Cells");

		o->labelfont (1);
	    }
	    y += 15;
	    {
		Fl_Text_Display *o = new Fl_Text_Display (250, y, 0, 0, "Nuclear ");

		o->labelfont (1);
	    }
	    y += 15;
	    {
		Fl_Text_Display *o = new Fl_Text_Display (95, y + 20, 0, 0, "Ma");

		o->labelfont (1);
	    }
	    {
		Fl_Input *o = arrows->Aa = new Fl_Input (110, y, 70, 20, "Na");

		o->tooltip ("a(nuclear) component of a(magnetic)");
		o->labelfont (1);
		o->align (FL_ALIGN_TOP);
	    }
	    {
		Fl_Input *o = arrows->Ab = new Fl_Input (210, y, 70, 20, "Nb");

		o->tooltip ("b(nuclear) component of a(magnetic)");
		o->labelfont (1);
		o->align (FL_ALIGN_TOP);
	    }
	    {
		Fl_Input *o = arrows->Ac = new Fl_Input (310, y, 70, 20, "Nc");

		o->tooltip ("c(nuclear) component of a(magnetic)");
		o->labelfont (1);
		o->align (FL_ALIGN_TOP);
	    }
	    yy = y;
	    {
		Fl_Text_Display *o = new Fl_Text_Display (75, yy, 0, 0, "M");

		o->labelfont (1);
	    }
	    {
		Fl_Text_Display *o = new Fl_Text_Display (75, yy + 12, 0, 0, "a");

		o->labelfont (1);
	    }
	    {
		Fl_Text_Display *o = new Fl_Text_Display (75, yy + 24, 0, 0, "g");

		o->labelfont (1);
	    }
	    {
		Fl_Text_Display *o = new Fl_Text_Display (75, yy + 36, 0, 0, "n");

		o->labelfont (1);
	    }
	    {
		Fl_Text_Display *o = new Fl_Text_Display (75, yy + 48, 0, 0, "e");

		o->labelfont (1);
	    }
	    {
		Fl_Text_Display *o = new Fl_Text_Display (75, yy + 60, 0, 0, "t");

		o->labelfont (1);
	    }
	    {
		Fl_Text_Display *o = new Fl_Text_Display (75, yy + 72, 0, 0, "i");

		o->labelfont (1);
	    }
	    {
		Fl_Text_Display *o = new Fl_Text_Display (75, yy + 84, 0, 0, "c");

		o->labelfont (1);
	    }
	    y += 20;
	    {
		Fl_Input *o = arrows->Ba = new Fl_Input (110, y, 70, 20, "Mb");

		o->tooltip ("a(nuclear) component of b(magnetic)");
		o->labelfont (1);
	    }
	    {
		Fl_Input *o = arrows->Bb = new Fl_Input (210, y, 70, 20);

		o->tooltip ("b(nuclear) component of b(magnetic)");
	    }
	    {
		Fl_Input *o = arrows->Bc = new Fl_Input (310, y, 70, 20);

		o->tooltip ("c(nuclear) component of b(magnetic)");
	    }
	    y += 20;
	    {
		Fl_Input *o = arrows->Ca = new Fl_Input (110, y, 70, 20, "Mc");

		o->tooltip ("a(nuclear) component of c(magnetic)");
		o->labelfont (1);
	    }
	    {
		Fl_Input *o = arrows->Cb = new Fl_Input (210, y, 70, 20);

		o->tooltip ("b(nuclear) component of c(magnetic)");
	    }
	    {
		Fl_Input *o = arrows->Cc = new Fl_Input (310, y, 70, 20);

		o->tooltip ("c(nuclear) component of c(magnetic)");
	    }
	    y += 30;
	    {
		Fl_Button *o = new Fl_Button (110, y, 70, 25, "Close");

		o->down_box (FL_DOWN_BOX);
		o->labelfont (1);
		o->tooltip ("Close this window and discard all changes.");
		o->callback ((Fl_Callback *) Edit_Arrow_Close_cb);
	    }
	    {
		Fl_Button *o = new Fl_Button (210, y, 70, 25, "Apply");

		o->down_box (FL_DOWN_BOX);
		o->tooltip
		    ("Apply current contents of top box to drawing, but leave this window open.");
		o->labelfont (1);
		o->callback ((Fl_Callback *) Edit_Arrow_Save_cb, &zero);
	    }
	    {
		Fl_Button *o = new Fl_Button (310, y, 70, 25, "Save");

		o->down_box (FL_DOWN_BOX);
		o->labelfont (1);
		o->tooltip
		    ("Apply current contents of top box to drawing, then close this window.");
		o->callback ((Fl_Callback *) Edit_Arrow_Save_cb, &one);
	    }
	}
	arrows->ArrowWindow->end ();
    }
    Arrow_Frame_Combo_cb (NULL, NULL);
    arrows->ArrowInstr1->show ();
    if (drvui->nmag) {		// populate widgets
	char string[128];

	sprintf (string, "%6.3f", drvui->mag_matrix[0][0]);
	arrows->Aa->value (string);
	sprintf (string, "%6.3f", drvui->mag_matrix[0][1]);
	arrows->Ab->value (string);
	sprintf (string, "%6.3f", drvui->mag_matrix[0][2]);
	arrows->Ac->value (string);
	sprintf (string, "%6.3f", drvui->mag_matrix[1][0]);
	arrows->Ba->value (string);
	sprintf (string, "%6.3f", drvui->mag_matrix[1][1]);
	arrows->Bb->value (string);
	sprintf (string, "%6.3f", drvui->mag_matrix[1][2]);
	arrows->Bc->value (string);
	sprintf (string, "%6.3f", drvui->mag_matrix[2][0]);
	arrows->Ca->value (string);
	sprintf (string, "%6.3f", drvui->mag_matrix[2][1]);
	arrows->Cb->value (string);
	sprintf (string, "%6.3f", drvui->mag_matrix[2][2]);
	arrows->Cc->value (string);
    }
    arrows->ArrowWindow->show ();
}

void
Edit_Arrow_Close_cb (class Fl_Button *, void *)
{
    arrows->ArrowWindow->hide ();
    Restore_Working_Copy ();	// undo any changes
    Generate_Drawing (1);	// regenerate
}

void
Edit_Arrow_Save_cb (Fl_Button *, int *save)
{
    int i = 0;

    char string[16382];

    float Px[3], Cx[3], length, diam;

    char color[40];

    int Frame_No = 1;

    char *selection;

    if (drvui->max_frame > 1)
	Frame_No = atoi (arrows->Frame_No->value ());

    memset (color, 0, 40);
    drvui->Str_File_Changed = 1;
    selection = arrows->ArrowBuffer->text ();
    strcpy (string, selection);
    free (selection);
    if (strlen (string) < 10) {
//      drvui->nmag = 0;
	strcpy (string, "");
	arrows->ArrowBuffer->text (string);
    } else {
	while (strlen (string) > 10) {
	    sscanf (string, " %f %f %f  %f %f %f %f %f %39c", &Px[0], &Px[1],
		    &Px[2], &Cx[0], &Cx[1], &Cx[2], &length, &diam, color);
	    int j, k;

	    for (j = 0; j < 3; j++) {
		drvui->arrows[i].mag_xp[j] = Px[j];
		drvui->arrows[i].mag_xc[j] = Cx[j];
	    }
	    drvui->arrows[i].arrow_length = length;
	    drvui->arrows[i].arrow_diam = diam;
	    drvui->arrows[i].arrow_fn = Frame_No;
	    trim_string (color, 40);
	    if (!strlen (color))
		strcpy (color, "Gray20");
	    strcpy (drvui->arrows[i].col_arrow, color);
	    i++;
	    for (j = 0; j < (int) strlen (string); j++) {
		if (string[j] == '\n')
		    break;
	    }
	    for (j++, k = 0; j < (int) strlen (string); j++)
		string[k++] = string[j];
	    string[k] = 0;
	}
    }
    drvui->nmag = i;
    drvui->mag_matrix[0][0] = (float) atof (arrows->Aa->value ());
    drvui->mag_matrix[0][1] = (float) atof (arrows->Ab->value ());
    drvui->mag_matrix[0][2] = (float) atof (arrows->Ac->value ());
    drvui->mag_matrix[1][0] = (float) atof (arrows->Ba->value ());
    drvui->mag_matrix[1][1] = (float) atof (arrows->Bb->value ());
    drvui->mag_matrix[1][2] = (float) atof (arrows->Bc->value ());
    drvui->mag_matrix[2][0] = (float) atof (arrows->Ca->value ());
    drvui->mag_matrix[2][1] = (float) atof (arrows->Cb->value ());
    drvui->mag_matrix[2][2] = (float) atof (arrows->Cc->value ());
    Update_Str (0);
    Generate_Drawing (0);
    if (*save != 3) {
	Save_Working_Copy ();	// commit changes except from Add_New_Arrow
	arrows->ArrowInstr2->hide ();
    }
    if (*save == 1) {
	arrows->ArrowWindow->hide ();
    }
    Fl::redraw ();
}

void
Edit_Bond_cb (void)
{
// routine to create the edit bonds screen
    char string[100];

    static int one = 1;

    int i;

    int newBonds = 0;

    int y = 25;

    if (!strlen (drvui->CurFile->value ())) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }
/*
    if (!natom) {         // No atoms, no bonds...
        Error_Box("This structure file does not contain any atoms.");
        return;
    }
*/
    Save_Working_Copy ();
    if (!Bonds) {
	Bonds = new BondParam;	// new instance of the bond parameters
	newBonds = 1;
	Bonds->Bond_Edit_Window =
	    new Fl_Window (50, 50, 560, 490, "Edit Bond Parameters");
	Bonds->Bond_Edit_Window->callback ((Fl_Callback *) Edit_Bond_Close_cb);
	if (drvui->max_frame > 1) {
	    Flu_Combo_List *o = Bonds->Frame_No =
		new Flu_Combo_List (220, y, 75, 25, "Frame No.");
	    o->align (FL_ALIGN_TOP);
	    o->callback (Bond_Frame_Combo_cb);
	    o->labelfont (1);
	    for (i = 1; i <= drvui->max_frame; i++) {
		sprintf (string, "%d", i);
		o->list.add (string);
	    }
	    o->pop_height (20 * drvui->max_frame);
	    o->value ("1");
	}
	y += 50;
	Bonds->Bond_Edit = new Fl_Text_Editor (25, y, 500, 95,
					       "  Type   From  To    Diameter  Min d    Max d    Color    ");
	Bonds->BondBuffer = new Fl_Text_Buffer;
	Bonds->BondBuffer->add_modify_callback ((Fl_Text_Modify_Cb) Modify_Bonds_cb,
						(void *) NULL);
	Bonds->Bond_Edit->textfont (FL_COURIER);
	Bonds->Bond_Edit->textsize (12);
	Bonds->Bond_Edit->buffer (Bonds->BondBuffer);
	Bonds->Bond_Edit->labelfont (FL_COURIER_BOLD);
	y += 95;
	Bonds->BondInstr = new Fl_Output (110, y, 300, 0, "Press 'Add' to replace "
					  "selected line - 'Remove' to delete it");
	Bonds->BondInstr->hide ();
	Bonds->BondInstr->align (FL_ALIGN_BOTTOM);
	Bonds->BondInstr1 = new Fl_Output (110, y, 300, 0, "Highlight text above "
					   "or double click to edit line");
	Bonds->BondInstr1->hide ();
	Bonds->BondInstr1->align (FL_ALIGN_BOTTOM);
	y += 50;
	Fl_Input *oo = Bonds->New_Bond_From = new Fl_Input (20, y, 50, 25, "From");

	oo->align (FL_ALIGN_TOP);
	oo->labelfont (1);
	Fl_Input *op = Bonds->New_Bond_To = new Fl_Input (80, y, 50, 25, "To");

	op->align (FL_ALIGN_TOP);
	op->callback ((Fl_Callback *) New_Bond_Input_cb);
	op->labelfont (1);
	Fl_Input *oq = Bonds->New_Bond_Dia = new Fl_Input (140, y, 50, 25, "Diameter");

	oq->align (FL_ALIGN_TOP);
	oq->callback ((Fl_Callback *) New_Bond_Input_cb);
	oq->labelfont (1);
	Fl_Input *m = Bonds->New_Bond_Min = new Fl_Input (200, y, 50, 25, "Min d");

	m->align (FL_ALIGN_TOP);
	m->callback ((Fl_Callback *) New_Bond_Input_cb);
	m->labelfont (1);
	Fl_Input *os = Bonds->New_Bond_Max = new Fl_Input (260, y, 50, 25, "Max d");

	os->align (FL_ALIGN_TOP);
	os->callback ((Fl_Callback *) New_Bond_Input_cb);
	os->labelfont (1);
	Flu_Combo_List *ot = Bonds->New_Bond_Color =
	    new Flu_Combo_List (320, y, 160, 25, "Color");
	ot->align (FL_ALIGN_TOP);
	ot->callback ((Fl_Callback *) New_Bond_Input_cb);
	Load_Color_Combo (ot);
	ot->labelfont (1);
	Fl_Group *qq = new Fl_Group (490, y, 40, 25, "Dashed");

	qq->labelfont (1);
	qq->box (FL_THIN_UP_BOX);
	Fl_Radio_Button *ou = Bonds->New_Bond_Style =
	    new Fl_Radio_Button (490, y, 10, 25, "");
	ou->callback ((Fl_Callback *) New_Bond_Input_cb);
	ou->labelfont (1);
	Bonds->New_Bond_Style->type (1);
	Bonds->New_Bond_Style->selection_color ((Fl_Color) 1);
	Bonds->New_Bond_Style->align (FL_ALIGN_RIGHT);
	Fl_Input *od = Bonds->New_Bond_Dashes = new Fl_Input (500, y, 30, 25, "");
	qq->end ();


	od->align (FL_ALIGN_TOP);
	od->callback ((Fl_Callback *) New_Bond_Input_cb);
	od->labelfont (1);

	y += 30;
	Fl_Button *om = Bonds->New_Bond_Add = new Fl_Button (180, y, 70, 25, "Add");

	om->callback ((Fl_Callback *) New_Bond_Add_cb, &zero);
	om->deactivate ();
	Fl_Button *mm = Bonds->New_Bond_Remove = new Fl_Button (270, y, 70, 25, "Remove");

	mm->callback ((Fl_Callback *) New_Bond_Add_cb, &one);
	mm->deactivate ();
	y += 30;
	Bonds->BondInstr2 = new Fl_Output (25, y, 470, 0,
					   "Changes are temporary until \"Apply\" or \"Save\" is pressed.");
	Bonds->BondInstr2->hide ();
	Bonds->BondInstr2->align (FL_ALIGN_BOTTOM);
	y += 40;
	Fl_Text_Editor *p =
	    new Fl_Text_Editor (250, y, 200, 100, "'To' Atoms & Distances");
	p->textfont (FL_COURIER);
	p->textsize (14);
	p->align (FL_ALIGN_TOP);
	p->labelfont (1);
	Bonds->Bond_Output_Buffer = new Fl_Text_Buffer;
	p->buffer (Bonds->Bond_Output_Buffer);
	Bonds->Bond_Output_Buffer->
	    add_modify_callback ((Fl_Text_Modify_Cb) Modify_Bonds_Distance_cb,
				 (void *) NULL);
#if !defined (WIN32) && !defined (__APPLE__)
	Bonds->Bond_Edit_Window->icon ((char *) drvui->icon);
#endif
	Flu_Combo_List *o = Bonds->Bond_Combo =
	    new Flu_Combo_List (100, y, 100, 25, "'From' Atom");
	o->align (FL_ALIGN_TOP);
	o->labelfont (1);
	o->callback (Bond_Combo_cb);
	y += 125;
	Fl_Button *r = new Fl_Button (125, y, 70, 25, "Close");

	r->tooltip ("Close this window and discard all changes.");
	r->callback ((Fl_Callback *) Edit_Bond_Close_cb);
	Fl_Button *s = new Fl_Button (325, y, 70, 25, "Save");

	s->tooltip
	    ("Apply current contents of top box to drawing, then close this window.");
	s->callback ((Fl_Callback *) Edit_Bond_Save_cb, &one);
	Fl_Button *a = new Fl_Button (225, y, 70, 25, "Apply");

	a->tooltip
	    ("Apply current contents of top box to drawing, but leave this window open.");
	a->callback ((Fl_Callback *) Edit_Bond_Save_cb, &zero);
	Bonds->Bond_Edit_Window->end ();

    } else {
	Bonds->BondInstr->hide ();
	Bonds->BondInstr1->hide ();
    }
    Bond_Frame_Combo_cb (NULL, NULL);
    Bonds->Bond_Edit_Window->show ();
}

void
Edit_Bond_Close_cb (void)
{
    Bonds->Bond_Edit_Window->hide ();
    Restore_Working_Copy ();	// undo any changes
    Generate_Drawing (1);	// regenerate
}

void
Edit_Bond_Save_cb (Fl_Button *, int *save)
{
// callback routine when 'save' or 'apply' button is pressed on the Edit Bonds screen
    char type[5], atom1[5], atom2[5], color[40], widget[16382];

    int i = 0;

    unsigned int j, k;

    float dia, d1, d2;

    int Frame_No = 1;

    int numd;

    char *selection;

    if (drvui->max_frame > 1) {
	Frame_No = atoi (Bonds->Frame_No->value ());
	for (i = 1, j = 1; i < drvui->nbond; i++) {	// copy parameters for other frames to
	    if (drvui->bonds[i].bond_fn != Frame_No) {	//   start of list
		if ((int) j != i) {
		    drvui->bonds[j].bond_fn = drvui->bonds[i].bond_fn;
		    strncpy (drvui->bonds[j].col_bond, drvui->bonds[i].col_bond, 25);
		    strncpy (drvui->bonds[j].bond_l1, drvui->bonds[i].bond_l1, 4);
		    strncpy (drvui->bonds[j].bond_l2, drvui->bonds[i].bond_l2, 4);
		    drvui->bonds[j].bond_l1[4] = '\0';
		    drvui->bonds[j].bond_l2[4] = '\0';
		    drvui->bonds[j].bond_max = drvui->bonds[i].bond_max;
		    drvui->bonds[j].bond_min = drvui->bonds[i].bond_min;
		    drvui->bonds[j].bond_size = drvui->bonds[i].bond_size;
		    drvui->bonds[j].bond_style = drvui->bonds[i].bond_style;
		}
		j++;
	    }
	}
	i = j - 1;
    }
    drvui->Str_File_Changed = 1;
    selection = Bonds->BondBuffer->text ();
    strcpy (widget, selection);
    free (selection);
    if (strlen (widget) < 10) {
	drvui->nbond = 0;
	Bonds->BondBuffer->text ("");
    } else {
	while (strlen (widget) > 10) {
	    i++;
	    sscanf (widget, " %s %s %s %f %f %f %s", type, atom1, atom2, &dia, &d1, &d2,
		    color);
	    if (!strcmp (type, "dash")) {
   		j = sscanf (widget, " %s %d %s %s %f %f %f %s", type, &numd, atom1, atom2, &dia, &d1, &d2,
		    color);
		if (j < 4) {
		    numd = 5;
		    sscanf (widget, " %s %s %s %f %f %f %s", type, atom1, atom2, &dia, &d1, &d2,
		    color);
		}
		drvui->bonds[i].bond_style = numd;
	    } else
		drvui->bonds[i].bond_style = 0;
	    strcpy (drvui->bonds[i].bond_l1, atom1);
	    strcpy (drvui->bonds[i].bond_l2, atom2);
	    drvui->bonds[i].bond_size = dia;
	    drvui->bonds[i].bond_min = d1;
	    drvui->bonds[i].bond_max = d2;
	    trim_string (color, 40);
	    if (!strlen (color))
		strcpy (color, "Gray20");
	    strncpy (drvui->bonds[i].col_bond, color, 40);
	    drvui->bonds[i].bond_fn = Frame_No;
	    drvui->nbond = i + 1;
	    for (j = 0; j < strlen (widget); j++) {
		if (widget[j] == '\n')
		    break;
	    }
	    for (j++, k = 0; j < strlen (widget); j++)
		widget[k++] = widget[j];
	    widget[k] = 0;
	}
    }
    Update_Str (0);
    Generate_Drawing (0);
    if (*save != 3) {
	Save_Working_Copy ();
	Bonds->BondInstr2->hide ();
    }
    if (*save == 1) {
	Bonds->Bond_Edit_Window->hide ();
    }
    Fl::redraw ();
}

void
Edit_Changed_cb (int, int nInserted, int nDeleted, int, const char *, void *v)
{
    char title[128];

    Fl_Window *w = (Fl_Window *) v;

    if ((nInserted || nDeleted) && !Edit_loading) {
	Edit_changed = 1;
	strcpy (title, Edit_title);
	strcat (title, " (modified)");
	w->label (title);
    }
}

void
Edit_Ellipsoid_cb (void)
{
// Callback routine to show the Ellipsoid Parameter screen and load the widgets on that page
    char string[100];

    int i, j;

    char atom[5];

    char color[40];

    char widget[16382];

    if (!strlen (drvui->CurFile->value ())) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }
/*
    if (!natom) {         // No atoms, no ellipsoids...
        Error_Box("This structure file does not contain any atoms.");
        return;
    }
*/
    Save_Working_Copy ();
    if (!ellipsoids) {
	ellipsoids = new Ellipsoids;	// Create the screen
#if !defined (WIN32) && !defined (__APPLE__)
	ellipsoids->Ellips_Window->icon ((char *) drvui->icon);
#endif
	ellipsoids->Color_Combo->list.add ("");
	ellipsoids->Color_Combo->value ("");
	Load_Color_Combo (ellipsoids->Color_Combo);
	ellipsoids->show ();	// show it
	drvui->frame_no = min (drvui->frame_no, drvui->max_frame);
	set_tf_status ();	// set the flags
	sprintf (string, "%5.2f", drvui->Ellipsoid_Prob);	// load the ellipsoid probability widget
	ellipsoids->Probability->value (string);
	ellipsoids->Probability->take_focus ();	// set this widget with focus
	ellipsoids->Probability->position (0);
	ellipsoids->Probability->mark (strlen (ellipsoids->Probability->value ()));
	if (drvui->do_ellipsoids)
	    ellipsoids->Show_Ellipsoids->set ();	// Show ellipsoids checkbox set
	if (drvui->El_Cutout) {
	    ellipsoids->Use_Cutouts->set ();	// cutout check box
	    if (!strlen (drvui->Cutout_color))
		strcpy (drvui->Cutout_color, "Gray20");
	    ellipsoids->Cutout_Color->value (drvui->Cutout_color);	// Cutout Color widget
	} else {
	    ellipsoids->Cutout_Color->deactivate ();
	}
	ellipsoids->Axis_Color->value (drvui->Ellipaxis_color);	// Load ellipsoid axis color
	sprintf (string, "%6.2f", drvui->Ellipaxis_width);
	ellipsoids->Axis_Width->value (string);	// and width
	ellipsoids->Atom_Combo->list.add ("");
	ellipsoids->Atom_Combo->value ("");
	for (i = 0; i < natom; i++) {
	    if (drvui->atoms[i].atom_fn != drvui->frame_no)
		continue;
	    if (drvui->atoms[i].TF_status > 0) {
		strcpy (string, drvui->atoms[i].atom_l);
		while (strlen (string) < 4)
		    strcat (string, " ");
		sprintf (widget, "%4s%2d", string, drvui->atoms[i].sv_atom_n);
		ellipsoids->Atom_Combo->list.add (widget);
	    }
	}
	for (i = 0; i < natom; i++) {
	    if (drvui->atoms[i].atom_fn != drvui->frame_no)
		continue;
	    if (drvui->atoms[i].TF_status > 0) {
		int j, haveit;

		haveit = 0;
		for (j = 0; j < i; j++)
		    if (drvui->atoms[j].TF_status &&
			check_atom_name (drvui->atoms[i].atom_l, drvui->atoms[j].atom_l))
			haveit = 1;
		if (haveit == 1)
		    continue;
		strcpy (string, drvui->atoms[i].atom_l);
		while (strlen (string) < 4)
		    strcat (string, " ");
		sprintf (widget, "%4s *", string);
		ellipsoids->Atom_Combo->list.add (widget);
	    }
	}
	widget[0] = '\0';
	for (i = 1; i < drvui->n_ellips; i++) {
	    if (drvui->ellips[i].ell_type < 1000)
		continue;	// skip for npd or not displayed
	    memset (atom, 0, 5);
	    strncpy (atom, drvui->ellips[i].ellips_l, 4);
	    atom[4] = 0;
	    for (j = 3; j >= 0; --j) {
		if (atom[j] == ' ')
		    atom[j] = 0;
	    }
	    j = drvui->ellips[i].save_el_number;
	    memset (color, 0, 40);
	    memset (string, 0, 100);
	    strncpy (color, drvui->ellips[i].ellips_col, 39);
	    trim_string (color, 40);
	    if (!strlen (color))
		strcpy (color, "Gray20");
	    if (j == -1) {
		int haveit = 0;

		for (j = 1; j < i; j++) {
		    if (check_atom_name
			(drvui->ellips[i].ellips_l, drvui->ellips[j].ellips_l))
			haveit = 1;
		}
		if (haveit == 0)
		    sprintf (string, "%4s   *   %s\n", atom, color);

	    } else
		sprintf (string, "%4s %3d   %s\n", atom, j, color);
	    strcat (widget, string);
	}
	ellipsoids->ColorInputBuf->text (widget);
    }
    ellipsoids->Ellips_Window->show ();
}

void
Edit_Ellipsoid_Close_cb (Fl_Button *, void *)
{
// callback from close button on ellipsoids screen
    drvui->destroy |= ELLIPSOID;
    ellipsoids->Ellips_Window->hide ();
    Restore_Working_Copy ();	// undo any changes
    Generate_Drawing (1);	// regenerate
}

void
Edit_Ellipsoid_Save_cb (Fl_Button *, int *save)
{
// callback routine when 'save' or 'apply' button is pressed on the Ellipsoids Parameter screen
    char atom[5], color[40], widget[16382], num[4];

    int i;

    int j, k;

    char string[128];

    char *selection;

    drvui->Str_File_Changed = 1;
    for (i = 1; i < drvui->n_ellips; i++)
	drvui->ellips[i].ell_type = 1;
    drvui->Ellipsoid_Prob = (float) atof (ellipsoids->Probability->value ());
    selection = ellipsoids->ColorInputBuf->text ();
    strcpy (widget, selection);
    free (selection);
    while (strlen (widget) > 10) {
	sscanf (widget, "%s %s %20s", atom, num, color);
	while (strlen (atom) < 4)
	    strcat (atom, " ");
	if (strstr (num, "*"))
	    j = -1;
	else
	    j = atoi (num);
	for (i = 1; i < drvui->n_ellips; i++) {
	    if (check_atom_name (drvui->ellips[i].ellips_l, atom) &&
		(j == drvui->ellips[i].ellips_n || j == -1)) {
		trim_string (color, 40);
		if (!strlen (color))
		    strcpy (color, "Gray20");
		strcpy (drvui->ellips[i].ellips_col, color);	/* copy color */
		drvui->ellips[i].ell_type = 1001;
		drvui->ellips[i].save_el_number = j;
	    }
	}
	for (j = 0; j < (int) strlen (widget); j++) {
	    if (widget[j] == '\n')
		break;
	}
	for (j++, k = 0; j < (int) strlen (widget); j++)
	    widget[k++] = widget[j];
	widget[k] = 0;
    }
    if (ellipsoids->Use_Cutouts->value () != 0) {
	drvui->El_Cutout = 1;
	strcpy (string, ellipsoids->Cutout_Color->value ());
	if (!strlen (string))
	    strcpy (string, "Gray20");
	trim_string (string, 40);
	strcpy (drvui->Cutout_color, string);
	strcpy (string, ellipsoids->Axis_Color->value ());	// Get ellipsoid axis color
	trim_string (string, 40);
	if (!strlen (string))
	    strcpy (string, "Gray20");
	strcpy (drvui->Ellipaxis_color, string);
    } else {
	strcpy (drvui->Cutout_color, "");
	drvui->El_Cutout = 0;
    }
    if (ellipsoids->Show_Ellipsoids->value () != 0)
	drvui->do_ellipsoids = 1;
    else
	drvui->do_ellipsoids = 0;
    drvui->Ellipaxis_width = (float) atof (ellipsoids->Axis_Width->value ());
    Update_Str (0);
    Generate_Drawing (0);
    if (*save != 3) {
	Save_Working_Copy ();
	ellipsoids->Instr2->hide ();
    }
    if (*save == 1) {
	ellipsoids->Ellips_Window->hide ();
    }
    Fl::redraw ();		// update the screen
}

void
Edit_LonePair_cb (void)
{
// Callback routine to show the Lone-Pair Parameter screen and load the widgets on that page
    char string[100];

    static int one = 1;

    static int zero = 0;

    int i;

    int y = 50;

    if (!strlen (drvui->CurFile->value ())) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }
/*
    if (!natom) {         // No atoms, no lone-pairs...
        Error_Box("This structure file does not contain any atoms.");
        return;
    }
*/
    Save_Working_Copy ();
    if (!LonePairs) {
	LonePairs = new LonePairParam;	// new instance of the lone-pair parameters
	LonePairs->LonePair_Edit_Window =
	    new Fl_Window (50, y, 520, 350, "Edit Lone-Pair Parameters");
	LonePairs->LonePair_Edit_Window->
	    callback ((Fl_Callback *) Edit_Lone_Pair_Close_cb);
	if (drvui->max_frame > 1) {
	    Flu_Combo_List *o = LonePairs->Frame_No = new Flu_Combo_List (220, y - 25, 75,
									  25,
									  "Frame No.");
	    o->align (FL_ALIGN_TOP);
	    o->callback (LonePair_Frame_Combo_cb);
	    o->labelfont (1);
	    for (i = 1; i <= drvui->max_frame; i++) {
		sprintf (string, "%d", i);
		o->list.add (string);
	    }
	    o->pop_height (20 * drvui->max_frame);
	    o->value ("1");
	}
	y += 25;
	LonePairs->LonePair_Edit = new Fl_Text_Editor (25, y, 460, 75,
						       "From  #e  Height   Rad 1   Rad 2    Color        ");
	LonePairs->LonePairBuffer = new Fl_Text_Buffer;
	LonePairs->LonePairBuffer->
	    add_modify_callback ((Fl_Text_Modify_Cb) Modify_LonePair_cb, (void *) NULL);
	LonePairs->LonePair_Edit->textfont (FL_COURIER);
	LonePairs->LonePair_Edit->textsize (12);
	LonePairs->LonePair_Edit->buffer (LonePairs->LonePairBuffer);
	LonePairs->LonePair_Edit->labelfont (FL_COURIER_BOLD);
	y += 85;
	LonePairs->LonePairInst = new Fl_Output (135, y, 250, 0, "Press 'Add' to replace "
						 "selected line - 'Remove' to delete it");
	LonePairs->LonePairInst->hide ();
	LonePairs->LonePairInst->align (FL_ALIGN_BOTTOM);
	LonePairs->LonePairInst1 = new Fl_Output (135, y, 250, 0, "Highlight text above "
						  "or double click to edit line");
	LonePairs->LonePairInst1->hide ();
	LonePairs->LonePairInst1->align (FL_ALIGN_BOTTOM);
	y += 40;
	Flu_Combo_List *o = LonePairs->LonePair_Combo =
	    new Flu_Combo_List (20, y, 50, 25, "From");
	o->align (FL_ALIGN_TOP);
	o->labelfont (1);
	o->callback (Lone_Pair_Combo_cb);
	Fl_Input *op = LonePairs->Number = new Fl_Input (80, y, 50, 25, "No. Elec.");

	op->align (FL_ALIGN_TOP);
	op->labelfont (1);
	op->callback (Lone_Pair_Combo_cb);
	Fl_Input *oq = LonePairs->Height = new Fl_Input (140, y, 50, 25, "Height");

	oq->align (FL_ALIGN_TOP);
	oq->labelfont (1);
	oq->callback (Lone_Pair_Combo_cb);
	Fl_Input *m = LonePairs->Radius1 = new Fl_Input (200, y, 50, 25, "Rad. 1");

	m->align (FL_ALIGN_TOP);
	m->labelfont (1);
	Fl_Input *os = LonePairs->Radius2 = new Fl_Input (260, y, 50, 25, "Rad. 2");

	os->align (FL_ALIGN_TOP);
	os->labelfont (1);
	os->callback (Lone_Pair_Combo_cb);
	Flu_Combo_List *ot = LonePairs->LonePair_Color =
	    new Flu_Combo_List (320, y, 160, 25, "Color");
	ot->align (FL_ALIGN_TOP);
	Load_Color_Combo (ot);
	ot->labelfont (1);
	ot->callback (Lone_Pair_Combo_cb);
	y += 30;
	Fl_Button *om = LonePairs->LonePair_Add = new Fl_Button (180, y, 70, 25, "Add");

	om->callback ((Fl_Callback *) New_Lone_Pair_Add_cb, &one);
	om->deactivate ();
	Fl_Button *pm = LonePairs->LonePair_Remove =
	    new Fl_Button (270, y, 70, 25, "Remove");
	pm->callback ((Fl_Callback *) New_Lone_Pair_Add_cb, &zero);
	pm->deactivate ();
	y += 30;
	LonePairs->LonePairInst2 = new Fl_Output (25, y, 470, 0,
						  "Changes are temporary until \"Apply\" or \"Save\" is pressed.");
	LonePairs->LonePairInst2->hide ();
	LonePairs->LonePairInst2->align (FL_ALIGN_BOTTOM);
	y += 30;
	Fl_Button *r = new Fl_Button (125, y, 70, 25, "Close");

	r->tooltip ("Close this window and discard all changes.");
	r->callback ((Fl_Callback *) Edit_Lone_Pair_Close_cb);
	Fl_Button *s = new Fl_Button (325, y, 70, 25, "Save");

	s->tooltip
	    ("Apply current contents of top box to drawing, then close this window.");
	s->callback ((Fl_Callback *) Edit_Lone_Pair_Save_cb, &one);
	Fl_Button *a = new Fl_Button (225, y, 70, 25, "Apply");

	a->tooltip
	    ("Apply current contents of top box to drawing, but leave this window open.");
	a->callback ((Fl_Callback *) Edit_Lone_Pair_Save_cb, &zero);
#if !defined (WIN32) && !defined (__APPLE__)
	LonePairs->LonePair_Edit_Window->icon ((char *) drvui->icon);
#endif
	LonePairs->LonePair_Edit_Window->end ();
    }
    LonePair_Frame_Combo_cb (NULL, NULL);
    LonePairs->LonePair_Edit_Window->show ();
}

void
Edit_Lone_Pair_Close_cb (void)
{
    drvui->destroy |= LONEPAIR;
    LonePairs->LonePair_Edit_Window->hide ();
    Restore_Working_Copy ();
    Generate_Drawing (1);
}

void
Edit_Lone_Pair_Save_cb (Fl_Button *, int *save)
{
// callback routine when 'save' or 'apply' button is pressed on the Edit Lone Pairs screen
    char atom[4], color[40], widget[16382];

    int i = 0;

    int no = 0;

    unsigned int j, k;

    float height, d1, d2;

    int Frame_No = 1;

    char *selection;

    if (drvui->max_frame > 1) {
	Frame_No = atoi (LonePairs->Frame_No->value ());
	for (i = 1, j = 1; i < drvui->ncone; i++) {	// copy parameters for other frames to
	    if (drvui->cones[i].cone_fn != Frame_No) {	//   start of list
		if ((int) j != i) {
		    drvui->cones[j].cone_height = drvui->cones[i].cone_height;
		    drvui->cones[j].cone_fn = drvui->cones[i].cone_fn;
		    strncpy (drvui->cones[j].col_cone, drvui->cones[i].col_cone, 25);
		    strncpy (drvui->cones[j].cone_l1, drvui->cones[i].cone_l1, 4);
		    drvui->cones[j].cone_max = drvui->cones[i].cone_max;
		    drvui->cones[j].cone_min = drvui->cones[i].cone_min;
		}
		j++;
	    }
	}
	i = j - 1;
    }

    drvui->Str_File_Changed = 1;
    selection = LonePairs->LonePairBuffer->text ();
    strcpy (widget, selection);
    free (selection);
    if (strlen (widget) < 10) {
	i = 0;
	strcpy (widget, "");
	LonePairs->LonePairBuffer->text (widget);
    } else {
	while (strlen (widget) > 10) {
	    i++;
	    sscanf (widget, " %s %d %f %f %f %s", atom, &no, &height, &d1, &d2, color);
	    strncpy (drvui->cones[i].cone_l1, atom, 4);
	    drvui->cones[i].numlonepairs = no;
	    drvui->cones[i].cone_height = height;
	    drvui->cones[i].cone_min = d1;
	    drvui->cones[i].cone_max = d2;
	    trim_string (color, 40);
	    if (!strlen (color))
		strcpy (color, "Gray20");
	    strcpy (drvui->cones[i].col_cone, color);
	    drvui->cones[i].cone_fn = Frame_No;
	    for (j = 0; j < strlen (widget); j++) {
		if (widget[j] == '\n')
		    break;
	    }
	    for (j++, k = 0; j < strlen (widget); j++)
		widget[k++] = widget[j];
	    widget[k] = 0;
	}
    }
    drvui->ncone = i + 1;
    Update_Str (0);
    Generate_Drawing (0);
    if (*save != 3) {
	Save_Working_Copy ();
	LonePairs->LonePairInst2->hide ();
    }
    if (*save == 1) {
	LonePairs->LonePair_Edit_Window->hide ();
    }
    Fl::redraw ();
}

void
Edit_Maps_cb (void)
{
// Callback routine to show the Map Parameter screen and load the widgets on that page
    char string[128];
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
    if (!Maps) {
	char title[30];

	if (drvui->Fourier2d)
	    strcpy (title, "Edit 2D Map parameters");
	else
	    strcpy (title, "Edit 3D Map Parameters");
	Maps = new MapsParam;	// new instance of the map parameters
	
	wh = 550;
	if (drvui->modulated)
	    wh += 55;
	if (drvui->max_frame > 1) 
	    wh += 60;
	Maps->Maps_Edit_Window = new Fl_Window (50, 50, 500, wh, title);

	Maps->Maps_Edit_Window->callback ((Fl_Callback *) Edit_Maps_Close_cb);
	y = 30;
	if (drvui->max_frame > 1) {
	    Maps->Frame_No =
		new Flu_Combo_List (220, y, 75, 25, "Frame No.");
	    Maps->Frame_No->align (FL_ALIGN_TOP);
	    Maps->Frame_No->callback (Maps_Frame_Combo_cb);
	    Maps->Frame_No->labelfont (1);
	    for (i = 1; i <= drvui->max_frame; i++) {
		sprintf (string, "%d", i);
		Maps->Frame_No->list.add (string);
	    }
	    Maps->Frame_No->pop_height (20 * drvui->max_frame);
	    Maps->Frame_No->value ("1");
            y += 50;
	}
	Maps->Filename = new Fl_Output (10, y, 260, 25, "Map Filename");
	Maps->Filename->align (FL_ALIGN_TOP);
	Maps->Filename->color ((Fl_Color) 17);
	Maps->Filename->labelfont (FL_COURIER_BOLD);

	Maps->Map_Browse = new Fl_Button (280, y, 70, 25, "Browse");
	Maps->Map_Browse->tooltip ("Browse for new Fourier file.");
	Maps->Map_Browse->callback ((Fl_Callback *) Browse_Map_File_cb);

	Maps->MapType = new Flu_Combo_List (360, y, 130, 25, "File Type");
	Maps->MapType->align (FL_ALIGN_TOP);
	Maps->MapType->list.add ("GSAS - grd");
	Maps->MapType->list.add ("JANA - stf");
	Maps->MapType->list.add ("WIEN - w2k");
	Maps->MapType->list.add ("VASP - vsp");
	Maps->MapType->list.add ("FullProf - flp");
	Maps->MapType->list.add ("CIF FoFc - fcf");
	Maps->MapType->list.add ("O Format - dn6");
	Maps->MapType->list.add ("JANA FoFc - m80");
	Maps->MapType->list.add ("Exciting - exc");
	Maps->MapType->list.add ("JANA - m81");
	Maps->MapType->list.add ("XCrysDen - xsf");
	Maps->MapType->pop_height (120);
	Maps->MapType->labelfont (FL_COURIER_BOLD);
	Maps->MapType->color ((Fl_Color) 17);
	Maps->MapType->callback ((Fl_Callback *) MapType_cb);

	y += 40;
	Maps->Map_Info = new Fl_Button (20, y - 10, 110, 25, "Show Map Info");
	Maps->Map_Info->tooltip ("Display header information from Fourier file.");
	Maps->Map_Info->callback ((Fl_Callback *) Map_Info_cb);

	Maps->MapCalc = new Fl_Button (20, y + 25, 110, 25, "Save Calc Map");
	Maps->MapCalc->tooltip ("Save Calculated 'fcf' Map in 'grd' format");
	Maps->MapCalc->callback ((Fl_Callback *) Write_Map_cb);
	Maps->MapCalcType = new Flu_Combo_List (20, y + 70, 110, 25, "Calc Type");
	Maps->MapCalcType->align (FL_ALIGN_TOP);
	Maps->MapCalcType->list.add ("Fo");
	Maps->MapCalcType->list.add ("Fc");
	Maps->MapCalcType->list.add ("Fo-Fc");
	Maps->MapCalcType->list.add ("2Fo-Fc");
	Maps->MapCalcType->list.add ("Fo2");
	Maps->MapCalcType->pop_height (80);
	Maps->MapCalcType->color ((Fl_Color) 17);
	if (Map_Info.map_type == 0)
	    Maps->MapCalcType->value ("Fo");
	else if (Map_Info.map_type == 1)
	    Maps->MapCalcType->value ("Fc");
	else if (Map_Info.map_type == 2)
	    Maps->MapCalcType->value ("Fo-Fc");
	else if (Map_Info.map_type == 3)
	    Maps->MapCalcType->value ("2Fo-Fc");
	else if (Map_Info.map_type == 4)
	    Maps->MapCalcType->value ("Fo2");

	Maps->Resolution = new Fl_Input (20, y + 110, 110, 25, "Resolution");
	Maps->Resolution->align (FL_ALIGN_TOP);
	sprintf(widget, "%i", Map_Info.res);
	Maps->Resolution->value(widget);

	y += 15;
	strcpy (widget, "");
	Maps->MapsBuffer = new Fl_Text_Buffer;
	Maps->MapsBuffer->add_modify_callback ((Fl_Text_Modify_Cb) Modify_Maps_cb,
					       (void *) NULL);
	if (!drvui->Fourier2d) {
	    Fl_Text_Editor *o =
		new Fl_Text_Editor (150, y, 300, 120, "Level     Type     Color     ");
	    o->labelfont (FL_COURIER_BOLD);
	    o->textfont (FL_COURIER);
	    o->buffer (Maps->MapsBuffer);
	} else {
	    Fl_Text_Editor *o =
		new Fl_Text_Editor (150, y, 300, 120, "Lower    Step   Upper Color ");
	    o->labelfont (FL_COURIER_BOLD);
	    o->textfont (FL_COURIER);
	    o->buffer (Maps->MapsBuffer);
	}
	Maps->MapsInstr = new Fl_Output (20, y + 130, 450, 0, "Press 'Add' to replace "
					 "selected line - 'Remove' to delete it");
	Maps->MapsInstr->hide ();
	Maps->MapsInstr->align (FL_ALIGN_BOTTOM);
	Maps->MapsInstr1 = new Fl_Output (70, y + 130, 450, 0, "Highlight text above "
					  "or double click to edit line");
	Maps->MapsInstr1->align (FL_ALIGN_BOTTOM);
	y += 160;
	if (!drvui->Fourier2d) {
	    Maps->Level = new Fl_Input (80, y, 70, 25, "Contour Level");
	    Maps->Level->align (FL_ALIGN_BOTTOM);
	    Maps->Level->callback ((Fl_Callback *) New_Map_Input_cb);
	    Maps->Type = new Flu_Combo_List (170, y, 90, 25, "Type");
	    Maps->Type->align (FL_ALIGN_BOTTOM);
	    Maps->Type->callback ((Fl_Callback *) New_Map_Input_cb);
	    Maps->Type->list.add ("mesh");
	    Maps->Type->list.add ("solid");
	    Maps->Type->pop_height (40);
	} else {
	    Maps->Level = new Fl_Input (40, y, 70, 25, "Lower");
	    Maps->Level->align (FL_ALIGN_BOTTOM);
	    Maps->Level->callback ((Fl_Callback *) New_Map_Input_cb);
	    Maps->Step = new Fl_Input (120, y, 70, 25, "Step");
	    Maps->Step->align (FL_ALIGN_BOTTOM);
	    Maps->Step->callback ((Fl_Callback *) New_Map_Input_cb);
	    Maps->Top = new Fl_Input (200, y, 70, 25, "Upper");
	    Maps->Top->align (FL_ALIGN_BOTTOM);
	    Maps->Top->callback ((Fl_Callback *) New_Map_Input_cb);
	}
	Maps->Color = new Flu_Combo_List (280, y, 160, 25, "Color");
	Load_Color_Combo (Maps->Color);
	Maps->Color->align (FL_ALIGN_TOP);
	Maps->Color->callback ((Fl_Callback *) New_Map_Input_cb);
	Load_Color_Combo (Maps->Color);
	Maps->Color->align (FL_ALIGN_BOTTOM);
	y += 45;
	Fl_Button *om;

	Fl_Button *mm;

	om = Maps->Add_Button = new Fl_Button (180, y, 70, 25, "Add");
	om->callback ((Fl_Callback *) New_Map_Add_cb, &one);
	om->tooltip ("When active, press to transfer data in boxes to window above");
	om->deactivate ();
	mm = Maps->Remove_Button = new Fl_Button (270, y, 70, 25, "Remove");
	mm->callback ((Fl_Callback *) New_Map_Add_cb, &zero);
	mm->tooltip ("When active, press to remove highlighted line.");
	mm->deactivate ();
	y += 30;
	Maps->MapsInstr2 = new Fl_Output (25, y, 450, 0,
					  "Changes are temporary until \"Apply\" or \"Save\" is pressed.");
	Maps->MapsInstr2->hide ();
	Maps->MapsInstr2->align (FL_ALIGN_BOTTOM);
	y += 30;
	Maps->XMin = new Fl_Input (80, y, 100, 25, "XMin");
	Maps->XMin->align (FL_ALIGN_BOTTOM);
	Maps->XMin->color ((Fl_Color) 17);
	Maps->YMin = new Fl_Input (200, y, 100, 25, "YMin");
	Maps->YMin->align (FL_ALIGN_BOTTOM);
	Maps->YMin->color ((Fl_Color) 17);
	Maps->ZMin = new Fl_Input (320, y, 100, 25, "ZMin");
	Maps->ZMin->align (FL_ALIGN_BOTTOM);
	Maps->ZMin->color ((Fl_Color) 17);
	y += 55;
	Maps->XMax = new Fl_Input (80, y, 100, 25, "XMax");
	Maps->XMax->align (FL_ALIGN_BOTTOM);
	Maps->XMax->color ((Fl_Color) 17);
	Maps->YMax = new Fl_Input (200, y, 100, 25, "YMax");
	Maps->YMax->align (FL_ALIGN_BOTTOM);
	Maps->YMax->color ((Fl_Color) 17);
	Maps->ZMax = new Fl_Input (320, y, 100, 25, "ZMax");
	Maps->ZMax->align (FL_ALIGN_BOTTOM);
	Maps->ZMax->color ((Fl_Color) 17);
	y += 55;
	if (drvui->modulated) {
	    Maps->X4 = new Flu_Spinner (110, y, 50, 25, "x4");
	    Maps->X4->value (x4Val);
	    Maps->X4->minimum (Map_Info.x4lim[0]);
	    Maps->X4->maximum (Map_Info.x4lim[1]);
	    Maps->X4->step (x4step);
	    Maps->X4->callback ((Fl_Callback *) Edit_Maps_Save_cb, &zero);
	    if (x4step < 0.001)
		Maps->X4->deactivate ();
	    Maps->X5 = new Flu_Spinner (230, y, 50, 25, "x5");
	    Maps->X5->value (x5Val);
	    Maps->X5->minimum (Map_Info.x5lim[0]);
	    Maps->X5->maximum (Map_Info.x5lim[1]);
	    Maps->X5->step (x5step);
	    Maps->X5->callback ((Fl_Callback *) Edit_Maps_Save_cb, &zero);
	    if (x5step < 0.001)
		Maps->X5->deactivate ();
	    Maps->X6 = new Flu_Spinner (350, y, 50, 25, "x6");
	    Maps->X6->value (x6Val);
	    Maps->X6->minimum (Map_Info.x6lim[0]);
	    Maps->X6->maximum (Map_Info.x6lim[1]);
	    Maps->X6->step (x6step);
	    Maps->X6->callback ((Fl_Callback *) Edit_Maps_Save_cb, &zero);
	    if (x6step < 0.001)
		Maps->X6->deactivate ();
	    y += 55;
	}
	Fl_Button *q = new Fl_Button (195, y, 130, 25, "Edit 2D Map Slice");
	q->tooltip("Choose 2D slice through arbitrary pplane");
	q->callback((Fl_Callback *) Edit_Slice_cb);
	y += 45;

	Fl_Button *r = new Fl_Button (125, y, 70, 25, "Close");
	r->tooltip ("Close this window and discard all changes.");
	r->callback ((Fl_Callback *) Edit_Maps_Close_cb);
	Fl_Button *s = new Fl_Button (325, y, 70, 25, "Save");

	s->tooltip
	    ("Apply current contents of top box to drawing, then close this window.");
	s->callback ((Fl_Callback *) Edit_Maps_Save_cb, &one);
	Fl_Button *a = new Fl_Button (225, y, 70, 25, "Apply");

	a->tooltip
	    ("Apply current contents of top box to drawing, but leave this window open.");
	a->callback ((Fl_Callback *) Edit_Maps_Save_cb, &zero);
#if !defined (WIN32) && !defined (__APPLE__)
	Maps->Maps_Edit_Window->icon ((char *) drvui->icon);
#endif
    }
    strcpy (widget, "");
    for (i = 1; i <= drvui->numOfFourierContours; i++) {	// fill in the widgets
	if (!drvui->Fourier2d) {
	    if (drvui->fourier[i].FourierContourSolid)
		strcpy (type, "solid");
	    else
		strcpy (type, "mesh ");
	    sprintf (string, "%8.3f     %6s    %s\n",
		     drvui->fourier[i].FourierContourLevel, type,
		     drvui->fourier[i].FourierContourColor);
	} else {
	    sprintf (string, "%7.3f %7.3f %7.3f %s\n",
		     drvui->fourier[i].FourierContourLevel,
		     drvui->fourier[i].FourierContourStep,
		     drvui->fourier[i].FourierContourTop,
		     drvui->fourier[i].FourierContourColor);
	}
	strcat (widget, string);
    }
    Maps->MapsBuffer->text (widget);
    if (ReadFourMap) {
	if (FourierMapType == 1)
	    Maps->MapType->value ("GSAS - grd");
	if (FourierMapType == 2)
	    Maps->MapType->value ("JANA - stf");
	if (FourierMapType == 3)
	    Maps->MapType->value ("WIEN - w2k");
	if (FourierMapType == 4)
	    Maps->MapType->value ("VASP - vsp");
	if (FourierMapType == 5)
	    Maps->MapType->value ("FullProf - flp");
	if (FourierMapType == 6)
	    Maps->MapType->value ("CIF FoFc - fcf");
	if (FourierMapType == 7)
	    Maps->MapType->value ("O Format - dn6");
	if (FourierMapType == 8)
	    Maps->MapType->value ("JANA FoFc - m80");
	if (FourierMapType == 9)
	    Maps->MapType->value ("Exciting - exc");
	if (FourierMapType == 10)
	    Maps->MapType->value ("JANA - m81");
	if (FourierMapType == 11)
	    Maps->MapType->value ("XCrysDen - xsf");
	Maps->Filename->value (FourierFileName);
	Maps_Frame_Combo_cb (NULL, NULL);
    }
    Maps->MapsInstr1->show ();

    Maps->Maps_Edit_Window->end ();
    Maps->Maps_Edit_Window->show ();
    MapType_cb ();
}

void
Edit_Maps_Close_cb (void)
{
    Maps->Maps_Edit_Window->~Fl_Window ();	// this window needs to be deleted
    delete (Maps->Maps_Edit_Window);	// not just killed (2d/3d nature might change)
    delete (Maps->MapsBuffer);
    delete (Maps);
    Maps = NULL;
    Restore_Working_Copy ();	// undo any changes
    Generate_Drawing (1);	// regenerate
}

void
Edit_Maps_Save_cb (Fl_Button *, int *save)
{
// callback routine when 'save' or 'apply' button is pressed on the Edit Maps screen
    char type[5], color[40], widget[16382];

    unsigned int j, k;

    int i = 0;

    int Frame_No = 1;

    float level, step, top;

    char *selection = Maps->MapsBuffer->text ();

    if (drvui->max_frame > 1)
	Frame_No = atoi (Maps->Frame_No->value ());

    strcpy (widget, selection);
    free (selection);
    if (strlen (widget) < 10) {
	drvui->numOfFourierContours = 0;
	strcpy (widget, "");
	Maps->MapsBuffer->text (widget);
    } else {
	while (strlen (widget) > 10) {
	    i++;
	    if (drvui->Fourier2d)
		sscanf (widget, "%f %f %f %s", &level, &step, &top, color);
	    else
		sscanf (widget, "%f %s %s", &level, type, color);
	    trim_string (color, 40);
	    if (!strlen (color))
		strcpy (color, "Gray20");
	    drvui->fourier[i].FourierContourLevel = level;
	    if (!drvui->Fourier2d) {
		if (strncmp (type, "mesh", 4) == 0)
		    drvui->fourier[i].FourierContourSolid = 0;
		else
		    drvui->fourier[i].FourierContourSolid = 1;
	    } else {
		drvui->fourier[i].FourierContourStep = step;
		drvui->fourier[i].FourierContourTop = top;
	    }
	    strcpy (drvui->fourier[i].FourierContourColor, color);
	    drvui->numOfFourierContours = i;
	    for (j = 0; j < strlen (widget); j++) {
		if (widget[j] == '\n')
		    break;
	    }
	    for (j++, k = 0; j < strlen (widget); j++)
		widget[k++] = widget[j];
	    widget[k] = 0;
	}
    }
    strcpy (widget, Maps->Filename->value ());
    if (strcmp (widget, FourierFileName)) {
	strcpy (FourierFileName, widget);	// File name changed
	if (FourierPt) {
	    free (FourierPt);
	    FourierPt = NULL;
	}
    }
    if (!strcmp (Maps->MapType->value (), "GSAS - grd"))
	FourierMapType = 1;
    if (!strcmp (Maps->MapType->value (), "JANA - stf"))
	FourierMapType = 2;
    if (!strcmp (Maps->MapType->value (), "WIEN - w2k"))
	FourierMapType = 3;
    if (!strcmp (Maps->MapType->value (), "VASP - vsp"))
	FourierMapType = 4;
    if (!strcmp (Maps->MapType->value (), "FullProf - flp"))
	FourierMapType = 5;
    if (!strcmp (Maps->MapType->value (), "CIF FoFc - fcf"))
	FourierMapType = 6;
    if (!strcmp (Maps->MapType->value (), "O Format - dn6"))
	FourierMapType = 7;
    if (!strcmp (Maps->MapType->value (), "JANA FoFc - m80"))
	FourierMapType = 8;
    if (!strcmp (Maps->MapType->value (), "Exciting - exc"))
	FourierMapType = 9;
    if (!strcmp (Maps->MapType->value (), "JANA - m81"))
	FourierMapType = 10;
    if (!strcmp (Maps->MapType->value (), "XCrysDen - xsf"))
	FourierMapType = 11;
    sscanf (Maps->XMin->value (), "%f", &drvui->frames[Frame_No].map_lim[0]);
    sscanf (Maps->YMin->value (), "%f", &drvui->frames[Frame_No].map_lim[1]);
    sscanf (Maps->ZMin->value (), "%f", &drvui->frames[Frame_No].map_lim[2]);
    sscanf (Maps->XMax->value (), "%f", &drvui->frames[Frame_No].map_lim[3]);
    sscanf (Maps->YMax->value (), "%f", &drvui->frames[Frame_No].map_lim[4]);
    sscanf (Maps->ZMax->value (), "%f", &drvui->frames[Frame_No].map_lim[5]);
    if (FourierMapType == 6 || FourierMapType == 8) {
	if (const char *temp = Maps->MapCalcType->value ()) {
	    int i = 0;

	    if (!strcmp (temp, "Fo"))
		i = 0;
	    else if (!strcmp (temp, "Fc"))
		i = 1;
	    else if (!strcmp (temp, "Fo-Fc"))
		i = 2;
	    else if (!strcmp (temp, "2Fo-Fc"))
		i = 3;
	    if (i != Map_Info.map_type) {
		Map_Info.map_type = i;
		Map_Info.info_valid = 0;
		if (FourierPt) {
		    free (FourierPt);
		    FourierPt = NULL;
		}
	    }
	}
    }
    i = Map_Info.res;
    sscanf (Maps->Resolution->value (), "%i", &Map_Info.res);
    if (Map_Info.res <= 0) { 
	Map_Info.res = 1;
	Maps->Resolution->value("1");
    }
    if ((int)i != Map_Info.res) {
 	if (FourierPt) {
	    free (FourierPt);
	    FourierPt = NULL;
	}
    }
    if (drvui->modulated) {
	if (x4step > 0.0f)
	    drvui->frames[Frame_No].map_lim[6] = (float) Maps->X4->value ();
	if (x5step > 0.0f)
	    drvui->frames[Frame_No].map_lim[7] = (float) Maps->X5->value ();
	if (x6step > 0.0f)
	    drvui->frames[Frame_No].map_lim[8] = (float) Maps->X6->value ();
    }
    drvui->Str_File_Changed = 1;
    if (*save != 3) {
	Save_Working_Copy ();
	Maps->MapsInstr2->hide ();
    }
    if (*save == 1) {
	Fl::delete_widget (Maps->Maps_Edit_Window);
//      Maps->Maps_Edit_Window->~Fl_Window();         // this window needs to be deleted
//      delete(Maps->Maps_Edit_Window);               // not just killed (2d/3d nature might change)
//      delete(Maps->MapsBuffer);
	delete (Maps);
	Maps = NULL;
    }
    Update_Str (0);
    Generate_Drawing (0);
    Fl::redraw ();
}

void
Edit_Parmeters_cb (void)
{
// Callback routine to show the Edit Parameter screen and load the widgets on that page
    char string[40];

    if (!strlen (drvui->CurFile->value ())) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }
    if (!edtprm) {
	edtprm = new EditScreen;	// Create the screen
	if (!Omit->nomits) {
	    edtprm->ClearLastOmit->deactivate ();
	    edtprm->ClearOmit->deactivate ();
	}
#if !defined (WIN32) && !defined (__APPLE__)
	edtprm->editWindow->icon ((char *) drvui->icon);
#endif
	edtprm->show ();	// show it
    }
    sprintf (string, "%6.1f", printdist);	// load the bond limit distance widget
    edtprm->List->value (string);
    sprintf (string, "%5.2f", Magnification);	// Magnification widget
    edtprm->Magnification->value (string);
    if (drvui->noshadow)
	edtprm->NoShadow->set ();
    if (Labels == 0)
	edtprm->NoLabels->set ();
    if (domolcomp) {
	sprintf (string, "%5.2f", drvui->mol_d);
	edtprm->Mol_Comp_Dist->value (string);	// Molecule completion stuff
	edtprm->MolCompButton->set ();	// set check button
    }
    sprintf (string, "%6.3f", DepthCue);
    edtprm->DepthCue->value (string);	// Depth Cue widget
    Load_Color_Combo (edtprm->Cell_Edge_Color);
    edtprm->Cell_Edge_Color->value (drvui->col_cell);	// load cell axis line color
    sprintf (string, "%6.3f", rad_cell);
    edtprm->Cell_Edge_Width->value (string);	//   and width
    sprintf (string, "%6.2f", drvui->polylimit);	// load polyhedral limit
    edtprm->Poly_Limit->value (string);
    sprintf (string, "%7.2f", drvui->Phong_Value);	// load phong value
    edtprm->Phong_Refl->value (string);
    sprintf (string, "%7.0f", drvui->Phong_Size);	//  and size
    edtprm->Phong_Size->value (string);
    sprintf (string, "%6.2f", drvui->ambient);	// load default POV ambient,
    edtprm->Ambient_Finish->value (string);
    sprintf (string, "%6.2f", drvui->diffuse);	// diffuse,
    edtprm->Diffuse_Finish->value (string);
    sprintf (string, "%6.2f", drvui->specular);	// specular finish and
    edtprm->Specular_Finish->value (string);
    sprintf (string, "%6.2f", drvui->roughness);	// surface roughness
    edtprm->Finish_Roughness->value (string);
    strcpy (string, drvui->col_bg);	// load background color
    Load_Color_Combo (edtprm->Background_Color);
    edtprm->Background_Color->value (string);
    sprintf (string, "%.3f", drvui->label_scale);
    edtprm->Label_Scale->value (string);
    edtprm->editWindow->show ();
}

void
Edit_Parmeters_Close_cb (Fl_Button *, void *)
{
// callback routine when the 'close' button is pushed
    edtprm->editWindow->hide ();	// hide the window
    Fl::redraw ();		// update the screen
}

void
Edit_Parmeters_Save_cb (Fl_Button *, int *tosave)
{
// callback routine when 'save' or 'apply' button is pressed on the Edit Parameter screen

    Magnification = (float) atof (edtprm->Magnification->value ());	// extract data from widgets
    printdist = (float) atof (edtprm->List->value ());
    DepthCue = (float) atof (edtprm->DepthCue->value ());
    drvui->label_scale = (float) atof (edtprm->Label_Scale->value ());
    strcpy (drvui->col_cell, edtprm->Cell_Edge_Color->value ());
    trim_string (drvui->col_cell, 39);
    if (!strlen (drvui->col_cell))
	strcpy (drvui->col_cell, "White");
    rad_cell = (float) atof (edtprm->Cell_Edge_Width->value ());
    drvui->polylimit = (float) atof (edtprm->Poly_Limit->value ());
    drvui->Phong_Value = (float) atof (edtprm->Phong_Refl->value ());
    drvui->Phong_Size = (float) atof (edtprm->Phong_Size->value ());
    drvui->ambient = (float) atof (edtprm->Ambient_Finish->value ());
    drvui->diffuse = (float) atof (edtprm->Diffuse_Finish->value ());
    drvui->specular = (float) atof (edtprm->Specular_Finish->value ());
    drvui->roughness = (float) atof (edtprm->Finish_Roughness->value ());
    strcpy (drvui->col_bg, edtprm->Background_Color->value ());
    trim_string (drvui->col_bg, 39);
    if (!strlen (drvui->col_bg))
	strcpy (drvui->col_bg, "White");
    drvui->noshadow = 0;
    if (edtprm->NoShadow->value ())
	drvui->noshadow = 1;
    Labels = 1;
    if (edtprm->NoLabels->value ()) {
	int i, j = 1, k, l;

	Labels = 0;
	for (i = 1; i < drvui->nlabel; i++) {
	    if (!strcmp (drvui->labels[i].label_label, "a"))
		drvui->labels[i].label_fn = 0;
	    if (!strcmp (drvui->labels[i].label_label, "b"))
		drvui->labels[i].label_fn = 0;
	    if (!strcmp (drvui->labels[i].label_label, "c"))
		drvui->labels[i].label_fn = 0;
	    if (!strcmp (drvui->labels[i].label_label, "o"))
		drvui->labels[i].label_fn = 0;
	    if (!strcmp (drvui->labels[i].label_label, "triple_vect"))
		drvui->labels[i].label_fn = 0;
	}
	for (i = 1; i < drvui->nlabel; i++) {
	    if (drvui->labels[i].label_fn > 0)
		j++;
	    if (j != i) {
		l = j;
		for (k = i + 1; k < drvui->nlabel; k++) {
		    drvui->labels[l].label_fn = drvui->labels[k].label_fn;
		    strcpy (drvui->labels[l].label_label, drvui->labels[k].label_label);
		    drvui->labels[l].label_x[0] = drvui->labels[k].label_x[0];
		    drvui->labels[l].label_x[1] = drvui->labels[k].label_x[1];
		    drvui->labels[l++].label_x[2] = drvui->labels[k].label_x[2];
		}
	    }
	}
	drvui->nlabel = j;
    } else {
	label_cell ();
    }
    if (edtprm->MolCompButton->value () != 0) {
	float temp;

	temp = (float) atof (edtprm->Mol_Comp_Dist->value ());
	if (temp > 3.5f) {
	    char string[100];

	    sprintf (string, "A molecular completion distance of %5.2f is very\n"
		     "large. Are you sure you wish to continue?", temp);
	    if (!fl_choice ("%s", "No", "Yes", NULL, string))
		return;
	}
	drvui->mol_d = temp;
	domolcomp = 1;
    } else {
	drvui->mol_d = 0.0;
	domolcomp = 0;
    }
    if (edtprm->ClearOmit->value () != 0) {
	Omit->nomits = 0;	// dump the omit list if ClearOmit button checked
	edtprm->ClearOmit->value (0);    // and disarm the checkbutton afterwards
	edtprm->ClearLastOmit->deactivate ();  
	edtprm->ClearOmit->deactivate ();
    }
    drvui->Str_File_Changed = 1;
    Update_Str (0);		// update the 'str' file
    Generate_Drawing (0);	// regenerate the drawing
    if (*tosave) {
	edtprm->editWindow->hide ();	// hide the window if 'save'
    }
    Fl::redraw ();		// update the screen
}

void
Edit_Polyhedra_cb (void)
{
// callback routine to load polyhedra edit screen
    char string[100];

    static int one = 1;

    static int zero = 0;

    int i;

    int y = 5;

    if (!strlen (drvui->CurFile->value ())) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }
/*
    if (!natom) {         // No atoms, no polyhedra...
        Error_Box("This structure file does not contain any atoms.");
        return;
    }
*/
    Save_Working_Copy ();
    if (!Polyhedra) {
	Polyhedra = new PolyParam;	// new instance of the polyhedra parameters
	Polyhedra->Polyhedra_Edit_Window =
	    new Fl_Window (50, 0, 560, 560, "Edit Polyhedral/Plane Parameters");
	Polyhedra->Polyhedra_Edit_Window->
	    callback ((Fl_Callback *) Edit_Polyhedra_Close_cb);
	if (drvui->max_frame > 1) {
	    Flu_Combo_List *o = Polyhedra->Frame_No =
		new Flu_Combo_List (220, y, 75, 25, "Frame No.");
	    o->align (FL_ALIGN_TOP);
	    o->callback (Polyhedra_Frame_Combo_cb);
	    o->labelfont (1);
	    for (i = 1; i <= drvui->max_frame; i++) {
		sprintf (string, "%d", i);
		o->list.add (string);
	    }
	    o->pop_height (20 * drvui->max_frame);
	    o->value ("1");
	}
	y += 50;
	Polyhedra->Polyhedra_Edit = new Fl_Text_Editor (25, y, 450, 95,
							"     Type    From   To   Min d   Max d    Color          ");
	Polyhedra->PolyhedraBuffer = new Fl_Text_Buffer;
	Polyhedra->PolyhedraBuffer->
	    add_modify_callback ((Fl_Text_Modify_Cb) Modify_Polyhedra_cb, (void *) NULL);
	Polyhedra->Polyhedra_Edit->textfont (FL_COURIER);
	Polyhedra->Polyhedra_Edit->textsize (12);
	Polyhedra->Polyhedra_Edit->buffer (Polyhedra->PolyhedraBuffer);
	Polyhedra->Polyhedra_Edit->labelfont (FL_COURIER_BOLD);
	y += 105;
	Polyhedra->PolyInstr = new Fl_Output (135, y, 250, 0, "Press 'Add' to replace "
					      "selected line - 'Remove' to delete it");
	Polyhedra->PolyInstr->hide ();
	Polyhedra->PolyInstr->align (FL_ALIGN_BOTTOM);
	Polyhedra->PolyInstr1 = new Fl_Output (135, y, 250, 0, "Highlight text above "
					       "or double click to edit line");
	Polyhedra->PolyInstr1->hide ();
	Polyhedra->PolyInstr1->align (FL_ALIGN_BOTTOM);
	y += 40;
	Polyhedra->New_Polyhedra_From = new Fl_Input (20, y, 50, 25, "From");
	Polyhedra->New_Polyhedra_From->align (FL_ALIGN_TOP);
	Polyhedra->New_Polyhedra_From->labelfont (1);
	Fl_Input *op = Polyhedra->New_Polyhedra_To = new Fl_Input (80, y, 50, 25, "To");

	op->align (FL_ALIGN_TOP);
	op->callback ((Fl_Callback *) New_Polyhedra_Input_cb);
	op->labelfont (1);
	Fl_Input *m = Polyhedra->New_Polyhedra_Min =
	    new Fl_Input (140, y, 50, 25, "Min d");
	m->align (FL_ALIGN_TOP);
	m->callback ((Fl_Callback *) New_Polyhedra_Input_cb);
	m->labelfont (1);
	Fl_Input *os = Polyhedra->New_Polyhedra_Max =
	    new Fl_Input (200, y, 50, 25, "Max d");
	os->align (FL_ALIGN_TOP);
	os->callback ((Fl_Callback *) New_Polyhedra_Input_cb);
	os->labelfont (1);
	Flu_Combo_List *ot = Polyhedra->New_Polyhedra_Color =
	    new Flu_Combo_List (260, y, 160, 25, "Color");
	ot->align (FL_ALIGN_TOP);
	ot->callback ((Fl_Callback *) New_Polyhedra_Input_cb);
	ot->labelfont (1);
	Load_Color_Combo (ot);
	Fl_Input *ou = Polyhedra->New_Polyhedra_Transp =
	    new Fl_Input (430, y, 50, 25, "Transp.");
	ou->callback ((Fl_Callback *) New_Polyhedra_Input_cb);
	ou->align (FL_ALIGN_TOP);
	ou->labelfont (1);
	y += 40;
	Polyhedra->Edge_Radius = new Fl_Input (150, y, 100, 25, "Edge Radius");
	Polyhedra->Edge_Radius->align (FL_ALIGN_TOP);
	Polyhedra->Edge_Radius->labelfont (1);
	Polyhedra->Edge_Color = new Flu_Combo_List (270, y, 120, 25, "Edge Color");
	Polyhedra->Edge_Color->align (FL_ALIGN_TOP);
	Polyhedra->Edge_Color->labelfont (1);
	Load_Color_Combo (Polyhedra->Edge_Color);
	Polyhedra->Edge_Color->value (drvui->col_edge);
	sprintf (string, "%6.3f", drvui->rad_edge);
	Polyhedra->Edge_Radius->value (string);
	y += 40;
	Fl_Button *om = Polyhedra->New_Polyhedra_Add =
	    new Fl_Button (180, y, 70, 25, "Add");
	om->callback ((Fl_Callback *) New_Polyhedra_Add_cb, &one);
	om->tooltip ("When active, press to transfer data in boxes to window above");
	om->deactivate ();
	Fl_Button *mm = Polyhedra->New_Polyhedra_Remove =
	    new Fl_Button (270, y, 70, 25, "Remove");
	mm->callback ((Fl_Callback *) New_Polyhedra_Add_cb, &zero);
	mm->tooltip ("When active, press to remove highlighted line.");
	mm->deactivate ();
	y += 30;
	Polyhedra->PolyInstr2 = new Fl_Output (25, y, 470, 0,
					       "Changes are temporary until \"Apply\" or \"Save\" is pressed.");
	Polyhedra->PolyInstr2->hide ();
	Polyhedra->PolyInstr2->align (FL_ALIGN_BOTTOM);
	y += 40;
	Fl_Group *qq = new Fl_Group (30, y, 150, 90, "Polyhedra/Plane Type");

	qq->labelfont (1);
	qq->box (FL_THIN_UP_BOX);
	Polyhedra->Polysz = new Fl_Radio_Button (40, y + 10, 12, 12, "Polysz (PS) Style");
	Polyhedra->Polysz->type (102);
	Polyhedra->Polysz->selection_color ((Fl_Color) 1);
	Polyhedra->Polysz->align (FL_ALIGN_RIGHT);
	Polyhedra->Polysz->set ();
	Polyhedra->Polyvert =
	    new Fl_Radio_Button (40, y + 30, 12, 12, "Polyvert (PV) Style");
	Polyhedra->Polyvert->type (102);
	Polyhedra->Polyvert->selection_color ((Fl_Color) 1);
	Polyhedra->Polyvert->align (FL_ALIGN_RIGHT);
	Polyhedra->Polyshell =
	    new Fl_Radio_Button (40, y + 50, 12, 12, "Shell (SH) Style");
	Polyhedra->Polyshell->type (102);
	Polyhedra->Polyshell->selection_color ((Fl_Color) 1);
	Polyhedra->Polyshell->align (FL_ALIGN_RIGHT);
	Polyhedra->Plane = new Fl_Radio_Button (40, y + 70, 12, 12, "Plane (PL) Style");
	Polyhedra->Plane->type (102);
	Polyhedra->Plane->selection_color ((Fl_Color) 1);
	Polyhedra->Plane->align (FL_ALIGN_RIGHT);
	qq->end ();
	Flu_Combo_List *o = Polyhedra->Polyhedra_Combo =
	    new Flu_Combo_List (200, y, 100, 25, "'From' Atom");
	o->align (FL_ALIGN_TOP);
	o->labelfont (1);
	o->callback (Polyhedra_Combo_cb);
	Fl_Text_Editor *p =
	    new Fl_Text_Editor (310, y, 200, 100, "'To' Atoms & Distances");
	p->textfont (FL_COURIER);
	p->textsize (14);
	p->align (FL_ALIGN_TOP);
	p->labelfont (1);
	Polyhedra->Polyhedra_Output_Buffer = new Fl_Text_Buffer;
	Polyhedra->Polyhedra_Output_Buffer->
	    add_modify_callback ((Fl_Text_Modify_Cb) Modify_Polyhedra_Distance_cb,
				 (void *) NULL);
	p->buffer (Polyhedra->Polyhedra_Output_Buffer);
	y += 120;
	Polyhedra->Def_Edge_Radius = new Fl_Input (145, y, 100, 25, "Default Edge Radius");
	Polyhedra->Def_Edge_Radius->align (FL_ALIGN_TOP);
	Polyhedra->Def_Edge_Radius->labelfont (1);
	Polyhedra->Def_Edge_Color = new Flu_Combo_List (275, y, 120, 25, "Default Edge Color");
	Polyhedra->Def_Edge_Color->align (FL_ALIGN_TOP);
	Polyhedra->Def_Edge_Color->labelfont (1);
	Load_Color_Combo (Polyhedra->Def_Edge_Color);
	Polyhedra->Def_Edge_Color->value (drvui->col_edge);
	sprintf (string, "%6.3f", drvui->rad_edge);
	Polyhedra->Def_Edge_Radius->value (string);
#if !defined (WIN32) && !defined (__APPLE__)
	Polyhedra->Polyhedra_Edit_Window->icon ((char *) drvui->icon);
#endif
	y += 35;
	Fl_Button *r = new Fl_Button (125, y, 70, 25, "Close");

	r->tooltip ("Close this window and discard all changes.");
	r->callback ((Fl_Callback *) Edit_Polyhedra_Close_cb);
	Fl_Button *a = new Fl_Button (225, y, 70, 25, "Apply");

	a->tooltip
	    ("Apply current contents of top box to drawing, but leave this window open.");
	a->callback ((Fl_Callback *) Edit_Polyhedra_Save_cb, &zero);
	Fl_Button *s = new Fl_Button (325, y, 70, 25, "Save");

	s->tooltip
	    ("Apply current contents of top box to drawing, then close this window.");
	s->callback ((Fl_Callback *) Edit_Polyhedra_Save_cb, &one);
	Polyhedra->Polyhedra_Edit_Window->end ();
    }
    Polyhedra_Frame_Combo_cb (NULL, NULL);
    Polyhedra->Polyhedra_Edit_Window->show ();
}

void
Load_Color_Combo (Flu_Combo_List * ot)
{
    int j, k = strlen (Colors_Combo), l = 0;
	
    char string[20];

    FILE *colinc;

    if (ot->list.size () > 1)
	return;

    // look for private colors.inc file first 
    colinc = NULL;
    if (strlen (drvui->POV_Include) > 10) {
	colinc = fopen (drvui->POV_Include, "r");
    }
    if (!colinc) {
	char* a, b[30];
	char POV_incpath[255] = "\0";
	strncpy (b, drvui->POV_Options, 29);
	b[29] = '\0';
	a = strstr (b, "+L");
	if (a != NULL) {
	    sscanf (a + 2, "%s", POV_incpath);
	    if (POV_incpath != NULL) {
		strcat (POV_incpath, "/colors.inc");
		colinc = fopen (POV_incpath, "r");
		if (colinc) {
		    strcpy (drvui->POV_Include, POV_incpath);
		    WriteConfig ();
		}
	    }
	}
    }
    if (colinc != NULL) { // load color names from private colors.inc file 
	char line[80];
	char keyword[20], thecolor[20];
	float thered,thegreen,theblue;

	while (!feof (colinc)) {
	    if (fgets (line, 80, colinc) == NULL)
		continue;
	    if (!strstr (line, "#declare") || strstr(line,"Colors_Inc"))
		continue;
	    sscanf (line, "%s %s = color red %f green %f blue %f", keyword, thecolor,
		    &thered, &thegreen, &theblue);
	    ot->list.add(thecolor);
	}
	fclose (colinc);
    } else {  // colors.inc not given, use built-in color list corresponding to povray 3.6 

	for (j = 0; j < k; j++) {
	    if (Colors_Combo[j] != '\n') {
		string[l++] = Colors_Combo[j];
	    } else {
		string[l++] = 0;
		l = 0;
		ot->list.add (string);
	    }
	}
    }
    ot->pop_height (300);
}



void
Edit_Modparms_Close_cb ()
{
    Modparms->Mods_Edit_Window->hide ();

    drvui->modulated = Modparms->saved_avg;
    if (drvui->modulated == -1)
	Modparms->Mod_average->set ();
    else
	Modparms->Mod_average->clear ();
    for (int i = 0; i < 3; i++)
	drvui->phaseshift[i] = Modparms->saved_t[i];
    Modparms->Mod_t0->value (drvui->phaseshift[0]);
    Modparms->Mod_t1->value (drvui->phaseshift[1]);
    Modparms->Mod_t2->value (drvui->phaseshift[2]);
    if (Modparms->Mod_average->value () != 0) {
	Modparms->Mod_t0->deactivate ();
	Modparms->Mod_t1->deactivate ();
	Modparms->Mod_t2->deactivate ();
    } else {
	Modparms->Mod_t0->activate ();
	if (drvui->no_cell_vec > 1)
	    Modparms->Mod_t1->activate ();
	if (drvui->no_cell_vec > 2)
	    Modparms->Mod_t2->activate ();
    }

    drvui->Str_File_Changed = 1;
    Update_Str (0);		// update the 'str' file
    Generate_Drawing (1);	// regenerate
    Fl::redraw ();		// update the screen
}

void
Edit_Modparms_Change_cb ()
{
    if (Modparms->Mod_average->value () != 0) {
	drvui->modulated = -1;
	Modparms->Mod_t0->deactivate ();
	Modparms->Mod_t1->deactivate ();
	Modparms->Mod_t2->deactivate ();
    } else {
	drvui->modulated = 1;
	Modparms->Mod_t0->activate ();
	if (drvui->no_cell_vec > 1)
	    Modparms->Mod_t1->activate ();
	if (drvui->no_cell_vec > 2)
	    Modparms->Mod_t2->activate ();
    }

    drvui->phaseshift[0] = (float) Modparms->Mod_t0->value ();
    drvui->phaseshift[1] = (float) Modparms->Mod_t1->value ();
    drvui->phaseshift[2] = (float) Modparms->Mod_t2->value ();

    drvui->Str_File_Changed = 1;
    Update_Str (0);		// update the 'str' file
    Generate_Drawing (0);
    Fl::redraw ();
}

void
Edit_Modparms_Save_cb (Fl_Button *, int *save)
{
    int i, j, k, n;

    float d1, d2;

    char atom1[5], numb[5], widget[16382];

    if (Modparms->Mod_average->value () != 0) {
	drvui->modulated = -1;
	Modparms->Mod_t0->deactivate ();
	Modparms->Mod_t1->deactivate ();
	Modparms->Mod_t2->deactivate ();
    } else {
	drvui->modulated = 1;
	Modparms->Mod_t0->activate ();
	if (drvui->no_cell_vec > 1)
	    Modparms->Mod_t1->activate ();
	if (drvui->no_cell_vec > 2)
	    Modparms->Mod_t2->activate ();
    }
    drvui->phaseshift[0] = (float) Modparms->Mod_t0->value ();
    drvui->phaseshift[1] = (float) Modparms->Mod_t1->value ();
    drvui->phaseshift[2] = (float) Modparms->Mod_t2->value ();
    Modparms->saved_avg = drvui->modulated;
    for (i = 0; i < 3; i++)
	Modparms->saved_t[i] = drvui->phaseshift[i];

    for (j = 0; j < natom; j++) {
	drvui->atoms[j].occupancy = 1.;
	drvui->atoms[j].min_occ = 0.;
    }
    char *selection = Modparms->Occ_Buffer->text ();

    strcpy (widget, selection);
    free (selection);
    if (strlen (widget) < 10) {
	strcpy (widget, "");
	Modparms->Occ_Buffer->text (widget);
    } else {
	while (strlen (widget) > 10) {
	    memset (atom1, 0, 5);
	    (void) sscanf (widget, "%4c %s %f %f", atom1, numb, &d1, &d2);
	    if (strstr (numb, "*"))
		n = -1;
	    else
		(void) sscanf (numb, "%d", &n);
	    for (j = 0; j < natom; j++) {
		if (check_atom_name (atom1, drvui->atoms[j].atom_l) &&
		    (n == -1 || n == drvui->atoms[j].sv_atom_n)) {
		    drvui->atoms[j].occupancy = d1;
		    drvui->atoms[j].min_occ = d2;
		    break;
		}
	    }
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
    Update_Str (0);		// update the 'str' file
    Generate_Drawing (0);	// regenerate the drawing

    if (*save)
	Modparms->Mods_Edit_Window->hide ();	// hide the window if 'save'
    Fl::redraw ();		// update the screen
}

void
Edit_Modparms_cb (void)
{
// Callback routine to show the Modulation Options screen and load the widgets on that page
    static int one = 1;

    static int zero = 0;

    char string[10], buffer[40];

    int i, y;

    if (!strlen (drvui->CurFile->value ())) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }
    if (drvui->modulated == 0) {	// Make sure we have a modulated structure
	Error_Box ("This Structure does not appear to be modulated.");
	return;
    }
    if (!Modparms) {
	Modparms = new ModParam;	// new instance of the modulation parameters
	Modparms->Mods_Edit_Window =
	    new Fl_Window (50, 50, 360, 500, "Edit Modulation Parameters");
	Modparms->Mods_Edit_Window->callback ((Fl_Callback *) Edit_Modparms_Close_cb);
	y = 30;
	Modparms->Mod_average =
	    new Fl_Check_Button (70, y, 225, 25, "Display average structure only");
	if (drvui->modulated == -1)
	    Modparms->Mod_average->set ();
	Modparms->Mod_average->callback ((Fl_Callback *) Edit_Modparms_Change_cb);
	y += 35;
	Modparms->Mod_t0 =
	    new Flu_Spinner (220, y, 100, 25, "Phase shift t (1st direction)");
	Modparms->Mod_t0->value (drvui->phaseshift[0]);
	Modparms->Mod_t0->callback ((Fl_Callback *) Edit_Modparms_Change_cb);
	y += 35;
	Modparms->Mod_t1 =
	    new Flu_Spinner (220, y, 100, 25, "Phase shift t (2nd direction)");
	Modparms->Mod_t1->value (drvui->phaseshift[1]);
	Modparms->Mod_t1->callback ((Fl_Callback *) Edit_Modparms_Change_cb);
	y += 35;
	Modparms->Mod_t2 =
	    new Flu_Spinner (220, y, 100, 25, "Phase shift t (3rd direction)");
	Modparms->Mod_t2->value (drvui->phaseshift[2]);
	Modparms->Mod_t2->callback ((Fl_Callback *) Edit_Modparms_Change_cb);

	if (Modparms->Mod_average->value () != 0) {
	    Modparms->Mod_t0->deactivate ();
	    Modparms->Mod_t1->deactivate ();
	    Modparms->Mod_t2->deactivate ();
	}
	if (drvui->no_cell_vec < 2)
	    Modparms->Mod_t1->deactivate ();
	if (drvui->no_cell_vec < 3)
	    Modparms->Mod_t2->deactivate ();
	y += 65;

	Modparms->Occ_Edit = new Fl_Text_Editor (25, y, 310, 80,
						 "Atom     Occupancy     Minimum");
	Modparms->Occ_Buffer = new Fl_Text_Buffer;
	for (i = 0; i < natom; i++) {
	    if (drvui->atoms[i].min_occ > 0.) {
		sprintf (buffer, "%4s%2d %15.3f %15.3f\n", drvui->atoms[i].atom_l,
			 drvui->atoms[i].sv_atom_n, drvui->atoms[i].occupancy,
			 drvui->atoms[i].min_occ);
		Modparms->Occ_Buffer->append (buffer);
	    }
	}
	Modparms->Occ_Buffer->add_modify_callback ((Fl_Text_Modify_Cb) Modify_Occ_cb,
						   (void *) NULL);
	y += 85;
	Modparms->Occ_Instr = new Fl_Output (0, y, 400, 0, "Press 'Add' to replace "
					     "line - 'Remove' to delete");
	Modparms->Occ_Instr->hide ();
	Modparms->Occ_Instr->align (FL_ALIGN_BOTTOM);
	Modparms->Occ_Instr1 = new Fl_Output (0, y, 400, 0,
					      "Highlight or double-click text above to edit");
	Modparms->Occ_Instr1->hide ();
	Modparms->Occ_Instr1->align (FL_ALIGN_BOTTOM);
	Modparms->Occ_Edit->textfont (FL_COURIER);
	Modparms->Occ_Edit->textsize (12);
	Modparms->Occ_Edit->buffer (Modparms->Occ_Buffer);
	Modparms->Occ_Edit->labelfont (FL_COURIER_BOLD);
	y += 40;
	Flu_Combo_List *o = Modparms->Occ_Combo =
	    new Flu_Combo_List (35, y, 100, 25, "Atom");
	o->align (FL_ALIGN_TOP);
	o->callback (Occ_Combo_cb);
	o->labelfont (1);
	for (i = 0; i < natom; i++)
	    if (drvui->atoms[i].occ_ismod != 0) {
		memset (string, 0, 10);
		sprintf (string, "%4s%2d", drvui->atoms[i].atom_l,
			 drvui->atoms[i].sv_atom_n);
		Modparms->Occ_Combo->list.add (string);
	    }
	Fl_Input *oa = Modparms->New_Occ_Avg = new Fl_Input (155, y, 60, 25, "Avg. occ.");

	oa->align (FL_ALIGN_TOP);
	oa->callback ((Fl_Callback *) New_Occ_Input_cb);
	oa->labelfont (1);
	Fl_Input *om = Modparms->New_Occ_Min = new Fl_Input (235, y, 60, 25, "Min. occ.");

	om->align (FL_ALIGN_TOP);
	om->callback ((Fl_Callback *) New_Occ_Input_cb);
	om->labelfont (1);
	y += 30;
	Fl_Button *oba = Modparms->New_Occ_Add = new Fl_Button (105, y, 70, 25, "Add");

	oba->callback ((Fl_Callback *) New_Occ_Add_cb, &one);
	oba->tooltip ("When active, press to transfer data in boxes to window above");
	oba->deactivate ();
	Fl_Button *obr = Modparms->New_Occ_Remove =
	    new Fl_Button (195, y, 70, 25, "Remove");
	obr->callback ((Fl_Callback *) New_Occ_Add_cb, &zero);
	obr->tooltip ("When active, press to remove highlighted line.");
	obr->deactivate ();

	y += 45;
	Fl_Button *q = new Fl_Button (125, y, 130, 25, "Movie generation");
	q->tooltip("Choose parameters for POVray movie generation from t1 sampling");
	q->callback((Fl_Callback *) Automation_Edit_cb);

	y += 45;
	Fl_Button *r = new Fl_Button (50, y, 70, 25, "Close");

	r->tooltip ("Close this window and discard all changes.");
	r->callback ((Fl_Callback *) Edit_Modparms_Close_cb);
	Fl_Button *a = new Fl_Button (150, y, 70, 25, "Apply");

	a->tooltip ("Apply current settings to drawing, but leave this window open.");
	a->callback ((Fl_Callback *) Edit_Modparms_Save_cb, &zero);
	Fl_Button *s = new Fl_Button (250, y, 70, 25, "Save");

	s->tooltip ("Apply current settings to drawing, then close this window.");
	s->callback ((Fl_Callback *) Edit_Modparms_Save_cb, &one);

#if !defined (WIN32) && !defined (__APPLE__)
	Modparms->Mods_Edit_Window->icon ((char *) drvui->icon);
#endif

	Modparms->Mods_Edit_Window->end ();
    }
    Modparms->saved_avg = drvui->modulated;
    for (i = 0; i < 3; i++)
	Modparms->saved_t[i] = drvui->phaseshift[i];
    Modparms->Mods_Edit_Window->show ();
}

void
Edit_Slice_cb (void)
{
// Callback routine to show the Map Slice screen and load the widgets on that page
    static int one = 1;
    static int zero = 0;
    int y = 30, x;

    if (!strlen (drvui->CurFile->value ())) {	// Make sure rendering enabled
	Error_Box ("A Structure File must be selected first.");
	return;
    }
    Save_Working_Copy ();
//    Maps->Maps_Edit_Window->hide();
    if (!Slice) {
	Slice = new SliceParam;	// new instance of the slice parameters
	Slice->Slice_Edit_Window =
	    new Fl_Window (100, 100, 400, 340, "Edit Map Slice Parameters");
#if !defined (WIN32) && !defined (__APPLE__)
	Slice->Slice_Edit_Window->icon ((char *) drvui->icon);
#endif
	if (drvui->max_frame > 1) {
	    int i;
	    char string[128];

	    Flu_Combo_List *o = Slice->Frame_No =
			        new Flu_Combo_List (65, y, 75, 25, "Frame No.");
	    o->align (FL_ALIGN_TOP);
	    o->callback (Slice_Frame_Combo_cb);
	    o->labelfont (0);
	    for (i = 1; i <= drvui->max_frame; i++) {
		sprintf (string, "%d", i);
		o->list.add (string);
	    }
	    o->pop_height (20 * drvui->max_frame);
	    o->value ("1");
	    x = 265;
	} else {
	    x = 165;
	}
	Fl_Text_Display *mh = new Fl_Text_Display (x, y, 70, 0,
						  "Contour Type for Slice");
	mh->labelfont(0);
	Fl_Input *mi = Slice->New_type =
	    new Fl_Input (x, y, 70, 25, "");
	mi->align (FL_ALIGN_BOTTOM);
	mi->tooltip("Choices for type:  1 - Contoured slice,"
		    "    2 - Solid slice from blue to red,"
		    "    3 - Solid slice from black to white");
	y += 70;
	Fl_Text_Display *m = new Fl_Text_Display (0, y, 400, 0,
						  "Point in Plane of Slice");
	m->labelfont(0);
	Fl_Input *ma = Slice->New_x =
	    new Fl_Input (65, y + 5, 70, 25, "x");
	ma->align (FL_ALIGN_BOTTOM);
	ma->tooltip ("Point (in fractional coordinates) through which slice passes");
	Fl_Input *mb = Slice->New_y =
	    new Fl_Input (165, y + 5, 70, 25, "y");
	mb->align (FL_ALIGN_BOTTOM);
	mb->tooltip ("Point (in fractional coordinates) through which slice passes");
	Fl_Input *mc = Slice->New_z =
	    new Fl_Input (265, y + 5, 70, 25, "z");
	mc->align (FL_ALIGN_BOTTOM);
	mc->tooltip ("Point (in fractional coordinates) through which slice passes");
	y += 75;
	Fl_Text_Display *md = new Fl_Text_Display (0, y, 400, 0,
						  "Normal to Plane of Slice");
	md->labelfont(0);
	Fl_Input *me = Slice->New_nx =
	    new Fl_Input (65, y + 5, 70, 25, "nx");
	me->align (FL_ALIGN_BOTTOM);
	Fl_Input *mf = Slice->New_ny =
	    new Fl_Input (165, y + 5, 70, 25, "ny");
	mf->align (FL_ALIGN_BOTTOM);
	Fl_Input *mg = Slice->New_nz =
	    new Fl_Input (265, y + 5, 70, 25, "nz");
	mg->align (FL_ALIGN_BOTTOM);

	y += 65;
	Fl_Check_Button *sc = Slice->Legend = new Fl_Check_Button (140,y,35,25, "Display Legend");
	sc->align (FL_ALIGN_RIGHT);

	y += 45;
	Fl_Button *r = new Fl_Button (65, y, 70, 25, "Close");

	r->tooltip ("Close this window and discard all changes.");
	r->callback ((Fl_Callback *) Edit_Slice_Close_cb);
	Fl_Button *a = new Fl_Button (165, y, 70, 25, "Apply");
	a->tooltip
	    ("Apply current contents of top box to drawing, but leave this window open.");
	a->callback ((Fl_Callback *) Edit_Slice_Save_cb, &zero);
	Fl_Button *s = new Fl_Button (265, y, 70, 25, "Save");
	s->tooltip
	    ("Apply current contents of top box to drawing, then close this window.");
	s->callback ((Fl_Callback *) Edit_Slice_Save_cb, &one);
	Slice->Slice_Edit_Window->end ();
    }
    Slice_Frame_Combo_cb(NULL, NULL);
    Slice->Slice_Edit_Window->show ();
}

void
Edit_Slice_Close_cb (void)
{
// callback routine when the 'close' button is pushed
    Slice->Slice_Edit_Window->hide ();	// hide the window
    Fl::redraw ();		// update the screen
}

void
Edit_Slice_Save_cb (Fl_Button *, int *save)
{
// callback routine when 'save' or 'apply' button is pressed on the Edit Slice screen
    int Frame_No = 1;

    if (drvui->max_frame > 1)
	Frame_No = atoi (Slice->Frame_No->value ());

    drvui->Str_File_Changed = 1;
    drvui->frames[Frame_No].mapslice[0] = (float) atof (Slice->New_x->value ());
    drvui->frames[Frame_No].mapslice[1] = (float) atof (Slice->New_y->value ());
    drvui->frames[Frame_No].mapslice[2] = (float) atof (Slice->New_z->value ());
    drvui->frames[Frame_No].mapnorm[0] = (float) atof (Slice->New_nx->value ());
    drvui->frames[Frame_No].mapnorm[1] = (float) atof (Slice->New_ny->value ());
    drvui->frames[Frame_No].mapnorm[2] = (float) atof (Slice->New_nz->value ());
    drvui->frames[Frame_No].slice = atoi (Slice->New_type->value ());
    ShowMapLegend = Slice->Legend->value();
    Update_Str (0);
    Generate_Drawing (0);
    Save_Working_Copy ();
    if (*save == 1) {
	Slice->Slice_Edit_Window->hide ();
    }
    Fl::redraw ();
}

void
Modify_Occ_cb (Fl_Widget *, void *)
{
    char atom[10], value[8], string[40];

    float avg, min;

    int start, end, n;

    const char *selection;

    if (!Modparms->Occ_Buffer->selected ()) {
	Modparms->Occ_Combo->value ("");
	Modparms->New_Occ_Avg->value ("");
	Modparms->New_Occ_Min->value ("");
	Modparms->New_Occ_Add->deactivate ();
	return;
    }
    memset (atom, 0, 10);
    Modparms->Occ_Buffer->selection_position (&start, &end);
    selection = Modparms->Occ_Buffer->line_text (start);
    Modparms->Occ_Instr1->hide ();
    Modparms->Occ_Instr->show ();
    if (strlen (selection) == 0)
	return;			// user clicked on empty line by mistake
    sscanf (selection, "%s %s %f %f", atom, string, &avg, &min);
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
	Modparms->Occ_Combo->value (string);
	sprintf (value, "%2.3f", avg);
	Modparms->New_Occ_Avg->value (value);
	sprintf (value, "%2.3f", min);
	Modparms->New_Occ_Min->value (value);
	Modparms->New_Occ_Add->activate ();
	Modparms->New_Occ_Remove->activate ();

	New_Occ_Input_cb (NULL, NULL);
    }
    return;
}

void
Occ_Combo_cb (Fl_Widget *, void *)
{
    Modparms->New_Occ_Avg->value ("1.0");
    Modparms->New_Occ_Add->activate ();
}

void
New_Occ_Input_cb (Fl_Widget *, void *)
{
    int n;

    char atom[5], number[5];

// callback routine to make 'Add' new occupancy limit button active whenever all 
//  required fields are non-blank
    if (atof (Modparms->New_Occ_Min->value ()) == 0.0)
	return;
    Modparms->New_Occ_Add->activate ();
    (void) sscanf (Modparms->Occ_Combo->value (), "%s %s", atom, number);
    if (strstr (number, "*")) {
	n = -1;
    } else {
	(void) sscanf (number, "%d", &n);
    }
    while (strlen (atom) < 4)
	strcat (atom, " ");
}

void
New_Occ_Add_cb (Fl_Widget *, int *action)
{
    char string[100];

    int start, end;

    float avg, min;

    char atom1[10];

    char *selection = NULL;

    avg = (float) atof (Modparms->New_Occ_Avg->value ());
    min = (float) atof (Modparms->New_Occ_Min->value ());

//    Modparms->OccInstr2->show();

    strcpy (atom1, Modparms->Occ_Combo->value ());
    if (Modparms->Occ_Buffer->selected ()) {
	Modparms->Occ_Buffer->selection_position (&start, &end);
	selection = Modparms->Occ_Buffer->line_text (start);
	if (*action != 1)
	    Modparms->Occ_Buffer->remove (Modparms->Occ_Buffer->line_start (start),
					  Modparms->Occ_Buffer->line_end (end) + 1);
    }
    if (*action == 1) {
	sprintf (string, "%6s %15.3f %15.3f\n", atom1, avg, min);
	if (selection)
	    Modparms->Occ_Buffer->replace (Modparms->Occ_Buffer->line_start (start),
					   Modparms->Occ_Buffer->line_end (end) + 1,
					   string);
	else
	    Modparms->Occ_Buffer->append (string);
    }
    free ((char *) selection);
    Modparms->New_Occ_Avg->value ("");
    Modparms->New_Occ_Min->value ("");
    Modparms->Occ_Combo->value ("");
    Modparms->New_Occ_Add->deactivate ();
    Modparms->New_Occ_Remove->deactivate ();
    Modparms->Occ_Instr->hide ();
    Modparms->Occ_Instr1->show ();
}

void
Maps_Frame_Combo_cb (Fl_Widget *, void *)
{
// update fourier map data when the frame number in the combo box is changed
    char string[10];

    int Frame_No = 1;

    if (drvui->max_frame > 1)
	Frame_No = atoi (Maps->Frame_No->value ());

	sprintf (string, "%.3f", drvui->frames[Frame_No].map_lim[0]);
	Maps->XMin->value (string);
	sprintf (string, "%.3f", drvui->frames[Frame_No].map_lim[3]);
	Maps->XMax->value (string);
	sprintf (string, "%.3f", drvui->frames[Frame_No].map_lim[1]);
	Maps->YMin->value (string);
	sprintf (string, "%.3f", drvui->frames[Frame_No].map_lim[4]);
	Maps->YMax->value (string);
	sprintf (string, "%.3f", drvui->frames[Frame_No].map_lim[2]);
	Maps->ZMin->value (string);
	sprintf (string, "%.3f", drvui->frames[Frame_No].map_lim[5]);
	Maps->ZMax->value (string);
}

void
Slice_Frame_Combo_cb (Fl_Widget *, void *)
{
// routine called when the frame number is changed on the Slice Edit Screen
    char string[10];
    int Frame_No = 1;

    if (drvui->max_frame > 1)
	Frame_No = atoi (Slice->Frame_No->value ());
    sprintf(string, "%.5f", drvui->frames[Frame_No].mapslice[0]);
    Slice->New_x->value(string);
    sprintf(string, "%.5f", drvui->frames[Frame_No].mapslice[1]);
    Slice->New_y->value(string);
    sprintf(string, "%.5f", drvui->frames[Frame_No].mapslice[2]);
    Slice->New_z->value(string);
    sprintf(string, "%.5f", drvui->frames[Frame_No].mapnorm[0]);
    Slice->New_nx->value(string);
    sprintf(string, "%.5f", drvui->frames[Frame_No].mapnorm[1]);
    Slice->New_ny->value(string);
    sprintf(string, "%.5f", drvui->frames[Frame_No].mapnorm[2]);
    Slice->New_nz->value(string);
    sprintf(string, "%i", drvui->frames[Frame_No].slice);
    Slice->New_type->value(string);
    Slice->Legend->value(ShowMapLegend);
}

