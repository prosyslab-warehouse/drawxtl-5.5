// $Id: Help.cxx 1107 2011-01-19 23:53:52Z martin $
//
// Help.cxx - routine for DRAWxtl V5.5 - the GUI version
//
// Coded using the FLTK 1.1.6 widget set
//
//     Larry W. Finger, Martin Kroeker and Brian Toby
//
// This module includes the help screens for the GUI
//
// routines contained within this file:
//
//  About_Help_cb - callback routine from About menu button
//  Color_Help_cb - callback routine to display table of RGB values
//  Graphics_Help_cb - callback routine to display graphics help
//  Input_Help_cb - callback routine for help on input file
//  Spacegroup_Help_cb - callback routine for help on spacegroups with alternate origins
//  View_Help_Close_cb - callback routine to hide help windows

#include "CrystalView.h"
#include "DRAWxtlViewUI.h"
#include "draw_ext.h"

#include "DRAWxtl_proto.h"

void
About_Help_cb (void)
{
    static int two = 2;

    char string[256];

    static char outstring[256];

    char tstring[30];

    char vstring[30];

// callback routine for "About Help" menu button
    if (helpwindow2) {
	helpwindow2->show ();
	return;
    }
    helpwindow2 = new Fl_Window (100, 100, 500, 430, "About DRAWxtl");
    helpwindow2->begin ();
    helpwindow2->callback ((Fl_Callback *) View_Help_Close_cb, &two);

    int y = 40;

    Fl_Text_Display *box1 = new Fl_Text_Display (20, y, 460, 0, "DRAWxtl V5.5");

    box1->box (FL_NO_BOX);
    box1->labelsize (24);
    box1->labelcolor ((Fl_Color) 1);

    y += 30;
    Fl_Text_Display *box2 =
	new Fl_Text_Display (20, y, 460, 0,
			     "Copyright C 1996-2011 by Larry Finger, Martin Kroeker and Brian Toby");
    box2->box (FL_NO_BOX);
    box2->labelsize (14);
    box2->labelcolor ((Fl_Color) 186);

    y += 20;
    strcpy (string, "$Id: Help.cxx 1107 2011-01-19 23:53:52Z martin $");
    sscanf (string, "%*s %*s %s %s", vstring, tstring);
    sprintf (outstring, "Modified  %s, Build %s", tstring, vstring);
    Fl_Text_Display *box3 = new Fl_Text_Display (20, y, 460, 20, outstring);

    box3->box (FL_NO_BOX);
    box3->labelsize (14);
    box3->labelcolor ((Fl_Color) 186);

    y += 30;
    Fl_Text_Display *box4 =
	new Fl_Text_Display (20, y, 460, 20,
			     "DRAWxtl uses the FLTK 1 Widget" " set by ");
    box4->box (FL_NO_BOX);
    box4->labelsize (14);
    box4->labelcolor ((Fl_Color) 186);

    y += 20;
    Fl_Text_Display *box5 =
	new Fl_Text_Display (20, y, 460, 20,
			     "Bill Spitzak and others, http://www.fltk.org");
    box5->box (FL_NO_BOX);
    box5->labelsize (14);
    box5->labelcolor ((Fl_Color) 186);

    y += 30;
    Fl_Text_Display *box7 =
	new Fl_Text_Display (20, y, 460, 20,
			     "Some classes from FLU 2.14 by Jason Bryan are also employed.");
    box7->box (FL_NO_BOX);
    box7->labelsize (14);
    box7->labelcolor ((Fl_Color) 186);

    y += 20;
    Fl_Text_Display *box8 =
	new Fl_Text_Display (20, y, 460, 20, "http://www.osc.edu/~jbryan/FLU");
    box8->box (FL_NO_BOX);
    box8->labelsize (14);
    box8->labelcolor ((Fl_Color) 186);

    y += 30;
    Fl_Text_Display *box9 =
	new Fl_Text_Display (20, y, 460, 20,
			     "Virtual Trackball by Gavin Bell, Thant Tessman, and Tom Holroyd");
    box9->box (FL_NO_BOX);
    box9->labelsize (14);
    box9->labelcolor ((Fl_Color) 186);

    y += 20;
    Fl_Text_Display *box10 =
	new Fl_Text_Display (20, y, 460, 20,
			     "with ideas from the August 1988 issue  of SigGraph's");
    box10->box (FL_NO_BOX);
    box10->labelsize (14);
    box10->labelcolor ((Fl_Color) 186);

    y += 20;
    Fl_Text_Display *box11 =
	new Fl_Text_Display (20, y, 460, 20, "\"Computer Graphics,\" pp. 121-129.");
    box11->box (FL_NO_BOX);
    box11->labelsize (14);
    box11->labelcolor ((Fl_Color) 186);

    y += 30;
    Fl_Text_Display *box12 =
	new Fl_Text_Display (20, y, 460, 20,
			     "Direct Postscript output uses GL2PS by Christophe Geuzaine, ");
    box12->box (FL_NO_BOX);
    box12->labelsize (14);
    box12->labelcolor ((Fl_Color) 186);

    y += 20;
    Fl_Text_Display *box13 =
	new Fl_Text_Display (20, y, 460, 20, "http://www.geuz.org/gl2ps");
    box13->box (FL_NO_BOX);
    box13->labelsize (14);
    box13->labelcolor ((Fl_Color) 186);

    y += 30;
    Fl_Text_Display *box14 =
	new Fl_Text_Display (20, y, 460, 20,
			     "Marching Cubes implementation" " of Fourier contours from");
    box14->box (FL_NO_BOX);
    box14->labelsize (14);
    box14->labelcolor ((Fl_Color) 186);

    y += 20;
    Fl_Text_Display *box15 = new Fl_Text_Display (20, y, 460, 20, "Michael Y. Polyakov, "
						  "http://www.angelfire.com/linux/myp");
    box15->box (FL_NO_BOX);
    box15->labelsize (14);
    box15->labelcolor ((Fl_Color) 186);

    y += 20;
    Fl_Button *o = new Fl_Button (210, y, 80, 30, "Close");

    o->callback ((Fl_Callback *) View_Help_Close_cb, &two);
#if !defined (WIN32) && !defined (__APPLE__)
    helpwindow2->icon ((char *) drvui->icon);
#endif
    helpwindow2->end ();
    helpwindow2->show ();
}

void
Color_Help_cb (void)
{
    static int one = 1;

// Callback routine to show RBG values
    const char *colortext = {
	"This default color table is taken from the POV file 'colors.inc'\n"
	    "and is sorted by RGB Values - Red first, etc.. If you wish to use\n"
	    "custom colors, edit the above file and use the 'Configure' screen\n"
	    "to include the path to the file."
	    "\n"
	    "\n"
	    "Color Name        R             G               B\n"
	    "\n"
	    "Black		  0		0		0\n"
	    "Mica		  0		0		0\n"
	    "NewMidnightBlue	  0		0		0.61\n"
	    "Blue		  0		0		1\n"
	    "SlateBlue	  0		0.498039	1\n"
	    "Green		  0		1		0\n"
	    "SpringGreen	  0		1		0.498039\n"
	    "Cyan		  0		1		1\n"
	    "Gray05		  0.05		0.05		0.05\n"
	    "Gray10		  0.1		0.1		0.1\n"
	    "HuntersGreen	  0.13		0.37		0.31\n"
	    "Navy		  0.137255	0.137255	0.556863\n"
	    "NavyBlue	  0.137255	0.137255	0.556863\n"
	    "SteelBlue	  0.137255	0.419608	0.556863\n"
	    "ForestGreen	  0.137255	0.556863	0.137255\n"
	    "SeaGreen	  0.137255	0.556863	0.419608\n"
	    "Gray15		  0.15		0.15		0.15\n"
	    "MidnightBlue	  0.184314	0.184314	0.309804\n"
	    "DarkGreen	  0.184314	0.309804	0.184314\n"
	    "DarkSlateGray	  0.184314	0.309804	0.309804\n"
	    "DarkSlateGrey	  0.184314	0.309804	0.309804\n"
	    "MediumBlue	  0.196078	0.196078	0.8\n"
	    "SkyBlue		  0.196078	0.6		0.8\n"
	    "LimeGreen	  0.196078	0.8		0.196078\n"
	    "MediumAquamarine  0.196078	0.8		0.6\n"
	    "Gray20		  0.2		0.2		0.2\n"
	    "SummerSky	  0.22		0.69		0.87\n"
	    "Gray25		  0.25		0.25		0.25\n"
	    "CornflowerBlue	  0.258824	0.258824	0.435294\n"
	    "MediumSeaGreen	  0.258824	0.435294	0.258824\n"
	    "DkGreenCopper	  0.29		0.46		0.43\n"
	    "Gray30		  0.3		0.3		0.3\n"
	    "NeonBlue	  0.3		0.3		1\n"
	    "IndianRed	  0.309804	0.184314	0.184314\n"
	    "Violet		  0.309804	0.184314	0.309804\n"
	    "DarkOliveGreen	  0.309804	0.309804	0.184314\n"
	    "GreenCopper	  0.32		0.49		0.46\n"
	    "DimGray		  0.329412	0.329412	0.329412\n"
	    "DimGrey		  0.329412	0.329412	0.329412\n"
	    "VeryDarkBrown	  0.35		0.16		0.14\n"
	    "Gray35		  0.35		0.35		0.35\n"
	    "RichBlue          0.35          0.35            0.67\n"
	    "BakersChoc	  0.36		0.2		0.09\n"
	    "DarkBrown	  0.36		0.25		0.2\n"
	    "CadetBlue	  0.372549	0.623529	0.623529\n"
	    "Gray40		  0.4		0.4		0.4\n"
	    "DarkSlateBlue	  0.119608	0.137255	0.556863\n"
	    "MediumForestGreen 0.419608	0.556863	0.137255\n"
	    "SemiSweetChoc	  0.42		0.26		0.15\n"
	    "Salmon		  0.435294	0.258824	0.258824\n"
	    "DarkTurquoise	  0.439216	0.576471	0.858824\n"
	    "Aquamarine	  0.439216	0.858824	0.576471\n"
	    "MediumTurquoise	  0.439216	0.858824	0.858824\n"
	    "Gray45		  0.45		0.45		0.45\n"
	    "MediumSlateBlue	  0.498039	0		1\n"
	    "MediumSpringGreen 0.498039	1		0		\n"
	    "Gray50		  0.5		0.5		0.5\n"
	    "DarkWood	  0.52		0.37		0.26\n"
	    "DustyRose	  0.52		0.39		0.39\n"
	    "DarkPurple	  0.53		0.12		0.47\n"
	    "Scarlet		  0.55		0.09		0.09\n"
	    "Bronze		  0.55		0.47		0.14\n"
	    "Gray55		  0.55		0.55		0.55\n"
	    "Firebrick	  0.556863	0.137255	0.137255\n"
	    "Maroon		  0.556863	0.137255	0.419608\n"
	    "Sienna		  0.556863	0.419608	0.137255\n"
	    "LightSteelBlue	  0.560784	0.560784	0.737255\n"
	    "PaleGreen	  0.560784	0.737255	0.560784\n"
	    "MediumOrchid	  0.576471	0.439216	0.858824\n"
	    "GreenYellow	  0.576471	0.858824	0.439216\n"
	    "DarkTan		  0.59		0.41		0.31\n"
	    "DarkOrchid	  0.6		0.196078	0.8\n"
	    "Gray60		  0.6		0.6		0.6\n"
	    "YellowGreen	  0.6		0.8		0.196078\n"
	    "BlueViolet	  0.62352	0.372549	0.623529\n"
	    "\n"
	    "Khaki		  0.623529	0.623529	0.372549\n"
	    "Brown		  0.647059	0.164706	0.164706\n"
	    "Bronze2		  0.65		0.49		0.24\n"
	    "MediumWood	  0.65		0.5		0.39\n"
	    "Gray65		  0.65		0.65		0.65\n"
	    "LightGray	  0.658824	0.658824	0.658824\n"
	    "LightGrey	  0.658824	0.658824	0.658824\n"
	    "Turquoise	  0.678431	0.917647	0.917647\n"
	    "Gray70		  0.7		0.7		0.7\n"
	    "Brass		  0.71		0.65		0.26\n"
	    "Copper		  0.72		0.45		0.2\n"
	    "Pink		  0.737255	0.560784	0.560784\n"
	    "LightBlue	  0.74902	0.847059	0.847059\n"
	    "Gray75		  0.75		0.75		0.75\n"
	    "Gray		  0.752941	0.752941	0.752941\n"
	    "Grey		  0.752941	0.752941	0.752941\n"
	    "VioletRed	  0.8		0.196078	0.6\n"
	    "Gold		  0.8		0.498039	0.196078\n"
	    "Gray80		  0.8		0.8		0.8\n"
	    "VLightGrey	  0.8		0.8		0.8\n"
	    "OldGold		  0.81		0.71		0.23\n"
	    "Feldspar	  0.82		0.57		0.46\n"
	    "Thistle		  0.847059	0.74902		0.847059\n"
	    "Wheat		  0.847059	0.847059	0.74902\n"
	    "CoolCopper	  0.85		0.53		0.1\n"
	    "BrightGold	  0.85		0.85		0.1\n"
	    "Gray85		  0.85		0.85		0.85\n"
	    "Quartz		  0.85		0.85		0.95\n"
	    "MediumVioletRed	  0.858824	0.439216	0.576471\n"
	    "Orchid		  0.858824	0.439216	0.858824\n"
	    "Tan		  0.858824	0.576471	0.439216\n"
	    "Goldenrod	  0.858824	0.858824	0.439216\n"
	    "MandarinOrange	  0.89		0.47		0.2\n"
	    "Gray90		  0.9		0.9		0.9\n"
	    "Silver		  0.9		0.91		0.98\n"
	    "LightWood	  0.91		0.76		0.65\n"
	    "Plum		  0.917647	0.678431	0.917647\n"
	    "MediumGoldenrod	  0.917647	0.917647	0.678431\n"
	    "NewTan		  0.92		0.78		0.62\n"
	    "Gray95		  0.95		0.95		0.95\n"
	    "Flesh		  0.96		0.8		0.69\n"
	    "Red		  1		0		0\n"
	    "Magenta		  1		0		1\n"
	    "SpicyPink	  1		0.11		0.68\n"
	    "NeonPink	  1		0.43		0.78\n"
	    "Coral		  1		0.498039	0\n"
	    "OrangeRed	  1		0.498039	0\n"
	    "Orange		  1		0.5		0\n"
	    "Yellow		  1		1		0\n"
	    "Clear		  1		1		1\n"
	    "White		  1		1		1"
    };

    if (helpwindow1) {
	helpwindow1->show ();
	return;
    }
    helpwindow1 = new Fl_Window (50, 50, 680, 550, "Color RGB Values");
    helpwindow1->resizable (helpwindow1);
    helpwindow1->begin ();
    helpwindow1->callback ((Fl_Callback *) View_Help_Close_cb, &one);
    Fl_Text_Editor *display = new Fl_Text_Editor (0, 0, 680, 510);

    display->textfont (FL_COURIER);
    helpbuf1 = new Fl_Text_Buffer;
    display->buffer (helpbuf1);
    helpbuf1->text (colortext);
    Fl_Button *o = new Fl_Button (300, 515, 80, 30, "Close");

    o->callback ((Fl_Callback *) View_Help_Close_cb, &one);
#if !defined (WIN32) && !defined (__APPLE__)
    helpwindow1->icon ((char *) drvui->icon);
#endif
    helpwindow1->end ();
    helpwindow1->show ();
}

void
Graphics_Help_cb (void)
{
    static int four = 4;

// routine to display graphics help
    if (helpwindow4) {
	helpwindow4->show ();
	return;
    }
    const char Text[] =
	{ "\n                                   Graphics Keyboard Shortcuts and Help\n\n"
"    'C' or 'c' - turn the graphics cursor on. Each successive press reduces the size of the\n"
"                   steps. When the size is reduced below 0.01 A, the cursor is turned off.\n"
"                   When atoms are selected as described below, the various distances, angles,\n"
"                   and the torsion angle 1-2-3-4 are shown in the display area.\n\n"
"    'x', 'y', 'z' - move the cursor in the positive direction parallel to the x-, y- or z-axis.\n\n"
"    'X', 'Y', 'Z' - move the cursor in the negative direction parallel to the x-, y- or z-axis.\n\n"
"    'P' or 'p' - place the graphics cursor on the atom nearest the mouse position.\n\n"
"    'A' or 'a' - place the cursor at the position of the atom nearest the cursor.\n\n"
"    'M' or 'm' - move the cursor to the min (M) or max (m) in the electron-density.\n\n"
"    'S' or 's' - draw the electron density in a plane through the last three atoms.\n\n"
"    'L' or 'l' - label the atom at the cursor position.\n\n"
"    'B' or 'b' - Label the bond distance between atoms 1 and 2.\n\n"
"    Left mouse and drag - rotate the graphics object using a virtual trackball.\n\n"
#if !defined(__APPLE__)
	    "    Right mouse and drag - zoom in/out.\n\n"
	    "    Middle mouse (both on 2-button mouse) and drag - pan motion.\n"
#else
	    "    Command (Apple) key and drag - zoom in/out.\n\n"
	    "    Alt (option) key and drag - pan the graphics view.\n"
#endif
	"                       The arrow keys may also be used to move the object.\n\n"
	    "    HOME key - remove all zoom and pan motions.\n\n"
	    "    Shift/leftclick - drag labels or the triple vector to a desired position. If the triple\n"
	    "                        vector is dragged, all 3 of its labels are moved with it.\n\n"
	    "    Ctrl/leftclick - remove the object at the position of the mouse.\n"
    };
    helpwindow4 = new Fl_Window (200, 100, 600, 650, "DRAWxtl V5.5 Screen Graphics Help");
    helpwindow4->resizable (helpwindow4);
    helpwindow4->begin ();
    helpwindow4->callback ((Fl_Callback *) View_Help_Close_cb, &four);
    int y = 40;

    Fl_Multiline_Output *a = new Fl_Multiline_Output (0, 0, 600, 580);

    a->textsize (13);
    a->box (FL_FLAT_BOX);
    a->textfont (FL_HELVETICA_BOLD);
    a->value (Text);
    a->color (FL_WHITE);
    y = 600;
    Fl_Button *o = new Fl_Button (260, y, 80, 30, "Close");

    o->callback ((Fl_Callback *) View_Help_Close_cb, &four);
    o->take_focus ();
#if !defined (WIN32) && !defined (__APPLE__)
    helpwindow4->icon ((char *) drvui->icon);
#endif
    helpwindow4->end ();
    helpwindow4->show ();
}

void
Input_Help_cb (void)
{
// callback routine for help on input - from menu item
    const char *helptext = {
	"Input Instructions for DRAWxtl:\n"
	    "\n"
	    "In this program, color is represented in symbolic form, and must be one of the color names\n"
	    "from POV's 'colors.inc' file. The color names, sorted by RGB values, are available in the colors\n"
	    "entry under the Help menu. Any of these colors may be made transparent by appending the phrase\n"
	    "'filter xx' after it, where xx is a number from 0.0 to 1.0. (The larger the value of xx, the more\n"
	    "transparent will be the entity.\n"
	    "\n"
	    "Each line of the input file is preceded by a character sequence that describes the type of\n"
	    "information, as follows:\n"
	    "\n"
	    "\n"
	    " aimsurf name number filename style color\n"
	    "\n"
	    "causes the program to read from the given file a precalculated surface mesh to display at the\n"
	    "position of the specified atom. The file must be in the format used by the aim program (part of\n"
	    "the WIEN2k program suite) for calculating Bader surfaces of atoms according to the AIM concept.\n"
	    "(The calculated surface should cover the whole range of 0 to pi in theta and 0 to 2*pi in phi,\n"
	    "as no symmetry expansion is performed.) The rendering style can be 'dots', 'mesh' or 'solid'.\n"
	    "\n"
	    "\n"
	    " arrow       xp yp zp xc yc zc length diameter color\n"
	    "\n"
	    "defines the position in fractional coordinates (xp, yp, zp) of the triclinic, nuclear cell, the\n"
	    "components (xc, yc, zc) of the spin vector, and the length diameter and color of the arrow.. The\n"
	    "reference direction for xc is parallel to direct space axis a, yc is parallel (a x b) x a,\n"
	    "and the reference direction for zc is perpendicular to xc and yc. The only space-group symmetry\n"
	    "elements used in placing arrows are the translations described by the mag_trans command below.\n"
	    "\n"
	    "\n"
	    " atom           name number x y z\n"
	    "\n"
	    "defines the atoms. The 1- to 4-character name will be used on the commands that describe the\n"
	    "objects to be created, the number is to identify which atom of this type, and x y z are the\n"
	    "fractional coordinates in the unit cell. If the position is more easily represented as a fraction\n"
	    "such as 1/3, 1/4, 5/6, etc., it may be given in this form.\n"
	    "\n"
	    "\n"
	    "average\n"
	    "\n"
	    "\n"
	    "causes the program to draw the average structure of an incommensurately modulated crystal\n"
	    "even if information about positional or occupancy modulation is available in its CIF file.\n"
	    "\n"
	    "\n"
	    "axislines    width color\n"
	    "\n"
	    "defines the width and color of the lines that depict the principal axes of ellipsoids. The color\n"
	    "defaults to dark gray (Gray20). If this is not given, any ellipsoids will be drawn with principal\n"
	    "axes of 0.00015 times the overall scale factor, which should normally be appropriate.\n"
	    "\n"
	    "\n"
	    " background color\n"
	    "\n"
	    "sets the color of the background of the graphical views. The default color is white.\n"
	    "\n"
	    "\n"
	    " bestplane number name1 name2 ... nameN width height color\n"
	    "\n"
	    "causes the program to calculate the best fitting plane through a subset of atoms, where number\n"
	    "is the number of unique atom names (name and number, e.g. C8) name1 to nameN that follow.\n"
	    "The plane is drawn as a rectangle of dimensions width x height in the given color.\n"
	    "\n"
	    "\n"
	    " betaij         name number color\n"
	    "\n"
	    "defines the anisotropic thermal coefficients for an atom and color of the ellipsoid. The name and\n"
	    "number should correspond to the atom input described above. In the POV version of the program, the\n"
	    "principal ellipses are drawn in black. The ellipses will, of course, be invisible if the ellipsoid\n"
	    "is also black.\n"
	    "\n"
	    "\n"
	    " bij or Bij     name color\n"
	    "\n"
	    "defines the anisotropic thermal coefficients for an atom and the color of the ellipsoid to be drawn.\n"
	    "The name and number should correspond to the atom input described above.\n"
	    "\n"
	    "\n"
	    " bond           name1 name2 radius min max color\n"
	    "\n"
	    "where name1 and name2 indicate the types of atoms to be connected by a bond, radius is the radius\n"
	    "of the resulting cylinder, and the minimum and maximum lengths are given in the same units as the\n"
	    "unit cell.\n"
	    "\n"
	    "\n"
	    " box            radius color\n"
	    "\n"
	    "defines the radius and color of the cylinders that form the unit cell boundary. If radius is 0.0,\n"
	    "plotting of the unit cell is suppressed. The radius of the cylinders will be scaled with the size\n"
	    "of the drawing. The default size is 0.02.\n"
	    "\n"
	    "\n"
	    " cell           a b c alpha beta gamma\n"
	    "\n"
	    "unit-cell lengths and angles. If no angles are listed, they are assumed to be the fixed values for\n"
	    "the symmetry class.\n"
	    "\n"
	    "\n"
	    " clip          xmin xmax ymin ymax zmin zmax\n"
	    "\n"
	    "defines a-,b-,c- clipping range in fractions of the axes. Any bonds extending beyond these limits\n"
	    "will be cut off at half-length. This command is to be used in conjunction with the pack keyword to\n"
	    "produce 'dangling' bonds in the display of framework structures.\n"
	    "\n"
	    "\n"
	    " cutout         color     (used only for POV and openGL)\n"
	    "\n"
	    "sets the POV generation of thermal ellipsoids to have one octant removed, as in the program ORTEP.\n"
	    "If this command is not given, all ellipsoids will be complete. The color is for the planes that\n"
	    "describe the edges of the cutout.\n"
	    "\n"
	    "\n"
	    " dash           name1 name2 radius min max color\n"
	    "\n"
	    "where name1 and name2 indicate the types of atoms to be connected by a dashed bond, radius is the\n"
	    "radius of the resulting cylinder, and the minimum and maximum lengths are given in the same units\n"
	    "as the unit cell.\n"
	    "\n"
	    "\n"
	    " depthcue       depth (POV and openGL only)\n"
	    "\n"
	    "defines the extent to which the size of polyhedral edges is increased as the edge is closer to the\n"
	    "viewer.\n"
	    "\n"
	    "\n"
	    " edges          radius color\n"
	    "\n"
	    "defines the thickness and color of cylinders along the edges of polyhedra that may be used to\n"
	    "emphasize the faces. The radius of these cylinders will also be scaled with the size of the drawing.\n"
	    "By default, black edges of size 0.02 will be drawn.\n"
	    "\n"
	    "\n"
	    " ellipcolor     name number color\n"
	    "\n"
	    "defines the color for ellipsoids when the thermal ellipsoid information has been read from a CIF,\n"
	    "GSAS, SCHAKAL or SHELX import or inline file. The name and number must match the identification\n"
	    "information in the input file. The parameter number may be an asterisk (*) to indicate all atoms\n"
	    "with that name. In addition, these input lines must be after the import or inline command.\n"
	    "\n"
	    "\n"
	    " ellipsoids     probability\n"
	    "\n"
	    "sets the size of the ellipsoid such that that fraction of the electron density is contained within\n"
	    "the bounding surface. Use either 0.50 or 50 to get the standard (default) 50% ellipsoids.\n"
	    "\n"
	    "\n"
	    " finish        ambient diffuse specular roughness\n"
	    "\n"
	    "defines parameters for the POV lighting functions that are applied to\n"
	    "all surfaces. Suggested values are 0.7 0.3 0.08 0.01 to reduce the harsh\n"
	    "contrasts that can result from the default material properties in POV.\n"
	    "\n"
	    "\n"
	    " frame        comment\n"
	    "\n"
	    "similar to 'end', marks the division between two complete sets of input that are to be superimposed\n"
	    "in a single output file. Use this where you would otherwise have to create separate datafiles (e.g.\n"
	    "one showing polyhedra, the other selected bonds) and join the resulting POV or VRML files.\n"
	    "For V4.1 and later, each frame may have different space groups, pack range, and object descriptions.\n"
	    "With this functionality, a single drawing can represent a cage-compound framework with an adsorbed\n"
	    "molecule that has lower symmetry than the framework. Another option is to draw ball-and-stick and\n"
	    "polyhedral pictures in side-by-side unit cells.\n"
	    "\n"
	    "\n"
	    " import         cif filename datablock\n"
	    " import         csd filename (or fdat filename)\n"
	    " import         gsas filename phasenumber\n"
	    " import         pcr filename phasenumber\n"
	    " import         schakal filename\n"
	    " import         shelx filename\n"
	    " import         wien2k filename\n"
	    " import         discus filename\n"
	    " import	 exciting filename\n"
	    " import	 elk filename\n"
	    "\n"
	    "causes the program to read information from the external file specified.\n"
	    "Import filters have been written for the CIF, FDAT (Cambridge Structure Database=CSD), FullProf (pcr),\n"
	    "GSAS, SCHAKAL, SHELX, DISCUS, WIEN2k and ELK (Exciting) formats.\n"
	    "For GSAS and FullProf format, the number of the phase should also be given. For CIF format, the\n"
	    "data block number should be given.\n"
	    "From these files, the atomic coordinates, thermal parameters, unit cell, and space group will be read.\n"
	    "To turn on ellipsoid output, and to set colors for the ellipsoids, use the ellipcolor command.\n"
	    "\n"
	    "\n"
	    " inline          csd (or fdat)\n"
	    " inline          schakal\n"
	    " inline          shelx\n"
	    " inline          wien2k\n"
	    "\n"
	    "is similar to import, except that the foreign input information is included in the DRAWxtl input file\n"
	    "in the lines immediately following this command. This form presently works for FDAT,\n"
	    "SCHAKAL, SHELX and WIEN2k data. To turn on ellipsoid output, and to set colors for the ellipsoids,\n"
	    "use the ellipcolor command.\n"
	    "\n"
	    "\n"
	    " labelscale size\n"
	    "\n"
	    "changes the relative size of 'labeltext' entries. The default is 1.0\n"
	    "\n"
	    "\n"
	    " labeltext x y z text\n"
	    "\n"
	    "places the text at position x,y,z\n"
	    "\n"
	    "\n"
	    " list maxdist\n"
	    "\n"
	    "causes the program to list bond distances up to 'maxdist' in the preliminary scan. If this command is\n"
	    "not given, 'maxdist' defaults to 3.5 of the input units. If pm are used, this input will definitely be needed\n"
	    "\n"
	    "\n"
	    " lonepair name number height radius1 radius2 color\n"
	    "\n"
	    "creates the specified 'number' (either 1 or 2) of cones representing free electron pairs extending from\n"
	    "atom 'name', where 'height' is the length of the cone, 'radius1' is the size of the tip, and 'radius2'\n"
	    "is the size of the spherical end cap. The position of the cone is derived from those of the neighboring\n"
	    "atoms within a search range governed by the 'list' keyword (or 3.0 units by default).\n"
	    "Electron pairs point at the vertices of a tetrahedron by default - use a negative value for\n"
	    "'number' to create one or two electron pairs on an atom in a planar environment.\n"
	    "\n"
	    "\n"
	    " lookat u1 u2 u3 v1 v2 v3\n"
	    "\n"
	    "causes the program to select an orientation such that vector u is towards the viewer, and\n"
	    "(u x v) * u is horizontal.\n"
	    "\n"
	    "\n"
	    " mag_trans      Aa Ab Ac Ba Bb Bc Ca Cb Cc\n"
	    "\n"
	    "describes the relationship between the magnetic and nuclear unit cells. In this notation, the\n"
	    "upper-case letter states which of the magnetic axes is being described, and the lower-case\n"
	    "letter corresponds to the nuclear cell axis. This matrix defaults to the identity.\n"
	    "\n"
	    "\n"
	    " magnification  factor\n"
	    "\n"
	    "sets the factor to modify the overall scaling in case the automatic value is not correct.\n"
	    "\n"
	    "\n"
	    " mapcalclimits xmin xmax ymin ymax zmin zmax\n"
	    "\n"
	    "describes the region of direct space (in fractional coordinates) for which the map has been calculated.\n"
	    "Map types that are self documenting such as FullProf's GFOURIER and JANA2000 do not need this line.\n"
	    " For other types, 0 to 1 in all three directions will be assumed.\n"
	    "\n"
	    "\n"
	    " mapcontour level style color\n"
	    "\n"
	    "defines a new contour at 'level'. The style can be either 'mesh' or 'solid', and the color is\n"
	    "set by 'color'\n"
	    "\n"
	    "\n"
	    " mapcontour2d lower step upper color\n"
	    "\n"
	    "defines a new set of 2d contours beginning at 'lower', with 'step' between contours. No contours\n"
	    "beyond 'top' will be drawn. The color is set by 'color'. See the 'mapregion' command to see how\n"
	    " to set 2d mode.\n"
	    "\n"
	    "\n"
	    " maplegend\n"
	    "\n"
	    "draws a color ramp legend in the top left corner of the image\n"
	    "\n"
	    "\n"
	    " mapread  maptype filename calctype resol\n"
	    "\n"
	    "reads a Fourier map of type 'maptype' from the file named 'filename'. At present, GSAS-style\n"
	    "(maptype = grd), JANA2000-style (maptype = stf), WIEN2k (maptype = w2k), VASP (maptype = vsp),\n"
	    "FullProf (GFOURIER output, maptype = flp), and O Format (maptype = dn6) electron density maps \n"
	    "are read, as are electron density and ELF files from the FP-LAPW program EXCITING (maptype=exc)\n"
	    "or its successor ELK and files in XCrysDen format (maptype=xsf).\n"
	    "If a Shelx/Cif-style Fo/Fc file (maptype = fcf) or a JANA-style M80 file (maptype = m80) is given,\n"
	    " the electron density is calculated during the initial read, which may take a few seconds. Both\n"
	    " A/B and Fo/phi data formats (Shelx commands LIST 3 and LIST 6) are supported. The calctype may be\n"
	    " 'Fo', 'Fc', 'Fo-Fc', or '2Fo-Fc' to indicate the type of map to calculate. If this parameter is not\n"
	    " given, an 'Fo' map is calculated. Parameter resol indicates the number of steps per cell unit.\n"
	    " If not given, the value is set to 4\n"
	    "\n"
	    "\n"
	    " mapregion xmin xmax ymin ymax zmin zmax\n"
	    "\n"
	    "describes the region of direct space (in fractional coordinates) that the map is to be displayed\n"
	    "in the output. If not entered, these values default to the values given under 'mapcalclimits'. If\n"
	    "the difference between the maximum and the minimum along one direction is zero, a 2d map will be\n"
	    "plotted for the other two coordinates.\n"
	    "\n"
	    "\n"
	    " mapslice px py pz nx ny nz type\n"
	    "\n"
	    "describes a planar slice through a fourier map at the location given by the coordinates px,py,pz\n"
	    "of a point in the plane and in the orientation determined by the plane normal nx,ny,nz.\n"
	    "Possible type values are 1 (contoured slice), 2 (solid, color-coded blue to red) or 3 (solid, black\n"
	    "to white).\n"
	    "\n"
	    "\n"
	    " molcomp        dist\n"
	    "\n"
	    "causes any incomplete molecules in the display box to be completed. The value of dist defines the\n"
	    "maximum intramolecular distance. Caution: If this distance is greater than any intermolecular distance,\n"
	    "or if the material is not molecular, the display list will overflow.\n"
	    "\n" "\n"
/*
" molecule       atom_name atom_number  dist\n"
"\n"
"sets up a molecular completion command about the atom specified by the atom name and number.\n"
"This command is identical to a pack command about the position of the specified atom, followed\n"
"by a molcomp command. The caution stated in the molcomp command also applies here.\n"
"\n"
"\n"
*/
	" nolabels\n"
	    "\n"
	    "removes all axis labels from the output diagrams.\n"
	    "\n"
	    "\n"
	    " noshadow\n"
	    "\n"
	    "causes objects in the POV file not to cast shadows.\n"
	    "\n"
	    "\n"
	    " occupancy     name  average  minimum\n"
	    "\n"
	    "defines the occupancy of the named site in the average structure of a modulated system,\n"
	    "and the occupancy threshold for including individual copies in a plot of the modulated\n"
	    "structure.\n"
	    "Use a negative value for the sphere radius to scale atom sizes by their individual site\n"
	    "occupancies.\n"
	    "\n"
	    "\n"
	    " origin         xcenter ycenter zcenter\n"
	    "\n"
	    "defines center of view box in crystal coordinates (defaults to 0.5 0.5 0.5).\n"
	    "\n"
	    "\n"
	    " orthographic\n"
	    "\n"
	    "causes the camera to be changed from the normal perspective view to an orthographic view.\n"
	    "\n"
	    "\n"
	    " pack           xmin xmax ymin ymax zmin zmax\n"
	    "\n"
	    "defines a-,b-,c- plotting range in fractions of the axes, similar to the PLUTO (Motherwell & Clegg 1978))\n"
	    "PACK RANGE command (this is especially useful for highly oblique cells, where the orthorhombic\n"
	    "view box does not always give satisfactory results). When drawing multiple frames, each frame can\n"
	    "have distinct packing values. If new values are not defined, they will be derived from the previous\n"
	    "frame.\n"
	    "\n"
	    "\n"
	    " phaseshift     value1   value2   value3\n"
	    "\n"
	    "defines the initial phases  t_n  of the n'th modulation wave in a modulated structure.\n"
	    "\n"
	    "\n"
	    " phong          value size                (used only for POV and openGL)\n"
	    "\n"
	    "defines the amount of Phong highlighting on spheres and ellipsoids.  The value ranges between 0.0\n"
	    "and 1.0, where 0.0 gives no highlight, and 1.0 causes complete saturation at the center of the highlight.\n"
	    "The size ranges from 1.0 (very dull) to 250 highly) polished). The default quantities are 0.1 and 1.0,\n"
	    "which gives a large, dull highlight.  If value is 0.0, the image can be rendered much more quickly.\n"
	    "\n"
	    "\n"
	    " plane          name length color\n"
	    "\n"
	    "defines the center of a plane group, such as CO3 that is to be drawn in a structure, where name is the\n"
	    "name of the atom at the center and length is the maximum distance to coordinating anions.\n"
	    "\n"
	    "\n"
	    " polyedge name radius color\n"
	    "\n"
	    "defines the thickness and color of cylinders used to emphasize the faces along the edges of polyhedra\n"
	    "for atom 'name'. The radius of these cylinders will also be scaled with the size of the drawing.\n"
	    "\n"
	    "\n"
	    " polysz         name length color\n"
	    "\n"
	    "defines a polyhedron, where name is the name of an atom at the center of a polyhedron and length is\n"
	    "the maximum length of distances to atoms that are to be considered as the vertices of the polyhedron.\n"
	    "The polyhedra can be of any desired complexity. For polyhedra with both upper and lower limits\n"
	    "(which might be desirable for intermetallic compounds), use the 'shell' command. To control both\n"
	    "center and target atoms, use the 'polyvert' command\n"
	    "\n"
	    "\n"
	    " polytolerance      factor\n"
	    "\n"
	    "modifies the internal limit for the deviation of vertices from the common plane. While the default\n"
	    "value (0.1) will always generate correct drawings, it may sometimes be desirable to increase it to\n"
	    "create idealized views of nearly symmetrical polyhedra that would otherwise show creased surfaces.\n"
	    "\n"
	    "\n"
	    " polyvert       name1 name2 length color\n"
	    "\n"
	    "defines a polyhedron, where atoms of type name1 are at the center of the polyhedron, atoms of\n"
	    "type name2 are at the vertices and length is the maximum distance to be included.\n"
	    "\n"
	    "\n"
	    " qvector       value1 value2 value3\n"
	    "\n"
	    "defines the components of the wave vector q for a modulated structure\n"
	    "\n"
	    "\n"
	    " rem            text\n"
	    " REM            text\n"
	    "\n"
	    "Any line preceded by this command is ignored.\n"
	    "\n"
	    "\n"
	    " shell         name length1 length2 color\n"
	    "\n"
	    "defines a polyhedral hull, where name is the name of an atom at the center of a polyhedron and\n"
	    "length1 and length2 are the minimum and maximum lengths of distances to atoms that are to be\n"
	    "considered as the vertices of the polyhedron. The polyhedra can be of any desired complexity,\n"
	    "and can be stacked as desired.\n"
	    "\n"
	    "\n"
	    " slab a b c alpha beta gamma xoff yoff zoff xrot yrot zrot flag\n"
	    "\n"
	    "defines a (possibly oblique) cutout box of the specified axis lengths and\n"
	    "angles that is offset by xoff,yoff zoff from the origin of the structure\n"
	    "and rotated at angles xrot yrot zrot relative to it. If flag is set to 1,\n"
	    "any part of the structure outside the box is deleted. If flag is 2, the\n"
	    "outline of the box is overlaid on the unchanged image to allow accurate\n"
	    "placement of the cutout box.\n"
	    "\n"
	    "\n"
	    " spgp    symbol\n"
	    " spgr    symbol\n"
	    " sgrp    symbol\n"
	    "\n"
	    "Space Group name consisting of the Bravais lattice symbol (must be upper case) followed by a space,\n"
	    "the elements parallel to the first axis followed by a space, etc. Examples are I 41/a m d, P 21/n, I a 3 d,\n"
	    "P b n m, etc. The generators will always select the origin choice with a center of symmetry at the\n"
	    "origin. Furthermore, all monoclinic cells will have the unique axis parallel to the b axis, unless the\n"
	    "full symbol is used, i.e. P 1 1 21/n describes a monoclinic cell with c as the unique axis. N.B.:\n"
	    "Rhombohedral space groups must be represented in the hexagonal form.\n"
	    "\n"
	    "\n"
	    " sphere         name radius color\n"
	    " sphere         name number radius color\n"
	    "\n"
	    "where name is a one- or two-character symbol of the atom type, radius is the radius of the sphere in\n"
	    "the input units, and color is the color of the sphere to be drawn. If the first form is used, all atoms with\n"
	    "that name will be drawn. The second restricts the command to that name and number only\n"
	    "\n"
	    "\n"
	    " title, titl    text\n"
	    "\n"
	    "General description of the structure - this line may appear anywhere in the file, but is generally first.\n"
	    "\n"
	    "\n"
	    " uij, Uij       name number u color\n"
	    "\n"
	    "defines the anisotropic thermal coefficients for an atom and color of the ellipsoid. The name and\n"
	    "number should correspond to the atom input described above.\n"
	    "\n"
	    "\n"
	    " values         name * radius\n"
	    " values         name number radius\n"
	    "\n"
	    "defines additional atomic properties for the given element or individual atom. Currently the\n"
	    "only supported property is the van der Waals radius to be used in cavity calculations (see the\n"
	    "''voids'' keyword).\n"
	    "\n"
	    "\n"
	    " vectors     posx  posy  posz\n"
	    "\n"
	    "turns on the orientation vector triple at a corner of the diagram. The optional parameters pos[xyz]\n"
	    "specify the location of the origin. If not given, the program will guess at a location\n"
	    "\n"
	    "\n"
	    " view  xrot  yrot  zrot\n"
	    "\n"
	    "where xrot, yrot and zrot are view rotation angles in Cartesian space. These values correspond to a\n"
	    "rotation of xrot about the x axis, followed by a rotation of yrot about the new y axis, and, a rotation\n"
	    "of zrot about the new z axis.\n"
	    "\n"
	    "\n"
	    " voids method probe_radius gridx gridy gridz color\n"
	    "\n"
	    "causes cavities in the structure to be determined using the given method and probe radius (e.g.\n"
	    "1.4 for a water molecule). gridx, gridy and gridz determine the resolution of the grid used\n"
	    "for subdividing the unit cell. Currently supported methods are 1 for a slow but reliable, sequential\n"
	    "test at all gridpoints, 2 for using the MSMS program of Sanner, and 3 for a pseudorandom sampling\n"
	    "of points. Methods 2 and 3 should be considered experimental.\n"
	    "\n"
	    "\n"
	    " vrml1\n"
	    "\n"
	    "causes the output VRML file to have the VRML1 syntax (as opposed to the newer VRML97 standard)\n"
	    "\n"
	    "\n"
	    " x3d\n"
	    "\n"
	    "causes the output VRML file to have X3D 'Classic VRML' encoding and .x3dv extension\n"
	    "\n"
	    "\n"
	    " xyzoff         u1  u2 u3\n"
	    "\n"
	    "causes all atom coordinates to be shifted by -u.  This command is used whenever the origin defined for\n"
	    "a structure does not conform to the standard origin selected by the space-group generator.\n"
	    "\n"
	    "\n"
	    " end\n"
	    " END\n"
	    "\n"
	    "The last line of a file that is read. Any information past this point will be ignored.\n"
	    "\n"
	    "\n"
	    "The maximum line length is 255 characters - everything beyond this limit is ignored.\n"
	    "\n"
	    "Sample input file\n"
	    "\n"
	    "title Buckyball with balls and sticks\n"
	    "cell 14.16 14.16  14.16\n"
	    "pack -.3 .3 -.3 .3 -.3 .3\n"
	    "spgp F m 3\n"
	    "sphere c 0.4 Red filter 0.3\n"
	    "bond c c 0.1 1.2 1.5 Gray30\n"
	    "atom c 1 0.04908 0.00000 0.24510\n"
	    "atom c 2 0.10028 0.08284 0.21346\n"
	    "atom c 3 0.18313 0.05120 0.16226\n"
	    "origin 0 0 0\n" "magnification 0.7\n" "view -18 0 0\n" "end\n"
    };

    static int zero = 0;

    if (helpwindow) {
	helpwindow->show ();
	return;
    }
    helpwindow = new Fl_Window (50, 50, 680, 450, "Input File Commands");
    helpwindow->resizable (helpwindow);
    helpwindow->begin ();
    helpwindow->callback ((Fl_Callback *) View_Help_Close_cb, &zero);
    Fl_Text_Editor *display = new Fl_Text_Editor (0, 0, 680, 410);

    helpbuf = new Fl_Text_Buffer;
    display->buffer (helpbuf);
    display->textsize (13);
    display->textfont (FL_HELVETICA_BOLD);
    helpbuf->text (helptext);
    Fl_Button *o = new Fl_Button (300, 415, 80, 30, "Close");

    o->callback ((Fl_Callback *) View_Help_Close_cb, &zero);
#if !defined (WIN32) && !defined (__APPLE__)
    helpwindow->icon ((char *) drvui->icon);
#endif
    helpwindow->end ();
    helpwindow->show ();
}

void
Spacegroup_Help_cb (void)
{
// callback routine for help on alternate origin space groups
    int y;

    char string2[100];

    static int five = 5;

    if (helpwindow5) {
	helpwindow5->show ();
	return;
    }
    y = 160 + 100 * drvui->origin1_flag;
    helpwindow5 = new Fl_Window (50, 50, 440, y, "Space Groups with Alternate Origin");
    helpwindow5->begin ();
    helpwindow5->callback ((Fl_Callback *) View_Help_Close_cb, &five);
#if !defined (WIN32) && !defined (__APPLE__)
    helpwindow5->icon ((char *) drvui->icon);
#endif
    y = 40;
    Fl_Text_Display *box1 = new Fl_Text_Display (20, y, 380, 0, "Space Group Help");

    box1->box (FL_NO_BOX);
    box1->labelsize (24);
    box1->labelcolor ((Fl_Color) 1);
    y += 30;
    if (drvui->origin1_flag) {
	Fl_Text_Display *box2 = new Fl_Text_Display (20, y, 380, 0,
						     "The current space group has multiple origins.");

	box2->box (FL_NO_BOX);
	box2->labelsize (16);
	box2->labelcolor ((Fl_Color) 186);
	y += 20;
	Fl_Text_Display *box3 = new Fl_Text_Display (20, y, 380, 0,
						     "The generator used by DRAWxtl always selects");

	box3->box (FL_NO_BOX);
	box3->labelsize (16);
	box3->labelcolor ((Fl_Color) 186);
	y += 20;
	Fl_Text_Display *box4 = new Fl_Text_Display (20, y, 380, 0,
						     "the origin at a center of inversion. If your structure");

	box4->box (FL_NO_BOX);
	box4->labelsize (16);
	box4->labelcolor ((Fl_Color) 186);
	y += 20;
	Fl_Text_Display *box5 = new Fl_Text_Display (20, y, 380, 0,
						     "does not display correctly, remove the 'rem' from the");

	box5->box (FL_NO_BOX);
	box5->labelsize (16);
	box5->labelcolor ((Fl_Color) 186);
	y += 20;
	Fl_Text_Display *box6 = new Fl_Text_Display (20, y, 380, 0,
						     "following line of your input 'str' file:");

	box6->box (FL_NO_BOX);
	box6->labelsize (16);
	box6->labelcolor ((Fl_Color) 186);
	y += 15;
	sprintf (string2, "  rem xyzoff  %.3f %.3f %.3f", -drvui->origin_offset[0],
		 -drvui->origin_offset[1], -drvui->origin_offset[2]);
	Fl_Text_Display *box7 = new Fl_Text_Display (20, y, 380, 20);

	box7->box (FL_NO_BOX);
	Fl_Text_Buffer *buff = new Fl_Text_Buffer;

	box7->textsize (16);
	box7->textcolor ((Fl_Color) 186);
	box7->buffer (buff);
	buff->text (string2);
	y += 40;
    } else {
	Fl_Text_Display *box2 = new Fl_Text_Display (20, y, 380, 0,
						     "The current space group has a single choice of origin, or an");

	box2->box (FL_NO_BOX);
	box2->labelsize (16);
	box2->labelcolor ((Fl_Color) 186);
	y += 20;
	Fl_Text_Display *box3 = new Fl_Text_Display (20, y, 380, 0,
						     "'xyzoff' command is specified. No further action is required.");

	box3->box (FL_NO_BOX);
	box3->labelsize (16);
	box3->labelcolor ((Fl_Color) 186);
	y += 20;
    }
    Fl_Button *p = new Fl_Button (160, y, 80, 30, "Close");

    p->callback ((Fl_Callback *) View_Help_Close_cb, &five);
    helpwindow5->end ();
    helpwindow5->show ();
}

void
View_Help_Close_cb (Fl_Window *, int *arg)	// callback to destruct help window
{
// callback to close the STR help screens
    switch (*arg) {
    case 0:
	helpwindow->hide ();
	break;
    case 1:
	helpwindow1->hide ();
	break;
    case 2:
	helpwindow2->hide ();
	break;
    case 3:
	helpwindow3->hide ();
	break;
    case 4:
	helpwindow4->hide ();
	break;
    case 5:
	helpwindow5->hide ();
	break;
    case 6:
	helpwindow6->hide ();
	break;
    }
}
