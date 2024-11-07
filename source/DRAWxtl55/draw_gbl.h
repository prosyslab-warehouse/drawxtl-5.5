// $Id: draw_gbl.h 1107 2011-01-19 23:53:52Z martin $
//
/* draw_gbl.h - global definitions for DRAWxtl */

#include "DRAWxtlViewUI.h"
#include "EditView.h"
#include "Ellipsoids.h"

DRAWxtlViewUI *drvui;

int cur_show = 0;

CrystalView *crystal;

ArrowParam *arrows;

SlabParam *Slabs;

Fl_Window *textwindow;

Fl_Window *listwindow1;

Fl_Window *listwindow2;

Fl_Window *listwindow3;

Fl_Window *listwindow4;

Fl_Window *listwindow5;

Fl_Window *helpwindow;

Fl_Window *helpwindow1;

Fl_Window *helpwindow2;

Fl_Window *helpwindow3;

Fl_Window *helpwindow4;

Fl_Window *helpwindow5;

Fl_Window *helpwindow6;

Fl_Text_Buffer *helpbuf = NULL;

Fl_Text_Buffer *helpbuf1 = NULL;

AutomationParam *Automate = NULL;

BondParam *Bonds = NULL;

PolyParam *Polyhedra = NULL;

SphereParam *Spheres = NULL;

ConfigParm *Configure;

ConfigMiscParm *MiscConfigure;

ConfigMSMSParm *MSMSConfigure;

OmitParam *Omit = NULL;

EditScreen *edtprm = NULL;

Ellipsoids *ellipsoids = NULL;

LonePairParam *LonePairs = NULL;

MapsParam *Maps = NULL;

SliceParam *Slice = NULL;

SurfParam *Surf = NULL;

ModParam *Modparms = NULL;

GLuint selectBuf[BUFSIZE];

Fl_Text_Buffer *textbuf = NULL;

Fl_Text_Buffer *textbuf1 = NULL;

Fl_Text_Buffer *textbuf2 = NULL;

Fl_Text_Buffer *textbuf3 = NULL;

Fl_Text_Buffer *textbuf4 = NULL;

//float gl_pos_x = 0.0;
//float gl_pos_y = 0.0;
//float gl_pos_z = 0.0;
float gl_size = 0.;

void *Edit_Str_Type;

QUAT Rotq;

double Xrot, Yrot, Zrot;	// Rotation angles

double xmin, ymin, zmin;	// Minimums

double xmax, ymax, zmax;	// Maximums

Fl_Window *errorbox;

int Edit_changed;

int Edit_loading;

/* global variables for atom lists */

int natom;			/* number of different types of atoms */

int ncell;			/* number of atoms in asymmetric unit */

float *xypos = NULL;		/* master list of fractional coordinates (modulated) */

float *xypos_nm = NULL;		/* master list of fractional coordinates (not modulated) */

float *o_vert = NULL;		/* crystal coordinates of saved vertices (modulated) */

float *o_vert_nm = NULL;	/* crystal coordinates of saved vertices (not modulated) */

float *s_vert = NULL;		/* saved coordinates of vertices (modulated) */

int *vert_sym_no = NULL;	/* saved symmetry operator number */

int *vert_sym_nos = NULL;	/* saved symmetry operator number including sign */

int *poly_list = NULL;		/* storage for polygon corner pointers */

int nvert;			/* number of vertices in list */

int NvertM;			/* number of vertices in master atom list */

/* global variables for polyhedra, planes and bonds */

int draw_list;			// number of items in polygon draw list

int *vertex_list = NULL;	/* storage for vertices about a cation */

int numb_list;			/* number in polygon list */

int domolcomp;			/* non-zero if molecule completion requested */

/* global variables that control size and orientation of output object */

GLdouble modelMatrix[16];

GLdouble projMatrix[16];

GLint viewport[4];

float DepthCue;			/* scale for Z-dependent thickness of edges */

float boxlim[3];		/* half limits of plotting box */

float origin[3];		/* position of plotting origin */

double G_Rot[3][3];		/* Grand rotation matrix - transforms from Cartesian to picture */

float POV_Max[3];		/* Maximum limits of POV after rotation */

float POV_Min[3];		/* Minimum limits of POV after rotation */

int boxflag;			/* Non-zero if 'bounds' command given */

int packflag;			/* Non-zero if 'pack' command given */

int clipflag;			/* Non-zero if 'clip' command given */

int docell;			// non-zero if unit-cell edges to be drawn

int Display_axes;		/* non-zero if axial triple to be drawn */

int Color_Warning;		/* Warning flag for non-standard color */

float Magnification;		// Magnification factor for image

float Scale;			/* Scale of diagram */

float offset[3];		/* offset for vector triple */

float Text_Size;		/* Size to make axis label text */

int Options;			/* Place to save command-line options flags */

float rad_cell;			/* radius of unit cell framebars */

int edges;			/* draw thin lines around the edges of polyhedra */

float xrot, yrot, zrot;		/* view rotation angles */

int Unit_Cell;			/* non-zero if unit cell to be drawn */

int no_comment;			/* non-zero to inhibit comment lines in VRML */

int M_cameras;			/* non-zero to inhibit multiple cameras */

int Vrml2;			/* non-zero to generate VRML97 (VRML2) output */

int X3D;			/* non-zero to write VRML2 in X3D-compatible syntax*/

int doVrml;			/* non-zero to allow VRML output */

int doPOV;			/* non-zero to allow povray output */

int doAsy;			/* non-zero to allow asymptote output */

float printdist;		/* Distance limit for tabulated output */

int Labels;			/* True if labels should be output */

int slabmode;			/* cutout: slab vertices */

float slabx1, slaby1, slabz1, slabx2, slaby2, slabz2;

float slabx3, slaby3, slabz3, slabx4, slaby4, slabz4;

float slabv[24];

float cur_cen[3];		/* location of cursor (fractional coordinates) */

int cur_atom[4];		/* sequence numbers of last four atoms under cursor */

char cur_name[4][10];		/* name of atom under cursor */

float dist12 = 0.0;

float dist23 = 0.0;

float dist34 = 0.0;

float ang123 = 0.0;

float ang234 = 0.0;

float torsion_ang = 0.0;

/* Fourier map stuff */

#include "drawmap.h"
int ReadFourMap = 0;		/* true if a Fourier map has been read */

float *FourierPt = NULL;	/* pointer to Fourier map contents */

float map_a, map_b, map_c, map_alpha, map_beta, map_gamma;	/* cell dimensions for map */

int mapstep_a, mapstep_b, mapstep_c;	/* map steps across unit cell */

float xMin, xMax, yMin, yMax, zMin, zMax;	/* bounding box for map */

float x4Val, x5Val, x6Val;	/* coordinates of map intersection with superspace */

float x4step, x5step, x6step;

char FourierFileName[256];

int FourierMapType;

int ShowMapLegend;

struct MAP_INFO Map_Info;
