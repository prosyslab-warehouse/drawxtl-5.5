// $Id: draw_ext.h 1107 2011-01-19 23:53:52Z martin $
//
/* draw_ext.h - external data definitions for DRAWxtl */

#include "DRAWxtlViewUI.h"
#include "EditView.h"
#include "Ellipsoids.h"

extern DRAWxtlViewUI *drvui;

extern int cur_show;

extern CrystalView *crystal;

extern ArrowParam *arrows;

extern SlabParam *Slabs;

extern Fl_Window *textwindow;

extern Fl_Window *listwindow1;

extern Fl_Window *listwindow2;

extern Fl_Window *listwindow3;

extern Fl_Window *listwindow4;

extern Fl_Window *listwindow5;

extern Fl_Window *helpwindow;

extern Fl_Window *helpwindow1;

extern Fl_Window *helpwindow2;

extern Fl_Window *helpwindow3;

extern Fl_Window *helpwindow4;

extern Fl_Window *helpwindow5;

extern Fl_Window *helpwindow6;

extern Fl_Text_Buffer *helpbuf;

extern Fl_Text_Buffer *helpbuf1;

extern AutomationParam *Automate;

extern BondParam *Bonds;

extern PolyParam *Polyhedra;

extern SphereParam *Spheres;

extern ConfigParm *Configure;

extern ConfigMiscParm *MiscConfigure;

extern ConfigMSMSParm *MSMSConfigure;

extern OmitParam *Omit;

extern EditScreen *edtprm;

extern Ellipsoids *ellipsoids;

extern LonePairParam *LonePairs;

extern MapsParam *Maps;

extern SliceParam *Slice;

extern SurfParam *Surf;

extern ModParam *Modparms;

extern GLuint selectBuf[BUFSIZE];

extern void *Edit_Str_Type;

extern char Edit_title[128];

extern int Edit_changed;

extern int Edit_loading;

extern Fl_Text_Buffer *textbuf;

extern Fl_Text_Buffer *textbuf1;

extern Fl_Text_Buffer *textbuf2;

extern Fl_Text_Buffer *textbuf3;

extern Fl_Text_Buffer *textbuf4;

//extern float gl_pos_x;
//extern float gl_pos_y;
//extern float gl_pos_z;
extern float gl_size;

extern void *Edit_Str_Type;

extern QUAT Rotq;

extern double Xrot, Yrot, Zrot;	// Rotation angles

extern double xmin, ymin, zmin;	// Minimums

extern double xmax, ymax, zmax;	// Maximums

extern Fl_Window *errorbox;

extern int Edit_changed;

extern int Edit_loading;

/* global variables for atom lists */

extern int natom;		/* number of different types of atoms */

extern int ncell;		/* number of atoms in asymmetric unit */

extern float *xypos;		/* master list fractional coordinates (modulated) */

extern float *xypos_nm;		/* master list fractional coordinates (not modulated) */

extern float *o_vert;		/* crystal coordinates of saved vertices (modulated) */

extern float *o_vert_nm;	/* crystal coordinates of saved vertices (not modulated) */

extern float *s_vert;		/* saved coordinates of vertices (modulated) */

extern int *vert_sym_no;	/* saved symmetry operator number */

extern int *vert_sym_nos;	/* saved symmetry operator number including sign */

extern int *poly_list;		/* storage for polygon corner pointers */

extern int nvert;		/* number of vertices in list */

extern int NvertM;		/* number of vertices in master atom list */

/* global variables for polyhedra, planes and bonds */

extern int draw_list;		/* number of items in polygon draw list */

extern int *vertex_list;	/* storage for vertices about a cation */

extern int numb_list;		/* number in polygon list */

extern int domolcomp;		/* non-zero if molecule completion requested */

/* global variables that control size and orientation of output object */

extern GLdouble modelMatrix[16];

extern GLdouble projMatrix[16];

extern GLint viewport[4];

extern float DepthCue;		/* scale for Z-dependent thickness of edges */

extern float boxlim[3];		/* half limits of plotting box */

extern float origin[3];		/* position of plotting origin */

extern double G_Rot[3][3];	/* Grand rotation matrix - transforms from Cartesian to picture */

extern float POV_Max[3];	/* Maximum limits of POV after rotation */

extern float POV_Min[3];	/* Minimum limits of POV after rotation */

extern int boxflag;		/* Non-zero if 'bounds' command given */

extern int packflag;		/* Non-zero if 'pack' command given */

extern int clipflag;		/* Non-zero if 'clip' command given */

extern int docell;		/* non-zero if unit-cell edges to be drawn */

extern int Display_axes;	/* non-zero if axial triple to be drawn */

extern int Color_Warning;	/* Warning flag for non-standard color */

extern float Magnification;	/* Magnification factor for image */

extern float Scale;		/* Scale of diagram */

extern float offset[3];		/* offset for vector triple */

extern float Text_Size;		/* Size to make axis label text */

extern int Options;		/* Place to save command-line options flags */

extern float rad_cell;		/* radius of unit cell framebars */

extern int edges;		/* draw thin black lines around the edges of polyhedra */

extern float xrot, yrot, zrot;	/* view rotation angles */

extern int Unit_Cell;		/* non-zero if unit cell to be drawn */

extern int no_comment;		/* non-zero to inhibit comment lines in VRML */

extern int M_cameras;		/* non-zero to inhibit multiple cameras */

extern int Vrml2;		/* non-zero to generate VRML97 (VRML2) output */

extern int X3D;			/* non-zero to write VRML2 in X3D-compatible syntax*/

extern int doVrml;		/* non-zero to allow VRML output */

extern int doPOV;		/* non-zero to allow povray output */

extern int doAsy;		/* non-zero to allow asymptote output */

extern int clipflag;

extern float printdist;		/* Distance limit for tabulated output */

extern int Labels;		/* True if labels should be output */

extern int slabmode;		/* cutout: slab vertices */

extern float slabx1, slaby1, slabz1, slabx2, slaby2, slabz2;

extern float slabx3, slaby3, slabz3, slabx4, slaby4, slabz4;

extern float slabv[24];

/* openGL cursor stuff */

extern float cur_cen[3];	/* location of cursor (fractional coordinates) */

extern int cur_atom[4];

extern char cur_name[4][10];

extern float dist12;

extern float dist23;

extern float dist34;

extern float ang123;

extern float ang234;

extern float torsion_ang;

/* Fourier map stuff */

#include "drawmap.h"
extern int ReadFourMap;		/* true if a Fourier map has been read */

extern float *FourierPt;	/* pointer to Fourier map contents */

extern float map_a, map_b, map_c, map_alpha, map_beta, map_gamma;	/* cell dimensions for map */

extern int mapstep_a, mapstep_b, mapstep_c;	/* map steps across unit cell */

extern float xMin, xMax, yMin, yMax, zMin, zMax;	/* bounding box for map */

extern float x4Val, x5Val, x6Val;	/* coordinates of map intersections with superspace */

extern float x4step, x5step, x6step;

extern char FourierFileName[256];

extern int FourierMapType;

extern int ShowMapLegend;

extern struct MAP_INFO Map_Info;

