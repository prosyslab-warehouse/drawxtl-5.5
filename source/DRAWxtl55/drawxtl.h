// $Id: drawxtl.h 1084 2010-11-17 22:58:30Z martin $
//
#ifndef DRAWxtl_h
#define DRAWxtl_h
#define MAX_VERTS 20000
#define VERSION "V5.5"
#define RAD 57.2957795		/* Conversion from radians to degrees */
#ifndef PI
#define PI 3.1415926
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define O_OPT 1
#define B_OPT 2
#define P_OPT 4
#define V_OPT 8
#define M_OPT 16
#define L_OPT 32
#define IDI_ICON1  101
#define BUFSIZE 512

#include <stdlib.h>
#include <string.h>
static inline void *
zalloc (size_t size)
{
    void *ret = malloc (size);

    if (ret)
	memset (ret, 0, size);
    return ret;
}

struct MAP_INFO
{
    int info_valid;
    char title[82];
    float lat_con[6];
    float rhomx;
    float rhomn;
    int map_int[3];
    float xlim[2];
    float ylim[2];
    float zlim[2];
    float x4lim[2];
    float x5lim[2];
    float x6lim[2];
    int map_type;
    int res;	/* points per A for calc maps */
};

struct arrow_struct
{
    int arrow_fn;		/* frame number for this arrow */
    char col_arrow[40];		/* the color of this one */
    float arrow_length;		/* the length of this arrow */
    float arrow_diam;		/*    and its diameter */
    float mag_xp[3];		/* position of arrow in real-space cell */
    float mag_xc[3];		/* components of arrow wrt XYZ */
};

struct cone_struct
{
    int cone_fn;		// frame number for the cone commands
    char cone_l1[5];		// label for atom end of cone
    int numlonepairs;		// number of lone pairs to draw at atom
    float cone_height;		// height of cone
    float cone_min;		// cone tip radius
    float cone_max;		// cone base radius
    char col_cone[40];		// color and properties of cones for POV
};

struct plane_struct
{
    int plane_fn;		// frame number for plane commands
    char plane_col[40];		// plane colors for POV 
    char plane_l[5];		// plane center label
    float plane_size;		// size of plane coordination
};

struct poly_struct
{
    int poly_fn;		// frame number for polyhedra commands
    char poly_l[5];		// polyhedral center label
    char poly_t[5];		// polyhedral vertex label
    float poly_size;		// size of polyhedron
    float poly_min;		// lower limit for shell-like polyhedron
    char poly_col[40];		// polyhedron colors for POV
    float poly_rad_edge;	// individual edge radius for polyhedron
    char poly_col_edge[40];	// color of individual poly edges for POV
};

struct edge_struct
{
    char name[5];
    char color[40];
    float radius;
};

struct bplane_struct
{
    char bplane_t[10][5];	// atom labels in plane list
    int bplane_n[10];		// atom numbers in plane list
    int nbatoms;		// number of items in list
    float bplane_d1;		// width of plane segment to draw
    float bplane_d2;		// second dimension
    char bplane_col[40];	// color
};

struct bond_struct
{
    int bond_fn;		// frame number for bond commands
    char bond_l1[5];		// label for one end of bond
    char bond_l2[5];		// label for other end
    float bond_size;		// bond diameter
    float bond_min;		// minimum bond distance
    float bond_max;		// maximum bond distance
    char col_bond[40];		// color and properties of bonds for POV
    int bond_style;		// solid bond or dash
};

struct sphere_struct
{
    int sphere_fn;		// frame number for sphere commands
    char sphere_l[5];		// sphere label
    int sphere_n;		// atom number - -1 for old behavior
    float sphere_size;		// sphere size
    char sphere_col[40];	// sphere colors
};

struct label_struct
{
    int label_fn;		// frame number
    float label_x[3];		// position of textual labels
    char label_label[64];	// label text
};

struct map_struct
{
    int FourierContourSolid;
    float FourierContourLevel;	// Level for 3d, initial level for 2d
    float FourierContourStep;	// Height between contours for 2d
    float FourierContourTop;	// Upper level for 2d
    char FourierContourColor[40];
    char FourierBackColor[40];
};

struct frame_struct		// allocate at max_frame
{
    float clip_lim[6];		/* limits of clip box in cell coords */
    float cryst_lim[6];		/* limits of view box in cell coords */
    float map_lim[9];
    int map_lim_set;
    int slice;
    float mapslice[3];
    float mapnorm[3];
    float planeeq[4];
    char molecule_atom_name[5];
    int molecule_atom_no;
    float molecule_distance;
};

struct atom_struct		// allocate at max_atom
{
    float atom_xyz[3];		// atomic coordinates
    char atom_l[5];		// atom labels
    int atom_n;			// input atom number
    int atom_fn;		// frame number for atom
    int sv_atom_n;		// saved input atom number
    int atom_ismod;		// flag for modulated position
    int occ_ismod;		// flag for modulated occupancy
    float occupancy;		// average occupancy
    float min_occ;		// minimum occupancy required for drawing the atom
    float saved_xyz[3];		// storage to save the atom position
    int TF_status;		// Temperature factor status (-1 no Uij, 0 sphere, 1 ellipsoid)
    float radius;		// vdW Radius for solvent calculations etc.
};

struct ellips_struct		// allocate at max_el
{
    char ellips_col[40];	// ellipsoid colors for POV
    float ellips_EV[3][3];	// Eigenvectors
    float ellips_RMS[3];	// RMS amplitudes
    char ellips_l[5];		// label for anisotropic atom
    int ellips_n;		// atom number
    int save_el_number;		// save original atom number
    int ell_type;		// form of ellipsoid input
    float ellips[6];		// thermal ellipsoid coefficients
    int ellips_ismod;		// flag for modulated Uij
};

struct mod_gbl_struct		// allocate at max_modulation
{
    int vector_mult[3];		// coefficients ijk for modvector = i*cell_vec_1 + j*cell_vec_2 ..
    float modvector[3];		// atom Fourier modulation vectors
};

struct mod_x_struct		// allocate at max_atom * max_modulation
{
    float atom_mod_sawtooth[5];	// amplitudes and interval for sawtooth modulation
    int ellips_modpar_atom;	// atom number for this item
    float atom_occpar[2];	// cosine and sine multiplier of occ fourier series
    int atom_occpar_atom;	// atom number for atom_occpar
    int atom_occpar_id;		// axis number (0-2) for atom_occpar
    float atom_occ_crenel[2];	// crenel lower and upper limit for occupancy
};

struct mod_3x_struct		// allocate at 3 * max_atom * max_modulation
{
    float atom_modpar[2];	// cosine and sine multiplier for an item
    int atom_modpar_atom;	// atom number for atom_modpar
    int atom_modpar_axis;	// axis number (0-2) for atom_modpar
    int atom_modpar_id;		// entry in vector_mult for this term (1 - no_mod_vectors)
};

struct mod_3t_struct		// allocate at 3 * max_el * max_modulation
{
    float ellips_modpar[2];	// wave vector id, cosine and sine multip
    int ellips_modpar_id;	// entry for this term (1 - no_mod_vectors)
    int ellips_modpar_term;	// identity of this term (1-6);
};

struct atprop_struct		// atomic properties
{
    int atprop_fn;		// frame number
    char atprop_l[5];		// label
    int atprop_n;		// atom number - -1 for old behavior
    float radius;		// van der waals size
};

#endif
