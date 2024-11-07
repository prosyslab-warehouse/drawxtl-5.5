// $Id: read_flp.h 900 2009-08-13 20:00:45Z larry $
//
// read_flp.h - header file to match the Fortran 90 TYPE statements used in FullProf
//
// part of DRAWxtl, Copyright 2005 by L. W. Finger, M. Kroeker, and B. Topy
//
#ifndef READ_FLP
#define READ_FLP
struct crystal_cell_type
{
    float cell[3];
    float ang[3];
    float cell_std[3];
    float ang_std[3];
    float rcell[3];
    float rang[3];
    float GD[3][3];
    float GR[3][3];
    float Gr_Orth_cel[3][3];
    float Orth_Gr_cel[3][3];
    float CellVol;
    float RCellVol;
    char CartType;
};

struct Sym_Oper_Type
{
    int Rot[3][3];
    float Tr[3];
};

struct wyck_pos_type
{
    int multp;
    char site[6];
    int norb;
    char str_orbit[48][40];
    char extra_orbit[MXSYM][40];
};

struct wyckoff_type
{
    int num_orbit;
    struct wyck_pos_type orbit[26];
};

struct space_group_type
{
    int NumSpg;
    char SPG_Symb[20];
    char Nall[16];
    char CrystalSys[12];
    char Laue[5];
    char PG[5];
    char Info[5];
    char SG_Setting[80];
    int Hexa;
    char SPG_lat;
    char SPG_latsy[2];
    int NumLat;
    float Latt_Trans[12][3];
    char Bravais[51];
    char Centre[80];
    int Centred;
    float Centre_coord[3];
    int NumOps;
    int Multip;
    int Num_gen;
    struct Sym_Oper_Type Centre_Coord[MXSYM];
    char SymopSymb[MXSYM][40];
    struct wyckoff_type Wyckoff;
    float R_Asym_Unit[2][3];
};

struct space_group
{
    char spgr[20];
    struct space_group_type grp_espacial;
};

struct Atom_Type
{
    char Lab[10];
    char ChemSymb[2];
    char SfacSymb[4];
    int active;
    int Z;
    int mult;
    float x[3];
    float x_std[3];
    int lx[3];
    float occ;
    float occ_std;
    float mOcc;
    int lOcc;
    float Biso;
    float Biso_std;
    float mBiso;
    int lBiso;
    char utype[4];
    char thtype[5];
    float U[6];
    float U_std[6];
    float Ueq;
    float mU[6];
    float lU[6];
    float Charge;
    float Moment;
    int Ind[5];
    int Nvar;
    float VarF[10];
};

struct FFT_param
{
    int ngrid[3];
    float denmin;
    float denmax;
    float xlim[2];
    float ylim[2];
    float zlim[2];
    float xinx;
    float yinc;
    float zinc;
};
#endif
