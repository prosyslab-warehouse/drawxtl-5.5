// $Id: DRAWxtl_proto.h 1107 2011-01-19 23:53:52Z martin $
//
// DRAWxtl_proto.h - routine prototypes for for DRAWxtl V5.5 - the GUI version 
//
// Coded using the FLTK 1.1.6 widget set
//
//    Larry W. Finger, Martin Kroeker and Brian Toby
//
#ifndef DRAWxtl_proto_h
#define DRAWxtl_proto_h
void About_Help_cb (void);

void Add_Frame_Main (void);

void add_to_list (float *xp, int no, int nos);

void add_vert (float vert[3], int no, int check, int modulate, int nos);

void add_vert_nc (float *vert);

void analyze_bonds (void);

double AngleInRange (double angle);

void Arrow_Frame_Combo_cb (Fl_Widget * a, void *b);

void Automation_Edit_cb(void);

void Automation_Go_cb(void);

void Automation_Abort_cb(void);

void Automation_Close_cb(void);

void axeqb (double a1[3][3], double x[3], double b[3]);

void Bond_Combo_cb (Fl_Widget * a, void *b);

void Bond_Frame_Combo_cb (Fl_Widget * a, void *b);

void Blank_Strip (char input[]);

void Browse_Map_File_cb (void);

void build_box_contents (void);

void build_poly_list (int numb_list, int polyno, int l0, int l1, float *);

void Calc_Rot (float v1[3], float v2[3]);

void Check_Box_cb (void);

int check_atom_name (char *, char *);

void check_dynamic_storage (void);

int check_vert_alloc (int number, int alloc_ok);

void Clean_Up_Files (void);

void Clear_Last_Omit_cb (Fl_Button *, void *);

void Color_Help_cb (void);

void Configure_cb (void);

void Configure_Close_cb (void);

void Configure_Save_cb (void);

void Configure_Misc_cb (void);

void Configure_Misc_Save_cb (void);

void Configure_Misc_Close_cb (void);

void Configure_MSMS_cb (void);

void Configure_MSMS_loc_cb (void);

void Configure_Close_MSMS_cb (void);

void Configure_Save_MSMS_cb (void);

void ConfigurePOV_cb (void);

void ConfigurePOVOptions_cb (void);

void ConfigureMencoder_cb (void);

void ConfigureFFmpeg_cb (void);

static inline void Convert_Cart_Cryst(float bmat_inv[3][3], float inp[3], float oup[3], float origin[3])
{
	int i, j;

	for (i = 0; i < 3; i++) {
		oup[i] = origin[i];
		for (j = 0; j < 3; j++)
			oup[i] += bmat_inv[i][j] * inp[j];
	}
}

static inline void Convert_Cryst_Cart(double b_mat[3][3], float inp[3], float oup[3], float origin[3])
{
	int i, j;

	for (i = 0; i < 3; i++) {
		oup[i] = 0.0f;
		for (j = 0; j < 3; j++)
			oup[i] += (float)b_mat[i][j] * (inp[j] - origin[j]);
	}
}

void convert_ellipsoid (void);

float convert_pos (char *);

void Conv_Sym_Mat (void);

void Destroy_Open_Windows (void);

float determinant (double rot[3][3]);

void display_cursor_text (void);

float dist (int j, int k);

void do_drawing (void);

float dot0_3d (float x0, float y0, float z0, float x1, float y1, float z1,
	       float x2, float y2, float z2);
void draw_cell (int docell);

void draw_GL_triple (void);

void draw_cursor (void);

void Dump_View_cb (void);

void Edit_Arrow_cb (void);

void Edit_Arrow_Close_cb (class Fl_Button *, void *);

void Edit_Arrow_Save_cb (Fl_Button *, int *save);

void Edit_Bond_cb (void);

void Edit_Bond_Close_cb (void);

void Edit_Bond_Save_cb (Fl_Button *, int *save);

void Edit_Changed_cb (int, int nInserted, int nDeleted, int, const char *, void *v);

void Edit_Ellipsoid_cb (void);

void Edit_Ellipsoid_Save_cb (Fl_Button *, int *save);

void Edit_LonePair_cb (void);

void Edit_Lone_Pair_Close_cb (void);

void Edit_Lone_Pair_Save_cb (Fl_Button *, int *save);

void Edit_Maps_cb (void);

void Edit_Maps_Close_cb (void);

void Edit_Maps_Save_cb (Fl_Button *, int *save);

void Edit_Slice_cb (void);

void Edit_Slice_Close_cb (void);

void Edit_Slice_Save_cb (Fl_Button *, int *save);

void Edit_Surfaces_cb (void);

void Edit_Surfaces_Close_cb (void);

void Edit_Surfaces_Save_cb (Fl_Button *, int *save);

void Surface_Frame_Combo_cb (Fl_Widget * a, void *b);

void Edit_Modparms_cb (void);

void Edit_Modparms_Close_cb (void);

void Edit_Modparms_Save_cb (Fl_Button *, int *save);

void Edit_Parmeters_cb (void);

void Edit_Parmeters_Close_cb (Fl_Button *, void *);

void Edit_Parmeters_Save_cb (Fl_Button *, int *tosave);

void Edit_Polyhedra_cb (void);

void Edit_Polyhedra_Close_cb (void);

void Edit_Polyhedra_Save_cb (Fl_Button *, int *save);

void Edit_Slab_cb (void);

void Edit_Slab_Close_cb (void);

void Edit_Slab_Save_cb (Fl_Button *, int *save);

void Edit_Spheres_cb (void);

void Edit_Spheres_Close_cb (void);

void Edit_Spheres_Save_cb (Fl_Button *, int *save);

void Edit_STR_cb (Fl_Menu_ *, void *arg);

void Edit_STR_Close_cb (void);

void Edit_STR_Save_cb (Fl_Button *, int *action);

int eigen (float *biso, float beta[6], float valu[3], float vect[3][3]);

int end_flip (int);

float end_flip_real (float);

void Error_Box (const char *message);

void Error_Box_cb (Fl_Widget * w, void *d);

void Exit_cb (void);

void expand_atom (int natom);

int fillcube (uchar * cube, FILE * mapin);

void findsys (void);

void find_all_in_box (int i);

int find_atom (void);

void find_lattice_type (void);

int find_proj_atom (int, int);

void generate_arrows (void);

void generate_bonds (void);

void generate_cones (void);

void Generate_Drawing (int);

void generate_ellipsoids (void);

void generate_gl_texts (void);

void generate_planes (void);

void generate_poly (void);

void generate_slab (void);

void generate_spheres (void);

void generate_texts (void);

void generate_triple (void);

void generate_lsq_planes (void);

void generate_aimsurf (void);

void generate_voids (void);

void getsym (char *text, int num, int kk);

void get_atom_id (void);

void get_input (int Quick);

int get_next_token (char *p, int max_len, FILE * fpin);

int Get_Unique_Atoms (char atoms[100][5], int Frame_No);

void get_label (char input[], char *c1, char *c2, char *c3, char *c4, int strip);

void Graphics_Help_cb (void);

void ImportDataFile_cb (Fl_Widget *, void *);

void import_cif (char *, int, int, int *, int);

void import_discus (char *, int, int, int);

void import_fdat (char *, int, int, int);

void import_gsas (char *, int, int, int);

void import_pcr (char *, int, int, int);

void import_schakal (char *, int, int, int);

void import_shelx (char *, int, int, int);

void import_wien (char *, int, int, int);

void import_exc (char *, int, int, int);

void Include_Cutouts_cb (void);

void init_dynamic_storage ();

void Init_DRAWxtl (void);

void Input_Help_cb (void);

float InterpolateMap (float x, float y, float z);

void label_cell (void);

void Load_Bond_Data (const char *atom, char table[20480]);

void Load_Color_Combo (Flu_Combo_List * ot);

void LoadConfig (bool full_load);

void Locate_Triple (void);

void Lone_Pair_Combo_cb (Fl_Widget * a, void *b);

void LonePair_Frame_Combo_cb (Fl_Widget * a, void *b);

int LookupMap (int ix, int iy, int iz);

void make_bmat (int sys, float lat_con[6], double b_mat[3][3], float ginv[3][3],
		float rec_lat_con[6]);
void Main_Frame_Combo_cb (Fl_Widget *, void *);

void Maps_Frame_Combo_cb (Fl_Widget *, void *);

void Map_Info_cb (void);

void MapType_cb (void);

float matinv (float a[3][3]);

void matmul (float a[3][3], float b[3][3], float c[3][3]);

void Max_Min_cb (void);

int Maximize_rho (int);

void Modify_Arrow_cb (Fl_Widget * a, void *b);

void Modify_Bonds_cb (Fl_Widget * a, void *b);

void Modify_Bonds_Distance_cb (Fl_Widget * a, void *b);

void Modify_LonePair_cb (Fl_Widget * a, void *b);

void Modify_Maps_cb (Fl_Widget * a, void *b);

void Modify_AimSurfaces_cb (Fl_Widget * a, void *b);

void Modify_Surfaces_cb (Fl_Widget * a, void *b);

void Modify_Polyhedra_cb (Fl_Widget * a, void *b);

void Modify_Polyhedra_Distance_cb (Fl_Widget * a, void *b);

void Modify_Spheres_cb (Fl_Widget * a, void *b);

void Modify_Occ_cb (Fl_Widget * a, void *b);

void Occ_Combo_cb (Fl_Widget * a, void *b);

void New_Occ_Add_cb (class Fl_Widget *, int *action);

void New_Occ_Input_cb (class Fl_Widget *, void *);

void Modify_Surfaces_cb (Fl_Widget * a, void *b);

void modulate_parameters (float vert[3], double *occ, int sym_no, int atom_no);

int modulate_uij (float vert[3], int ellips_no, int atom_no, int sym_no, float uij[6]);

void move_cursor (int axis, float inc_amt);

void update_cursor_window (void);

//void moveto_atom (int, int, int, int);

int pick_label (int, int, int, int);

void New_Arrow_Add_cb (class Fl_Widget *, int *action);

void New_Arrow_Input_cb (class Fl_Widget *, void *);

void New_Bond_Add_cb (class Fl_Widget *, int *action);

void New_Bond_Input_cb (class Fl_Widget *, void *);

void New_Ellipsoid_Input_cb (class Fl_Widget *, void *);

void New_Map_Add_cb (class Fl_Widget *, int *action);

void New_Map_Input_cb (Fl_Widget *, void *);

void New_AimSurf_Add_cb (class Fl_Widget *, int *action);

void New_AimSurf_Input_cb (Fl_Widget *, void *);

void AimSurf_Combo_cb (Fl_Widget * a, void *b);

void Surf_Combo_cb (Fl_Widget * a, void *b);

void New_Radius_Add_cb (class Fl_Widget *, int *action);

void New_Surf_Input_cb (Fl_Widget *, void *);

void New_Polyhedra_Add_cb (class Fl_Widget *, int *action);

void New_Lone_Pair_Add_cb (class Fl_Widget *, int *action);

void New_Polyhedra_Input_cb (class Fl_Widget *, void *);

void New_Sphere_Add_cb (class Fl_Widget *, int *action);

void New_Sphere_Input_cb (class Fl_Widget *, void *);

void next_focus (void);

int not_in_slab (float x, float y, float z);

void Offset_cb (void);

void Output_Spheres (float *radii, int i);

float P_to_C (float prob);

int pick_box (int x, int y, int w, int h);

void plot_vrml_poly (int polyno);

void Polyhedra_Combo_cb (Fl_Widget * a, void *b);

void Polyhedra_Frame_Combo_cb (Fl_Widget * a, void *b);

void print_sym (void);

void process_hits (int hits, GLuint buffer[]);

void Process_Inp (int i);

void Progress_Window (int, const char *, float);

void push_cylinder (float df[3], float at[3], float radius, char *color);

void read_aim (char *infile, int Quick);

void read_dn6 (char *infile, int Quick);

void read_exc (char *infile, int Quick);

void read_xsf (char *infile, int Quick);

void read_fcf (char *infile, int Quick);

void read_flp (char *infile, int Quick);

void read_grd (char *infile, int Quick);

void read_inp (int Quick);

void read_m80 (char *infile, int Quick);

void read_m81 (char *infile, int Quick);

void read_stf (char *infile, int Quick);

void read_w2k (char *infile, int Quick);

void read_vasp (char *infile, int Quick);

void Restore_Working_Copy (void);

void Rotation_cb (void);

void Save_Current_cb (void);

void Save_Working_Copy (void);

void SelectDataFile_cb (void);

void set_tf_status (void);

void show_slab (void);

void show_slab_ovl (void);

void skip_blocks (int i, FILE * in);

void Slice_Frame_Combo_cb (Fl_Widget *, void *);

void Spacegroup_Help_cb (void);

void Sphere_Combo_cb (Fl_Widget * a, void *b);

void Sphere_Frame_Combo_cb (Fl_Widget * a, void *b);

void start_picking (int x, int y, int w, int h);

void sub_add_vert (float *vert, int no_cell);

void sub_add_vert_nc (float *vert, int no_cell);

void symop (char *input);

void Token_Strip (char string[], int no);

void Transform_POV_Color (char *);

void Transform_VRML_Color (char *);

void trim_string (char string[], int len);

void Update_Objects (int Frame_No, FILE * out);

int update_box (int hits, GLuint buffer[]);

void Update_Cursor_List (int);

void Update_Str (int overwrite);

int Unique_Atom (void);

void View_Console_cb (void);

void View_Cursor_cb (void);

void Cursor_Reset_Combo_cb (Fl_Widget * , void *);

void View_File_cb (void);

void View_Help_Close_cb (Fl_Window *, int *);

void View_Listing_cb (void);

void View_Listing_Close_cb (Fl_Window *, int *arg);

void View_POV_cb (void);

int vec_dif (int n1, float v1[3], int n2, float v2[3], int n3, float v3[3], float v[3]);

void WriteConfig (void);

void Write_Map_cb (void);

QUAT XYZ_Rot_to_Q (double Xrot, double Yrot, double Zrot);

void calculate_voids (void);

void calc_simplevoids (void);

void calculate_sas (void);

void calculate_msms (void);

void generate_msms (void);

void generate_simplevoids (void);

void dump_gif (void);

void Add_mapslice (int);
#endif
