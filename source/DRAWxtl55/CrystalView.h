// $Id: CrystalView.h 1093 2010-12-20 21:48:21Z martin $
//
// main header for DRAWxtl V5.5 - the GUI version
// Coded using the FLTK 1.1.6 widget set
//
//     Larry W. Finger, Martin Kroeker and Brian Toby
//
// defines the main mindow class and the various data classes

#ifndef CRYSTALVIEW_H
#define CRYSTALVIEW_H 1
#include <FL/Fl.H>
#include <FL/glut.H>
#include "Tb_Window.h"
#include <FL/Fl_Window.H>
#include <FL/gl.h>
#if defined(__APPLE__)
#  include <openGL/glu.h>
#else
#  include <GL/glu.h>
#endif

#include <stdlib.h>
#include <FL/glut.H>
#include <FL/Fl_Text_Editor.H>
#include <FL/Fl_File_Chooser.H>
#include "Flu_Spinner.h"
#include "Flu_Combo_List.h"
#include <FL/Fl_Multiline_Output.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Radio_Button.H>
#ifndef max
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#endif

class CrystalView:public Tb_Window
{
  public:
    CrystalView (int x, int y, int w, int h, const char *t);
    void draw (void);
    void draw_overlay (void);
     ~CrystalView ();
};

class BondParam
{				// class to hold bond screen parameters
  public:
    Fl_Window * Bond_Edit_Window;
    Fl_Text_Editor *Bond_Edit;
    Fl_Text_Buffer *BondBuffer;
    Fl_Text_Buffer *Bond_Output_Buffer;
    Flu_Combo_List *Bond_Combo;
    Fl_Input *New_Bond_From;
    Fl_Input *New_Bond_To;
    Fl_Input *New_Bond_Dia;
    Fl_Input *New_Bond_Min;
    Fl_Input *New_Bond_Max;
    Fl_Input *New_Bond_Dashes;
    Fl_Output *BondInstr;
    Fl_Output *BondInstr1;
    Fl_Output *BondInstr2;
    Flu_Combo_List *New_Bond_Color;
    Fl_Radio_Button *New_Bond_Style;
    Fl_Button *New_Bond_Add;
    Fl_Button *New_Bond_Remove;
    Flu_Combo_List *Frame_No;
};

class LonePairParam
{				// class to hold lone-pair screen parameters
  public:
    Fl_Window * LonePair_Edit_Window;
    Fl_Text_Editor *LonePair_Edit;
    Fl_Text_Buffer *LonePairBuffer;
    Flu_Combo_List *LonePair_Combo;
    Fl_Input *Number;
    Fl_Input *Height;
    Fl_Input *Radius1;
    Fl_Input *Radius2;
    Fl_Output *LonePairInst;
    Fl_Output *LonePairInst1;
    Fl_Output *LonePairInst2;
    Flu_Combo_List *LonePair_Color;
    Fl_Button *LonePair_Add;
    Fl_Button *LonePair_Remove;
    Flu_Combo_List *Frame_No;
};

class MapsParam
{				// class to hold map parameters
  public:
    Fl_Window * Maps_Edit_Window;
    Fl_Text_Buffer *MapsBuffer;
    Fl_Button *Map_Browse;
    Fl_Button *Map_Info;
    Fl_Output *MapsInstr;
    Fl_Output *MapsInstr1;
    Fl_Output *MapsInstr2;
    Fl_Input *Level;
    Fl_Input *Step;
    Fl_Input *Top;
    Flu_Combo_List *Type;
    Flu_Combo_List *Color;
    Fl_Output *Filename;
    Fl_Button *MapCalc;
    Flu_Combo_List *MapCalcType;
    Flu_Combo_List *MapType;
    Fl_Input *XMin;
    Fl_Input *YMin;
    Fl_Input *ZMin;
    Fl_Input *XMax;
    Fl_Input *YMax;
    Fl_Input *ZMax;
    Fl_Input *Resolution;
    Flu_Spinner *X4;
    Flu_Spinner *X5;
    Flu_Spinner *X6;
    Fl_Button *New_Map_Add;
    Fl_Button *New_Map_Remove;
    Fl_Button *Add_Button;
    Fl_Button *Remove_Button;
    Flu_Combo_List *Frame_No;
};

class SliceParam
{				// class to hold map slice parameters
  public:
    Fl_Window * Slice_Edit_Window;
    Fl_Input *New_x;
    Fl_Input *New_y;
    Fl_Input *New_z;
    Fl_Input *New_nx;
    Fl_Input *New_ny;
    Fl_Input *New_nz;
    Fl_Input *New_type;
    Fl_Check_Button *Legend;
    Flu_Combo_List *Frame_No;
};

class AutomationParam
{				// class to hold automation parameters
  public:
    Fl_Window * Automation_Edit_Window;
    Fl_Input *t_start;
    Fl_Input *t_end;
    Fl_Input *t_step;
    Fl_Input *width;
    Fl_Input *height;
    Fl_Input *fps;
    Fl_Input *POV_Filename;
    Fl_Radio_Button *NoMovie;
    Fl_Radio_Button *Mencoder;
    Fl_Radio_Button *Ffmpeg;
    Fl_Radio_Button *FfmpegG;
    Fl_Check_Button *keeptemps;
    Fl_Button *Go;
    Fl_Button *Abort;
    Fl_Button *Close;
};

class SurfParam
{				// class to hold map parameters
  public:
    Fl_Window * Surfaces_Edit_Window;
    Fl_Text_Buffer *AimSurfBuffer;
    Fl_Output *AimSurfInstr;
    Fl_Output *AimSurfInstr1;
    Fl_Output *AimSurfInstr2;
    Flu_Combo_List *AimSurfType;
    Flu_Combo_List *AimSurfColor;
    Flu_Combo_List *AimSurf_Combo;
    Fl_Input *AimFile;
    Flu_Combo_List *SurfType;
    Flu_Combo_List *SurfColor;
    Fl_Input *Probe;
    Fl_Input *GridX;
    Fl_Input *GridY;
    Fl_Input *GridZ;
    Fl_Text_Buffer *SurfBuffer;
    Fl_Output *SurfInstr;
    Fl_Output *SurfInstr1;
    Fl_Output *SurfInstr2;
    Flu_Combo_List *Surf_Combo;
    Fl_Input *Radius;
    Fl_Button *New_AimSurf_Add;
    Fl_Button *New_AimSurf_Remove;
    Fl_Button *Add_Button;
    Fl_Button *Remove_Button;
    Fl_Button *Add_Button2;
    Fl_Button *Remove_Button2;
    Flu_Combo_List *Frame_No;
};

class ModParam
{				// class to hold modulation parameters
  public:
    Fl_Window * Mods_Edit_Window;
    Flu_Spinner *Mod_t0;
    Flu_Spinner *Mod_t1;
    Flu_Spinner *Mod_t2;
    Fl_Check_Button *Mod_average;
    Fl_Text_Editor *Occ_Edit;
    Fl_Text_Buffer *Occ_Buffer;
    Fl_Output *Occ_Instr;
    Fl_Output *Occ_Instr1;
    Flu_Combo_List *Occ_Combo;
    Fl_Input *New_Occ_Avg;
    Fl_Input *New_Occ_Min;
    Fl_Button *New_Occ_Add;
    Fl_Button *New_Occ_Remove;
    int saved_avg;
    float saved_t[3];
};

class PolyParam
{				// class to hold polyhedra screen parameters
  public:
    Fl_Window * Polyhedra_Edit_Window;
    Fl_Text_Editor *Polyhedra_Edit;
    Fl_Text_Buffer *PolyhedraBuffer;
    Fl_Text_Buffer *Polyhedra_Output_Buffer;
    Flu_Combo_List *Polyhedra_Combo;
    Fl_Input *New_Polyhedra_From;
    Fl_Input *New_Polyhedra_To;
    Fl_Input *New_Polyhedra_Min;
    Fl_Input *New_Polyhedra_Max;
    Fl_Input *New_Polyhedra_Transp;
    Flu_Combo_List *New_Polyhedra_Color;
    Fl_Input *Edge_Radius;
    Flu_Combo_List *Edge_Color;
    Fl_Input *Def_Edge_Radius;
    Flu_Combo_List *Def_Edge_Color;
    Fl_Button *New_Polyhedra_Add;
    Fl_Button *New_Polyhedra_Remove;
    Fl_Radio_Button *Polysz;
    Fl_Radio_Button *Polyvert;
    Fl_Radio_Button *Polyshell;
    Fl_Output *PolyInstr;
    Fl_Output *PolyInstr1;
    Fl_Output *PolyInstr2;
    Fl_Radio_Button *Plane;
    Flu_Combo_List *Frame_No;
};

class SphereParam
{				// class to hold sphere screen parameters
  public:
    Fl_Window * Sphere_Edit_Window;
    Fl_Text_Editor *Sphere_Edit;
    Fl_Text_Buffer *SphereBuffer;
    Fl_Output *SphereInstr;
    Fl_Output *SphereInstr1;
    Fl_Output *SphereInstr2;
    Fl_Text_Buffer *Sphere_Output_Buffer;
    Flu_Combo_List *Sphere_Combo;
    Fl_Input *New_Sphere_Size;
    Flu_Combo_List *New_Sphere_Color;
    Fl_Button *New_Sphere_Add;
    Fl_Button *New_Sphere_Remove;
    Fl_Button *New_Sphere_Convert;
    Flu_Combo_List *Frame_No;
};

class ConfigParm
{				// class to hold POV configuration parameters
  public:
    Fl_Window * ConfigWindow;
    Fl_Input *POVOptions;
    Fl_Input *POVPath;
    Fl_Input *POVIncludePath;
    Fl_Input *POVDefaultFinish;
    Fl_Check_Button *Stereo;
    Fl_Check_Button *StereoMesh;
    Fl_Check_Button *CrossEyed;
    Fl_Input *Stereo_Base;
    Fl_Input *MencoderPath;
    Fl_Input *FFmpegPath;
};

class ConfigMiscParm
{				// class to hold miscellaneous configuration parameters
  public:
    Fl_Window * MiscConfigWindow;
    Fl_Check_Button *LoadLast;
    Fl_Check_Button *AutoLabel;
    Fl_Check_Button *doVrml;
    Fl_Check_Button *doPOV;
    Fl_Check_Button *doAsy;
};

class ConfigMSMSParm
{				// class to hold MSMS configuration
  public:
    Fl_Window * MSMSConfigWindow;
    Fl_Input *MSMSPath;
};

class OmitParam
{				// class to hold omit parameters
  public:
    int nomits;			// number of omits
    int omit1[1000];
    int omit2[1000];
};

class SlabParam
{				// class to hold slab edit parameters
  public:
    Fl_Window * SlabWindow;
    Fl_Input *Slab_A;
    Fl_Input *Slab_B;
    Fl_Input *Slab_C;
    Fl_Input *Slab_Alpha;
    Fl_Input *Slab_Beta;
    Fl_Input *Slab_Gamma;
    Fl_Input *Slab_Off_X;
    Fl_Input *Slab_Off_Y;
    Fl_Input *Slab_Off_Z;
    Fl_Input *Slab_Rot_X;
    Fl_Input *Slab_Rot_Y;
    Fl_Input *Slab_Rot_Z;
    Fl_Choice *Slab_Mode;
};

class ArrowParam
{				// class to hold arrow parameters
  public:
    Fl_Window * ArrowWindow;
    Fl_Input *Px;
    Fl_Input *Py;
    Fl_Input *Pz;
    Fl_Input *Cx;
    Fl_Input *Cy;
    Fl_Input *Cz;
    Fl_Input *Aa;
    Fl_Input *Ba;
    Fl_Input *Ca;
    Fl_Input *Ab;
    Fl_Input *Bb;
    Fl_Input *Cb;
    Fl_Input *Ac;
    Fl_Input *Bc;
    Fl_Input *Cc;
    Fl_Input *Length;
    Fl_Input *Diameter;
    Fl_Output *ArrowInstr;
    Fl_Output *ArrowInstr1;
    Fl_Output *ArrowInstr2;
    Flu_Combo_List *Color;
    Fl_Text_Buffer *ArrowBuffer;
    Fl_Button *AddButton;
    Fl_Button *RemoveButton;
    Flu_Combo_List *Frame_No;
};


#endif
