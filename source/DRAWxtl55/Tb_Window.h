// $Id: Tb_Window.h 900 2009-08-13 20:00:45Z larry $
//
/* A GL window with a trackball interface.  This is a subclass of Fl_Gl_Window
that only implements the handle() method to keep track of mouse motions.
You subclass Tb_Window and implement draw(), and call the transform() method
to orient the scene. */

#ifndef Tb_Window_h
#define Tb_Window_h

#include <FL/Fl.H>
#include <FL/glut.H>
#include <FL/Fl_Window.H>
#include <FL/gl.h>

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/* Handy typedefs for declaring storage of 2- & 3-vectors. */

typedef float XY[2];

typedef float XYZ[3];


/* Routines for the quaternion representation of rotations. [quat.c] */

typedef struct
{
    float s;
    float v[3];
} QUAT;


float *vnew (float x, float y, float z);

void vset (float *v, float x, float y, float z);

void vcopy (float *vsrc, float *vdst);

void vprint (float *v);

void vzero (float *v);

void vnormalize (float *v);

float vlength (float *v);

void vscale (float *v, float f);

void vmult (float *src1, float *src2, float *dst);

void vadd (float *src1, float *src2, float *dst);

void vsub (float *src1, float *src2, float *dst);

float vdot (float *v1, float *v2);

void vcross (float *v1, float *v2, float *cross);

void axis_to_quaternion (float *axis, float theta, QUAT * quat);

void qmult (QUAT * q1, QUAT * q2, QUAT * dest);

void qnormalize (QUAT * q);

void quaternion_to_rotmatrix (QUAT * q, float *m);


class Tb_Window:public Fl_Gl_Window
{
  public:
    Tb_Window (int x, int y, int w, int h, const char *l = 0)
  :	Fl_Gl_Window (x, y, w, h, l) {
	glutInitWindowSize (w, h);
	glutInitWindowPosition (x, y);
	init ();
    }
    Tb_Window (int w, int h, const char *l = 0)
  :	Fl_Gl_Window (w, h, l) {
	glutInitWindowSize (w, h);
	init ();
    }

    ~Tb_Window ();

    /* Set the trackball size and mouse scale. */

    void tbsize (float f)
    {
	Tbsize = f;
    }
    void mscale (float f)
    {
	Mscale = 30.f * f;
    }

    // Set the initial trackball origin, translation and rotation.

    void origin (float x, float y, float z);

    void translate (float *v);

    void rotate (float *axis, float theta);

    // Call this from draw() to calculate the view transformation.

    void calculate (float a[16]);

    /* Call this when the geometry changes to force a redraw(),
       or just call redraw(). */

    void set_changed ()
    {
	changed = TRUE;
    }

    /* Call idle_redraw(TRUE) to force a redraw() every time idle() is
       called.  It is initially FALSE. */

    void idle_redraw (int b)
    {
	idle_redraw_ = b;
    }
//      XYZ Trans;              // total translation

  private:
    void init ();

    int handle (int);

    void idle ();

    QUAT Spin;			// quaternion incremental rotation

    float Tbsize;		// fraction of window to fill with trackball

    float Mscale;		// mouse movement multiplier

    ulong Event;		// mouse button and modifier keys

    int Mx;			// place where mouse is now

    int My;

    int Oldx;			// place where mouse used to be

    int Oldy;

    int spinning;		// if the user did a SPIN

    int changed;		// if the view changed

    int idle_redraw_;		// call redraw() in any case

    void trackball ();

    void window_to_tb (float mx, float my, float *x, float *y);

    float tb_project_to_sphere (float, float, float);
};

#endif
