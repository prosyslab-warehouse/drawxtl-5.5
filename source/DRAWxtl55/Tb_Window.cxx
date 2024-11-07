// $Id: Tb_Window.cxx 1079 2010-11-10 23:03:23Z martin $
//
// Coded using the FLTK 1.1.6 widget set
//
//     Larry W. Finger, Martin Kroeker and Brian Toby
//
// An openGL window with a trackball interface.
//
// routines contained within this file:
//
// Tb_Window::translate - set the initial translation
// Tb_Window::rotate - set the initial rotation
// Tb_Window::calculate - performs the translation and calculates the rotation
// Tb_Window::init - initialize the trackball
// Tb_Window::~Tb_Window - destructor
// Tb_Window::handle - keyboard and mouse event handler
// Tb_Window::idle - routine to perform actual work detected by handle method
// Tb_Window::trackball - simulates a trackball using the mouse
// Tb_Window::window_to_tb - maps mouse mx, my to -1,1 range
// Tb_Window::tb_project_to_sphere - projects x,y onto sphere or hyperboloid
// axis_to_quaternion - converts axis and angle to quaternion rotation matrix
// qmult - multiplies two quaternions
// qnormalize - normalizes greatest component of quaternion
// quaternion_to_rotmatrix - convert quaternion to a rotation matrix
// vnew - returns a vector of length 3 from components
// vset - loads components in a vector of length 2
// vcopy - copies a vector of length 3
// vzero - clears a vector of length 3
// vnormalize - normalizes (sets to length = 1) a 3-vector
// vlength - returns the length of a 3-vector
// vscale - multiplies a 3-vector by a constant
// vmult - multiplies the components of two 3-vectors
// vadd - adds two 3-vectors
// vsub - subtracts one 3-vector from another
// vdot - returns dot product of two 3-vectors
// vcross - calculates cross product of two 3-vectors
// Update_Cursor_List - update atom list used to calculate bond, etc.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Tb_Window.h"
#include "DRAWxtlViewUI.h"
#include "draw_ext.h"

static int c_pick = 0;

#include "DRAWxtl_proto.h"

#define SPIN (FL_BUTTON1)
#define EPS 1e-7
#ifdef __APPLE__
#define ZOOM (FL_META | FL_BUTTON1)	//  Command + MOUSE (Mac)
#define PAN (FL_ALT | FL_BUTTON1)	//  Alt + MOUSE (Mac)
#else
#define ZOOM (FL_BUTTON3)	// just right mouse on other platforms
#define PAN (FL_BUTTON2)	// middle button
#endif
#define BUTTON_MASK (PAN | ZOOM | SPIN )

/* Set the trackball's initial center of rotation -- this can be changed
using the Shift- PAN or ZOOM buttons. */

/* Set the initial translation -- this is changed by PAN or ZOOM. */

void
Tb_Window::translate (float *v)
{
    vadd (drvui->Trans, v, drvui->Trans);
}

/* Set the initial rotation of the trackball (and the scene) to theta
degrees around the given axis -- this is changed by the SPIN button. */

void
Tb_Window::rotate (float *axis, float theta)
{
    QUAT q;

    axis_to_quaternion (axis, theta * (float) RAD, &q);
    qmult (&q, &Rotq, &Rotq);
}

/* Call this in the draw() method just before rendering the scene.  This
performs the translation and calculates the rotation specified by the trackball.  It is
the caller's responsibility to push and pop the gl matrix. */

void
Tb_Window::calculate (float m[16])
{

    /* Translate (rotated coordinates). */

    glTranslatef (drvui->Trans[0], drvui->Trans[1], drvui->Trans[2]);

    /* Calculate Rotation. */

    quaternion_to_rotmatrix (&Rotq, m);

}


/* Initialize the trackball interface. */

void
Tb_Window::init (void)
{

    /* Set the default trackball size and mouse scale factor. */

    tbsize (0.5f);
    mscale (0.25f);

    /* Default is no rotation or translation. */

    Spin.s = 1.0;
    vzero (Spin.v);
    spinning = FALSE;
    Rotq.s = 1.0;
    vzero (Rotq.v);
//      vzero(Trans);                    // Trans moved to drvui, which is not defined when Tb_Window

    /* Initialize the other instance variables.  The rest are set the
       first time the user clicks the mouse. */

    Event = 0;
    idle_redraw_ = FALSE;
    changed = TRUE;
}

Tb_Window::~Tb_Window ()
{
}

/* The handle() method just records the events, and the idle callback
does the actual work. */

int
Tb_Window::handle (int e)
{
    int key;

    int state;

    int Sense = 1;

    int atomno;

    float dx, dy, v[3];

    switch (e) {
    case FL_ENTER:
	Tb_Window::focus (this);
	return (1);
    case FL_KEYDOWN:
	key = Fl::event_key ();
	state = Fl::event_state ();
	if (((key >= 'a') && (key <= 'z')) && (state & FL_SHIFT || state & FL_CAPS_LOCK))
	    key -= 'a' - 'A';	/* change to UC */
	switch (key) {
	case 'C':
	case 'c':
	    if (cur_show) {
		drvui->cur_step *= 0.2f;
		if (drvui->cur_step < 0.01) {
		    cur_show = 0;
		    drvui->cur_step = 0.5;
		    drvui->Cursor_pos->value ("");
		    cur_atom[0] = cur_atom[1] = cur_atom[2] = cur_atom[3] = 0;
		    strcpy (cur_name[0], "");
		}
	    } else {
		cur_show = 1;
		drvui->cur_step = 0.5f;
		cur_atom[0] = cur_atom[1] = cur_atom[2] = cur_atom[3] = 0;
		strcpy (cur_name[0], "");
	    }
	    Tb_Window::idle ();
	    break;
	case 'x':
	    move_cursor (0, 1.0f);
	    break;
	case 'y':
	    move_cursor (1, 1.0f);
	    break;
	case 'z':
	    move_cursor (2, 1.0f);
	    break;
	case 'X':
	    move_cursor (0, -1.0f);
	    break;
	case 'Y':
	    move_cursor (1, -1.0f);
	    break;
	case 'Z':
	    move_cursor (2, -1.0f);
	    break;
	case 'M':
	    Sense = -1;
	case 'm':
	    atomno = Maximize_rho (Sense);	/* find local maximum (minimum if Sense < 0) */
	    Update_Cursor_List (atomno);	// add this position to cursor list
	    break;
	case 's':
	    if (cur_atom[2] <1) break;
	    Add_mapslice(2);
	    break;
	case 'S':
	    if (cur_atom[2] <1) break;
	    Add_mapslice(3);
	    break;
	case 'L':
	case 'l':
	    unsigned int j, nn, theatom;

	    if (strlen (cur_name[0]) == 0)
		break;
	    if (cur_atom[0] < 1)
		break;
	    theatom = 0;

	    if (cur_atom[1] > 0)
		theatom = 1;
	    if (cur_atom[2] > 0)
		theatom = 2;
	    if (cur_atom[3] > 0)
		theatom = 3;
	    if (drvui->max_frame > 1)
		drvui->labels[drvui->nlabel].label_fn = drvui->frame_no - 1;
	    else
		drvui->labels[drvui->nlabel].label_fn = drvui->frame_no;
	    drvui->labels[drvui->nlabel].label_x[0] = cur_cen[0];
	    drvui->labels[drvui->nlabel].label_x[1] = cur_cen[1];
	    drvui->labels[drvui->nlabel].label_x[2] = cur_cen[2];
	    nn = 0;
	    for (j = 0; j < strlen (cur_name[theatom]); j++) {
		drvui->labels[drvui->nlabel].label_label[nn++] = cur_name[theatom][j];
	    }
	    drvui->labels[drvui->nlabel].label_label[nn] = '\0';
	    drvui->nlabel++;
	    check_dynamic_storage ();
	    Update_Str (0);
	    break;
	case 'B':
	case 'b':
	    int first, second;

	    double curdist;

	    if (dist12 < 0.01)
		break;

//fprintf(stderr,"atoms in cache : %d %d %d\n",cur_atom[0],cur_atom[1],cur_atom[2]);      
	    if (cur_atom[0] < 1)
		break;
	    if (cur_atom[1] < 1)
		break;
	    first = cur_atom[0];
	    second = cur_atom[1];
	    curdist = dist12;

	    if (cur_atom[2] > 0) {
		first = cur_atom[2];
		second = cur_atom[1];
		curdist = dist23;
	    }
	    if (cur_atom[3] > 0) {
		first = cur_atom[3];
		second = cur_atom[2];
		curdist = dist34;
	    }


	    if (drvui->max_frame > 1)
		drvui->labels[drvui->nlabel].label_fn = drvui->frame_no - 1;
	    else
		drvui->labels[drvui->nlabel].label_fn = drvui->frame_no;
	    drvui->labels[drvui->nlabel].label_x[0] =
		(o_vert[3 * second] + o_vert[3 * first]) / 2.0f;
	    drvui->labels[drvui->nlabel].label_x[1] =
		(o_vert[3 * second + 1] + o_vert[3 * first + 1]) / 2.0f;
	    drvui->labels[drvui->nlabel].label_x[2] =
		(o_vert[3 * second + 2] + o_vert[3 * first + 2]) / 2.0f;
	    sprintf (drvui->labels[drvui->nlabel].label_label, "%.3f", curdist);
	    drvui->nlabel++;
	    check_dynamic_storage ();
	    Update_Str (0);
	    break;
	case 'P':
	case 'p':
	    atomno = find_proj_atom (Fl::event_x (), Fl::event_y ());
	    if (atomno >= 0)
		Update_Cursor_List (atomno);	// add this position to cursor list
	    break;
	case 'a':
	case 'A':		// find nearest atom, calculate distance
	    atomno = find_atom ();
	    Update_Cursor_List (atomno);	// add this position to cursor list
	    break;
	case FL_Down:
	    dy = -10.0f * Mscale / (float) h ();
	    vset (v, 0., dy, 0.);
	    vadd (v, drvui->Trans, drvui->Trans);
	    break;
	case FL_Up:
	    dy = 10.0f * Mscale / (float) h ();
	    vset (v, 0., dy, 0.);
	    vadd (v, drvui->Trans, drvui->Trans);
	    break;
	case FL_Left:
	    dx = -10.0f * Mscale / (float) h ();
	    vset (v, dx, 0., 0.);
	    vadd (v, drvui->Trans, drvui->Trans);
	    break;
	case FL_Right:
	    dx = 10.0f * Mscale / (float) h ();
	    vset (v, dx, 0., 0.);
	    vadd (v, drvui->Trans, drvui->Trans);
	    break;
	case FL_Home:
	    vset (drvui->Trans, 0., 0., 0.);
	    break;

	}
	Tb_Window::idle ();
	Fl::redraw ();
	return (1);
    case FL_FOCUS:
	return (1);
    case FL_UNFOCUS:
	return (1);
    case FL_PUSH:
	Oldx = Mx = Fl::event_x ();
	Oldy = My = Fl::event_y ();
	Event = Fl::event_state ();
	Tb_Window::idle ();
	return (1);

    case FL_RELEASE:
	Event = Fl::event_state ();
	Tb_Window::idle ();
	return (1);

    case FL_DRAG:
	Mx = Fl::event_x ();
	My = Fl::event_y ();
	Event = Fl::event_state ();
	Tb_Window::idle ();
	return (1);
    case FL_SHOW:
	return Fl_Gl_Window::handle (e);
    }
    return (0);
}

void
Tb_Window::idle (void)
{
    int m;

    float *v, dx, dy;

    XYZ trans;

    float ratio;

    float factor;

    /* Shift- PAN or ZOOM translates the unrotated scene. */

    v = drvui->Trans;

    /* If the user drags the mouse, Spin or [O]Trans are updated.  If no
       mouse buttons are pressed, just keep on applying the previous spin
       rotation, over and over. */

    m = Event & BUTTON_MASK;
    if (m) {
	dx = (float) (Mx - Oldx) * Mscale / (float) w ();
	dy = (float) (Oldy - My) * Mscale / (float) h ();
	if (m == PAN) {
	    vset (trans, dx, dy, 0.);
	    vadd (trans, v, drvui->Trans);
	} else if (m == ZOOM) {
	    if (M_cameras == 1) {
		glLoadIdentity ();
		vset (trans, 0., 0., dx + dy);
		vadd (trans, v, drvui->Trans);
	    } else {
		gl_size -= (dx + dy) / 2.0f;
		glMatrixMode (GL_PROJECTION);
		ratio = 1.0f * w () / h ();
		if (w () <= h ())
		    glOrtho (-gl_size, gl_size, -gl_size / ratio, gl_size / ratio,
			     -10000., 10000.);
		else
		    glOrtho (-gl_size * ratio, gl_size * ratio, -gl_size, gl_size,
			     -10000., 10000.);
		glMatrixMode (GL_MODELVIEW);
	    }
	} else if (m == SPIN) {
	    if (Event & FL_CTRL) {
		start_picking (Mx, My, w (), h ());
		return;
	    } else if (Event & FL_SHIFT && slabmode == 2) {
		if (Mx == Oldx && My == Oldy) {
		    Fl_Gl_Window::redraw_overlay ();
		    c_pick = pick_box (Mx, My, w (), h ());
		} else {
		    switch (c_pick) {
		    case 0:
		    default:
			break;
		    case 1:
			drvui->slab_off[0] += dx;
			drvui->slab_off[1] += dy;
			Omit->nomits = 0;
			break;
		    case 2:
			drvui->slab_con[2] += dx;
			drvui->slab_con[2] += dy;
			Omit->nomits = 0;
			break;
		    case 3:
			drvui->slab_con[1] += dx;
			Omit->nomits = 0;
			drvui->slab_con[1] += dy;
			break;
		    case 4:
			drvui->slab_con[0] += dx;
			drvui->slab_con[0] += dy;
			Omit->nomits = 0;
			break;
		    }
		    Oldx = Mx;
		    Oldy = My;
		    generate_slab ();
		    Fl_Gl_Window::redraw_overlay ();
		    Fl::redraw ();
		}
		return;
	    } else if (Event & FL_SHIFT) {
		if (drvui->nlabel == 1)
		    return;
		if (Mx == Oldx && My == Oldy) {
		    c_pick = pick_label (Mx, My, w (), h ());
		} else {
		    GLdouble objx, objy, objz, objx1, objy1, objz1;

		    if (c_pick == 0)
			return;
		    factor = (float) w () * Scale * 0.01f;
		    if (M_cameras == 0)
			factor = gl_size * Scale * 0.01f;
		    gluUnProject (Oldx, Oldy, 0, modelMatrix, projMatrix, viewport,
				  &objx, &objy, &objz);
		    gluUnProject (Oldx + dx, Oldy + dy, 0, modelMatrix, projMatrix,
				  viewport, &objx1, &objy1, &objz1);
		    if (c_pick != drvui->triple[0]) {
			drvui->labels[c_pick].label_x[0] +=
			    factor * (float) (objx1 - objx);
			drvui->labels[c_pick].label_x[1] +=
			    factor * (float) (objy1 - objy);
			drvui->labels[c_pick].label_x[2] +=
			    factor * (float) (objz1 - objz);
		    } else {
			offset[0] += 5.0f * factor * (float) (objx1 - objx);
			offset[1] += 5.0f * factor * (float) (objy1 - objy);
			offset[2] += 5.0f * factor * (float) (objz1 - objz);
		    }
		    Oldx = Mx;
		    Oldy = My;
		    drvui->Str_File_Changed = 1;
		    Fl::redraw ();
		}
		return;
//            } else if (Event & FL_SHIFT && cur_show) {
//fprintf(stderr,"pick/moveto\n");
//                moveto_atom(Mx,My,w(),h(),drvui->Trans);
//                return;
	    } else {
/* trackball() updates Spin, which we then apply. */

		trackball ();
		qmult (&Spin, &Rotq, &Rotq);
	    }
	}
	Oldx = Mx;
	Oldy = My;
	changed = TRUE;
    } else {
	if (drvui->slab_con[0] > 0. && (Event & FL_SHIFT)) {	/* slab changed */
	    Update_Str (0);
	    drvui->crystalDL = glGenLists (1);
	    Generate_Drawing (0);	/* regenerate drawing */
	}
    }

    /* Tell the window to call draw(). */

    if (idle_redraw_ || spinning || changed) {
	redraw ();
    }
    changed = FALSE;
}

/* Implementation of a virtual trackball.  Original code by Gavin Bell, lots
of ideas from Thant Tessman and the August '88 issue of SigGraph's "Computer
Graphics," pp. 121-129. */

/* Ok, simulate a trackball.  Project the mouse positions onto the virtual
trackball, then figure out the axis of rotation, which is the cross product
of O P1 and O P2 (O is the center of the ball).  Note: This is a deformed
trackball -- it's a trackball in the center, but is deformed into a
hyperbolic solid of rotation away from the center. */

void
Tb_Window::trackball (void)
{
    float p1x, p1y, p2x, p2y;

    float theta, t;

    XYZ p1, p2, axis, d;

    QUAT *q = &Spin;

    if (Mx == Oldx && My == Oldy) {

	/* Zero rotation. */

	q->s = 1.;
	vzero (q->v);
	spinning = FALSE;
	return;
    }
    spinning = FALSE;		// This is changed by LWF from TRUE to FALSE
    // Scale the old and new mouse positions to (-1, 1).

    window_to_tb ((float) Oldx, (float) Oldy, &p1x, &p1y);
    window_to_tb ((float) Mx, (float) My, &p2x, &p2y);

    //First, figure out z-coordinates for projection of P1 and P2 to
    //the deformed sphere. 

    vset (p1, p1x, p1y, tb_project_to_sphere (Tbsize, p1x, p1y));
    vset (p2, p2x, p2y, tb_project_to_sphere (Tbsize, p2x, p2y));

    // Now the axis of rotation.

    vcross (p1, p2, axis);

    // If this is a Shift-SPIN, force it to be either the X or Y axis,
    // for better control

    if (Event & FL_SHIFT) {
	axis[2] = 0.;
	if (fabs (axis[0]) > fabs (axis[1])) {
	    axis[1] = 0.;
	} else {
	    axis[0] = 0.;
	}
    }
    // Figure out how much to rotate around that axis. 

    vsub (p1, p2, d);
    t = vlength (d) * Mscale / 20.0f;
    if (t > 1.0)
	t = 1.0;
    if (t < -1.0)
	t = -1.0;
    theta = (float) asin (t);

    // Return a rotation quaternion

    axis_to_quaternion (axis, theta, q);
}

/* Map mouse click mx, my to a more convenient (-1.0, 1.0) range, based on
window size. */

void
Tb_Window::window_to_tb (float mx, float my, float *x, float *y)
{
    *x = (2.0f * mx) / (float) w () - 1.f;
    *y = 1.0f - (2.0f * my) / (float) h ();
}

/* Project an x, y pair onto a sphere of radius r OR a hyperbolic sheet if we
are away from the center of the sphere.  */

float
Tb_Window::tb_project_to_sphere (float r, float x, float y)
{
    float d, t, z;

    d = (float) sqrt (x * x + y * y);
    if (d < r * 0.70710678118654752440) {	/* inside sphere */
	z = (float) sqrt (r * r - d * d);
    } else {			/* on hyperbola */
	t = r / (float) sqrt (2.0);
	z = t * t / d;
    }
    return (z);
}

/* Given an axis and an angle, compute a rotation quaternion. */

void
axis_to_quaternion (float *axis, float theta, QUAT * quat)
{
    quat->s = (float) cos (theta / 2.0);
    vcopy (axis, quat->v);
    vnormalize (quat->v);
    vscale (quat->v, (float) sin (theta / 2.0));
}

/* Given two rotations, q1 and q2, expressed as quaternions,
quaternion multiplication yields the equivalent product rotation.
This routine also normalizes the result every COUNT times it is
called, to keep error from creeping in. */

#define COUNT 100

void
qmult (QUAT * q1, QUAT * q2, QUAT * dest)
{
    XYZ v1, v2;

    QUAT q;

    static int count = 0;

    vcopy (q1->v, v1);
    vcopy (q2->v, v2);
    vscale (v1, q2->s);
    vscale (v2, q1->s);
    vcross (q1->v, q2->v, q.v);
    vadd (v1, q.v, q.v);
    vadd (v2, q.v, q.v);
    q.s = q1->s * q2->s - vdot (q1->v, q2->v);

    *dest = q;

    if (++count > COUNT) {
	count = 0;
	qnormalize (dest);
    }
}

/* Normalize greatest component, to avoid problems that occur when the
component we're normalizing gets close to zero (and the other components may
add up to more than 1.0 because of rounding error).  */

void
qnormalize (QUAT * q)
{
    int which, i;

    float gr, v[4];

    vcopy (q->v, v);
    v[3] = q->s;

    which = 0;
    gr = v[0];
    for (i = 1; i < 4; i++) {
	if (fabs (v[i]) > fabs (gr)) {
	    gr = v[i];
	    which = i;
	}
    }

    v[which] = 0.;
    v[which] = (float) sqrt (1.0 - (v[0] * v[0] + v[1] * v[1] +
				    v[2] * v[2] + v[3] * v[3]));

    /* Check to see if we need negative square root. */

    if (gr < 0.0) {
	v[which] = -v[which];
    }
    vcopy (v, q->v);
    q->s = v[3];
}

/* Build a rotation matrix from a rotation quaternion. */

void
quaternion_to_rotmatrix (QUAT * q, float *m)
{
    float w, x, y, z;

    w = q->s;
    x = q->v[0];
    y = q->v[1];
    z = q->v[2];

    /* Column major order, used to left multiply points, i.e., P'=M*P */

    m[0] = 1.f - 2.f * (y * y + z * z);
    m[1] = 2.f * (x * y + w * z);
    m[2] = 2.f * (x * z - w * y);
    m[3] = 0.;

    m[4] = 2.f * (x * y - w * z);
    m[5] = 1.f - 2.f * (x * x + z * z);
    m[6] = 2.f * (y * z + w * x);
    m[7] = 0.;

    m[8] = 2.f * (x * z + w * y);
    m[9] = 2.f * (y * z - w * x);
    m[10] = 1.f - 2.f * (x * x + y * y);
    m[11] = 0.;

    m[12] = 0.;
    m[13] = 0.;
    m[14] = 0.;
    m[15] = 1.;
}

float *
vnew (float x, float y, float z)
{
    float *v;

    v = (float *) malloc (3 * sizeof (float));	//(float, 3);
    v[0] = x;
    v[1] = y;
    v[2] = z;

    return (v);
}

void
vset (float *v, float x, float y, float z)
{
    v[0] = x;
    v[1] = y;
    v[2] = z;
}

void
vcopy (float *vsrc, float *vdst)
{
    vdst[0] = vsrc[0];
    vdst[1] = vsrc[1];
    vdst[2] = vsrc[2];
}

void
vzero (float *v)
{
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
}

void
vnormalize (float *v)
{
#if defined(_IEEE) || defined(_IEEE_FP)
    vscale (v, 1.0 / vlength (v));
#else
    float l;

    l = vlength (v);
    if (l > EPS) {
	vscale (v, 1.0f / l);
    }
#endif
}

float
vlength (float *v)
{
    return (float) (sqrt (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
}

void
vscale (float *v, float f)
{
    v[0] *= f;
    v[1] *= f;
    v[2] *= f;
}

void
vmult (float *src1, float *src2, float *dst)
{
    dst[0] = src1[0] * src2[0];
    dst[1] = src1[1] * src2[1];
    dst[2] = src1[2] * src2[2];
}

void
vadd (float *src1, float *src2, float *dst)
{
    dst[0] = src1[0] + src2[0];
    dst[1] = src1[1] + src2[1];
    dst[2] = src1[2] + src2[2];
}

void
vsub (float *src1, float *src2, float *dst)
{
    dst[0] = src1[0] - src2[0];
    dst[1] = src1[1] - src2[1];
    dst[2] = src1[2] - src2[2];
}

float
vdot (float *v1, float *v2)
{
    return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
}

void
vcross (float *v1, float *v2, float *cross)
{
    XYZ temp;

    /* x        y       z         z       y */
    temp[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
    temp[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
    temp[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
    vcopy (temp, cross);
}

void
Update_Cursor_List (int i)
{
    float dot;

    char atnum[20];

    int n;

    char cur_name_t[10];

    n = drvui->orig_atom_no[i];
    strcpy (cur_name_t, "");

    if (drvui->cur_reset > -1 ) {
	if (cur_atom[drvui->cur_reset] > 0) {
	    cur_atom[0] = cur_atom[1] = cur_atom[2] = cur_atom[3] = 0;
	    strcpy (cur_name[0], "");
	    strcpy (cur_name[1], "");
	    strcpy (cur_name[2], "");
	    strcpy (cur_name[3], "");
	}
    }

    if (cur_atom[1] > 0)
	dist12 = dist (cur_atom[0], cur_atom[1]);
    if (cur_atom[2] > 0) {
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
    }
    if (n < 0) {
	strcat (cur_name_t, "rho");
    } else {
	strcat (cur_name_t, drvui->atoms[n].atom_l);
	trim_string (cur_name_t, 5);
	sprintf (atnum, "%d", drvui->atoms[n].sv_atom_n);
	strcat (cur_name_t, atnum);
    }
    if (cur_atom[0] <= 0) {
	cur_atom[0] = i;
	strcpy (cur_name[0], cur_name_t);
    } else if (cur_atom[1] <= 0) {
	cur_atom[1] = i;
	strcpy (cur_name[1], cur_name_t);
	dist12 = dist (cur_atom[0], cur_atom[1]);
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

	if (cur_atom[3] > 0) {
	    cur_atom[0] = cur_atom[1];
	    cur_atom[1] = cur_atom[2];
	    cur_atom[2] = cur_atom[3];
	    strcpy (cur_name[0], cur_name[1]);
	    strcpy (cur_name[1], cur_name[2]);
	    strcpy (cur_name[2], cur_name[3]);
	    dist12 = dist23;
	    dist23 = dist34;
	}
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
	vcross (v0, v1, p1);	// perpendicular to plane 1-2-3
	vcross (v2, v1, p2);	// perpendicular to plane 2-3-4
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
	if (vdot (p1, v2) > 0.0f)	// get the sign
	    torsion_ang *= -1.0f;
    }
}
