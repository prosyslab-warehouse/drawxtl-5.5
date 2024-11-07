// Gmsh - Copyright (C) 1997-2009 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <gmsh@geuz.org>.

#ifndef _PIXEL_BUFFER_H_
#define _PIXEL_BUFFER_H_

#include <string.h>
#include <FL/gl.h>

#if defined(HAVE_OSMESA)
#include <GL/osmesa.h>
#endif

class PixelBuffer {
  private:
    int _width, _height, _numComp, _dataSize;
    GLenum _format, _type;
    unsigned char *_pixels;
  public:
     PixelBuffer(int width, int height, GLenum format, GLenum type)
    :_width(width), _height(height), _format(format), _type(type) {
	if (format == GL_RGB) {
	    _numComp = 3;
	} else if (format == GL_RGBA) {
	    _numComp = 4;
	} else {
	    Error_Box("Unknown pixel format: assuming RGB");
	    _format = GL_RGB;
	    _numComp = 3;
	}

	if (type == GL_UNSIGNED_BYTE) {
	    _dataSize = sizeof(unsigned char);
	} else if (type == GL_FLOAT) {
	    _dataSize = sizeof(float);
	} else {
	    Error_Box("Unknown pixel storage type: assuming unsigned byte");
	    _type = GL_UNSIGNED_BYTE;
	    _dataSize = sizeof(unsigned char);
	}
	int n = _numComp * _width * _height * _dataSize;
	_pixels = new unsigned char[n];
	for (int i = 0; i < n; i++)
	    _pixels[i] = 0;
    }
    ~PixelBuffer() {
	delete[]_pixels;
    }
    int getWidth() {
	return _width;
    }
    int getHeight() {
	return _height;
    }
    int getNumComp() {
	return _numComp;
    }
    int getDataSize() {
	return _dataSize;
    }
    GLenum getFormat() {
	return _format;
    }
    GLenum getType() {
	return _type;
    }
    void *getPixels() {
	return (void *) _pixels;
    }
    void copyPixels(int x, int y, PixelBuffer * buffer) {
	if (x + buffer->getWidth() > _width || y + buffer->getHeight() > _height) {
	    Error_Box("Destination pixel buffer too small for holding copy");
	    return;
	}
	if (buffer->getNumComp() != _numComp || buffer->getDataSize() != _dataSize ||
	    buffer->getFormat() != _format || buffer->getType() != _type) {
	    Error_Box("Pixel buffer type mismatch: impossible to copy");
	    return;
	}
	for (int i = 0; i < buffer->getWidth(); i++)
	    for (int j = 0; j < buffer->getHeight(); j++)
		memcpy(_pixels + ((j + y) * _width + (i + x)) * _dataSize * _numComp,
		       (unsigned char *) buffer->getPixels() + (j * buffer->getWidth() +
								i) * _dataSize * _numComp,
		       _dataSize * _numComp);
    }
    void fill(int offscreen) {
	float cpx, cpy, cpz;
	float m[16];
	GLfloat mat_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat light_ambient[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	GLfloat light_specular[] = { 0.1f, 0.1f, 0.1f, 1.0f };
	GLfloat light_diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat mat_shininess[] = { 0.8f };
	GLfloat light_position[] = { 0.0f, 1.0f, 1.0f, 0.0f };
	GLfloat light_direction[] = { 0.0f, -1.0f, 0.0f };

	if (!offscreen) {
	    glMatrixMode(GL_PROJECTION);
	    glLoadIdentity();
	    float ratio=1.0f*(float)_width/(float)_height;
	    if (M_cameras == 0) {
	        if (_width <=_height)
		glOrtho(-gl_size, gl_size, -gl_size/ratio, gl_size/ratio, -10000., 10000.);
		else
		glOrtho(-gl_size*ratio, gl_size*ratio, -gl_size, gl_size, -10000., 10000.);
	    } else {
		gluPerspective(17., ratio, 0.01, 1000.);
	    }
	    glMatrixMode(GL_MODELVIEW);
	    glLoadIdentity();
	    cpx = (POV_Max[0] + POV_Min[0]) / 2.0f;
	    cpy = (POV_Max[1] + POV_Min[1]) / 2.0f;
	    cpz = (POV_Max[2] + POV_Min[2]) / 2.0f;
//    glShadeModel (GL_SMOOTH);
	    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	    glEnable(GL_LIGHTING);
	    glEnable(GL_LIGHT0);
	    glEnable(GL_COLOR_MATERIAL);
	    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	    glEnable(GL_DEPTH_TEST);
	    glDepthFunc(GL_LESS);
	    glClearColor(drvui->glback[0], drvui->glback[1], drvui->glback[2], 0.0f);
	    glClear((GLbitfield) (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT));
	    gluLookAt(cpx, cpy, Scale * .50,	// camera position
		      cpx, cpy, -1.0,	// camera lookat point
		      0.0f, 1.0f, 0.0f);	// camera "up" vector
	    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	    glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, light_direction);
	    glColor3f(0.0f, 0.0f, 0.0f);
	    glPushMatrix();
	    glTranslatef(drvui->Trans[0], drvui->Trans[1], drvui->Trans[2]);
	    quaternion_to_rotmatrix(&Rotq, m);

	    glMultMatrixf(m);
	    glTranslatef(-cpx, -cpy, -cpz);
	    draw_cursor();
	    glCallList(drvui->crystalDL);
	    for (drvui->frame_no = 1; drvui->frame_no <= drvui->max_frame;
		 drvui->frame_no++)
		generate_gl_texts();
	    glPopMatrix();
	    drvui->frame_no = drvui->max_frame;

	    glFinish();
	    glPixelStorei(GL_PACK_ALIGNMENT, 1);
	    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	    glReadPixels(0, 0, _width, _height, _format, _type, (void *) _pixels);
	} else {
#if defined(HAVE_OSMESA)
	    if (_format != GL_RGB && _type != GL_UNSIGNED_BYTE) {
		Error_Box
		    ("Offscreen rendering only implemented for GL_RGB/GL_UNSIGNED_BYTE");
		return;
	    }
	    OSMesaContext ctx = OSMesaCreateContextExt(OSMESA_RGB, 16, 0, 0, NULL);
	    if (!ctx) {
		Error_Box("OSMesaCreateContext failed");
		return;
	    }
	    if (!OSMesaMakeCurrent
		(ctx, (void *) _pixels, GL_UNSIGNED_BYTE, _width, _height)) {
		Error_Box("OSMesaMakeCurrent failed");
	    }
	    DrawCurrentOpenglWindow(false);
	    glFinish();
	    OSMesaDestroyContext(ctx);
#else
	    Error_Box("Gmsh must be compiled with OSMesa to support offscreen rendering");
#endif
	}
    }
};

#endif
