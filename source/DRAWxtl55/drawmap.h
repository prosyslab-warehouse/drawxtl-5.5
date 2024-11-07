// $Id: drawmap.h 1034 2010-10-17 13:16:30Z martin $
//
#ifndef DRAWMAP_h
#define DRAWMAP_h
#include "mpVector.h"
#include "MC.h"

/* Header file to support Fourier map */

/* Structure for storing vectors */
/*
typedef struct {
	float x, y, z;
} mpVector;

typedef struct {
  float x, y, z;    // orthonormal coordinates
  float val;        // rho value
} mp4Vector;

// struct for storing triangle information - 3 vertices and 3 normal vectors for each vertex
typedef struct {
	mpVector p[3];
	mpVector norm[3];
} TRIANGLE;
*/

TRIANGLE *MC_c (int ncellsX, int ncellsY, int ncellsZ,
		float gradFactorX, float gradFactorY, float gradFactorZ,
		float minValue, mp4Vector * points, int &numTriangles);

void MCfree_c (TRIANGLE * trianglePt);	/* use this to free memory created in C++  MCfree_c (pTriangles);  */

void generate_map (float minValue, int Solid, char *Color, char *BackColor);

void generate_slice (void);

void colorramp (float, float*, float* , float*);

void MapLegend (void);

int ContourFacet (TRIANGLE tri, mpVector *p1, mpVector *p2);
#endif
