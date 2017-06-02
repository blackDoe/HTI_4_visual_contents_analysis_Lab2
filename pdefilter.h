/*
 *	pdefilter.h
 *	Auteur:		Nicolas ROUGON
 *			Institut Mines-Telecom / Telecom SudParis
 *      Date:		01/06/2015
 *
 *---------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "globals.h"  

//---- Functions
void Sobel(Matrix2D *, Matrix2D *);
void usage(char *);

extern void alloc_Matrix2D(Matrix2D *, int, int, int);
extern void free_Matrix2D(Matrix2D *);
extern void read_PGM_file(char *, Matrix2D *);
extern void write_PGM_file(char *, Matrix2D *);


