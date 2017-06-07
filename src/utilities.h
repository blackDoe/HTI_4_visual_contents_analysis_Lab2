/*
 *	utilities.h
 *	Auteur:		Nicolas ROUGON
 *			Institut Mines-Telecom / Telecom SudParis
 *      Date:		01/06/2015
 *
 *---------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "globals.h"

//---- Functions
void alloc_Matrix2D(Matrix2D *, int, int, int);
void free_Matrix2D(Matrix2D *);
void read_PGM_file(char *, Matrix2D *);
void write_PGM_file(char *, Matrix2D *);
