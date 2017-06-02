/*
 *	globals.h
 *	Auteur:		Nicolas ROUGON
 *			Institut Mines-Telecom / Telecom SudParis
 *      Date:		01/06/2015
 *
 *---------------------------------------------------------------------*/

#define UCHAR_DATA   1
#define DOUBLE_DATA  2

//---- Data types
typedef struct _Matrix2D {
  int           xsize;
  int           ysize;
  int           data_type;
  unsigned char **udata;  
  double        **ddata;  
} Matrix2D;
  
