/*	
 *	utitilies.c
 *	Auteur:		Nicolas ROUGON
 *			Institut Mines-Telecom / Telecom SudParis
 *      Date:		01/06/2015
 *
 *---------------------------------------------------------------------*/

#include "utilities.h"

/*---------------------------------------------------------------------*/
/* Description: Allocate 2D array                                      */

void alloc_Matrix2D(Matrix2D *mat, int xsize, int ysize, int data_type)
{
  register int j;  

  mat->xsize = xsize;
  mat->ysize = ysize;  
  mat->data_type = data_type;  

  switch(data_type) {
  case DOUBLE_DATA:
    mat->ddata = (double **)malloc(ysize*sizeof(double *));
    if (mat->ddata == NULL) {//
// Created by blackdoe on 6/2/17.
//


      fprintf(stderr, "Memoire insuffisante\n");//
// Created by blackdoe on 6/2/17.
//


      exit(1);
    }
    
    for (j=0; j<ysize; j++) {
      (mat->ddata)[j] = (double *)malloc(xsize*sizeof(double));
      if ((mat->ddata)[j] == NULL) {
	fprintf(stderr, "Memoire insuffisante\n");  
	exit(1);
      }
    }
    break;
  case UCHAR_DATA:
    mat->udata = (unsigned char **)malloc(ysize*sizeof(unsigned char *));
    if (mat->udata == NULL) {
      fprintf(stderr, "Memoire insuffisante\n");  
      exit(1);
    }
    
    for (j=0; j<ysize; j++) {
      (mat->udata)[j] = (unsigned char *)malloc(xsize*sizeof(unsigned char));
      if ((mat->udata)[j] == NULL) {
	fprintf(stderr, "Memoire insuffisante\n");  
	exit(1);
      }
    }
    break;    
  default:
    fprintf(stderr, "Type de donnees inconnu\n");  
    exit(1);
    break;    
  }
}

/*---------------------------------------------------------------------*/
/* Description: Free 2D array                                          */

void free_Matrix2D(Matrix2D *mat)
{
  register int j;

  switch(mat->data_type) {
  case DOUBLE_DATA:
    for (j=0; j<mat->ysize; j++)
      free((char *)(mat->ddata)[j]);    

    free((char *)mat->ddata);
    break;
  case UCHAR_DATA:
    for (j=0; j<mat->ysize; j++)
      free((char *)(mat->udata)[j]);    

    free((char *)mat->udata);
    break;    
  default:
    fprintf(stderr, "Type de donnees inconnu\n");  
    exit(1);
    break;    
  }
}

/*---------------------------------------------------------------------*/
/* Description: Read PGM raw format file into a Matrix2D array         */

void read_PGM_file(char *filename, Matrix2D *mat)
{   
  FILE         *fp;
  register int j;
  int          xsize, ysize, maxval;

  if ((fp = fopen(filename, "r")) == NULL ) { 
    fprintf(stderr, "Impossible d'ouvrir %s\n", filename);  
    exit(1); 
  }

  /* Lecture de l'ent�te */
  if ( getc(fp)!='P' || getc(fp)!='5' ) {
    fprintf (stderr, "%s n'est pas au format PGM binaire (P5)\n", filename);
    exit(1) ;
  }

  fscanf(fp, "%d %d\n", &xsize, &ysize);
  fscanf(fp, "%d\n", &maxval) ;

  /* Allocation m�moire pour les donn�es */
  alloc_Matrix2D(mat, xsize, ysize, UCHAR_DATA);

  /* Lecture des donn�es */
  for (j=0; j<ysize; j++)
    if (fread ((mat->udata)[j], sizeof(unsigned char), xsize, fp) != xsize ) { 
      fprintf(stderr, "Erreur de lecture du fichier %s\n", filename);  
      exit(1); 
    }

  fclose(fp);
}

/*---------------------------------------------------------------------*/
/* Description: Save image data to PGM raw format                      */

void write_PGM_file(char *filename, Matrix2D *mat) 
{ 
  FILE         *fp;
  register int j;

  if ((fp = fopen(filename, "w")) == NULL ) { 
    fprintf(stderr, "Impossible d'ouvrir %s\n", filename);  
    exit(1); 
  }
  
  /* Ecriture de l'ent�te */
  fprintf(fp, "P5\n");
  fprintf(fp, "%d %d\n", mat->xsize, mat->ysize);
  fprintf(fp, "255\n");
  
  /* Ecriture des donn�es */
  for (j=0; j<mat->ysize; j++)
    if (fwrite(mat->udata[j], sizeof(unsigned char), mat->xsize, fp) != mat->xsize ) { 
      fprintf(stderr, "Erreur d'ecriture du fichier %s\n", filename);  
      exit(1); 
    }

  fclose(fp);
}

/*---------------------------------------------------------------------*/
