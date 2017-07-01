/*	NOMS DU BINOME :
 *
 *	pdefilter.c
 *	Compilation:	make pdefilter
 *	Auteur:		Nicolas ROUGON
 *			Institut Mines-Telecom / Telecom SudParis
 *      Date:		01/06/2015
 *
 *---------------------------------------------------------------------*/

#include "pdefilter.h"

/*---------------------------------------------------------------------*/
/* Description: compute contrast map using Sobel filter                */

void Sobel(Matrix2D *mat_in, Matrix2D *mat_out)
{
    register int  i, j;
    int           xsize_1, ysize_1;
    int           n, s, w, e;
    int           jn, js, iw, ie;
    int           gx, gy;
    unsigned char **in, **out;

    in  = mat_in->udata;
    out = mat_out->udata;

    // Pre-compute constants
    xsize_1 = mat_in->xsize - 1;
    ysize_1 = mat_in->ysize - 1;

    for (j=0; j < mat_in->ysize; j++) {
        // Row index offsets (Neumann conditions)
        n = ((j==0) ? 1 : -1);
        s = ((j==ysize_1) ? -1 : 1);
        // Pre-compute row indices
        jn = j + n;
        js = j + s;
        for (i=0 ; i < mat_in->xsize ; i++) {
            // Column index offsets (Neumann)
            w = ((i==0) ? 1 : -1);
            e = ((i==xsize_1) ? -1 : 1);
            // Pre-compute column indices
            iw = i + w;
            ie = i + e;

            // Compute gradient components
            // - kernel is not normalized to save 1 multiplication
            gx = in[jn][ie] + 2*in[j][ie] + in[js][ie] -
                 in[jn][iw] - 2*in[j][iw] - in[js][iw];
            gy = in[js][iw] + 2*in[js][i] + in[js][ie] -
                 in[jn][iw] - 2*in[jn][i] - in[jn][ie];

            // Compute contrast
            // - kernel normalization is performed here
            out[j][i] = (int)(0.25*sqrt((double)(gx*gx + gy*gy)));
        }
    }
}

/*---------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    register int   i;
    char           filename_in[256];
    char           filename_out[256];
    Matrix2D in, out;

    sprintf(filename_out, "resultat.pgm");

    /* Decodage des arguments */
    if (argc > 1)
        for (i=1; i< argc; i++) {
            if (!strcmp(argv[i], "-i")) {       /* Input file name */
                if (++i > argc)
                    usage(argv[0]);
                sprintf(filename_in, "%s", argv[i]);
            }
            else if (!strcmp(argv[i], "-o")) {  /* Output file name */
                if (++i > argc)
                    usage(argv[0]);
                sprintf(filename_out, "%s", argv[i]);
            }
        } else {
        usage(argv[0]);
        exit(1);
    }

    /* Lecture et allocation m�moire de l'image d'entr�e */
    read_PGM_file(filename_in, &in);

    /* Allocation m�moire de l'image de sortie */
    alloc_Matrix2D(&out, in.xsize, in.ysize, UCHAR_DATA);

    /* Traitement */
    Sobel(&in, &out);    // A remplacer par ***

    /* Ecriture de l'image de sortie */
    write_PGM_file(filename_out, &out);

    /* Lib�ration de la m�moire */
    free_Matrix2D(&in);
    free_Matrix2D(&out);

    return 0;
}

/*---------------------------------------------------------------------*/
/* Description: Print out program syntax                               */

void usage(char *prgnam)
{

    (void) fprintf(stderr, "Usage: %s ", prgnam);
    fprintf(stderr, "-i input_filename [-o output_filename]\n");
    fprintf(stderr, "\n");
}


/*---------------------------------------------------------------------*/