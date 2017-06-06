/*	NOMS DU BINOME : WILL KOMBO
 *
 *	pdefilter.c
 *	Compilation:	make pdefilter
 *	Auteur:		Nicolas ROUGON
 *			Institut Mines-Telecom / Telecom SudParis
 *      Date:		01/06/2015
 *
 *---------------------------------------------------------------------*/

#include "pdefilter.h"
#include <math.h>

/*---------------------------------------------------------------------*/

double g1(int x, double K){
    return exp(-(x/K)*(x/K));
}

double g2(int x, double K){
    return 1/(1+(x/K)*(x/K));
}

/* Description: compute contrast map using Sobel filter                */

void Permal(Matrix2D *mat_in, Matrix2D *mat_out, double dt, int K, int g)
{
    register int  i, j;
    int           xsize_1, ysize_1;
    int           n, s, w, e;
    int           jn, js, iw, ie;
    double        g1e, g1w, g1n, g1s, sigma1;
    double        g2e, g2w, g2n, g2s, sigma2;
    double        **in, **out, **tmp;

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
            if (g == 1){
                g1e = g1(in[j][ie]-in[j][i]);
                g1w = g1(in[j][iw]-in[j][i]);
                g1n = g1(in[js][i]-in[j][i]);
                g1s = g1(in[jn][i]-in[j][i]);
                sigma1 = g1e + g1w + g1n + g1s;
                // Compute contrast
                // - kernel normalization is performed here
                out[j][i] = (int)(in[j][i]+dt*(g1e*in[j][ie] + g1w*in[j][iw] + g1n*in[jn][i] + g1s*in[js][i] - sigma1));
            }

            else if(g == 2){
                g2e = g2(in[j][ie]-in[j][i]);
                g2w = g2(in[j][iw]-in[j][i]);
                g2n = g2(in[js][i]-in[j][i]);
                g2s = g2(in[jn][i]-in[j][i]);
                sigma2 = g2e + g2w + g2n + g2s;
                // Compute contrast
                // - kernel normalization is performed here
                out[j][i] = (int)(in[j][i]+dt*(g2e*in[j][ie] + g2w*in[j][iw] + g2n*in[jn][i] + g2s*in[js][i] - sigma2));
            }
        }
    }
}

/*---------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    register int   i;
    char           filename_in[256];
    char           filename_out[256];
    double         dt, K;
    int            t, g;
    Matrix2D in, out, tmp;

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
    int i;
    for(i = 0; i < t; i++){
        Permal(&in, &out, dt, K, g);
        tmp = out;
        out = in;
        in = tmp;
    }

    /* Ecriture de l'image de sortie */
    write_PGM_file(filename_out, &out);

    /* Lib�ration de la m�moire */
    free_Matrix2D(&in);
    free_Matrix2D(&out);
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
