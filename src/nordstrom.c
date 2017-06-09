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
#include <stdlib.h>
#include <math.h>

/*---------------------------------------------------------------------*/

double g1(int x, double K){
    return exp(-(x/K)*(x/K));
}

double g2(int x, double K){
    return 1/(1+(x/K)*(x/K));
}

/* Description: compute contrast map using Nordstorm filter                */

void Nordstrom(Matrix2D *mat_ut, Matrix2D *mat_utt, Matrix2D *mat_in, Matrix2D *mat_out,
          Matrix2D *mat_ut_0, Matrix2D *mat_in_0, double dt, int t, int g, double K, double lambda)
{
    register int  i, j, k;
    int           xsize_1, ysize_1;
    int           n, s, w, e;
    int           jn, js, iw, ie;
    double        g1e, g1w, g1n, g1s, sigma1;
    double        g2e, g2w, g2n, g2s, sigma2;
    unsigned char **in, **in_0, **out;
    double        **ut, **ut_0, **utt, **tmp;



    in  = mat_in->udata;
    in_0  = mat_in_0->udata;
    out = mat_out->udata;
    ut  = mat_ut->ddata;
    ut_0  = mat_ut_0->ddata;
    utt = mat_utt->ddata;


    // Pre-compute constants
    xsize_1 = mat_in->xsize - 1;
    ysize_1 = mat_in->ysize - 1;

    /* remplissage de ut pour les calculs en double */

    for(j = 0; j < mat_in->ysize; j++) {
        for (i = 0; i < mat_in->xsize; i++) {
            ut[j][i] = (double) (in[j][i]);
            ut_0[j][i] = (double) (in_0[j][i]);
            //printf("avant datalink\n");
            //datalink[j][i] = 0;
        }
        //printf("ut[j][2] = %f \nut_0[j][2] = %f\n", ut[j][2], ut_0[j][2]);
    }


    /* boucle principale de Nordstrom */

    for(k = 0; k < t; k++) {
        for (j = 0; j < mat_in->ysize; j++) {
            // Row index offsets (Neumann conditions)
            n = ((j == 0) ? 1 : -1);
            s = ((j == ysize_1) ? -1 : 1);
            // Pre-compute row indices
            jn = j + n;
            js = j + s;
            for (i = 0; i < mat_in->xsize; i++) {
                // Column index offsets (Neumann)
                w = ((i == 0) ? 1 : -1);
                e = ((i == xsize_1) ? -1 : 1);
                // Pre-compute column indices
                iw = i + w;
                ie = i + e;

                //datalink[j][i] = lambda * (ut_0[j][i] - ut[j][i]);

                // Compute gradient components
                // - kernel is not normalized to save 1 multiplication
                if (g == 1) {
                    g1e = g1(ut[j][ie] - ut[j][i], K);
                    g1w = g1(ut[j][iw] - ut[j][i], K);
                    g1n = g1(ut[js][i] - ut[j][i], K);
                    g1s = g1(ut[jn][i] - ut[j][i], K);
                    sigma1 = g1e + g1w + g1n + g1s;
                    // Compute contrast
                    // - kernel normalization is performed here
                    utt[j][i] = ut[j][i] + dt * (g1e * ut[j][ie] + g1w * ut[j][iw] + g1n * ut[jn][i] + g1s * ut[js][i]
                                                 - sigma1*ut[j][i]);

                    //printf("utt[j][i] = %f \n", utt[j][i]);
                } else if (g == 2) {
                    g2e = g2(ut[j][ie] - ut[j][i], K);
                    g2w = g2(ut[j][iw] - ut[j][i], K);
                    g2n = g2(ut[js][i] - ut[j][i], K);
                    g2s = g2(ut[jn][i] - ut[j][i], K);
                    sigma2 = g2e + g2w + g2n + g2s;
                    // Compute contrast
                    // - kernel normalization is performed here
                    utt[j][i] = ut[j][i] + dt * (g2e * ut[j][ie] + g2w * ut[j][iw] + g2n * ut[jn][i] + g2s * ut[js][i]
                                                 - sigma2*ut[j][i]);

                }

            }
        }

        //printf("itération = %d, ut = %p , utt = %p \n", k, ut, utt);
        tmp = utt;
        utt = ut;
        ut = tmp;

    }

    /* réécriture de la sortie dans out pour l'affichage en image */

    for(j = 0; j < mat_in->ysize; j++)
        for(i = 0; i < mat_in->xsize; i++) {
            out[j][i] = (int) (ut[j][i] + lambda * (ut_0[j][i] - utt[j][i]));
            printf("out[j][i] = %d \n", out[j][i]);
        }

}

/*---------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    register int   i;
    char           filename_in[256];
    char           filename_in_0[256];
    char           filename_out[256];
    double         dt, K, lambda;
    int            t, g;
    Matrix2D       in, in_0, out;
    Matrix2D       ut, ut_0, utt;    // buffers qui contiennent les in et out des itérations successives

    //sprintf(filename_out, "resultat.pgm");



    // ./permal -i "nomfic.pgm" dt t g K

    // printf("argv[3] = %s, argv[4] = %s, argv[5] = %s, argv[6] = %s \n", argv[3], argv[4], argv[5], argv[6]);

    dt = atof(argv[4]);
    t = atoi(argv[5]);
    g = atoi(argv[6]);
    K = atof(argv[7]);
    lambda = atof(argv[8]);

    /* Decodage des arguments */
    if (argc == 9)
        for (i=1; i< argc; i++) {
            if (!strcmp(argv[i], "-i")) {       /* Input file name */
                if (++i > argc-1)
                    usage(argv[0]);
                sprintf(filename_in, "%s", argv[i]);
                sprintf(filename_in_0, "%s", argv[i+1]);
            }
            else if (!strcmp(argv[i], "-o")) {  /* Output file name */
                if (++i > argc)
                    usage(argv[0]);
                //sprintf(filename_out, "resultat_%f_%d_%d_%f.pgm", dt, t, g, K);
//   printf("filename_out = %s", filename_out);
            }
        } else {
        usage(argv[0]);
        exit(1);
    }

    sprintf(filename_out, "images_nord/resultat_%f_%d_%d_%f_%f.pgm", dt, t, g, K, lambda);

    printf("dt = %f, t = %d, g = %d, K = %f, lambda = %f \n", dt, t, g, K, lambda);


////////////////////////////////
    /* Lecture et allocation m�moire de l'image d'entr�e */
    read_PGM_file(filename_in, &in);
    read_PGM_file(filename_in_0, &in_0);


    /* Allocation m�moire de l'image de sortie et des buffer intermédiaires */
    alloc_Matrix2D(&out, in.xsize, in.ysize, UCHAR_DATA);
    alloc_Matrix2D(&ut, in.xsize, in.ysize, DOUBLE_DATA);
    alloc_Matrix2D(&ut_0, in.xsize, in.ysize, DOUBLE_DATA);
    alloc_Matrix2D(&utt, in.xsize, in.ysize, DOUBLE_DATA);


    /* Traitement */

    Nordstrom(&ut, &utt, &in, &out, &ut_0, &in_0, dt, t, g, K, lambda);

    /* Ecriture de l'image de sortie */
    write_PGM_file(filename_out, &out);

    /* Lib�ration de la m�moire */
    free_Matrix2D(&in);
    free_Matrix2D(&in_0);
    free_Matrix2D(&out);
    free_Matrix2D(&ut);
    free_Matrix2D(&ut_0);
    free_Matrix2D(&utt);

    return 0;
}

/*---------------------------------------------------------------------*/
/* Description: Print out program syntax                               */

void usage(char *prgnam)
{

    (void) fprintf(stderr, "Usage: %s ", prgnam);
    fprintf(stderr, "./permal -i \"nomfic.pgm\" \"nomfic0.pgm\" dt t g K lambda\n");
    fprintf(stderr, "\n");
}


/*---------------------------------------------------------------------*/
