/***************************************************************************
 *
 * This function is used to calculate differential coefficient of the
 *     staggered grid method.
 *
 *  'o' is the half stencil of difference.
 *  'n' is the n-th coefficent
 *
 * Authors: Luqian Jiang   Email: jlq0608@gmail.com
 *          Wei Zhang      Email: zhangwei@sustc.edu.cn
 *
 * Copyright (c) 2018 - 2019 zwlab
 *
 * History: 08/2019: Original version created by Luqian Jiang
 *
 ***************************************************************************/
#include "elastic2d_math.h"
#include "staggered_fd_coef.h"

double ssg_coef(int o, int n)
{
    double B[o], A[o][o+1], *C;
    int i,j;
    B[0] = 1.0;
    for (i = 1; i < o; i++) {
        B[i] = 0.0;
    }
    for (i = 0; i < o; i++) {
        for (j = 0; j < o; j++) {
            A[i][j] = pow( (2*(j+1)-1), (2*(i+1)-1) );
         //   printf("%f\n",A[i][j] );
        }
        A[i][o] = B[i];
    }

    C = DirectLU( (double *)A, o);

   //for(i = 0; i < o; i++)
   //    printf("%f  \n", C[i]);
   //printf("\n");

    return C[n];

}

double *DirectLU(double *u, int n)
{
    int i, r, k;
    double *x;
    x = (double*) malloc(sizeof(double)*n);
    for (r = 0; r <= n-1; r++) {
        for (i = r; i <= n; i++)
            for (k = 0; k <= r-1; k++)
                *(u + r*(n+1) + i) -= * (u + r*(n+1) + k) * ( *(u + k*(n+1) + i ) );

        for (i = r+1; i <= n-1; i++) {
            for (k = 0; k <= r-1; k++)
                *(u + i*(n+1) + r) -= * (u + i*(n+1) + k) * ( *(u + k*(n+1) + r) );
            *(u + i*(n+1) + r) /= *(u + r*(n+1) + r);
        }
    }

    for (i = n-1; i >= 0; i--) {
        for (r = n-1; r >= i+1; r--)
            *(u + i*(n+1) + n) -= * (u + i*(n+1) + r) * x[r];
        x[i] = * (u + i*(n+1) + n)/( *(u + i*(n+1) + i));
    }

    return x;
    free(x);
}

/********************************************************************
 * Optimal staggered-grid finite-difference schemes
 *
 *  ref: Liu, 2014, Optimal staggered-grid finite-difference schemes
 *              based on least-squares for wave equation modelling
 *
 ********************************************************************/

float esgfd_coef(int o, int n)
{
    float A[o];
    switch(o) {
        case 2 :
            A[0] = -0.1129042e+1;
            A[1] =  0.4301412e-1;
            break;
        case 3 :
            A[0] =  0.1185991e+1;
            A[1] = -0.7249965e-1;
            A[2] =  0.6301572e-2;
            break;
        case 4 :
            A[0] =  0.1217990e+1;
            A[1] = -0.9382142e-1;
            A[2] =  0.1507536e-1;
            A[3] = -0.1700324e-2;
            break;
        case 5 :
            A[0] =  0.1236607e+1;
            A[1] = -0.1082265e+0;
            A[2] =  0.2343440e-1;
            A[3] = -0.5033546e-2;
            A[4] =  0.6817483e-3;
            break;
        case 6 :
            A[0] =  0.1247662e+1;
            A[1] = -0.1175538e+0;
            A[2] =  0.2997970e-1;
            A[3] = -0.8719078e-2;
            A[4] =  0.2215897e-2;
            A[5] = -0.3462075e-3;
            break;
    }

    return A[n];

}
