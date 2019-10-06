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

float coef(int o, int n)
{
    float B[o], A[o][o+1], *C;
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

    C = DirectLU( (float *)A, o);

   //for(i = 0; i < o; i++)
   //    printf("%f  \n", C[i]);
   //printf("\n");

    return C[n];

}

float *DirectLU(float *u, int n)
{
    int i, r, k;
    float *x;
    x = (float*) malloc(sizeof(float)*n);
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
}
