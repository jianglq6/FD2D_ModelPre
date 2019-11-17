/***************************************************************************
 *
 * This function is used to add source in Lebedev or staggered grid scheme.
 *
 * Authors: Luqian Jiang   Email: jlq0608@gmail.com
 *          Wei Zhang      Email: zhangwei@sustc.edu.cn
 *
 * Copyright (c) 2018 - 2019 zwlab
 *
 * History: 09/2019: Original version created by Luqian Jiang
 *          10/2019: Add some source impulse method. Luqian
 *          11/2019: Add staggered source impulse method by simplifying
 *                   the method in Lebedev grid.
 *
 ***************************************************************************/
#include "elastic2d_src.h"
#include "elastic2d_stf.h"
#include "elastic2d_math.h"


/* Moment source */
int src_moment_lebedev(float *hTxx_1, float *hTxx_2,
                       float *hTxz_1, float *hTxz_2,
                       float *hTzz_1, float *hTzz_2,
                       float *xvec1, float *zvec1,
                       float *xvec2, float *zvec2,
                       int source_impulse_method, struct Src src,
                       float current_time, float dt,
                       float dx, float dz, int nx)
{
    int i_source, indx_xs, indx_zs, indx_source;
    float stf;

    FILE *fp1;
    if ((fp1=fopen(STF_FILE, "a+")) == NULL) {
        printf("%s cannot be opened\n",STF_FILE );
    }

    stf = cal_source_time_function(src.stf_type_id[0], current_time,
                        src.stf_timefactor[0], src.stf_freqfactor[0]);
    /* write source time function, for check */
    fprintf(fp1, "%12.6f %20.6f\n", current_time, stf);
    fclose(fp1);

    /* If there i_source no moment source */
    if (src.number_of_src <= 0)
        return;

    /* source is just the first set of grid */
    /* !!! Just for test */
    if (source_impulse_method == IMPULSE_DIRECT) {

        for (i_source = 0; i_source < src.number_of_src; i_source++) {
        /* source time function */
            stf = cal_source_time_function(src.stf_type_id[i_source], current_time,
                        src.stf_timefactor[i_source], src.stf_freqfactor[i_source]);
            indx_xs = round( (src.xs[i_source]-xvec1[0])/dx );
            indx_zs = round( (src.zs[i_source]-zvec1[0])/dz );
            indx_source = indx_zs*nx+indx_xs;
            hTxx_1[indx_source] -=  src.Mxx[i_source]/dx/dz * stf;
            hTzz_1[indx_source] -=  src.Mzz[i_source]/dx/dz * stf;

        }
    }

    /* integer grid weight 1, half-grid weight 1/4 */
    /* !!! Just for test */
    else if (source_impulse_method == IMPULSE_AVERAGE) {
        for (i_source = 0; i_source < src.number_of_src; i_source++) {

            stf = cal_source_time_function(src.stf_type_id[i_source], current_time,
                        src.stf_timefactor[i_source], src.stf_freqfactor[i_source]);

            indx_xs = round( (src.xs[i_source]-xvec1[0])/dx );
            indx_zs = round( (src.zs[i_source]-zvec1[0])/dz );
            indx_source = indx_zs*nx+indx_xs;
            hTxx_1[indx_source] -=  src.Mxx[i_source]/dx/dz * stf;
            hTzz_1[indx_source] -=  src.Mzz[i_source]/dx/dz * stf;

            hTxx_2[indx_source]      -=  1.0/4.0*src.Mxx[i_source]/dx/dz * stf;
            hTzz_2[indx_source]      -=  1.0/4.0*src.Mzz[i_source]/dx/dz * stf;
            hTxx_2[indx_source-1]    -=  1.0/4.0*src.Mxx[i_source]/dx/dz * stf;
            hTzz_2[indx_source-1]    -=  1.0/4.0*src.Mzz[i_source]/dx/dz * stf;

            hTxx_2[indx_source-nx]   -=  1.0/4.0*src.Mxx[i_source]/dx/dz * stf;
            hTzz_2[indx_source-nx]   -=  1.0/4.0*src.Mzz[i_source]/dx/dz * stf;
            hTxx_2[indx_source-nx-1] -=  1.0/4.0*src.Mxx[i_source]/dx/dz * stf;
            hTzz_2[indx_source-nx-1] -=  1.0/4.0*src.Mzz[i_source]/dx/dz * stf;

        }

    }

    /* arbitary source */
    else if (source_impulse_method == IMPULSE_SINC) {

        src_sinc_moment_lebedev(hTxx_1, hTxx_2,
                                hTxz_1, hTxz_2,
                                hTzz_1, hTzz_2,
                                xvec1, zvec1,
                                xvec2, zvec2,
                                src,
                                current_time, dt,
                                dx, dz, nx);
    }


}

/* Force source */
int src_force_lebedev(float *hVx_1, float *hVz_1,
                      float *hVx_2, float *hVz_2,
                      float *xvec1, float *zvec1,  /* for the integer grids  */
                      float *xvec2, float *zvec2,  /* for the half grids     */
                      int source_impulse_method, struct Src src,
                      float current_time, float dt,
                      float dx, float dz, int nx)
{


    /* arbitary source */
    if (source_impulse_method == IMPULSE_SINC) {

        src_sinc_force_lebedev(hVx_1, hVz_1,
                               hVx_2, hVz_2,
                               xvec1, zvec1,  /* for the integer grids  */
                               xvec2, zvec2,  /* for the half grids     */
                               src,
                               current_time, dt,
                               dx, dz, nx);
    }

}


/* Moment source */
int src_sinc_moment_lebedev(float *hTxx_1, float *hTxx_2,
                            float *hTxz_1, float *hTxz_2,
                            float *hTzz_1, float *hTzz_2,
                            float *xvec1, float *zvec1,
                            float *xvec2, float *zvec2,
                            struct Src src,
                            float current_time, float dt,
                            float dx, float dz, int nx)
{
    int   i_source, i, k, indx;  // i-th source
    // index of the source location in  grid
    float indx_xs_1, indx_zs_1, indx_xs_2, indx_zs_2;
    int   ixvec1_src[2*NSRCEXT], ixvec2_src[2*NSRCEXT], izvec1_src[2*NSRCEXT], izvec2_src[2*NSRCEXT];
    double DNx1[2*NSRCEXT], DNx2[2*NSRCEXT], DNz1[2*NSRCEXT], DNz2[2*NSRCEXT];
    float  V, Mxx, Mzz, Mxz, d, stf;
    double moment_src;

   //ixvec1_src = (int*)malloc(2*NSRCEXT*sizeof(int));
   //ixvec2_src = (int*)malloc(2*NSRCEXT*sizeof(int));
   //izvec1_src = (int*)malloc(2*NSRCEXT*sizeof(int));
   //izvec2_src = (int*)malloc(2*NSRCEXT*sizeof(int));

   ///* Damp */
   //DNx1 = (double*) malloc(2*NSRCEXT*sizeof(double));
   //DNz1 = (double*) malloc(2*NSRCEXT*sizeof(double));
   //DNx2 = (double*) malloc(2*NSRCEXT*sizeof(double));
   //DNz2 = (double*) malloc(2*NSRCEXT*sizeof(double));


    /* If there i_source no moment source */
    if (src.number_of_src <= 0)
        return;

    /* Loop each moment source */
    for (i_source = 0; i_source < src.number_of_src; i_source++) {

        /* the factor of the i_source-th source */
        stf = cal_source_time_function(src.stf_type_id[i_source], current_time,
                        src.stf_timefactor[i_source], src.stf_freqfactor[i_source]);

        Mxx = src.Mxx[i_source];
        Mzz = src.Mzz[i_source];
        Mxz = src.Mxz[i_source];

        indx_xs_1 = (src.xs[i_source] - xvec1[0])/dx;
        indx_zs_1 = (src.zs[i_source] - zvec1[0])/dz;

        indx_xs_2 = (src.xs[i_source] - xvec2[0])/dx;
        indx_zs_2 = (src.zs[i_source] - zvec2[0])/dz;

        /* calculate index of discrete source */
        for (i = 0; i < NSRCEXT; i++) {
            ixvec1_src[NSRCEXT-i-1] = floor(indx_xs_1) - i;
            ixvec1_src[NSRCEXT+i  ] =  ceil(indx_xs_1) + i;
            izvec1_src[NSRCEXT-i-1] = floor(indx_zs_1) - i;
            izvec1_src[NSRCEXT+i  ] =  ceil(indx_zs_1) + i;

            ixvec2_src[NSRCEXT-i-1] = floor(indx_xs_2) - i;
            ixvec2_src[NSRCEXT+i  ] =  ceil(indx_xs_2) + i;
            izvec2_src[NSRCEXT-i-1] = floor(indx_zs_2) - i;
            izvec2_src[NSRCEXT+i  ] =  ceil(indx_zs_2) + i;
        }

        /* calculate discrete delta function */
        for (i = 0; i < 2*NSRCEXT; i++) {
            DNx1[i] =  calculate_windowed_impulse(ixvec1_src[i]-indx_xs_1, KAISERB, NSRCEXT);
            DNz1[i] =  calculate_windowed_impulse(izvec1_src[i]-indx_zs_1, KAISERB, NSRCEXT);
            DNx2[i] =  calculate_windowed_impulse(ixvec2_src[i]-indx_xs_2, KAISERB, NSRCEXT);
            DNz2[i] =  calculate_windowed_impulse(izvec2_src[i]-indx_zs_2, KAISERB, NSRCEXT);
        }

        /* add source */
        for (i = 0; i < 2*NSRCEXT; i++) {
            for (k = 0; k < 2*NSRCEXT; k++ ) {

                V = dx * dz;

                //-- for 00
                indx = izvec1_src[k] * nx + ixvec1_src[i];
                d = DNx1[i] * DNz1[k]/4.0;        /* d: Hicks, 2002 */
                moment_src = d*stf/V;
               // printf("%f\n", d);
                hTxx_1[indx] -= Mxx * moment_src;
                hTzz_1[indx] -= Mzz * moment_src;
                hTxz_2[indx] -= Mxz * moment_src;

                //-- for 11
                indx = izvec2_src[k] * nx + ixvec2_src[i];
                d = DNx2[i] * DNz2[k]/4.0;
                moment_src = d*stf/V;
               // printf("%f\n", d);
                hTxx_2[indx] -= Mxx * moment_src;
                hTzz_2[indx] -= Mzz * moment_src;
                hTxz_1[indx] -= Mxz * moment_src;

            }
        }
    }


    //free(ixvec1_src); free(ixvec2_src);
    //free(izvec1_src); free(izvec2_src);
    //free(DNx1); free(DNx2); free(DNz1); free(DNz2);
    return 0;
}

/* Force source */
int src_sinc_force_lebedev(float *hVx_1, float *hVz_1,
                           float *hVx_2, float *hVz_2,
                           float *xvec1, float *zvec1,  /* for the integer grids  */
                           float *xvec2, float *zvec2,  /* for the half grids     */
                           struct Src src,
                           float current_time, float dt,
                           float dx, float dz, int nx)
{
    int   i_source, i, k, indx;  // i-th source
    // index of the source location in  grid
    float indx_xs_1, indx_zs_1, indx_xs_2, indx_zs_2;
    int   *ixvec1_src, *ixvec2_src, *izvec1_src, *izvec2_src;
    double *DNx1, *DNx2, *DNz1, *DNz2;
    float V, Fx, Fz, Mxz, d, stf;
    double force_src;

    ixvec1_src = (int*)malloc(2*NSRCEXT*sizeof(int));
    ixvec2_src = (int*)malloc(2*NSRCEXT*sizeof(int));
    izvec1_src = (int*)malloc(2*NSRCEXT*sizeof(int));
    izvec2_src = (int*)malloc(2*NSRCEXT*sizeof(int));

    DNx1 = (double*) malloc(2*NSRCEXT*sizeof(double));
    DNz1 = (double*) malloc(2*NSRCEXT*sizeof(double));
    DNx2 = (double*) malloc(2*NSRCEXT*sizeof(double));
    DNz2 = (double*) malloc(2*NSRCEXT*sizeof(double));


    /* If there i_source no moment source */
    if (src.number_of_src <= 0)
        return;

    /* Loop each moment source */
    for (i_source = 0; i_source < src.number_of_src; i_source++) {

        Fx = src.Fx[i_source];
        Fz = src.Fz[i_source];

        indx_xs_1 = (src.xs[i_source] - xvec1[0])/dx;
        indx_zs_1 = (src.zs[i_source] - zvec1[0])/dz;

        indx_xs_2 = (src.xs[i_source] - xvec2[0])/dx;
        indx_zs_2 = (src.zs[i_source] - zvec2[0])/dz;

        /* calculate index of discrete source */
        for (i = 0; i < NSRCEXT; i++) {
            ixvec1_src[NSRCEXT-i-1] = floor(indx_xs_1) - i;
            ixvec1_src[NSRCEXT+i  ] =  ceil(indx_xs_1) + i;
            izvec1_src[NSRCEXT-i-1] = floor(indx_zs_1) - i;
            izvec1_src[NSRCEXT+i  ] =  ceil(indx_zs_1) + i;

            ixvec2_src[NSRCEXT-i-1] = floor(indx_xs_2) - i;
            ixvec2_src[NSRCEXT+i  ] =  ceil(indx_xs_2) + i;
            izvec2_src[NSRCEXT-i-1] = floor(indx_zs_2) - i;
            izvec2_src[NSRCEXT+i  ] =  ceil(indx_zs_2) + i;
        }

        /* calculate discrete delta function */
        for (i = 0; i < 2*NSRCEXT; i++) {
            DNx1[i] =  calculate_windowed_impulse(ixvec1_src[i]-indx_xs_1, KAISERB, NSRCEXT);
            DNz1[i] =  calculate_windowed_impulse(izvec1_src[i]-indx_zs_1, KAISERB, NSRCEXT);
            DNx2[i] =  calculate_windowed_impulse(ixvec2_src[i]-indx_xs_2, KAISERB, NSRCEXT);
            DNz2[i] =  calculate_windowed_impulse(izvec2_src[i]-indx_zs_2, KAISERB, NSRCEXT);
        }

        /* the factor of the i_source-th source */
        stf = cal_source_time_function(src.stf_type_id[i_source], \
                    current_time, src.stf_timefactor[i_source], src.stf_freqfactor[i_source]);

        /* add source */
        for (i = 0; i < 2*NSRCEXT; i++) {
            for (k = 0; k < 2*NSRCEXT; k++ ) {

                V = dx * dz;

                //-- for points 01
                indx = izvec1_src[k] * nx + ixvec2_src[i];
                d = DNx2[i] * DNz1[k]/4.0;        /* d: Hicks, 2002 */
                force_src = d*stf/V/dt;
                //printf("%f%f %f %f\n",d, force_src, DNx2[i], DNz1[k]);
                hVx_1[indx] += Fx * force_src;
                hVz_2[indx] += Fz * force_src;

                //-- for points 10
                indx = izvec2_src[k] * nx + ixvec1_src[i];
                d = DNx1[i] * DNz2[k]/4.0;
                force_src = d*stf/V/dt;
                //printf("%f%f %f   %f\n",d, force_src, DNx1[i], DNz2[k]);
                hVx_2[indx] += Fx * force_src;
                hVz_1[indx] += Fz * force_src;

            }
        }
    }

    free(ixvec1_src); free(ixvec2_src);
    free(izvec1_src); free(izvec2_src);
    free(DNx1); free(DNx2); free(DNz1); free(DNz2);

    return 0;

}

float cal_source_time_function(int flag_stf_type, float t, float t0, float f0)
{
    float stf;

    // TODO: have not done
    if (flag_stf_type == SIG_STF_RICKER) {
        stf = fun_ricker(t,f0,t0);
    } else if (flag_stf_type == SIG_STF_GAUSS) {
        stf = fun_gauss(t,f0,t0);
    }

    return stf;

}


/* Dicrete source function                                               *
 *  ref: Hicks, 2002, Arbitrary source and receiver positioning in       *
 *        finite-difference schemes using Kai_sourceer windowwd sinc function  */
double calculate_windowed_impulse(float x, float b, float r)
{
    double y = 0.0;
    if (abs(x) < 1e-5) {
        y = 1.00;
    }
    else if (1.0-x*x/r*r > 0.0) {
    /* Kai_sourceer finite impulse response filters. */
        y = first_modified_Bessel( 0, b*sqrt(1-x*x/r*r) ) / first_modified_Bessel(0, b);
    /* sinc */
        y = y * sin(PI*x) / (PI*x);
    }
    else {
        y = 0.0;
    }

    return y;
}

/* Moment source */
int src_moment_ssg(float *hTxx, float *hTxz, float *hTzz,
                         float *xvec1, float *zvec1,
                         float *xvec2, float *zvec2,
                         int source_impulse_method, struct Src src,
                         float current_time, float dt,
                         float dx, float dz, int nx)
{
    int i_source, indx_xs, indx_zs, indx_source;
    float stf;

    FILE *fp1;
    if ((fp1=fopen(STF_FILE, "a+")) == NULL) {
        printf("%s cannot be opened\n",STF_FILE );
    }

    stf = cal_source_time_function(src.stf_type_id[0], current_time,
                        src.stf_timefactor[0], src.stf_freqfactor[0]);

    /* write source time function, for check */
    fprintf(fp1, "%12.6f %20.6f\n", current_time, stf);
    fclose(fp1);

    /* If there i_source no moment source */
    if (src.number_of_src <= 0)
        return;

    /* source is just the first set of grid */
    /* !!! Just for test */
    if (source_impulse_method == IMPULSE_DIRECT) {

        for (i_source = 0; i_source < src.number_of_src; i_source++) {
        /* source time function */
            stf = cal_source_time_function(src.stf_type_id[i_source], current_time,
                        src.stf_timefactor[i_source], src.stf_freqfactor[i_source]);
            indx_xs = round( (src.xs[i_source]-xvec1[0])/dx );
            indx_zs = round( (src.zs[i_source]-zvec1[0])/dz );
            indx_source = indx_zs*nx+indx_xs;
            hTxx[indx_source] -= src.Mxx[i_source]/dx/dz * stf;
            hTzz[indx_source] -= src.Mzz[i_source]/dx/dz * stf;
        }
    }

    /* In staggered grid, this is the same as first situation */
    /* !!! Just for test */
    else if (source_impulse_method == IMPULSE_AVERAGE) {
        for (i_source = 0; i_source < src.number_of_src; i_source++) {

            stf = cal_source_time_function(src.stf_type_id[i_source], current_time,
                        src.stf_timefactor[i_source], src.stf_freqfactor[i_source]);

            indx_xs = round( (src.xs[i_source]-xvec1[0])/dx );
            indx_zs = round( (src.zs[i_source]-zvec1[0])/dz );
            indx_source = indx_zs*nx+indx_xs;
            hTxx[indx_source] -=  src.Mxx[i_source]/dx/dz * stf;
            hTzz[indx_source] -=  src.Mzz[i_source]/dx/dz * stf;

        }

    }

    /* arbitary source */
    else if (source_impulse_method == IMPULSE_SINC) {

        src_sinc_moment_ssg(hTxx, hTxz, hTzz,
                            xvec1, zvec1,
                            xvec2, zvec2,
                            src,
                            current_time, dt,
                            dx, dz, nx);
    }


}

/* Force source */
int src_force_ssg(float *hVx, float *hVz,
                        float *xvec1, float *zvec1,  /* for the integer grids  */
                        float *xvec2, float *zvec2,  /* for the half grids     */
                        int source_impulse_method, struct Src src,
                        float current_time, float dt,
                        float dx, float dz, int nx)
{


    /* arbitary source */
    if (source_impulse_method == IMPULSE_SINC) {

        src_sinc_force_ssg(hVx, hVz,
                           xvec1, zvec1,  /* for the integer grids  */
                           xvec2, zvec2,  /* for the half grids     */
                           src,
                           current_time, dt,
                           dx, dz, nx);
    }

}


/* Moment source */
int src_sinc_moment_ssg(float *hTxx, float *hTxz, float *hTzz,
                        float *xvec1, float *zvec1,
                        float *xvec2, float *zvec2,
                        struct Src src,
                        float current_time, float dt,
                        float dx, float dz, int nx)
{
    int   i_source, i, k, indx;  // i-th source
    // index of the source location in  grid
    float indx_xs_1, indx_zs_1, indx_xs_2, indx_zs_2;
    int   ixvec1_src[2*NSRCEXT], ixvec2_src[2*NSRCEXT], izvec1_src[2*NSRCEXT], izvec2_src[2*NSRCEXT];
    double DNx1[2*NSRCEXT], DNx2[2*NSRCEXT], DNz1[2*NSRCEXT], DNz2[2*NSRCEXT];
    float  V, Mxx, Mzz, Mxz, d, stf;
    double moment_src;

   //ixvec1_src = (int*)malloc(2*NSRCEXT*sizeof(int));
   //ixvec2_src = (int*)malloc(2*NSRCEXT*sizeof(int));
   //izvec1_src = (int*)malloc(2*NSRCEXT*sizeof(int));
   //izvec2_src = (int*)malloc(2*NSRCEXT*sizeof(int));

   ///* Damp */
   //DNx1 = (double*) malloc(2*NSRCEXT*sizeof(double));
   //DNz1 = (double*) malloc(2*NSRCEXT*sizeof(double));
   //DNx2 = (double*) malloc(2*NSRCEXT*sizeof(double));
   //DNz2 = (double*) malloc(2*NSRCEXT*sizeof(double));


    /* If there i_source no moment source */
    if (src.number_of_src <= 0)
        return;

    /* Loop each moment source */
    for (i_source = 0; i_source < src.number_of_src; i_source++) {

        /* the factor of the i_source-th source */
        stf = cal_source_time_function(src.stf_type_id[i_source], current_time,
                        src.stf_timefactor[i_source], src.stf_freqfactor[i_source]);

        Mxx = src.Mxx[i_source];
        Mzz = src.Mzz[i_source];
        Mxz = src.Mxz[i_source];

        indx_xs_1 = (src.xs[i_source] - xvec1[0])/dx;
        indx_zs_1 = (src.zs[i_source] - zvec1[0])/dz;

        indx_xs_2 = (src.xs[i_source] - xvec2[0])/dx;
        indx_zs_2 = (src.zs[i_source] - zvec2[0])/dz;

        /* calculate index of discrete source */
        for (i = 0; i < NSRCEXT; i++) {
            ixvec1_src[NSRCEXT-i-1] = floor(indx_xs_1) - i;
            ixvec1_src[NSRCEXT+i  ] =  ceil(indx_xs_1) + i;
            izvec1_src[NSRCEXT-i-1] = floor(indx_zs_1) - i;
            izvec1_src[NSRCEXT+i  ] =  ceil(indx_zs_1) + i;

            ixvec2_src[NSRCEXT-i-1] = floor(indx_xs_2) - i;
            ixvec2_src[NSRCEXT+i  ] =  ceil(indx_xs_2) + i;
            izvec2_src[NSRCEXT-i-1] = floor(indx_zs_2) - i;
            izvec2_src[NSRCEXT+i  ] =  ceil(indx_zs_2) + i;
        }

        /* calculate discrete delta function */
        for (i = 0; i < 2*NSRCEXT; i++) {
            DNx1[i] =  calculate_windowed_impulse(ixvec1_src[i]-indx_xs_1, KAISERB, NSRCEXT);
            DNz1[i] =  calculate_windowed_impulse(izvec1_src[i]-indx_zs_1, KAISERB, NSRCEXT);
            DNx2[i] =  calculate_windowed_impulse(ixvec2_src[i]-indx_xs_2, KAISERB, NSRCEXT);
            DNz2[i] =  calculate_windowed_impulse(izvec2_src[i]-indx_zs_2, KAISERB, NSRCEXT);
        }

        /* add source */
        for (i = 0; i < 2*NSRCEXT; i++) {
            for (k = 0; k < 2*NSRCEXT; k++ ) {

                V = dx * dz;

                //-- for 00
                indx = izvec1_src[k] * nx + ixvec1_src[i];
                d = DNx1[i] * DNz1[k];        /* d: Hicks, 2002 */
                moment_src = d*stf/V;
               // printf("%f\n", d);
                hTxx[indx] -= Mxx * moment_src;
                hTzz[indx] -= Mzz * moment_src;

                //-- for 11
                indx = izvec2_src[k] * nx + ixvec2_src[i];
                d = DNx2[i] * DNz2[k];
                moment_src = d*stf/V;
               // printf("%f\n", d);
                hTxz[indx] -= Mxz * moment_src;

            }
        }
    }


    //free(ixvec1_src); free(ixvec2_src);
    //free(izvec1_src); free(izvec2_src);
    //free(DNx1); free(DNx2); free(DNz1); free(DNz2);
    return 0;
}

/* Force source */
int src_sinc_force_ssg(float *hVx, float *hVz,
                             float *xvec1, float *zvec1,  /* for the integer grids  */
                             float *xvec2, float *zvec2,  /* for the half grids     */
                             struct Src src,
                             float current_time, float dt,
                             float dx, float dz, int nx)
{
    int   i_source, i, k, indx;  // i-th source
    // index of the source location in  grid
    float indx_xs_1, indx_zs_1, indx_xs_2, indx_zs_2;
    int   *ixvec1_src, *ixvec2_src, *izvec1_src, *izvec2_src;
    double *DNx1, *DNx2, *DNz1, *DNz2;
    float V, Fx, Fz, Mxz, d, stf;
    double force_src;

    ixvec1_src = (int*)malloc(2*NSRCEXT*sizeof(int));
    ixvec2_src = (int*)malloc(2*NSRCEXT*sizeof(int));
    izvec1_src = (int*)malloc(2*NSRCEXT*sizeof(int));
    izvec2_src = (int*)malloc(2*NSRCEXT*sizeof(int));

    DNx1 = (double*) malloc(2*NSRCEXT*sizeof(double));
    DNz1 = (double*) malloc(2*NSRCEXT*sizeof(double));
    DNx2 = (double*) malloc(2*NSRCEXT*sizeof(double));
    DNz2 = (double*) malloc(2*NSRCEXT*sizeof(double));


    /* If there i_source no moment source */
    if (src.number_of_src <= 0)
        return;

    /* Loop each moment source */
    for (i_source = 0; i_source < src.number_of_src; i_source++) {

        Fx = src.Fx[i_source];
        Fz = src.Fz[i_source];

        indx_xs_1 = (src.xs[i_source] - xvec1[0])/dx;
        indx_zs_1 = (src.zs[i_source] - zvec1[0])/dz;

        indx_xs_2 = (src.xs[i_source] - xvec2[0])/dx;
        indx_zs_2 = (src.zs[i_source] - zvec2[0])/dz;

        /* calculate index of discrete source */
        for (i = 0; i < NSRCEXT; i++) {
            ixvec1_src[NSRCEXT-i-1] = floor(indx_xs_1) - i;
            ixvec1_src[NSRCEXT+i  ] =  ceil(indx_xs_1) + i;
            izvec1_src[NSRCEXT-i-1] = floor(indx_zs_1) - i;
            izvec1_src[NSRCEXT+i  ] =  ceil(indx_zs_1) + i;

            ixvec2_src[NSRCEXT-i-1] = floor(indx_xs_2) - i;
            ixvec2_src[NSRCEXT+i  ] =  ceil(indx_xs_2) + i;
            izvec2_src[NSRCEXT-i-1] = floor(indx_zs_2) - i;
            izvec2_src[NSRCEXT+i  ] =  ceil(indx_zs_2) + i;
        }

        /* calculate discrete delta function */
        for (i = 0; i < 2*NSRCEXT; i++) {
            DNx1[i] =  calculate_windowed_impulse(ixvec1_src[i]-indx_xs_1, KAISERB, NSRCEXT);
            DNz1[i] =  calculate_windowed_impulse(izvec1_src[i]-indx_zs_1, KAISERB, NSRCEXT);
            DNx2[i] =  calculate_windowed_impulse(ixvec2_src[i]-indx_xs_2, KAISERB, NSRCEXT);
            DNz2[i] =  calculate_windowed_impulse(izvec2_src[i]-indx_zs_2, KAISERB, NSRCEXT);
        }

        /* the factor of the i_source-th source */
        stf = cal_source_time_function(src.stf_type_id[i_source], \
                    current_time, src.stf_timefactor[i_source], src.stf_freqfactor[i_source]);

        /* add source */
        for (i = 0; i < 2*NSRCEXT; i++) {
            for (k = 0; k < 2*NSRCEXT; k++ ) {

                V = dx * dz;

                //-- for points 01
                indx = izvec1_src[k] * nx + ixvec2_src[i];
                d = DNx2[i] * DNz1[k];        /* d: Hicks, 2002 */
                force_src = d*stf/V/dt;
                //printf("%f%f %f %f\n",d, force_src, DNx2[i], DNz1[k]);
                hVx[indx] += Fx * force_src;

                //-- for points 10
                indx = izvec2_src[k] * nx + ixvec1_src[i];
                d = DNx1[i] * DNz2[k];
                force_src = d*stf/V/dt;
                //printf("%f%f %f   %f\n",d, force_src, DNx1[i], DNz2[k]);
                hVz[indx] += Fz * force_src;

            }
        }
    }

    free(ixvec1_src); free(ixvec2_src);
    free(izvec1_src); free(izvec2_src);
    free(DNx1); free(DNx2); free(DNz1); free(DNz2);

    return 0;

}
