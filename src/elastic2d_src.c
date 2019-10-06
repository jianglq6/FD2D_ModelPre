/***************************************************************************
 *
 * Thi_source function i_source used to add source in Lebedev grid scheme
 *
 * Authors: Luqian Jiang   Email: jlq0608@gmail.com
 *          Wei Zhang      Email: zhangwei@sustc.edu.cn
 *
 * Copyright (c) 2018 - 2019 zwlab
 *
 * History: 09/2019: Original version created by Luqian Jiang
 *          10/2019: Add some source impulse method
 *
 ***************************************************************************/
#include "elastic2d_src.h"
#include "elastic2d_stf.h"
#include "elastic2d_math.h"


/* Moment source */
int src_tti_moment(float *hTxx_1, float *hTxx_2,
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

    /* source is just the first set of grid */
    /* !!! Just for test */
    for (i_source = 0; i_source < src.number_of_src; i_source++)
    {
        /* source time function */
        stf = cal_source_time_function(src.stf_type_id[i_source], current_time,
                    src.stf_timefactor[i_source], src.stf_freqfactor[i_source]);

        if (source_impulse_method == IMPULSE_DIRECT) {
            indx_xs = round( (src.xs[i_source]-xvec1[0])/dx );
            indx_zs = round( (src.zs[i_source]-zvec1[0])/dz );
            indx_source = indx_zs*nx+indx_xs;
            hTxx_1[indx_source] -=  src.Mxx[i_source]/dx/dz * stf;
            hTzz_1[indx_source] -=  src.Mzz[i_source]/dx/dz * stf;

            if (i_source==0)
            /* write source time function */
                fprintf(fp1, "%12.6f %20.6f\n", current_time, stf);

        }
        else if (source_impulse_method == IMPULSE_AVERAGE) {

        }
        else if (source_impulse_method == IMPULSE_SINC) {

        }
    }

    fclose(fp1);

}

/* Force source */
int src_tti_force(float *hVx_1, float *hVz_1,
                  float *hVx_2, float *hVz_2,
                  float *xvec1, float *zvec1,  /* for the integer grids  */
                  float *xvec2, float *zvec2,  /* for the half grids     */
                  int source_impulse_method, struct Src src,
                  float current_time, float dt,
                  float dx, float dz, int nx,
                  float *B01, float *B10)
{

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
float calculate_windowed_impulse(float x, float b, float r)
{
    float y = 0.0;
    if (1.0-x*x/r*r > 0.0) {
    /* Kai_sourceer finite impulse response filters. */
        y = first_modified_Bessel( 0, b*sqrt(1-x*x/r*r) ) / first_modified_Bessel(0, b);
    /* sinc */
        y = y * sin(PI*x) / (PI*x);
    }
    else
        y = 0.0;

    return y;
}

/* Moment source */
//int src_tti_moment(float *Txx_1, float *Txx_2,
//                   float *Txz_1, float *Txz_2,
//                   float *Tzz_1, float *Tzz_2,
//                   float *xvec1, float *zvec1,
//                   float *xvec2, float *zvec2,
//                   int source_impulse_method, struct Src src,
//                   float current_time, float dt,
//                   float dx, float dz, int nx)
//{
//    int   i_source, i, k, indx;  // i-th source
//    // index of the source location in  grid
//    float indx_xs_1, indx_zs_1, indx_xs_2, indx_zs_2;
//    int   *ixvec1_src=NULL, *ixvec2_src=NULL, *izvec1_src=NULL, *izvec2_src=NULL;
//    float *DNx1=NULL, *DNx2=NULL, *DNz1=NULL, *DNz2=NULL;
//    float moment_src, V, Mxx, Mzz, Mxz, d;
//
//    ixvec1_src = (int*)malloc(2*NSRCEXT*sizeof(int));
//    ixvec2_src = (int*)malloc(2*NSRCEXT*sizeof(int));
//    izvec1_src = (int*)malloc(2*NSRCEXT*sizeof(int));
//    izvec2_src = (int*)malloc(2*NSRCEXT*sizeof(int));
//
//    /* Damp */
//    DNx1 = (float*) malloc(2*NSRCEXT*sizeof(float));
//    DNz1 = (float*) malloc(2*NSRCEXT*sizeof(float));
//    DNx2 = (float*) malloc(2*NSRCEXT*sizeof(float));
//    DNz2 = (float*) malloc(2*NSRCEXT*sizeof(float));
//
//
//    /* If there i_source no moment source */
//    if (src.number_of_src <= 0)
//        return;
//
//    /* Loop each moment source */
//    for (i_source = 0; i_source < src.number_of_src; i_source++) {
//
//        /* the factor of the i_source-th source */
//        moment_src = cal_source_time_function(src.stf_type_id[i_source], \
//                       current_time, src.stf_timefactor[i_source], src.stf_freqfactor[i_source]);
//        Mxx = src.Mxx[i_source];
//        Mzz = src.Mzz[i_source];
//        Mxz = src.Mxz[i_source];
//
//        indx_xs_1 = (src.xs[i_source] - xvec1[0])/dx;
//        indx_zs_1 = (src.zs[i_source] - zvec1[0])/dz;
//
//        indx_xs_2 = (src.xs[i_source] - xvec2[0])/dx;
//        indx_zs_2 = (src.zs[i_source] - zvec2[0])/dz;
//
//        /* calculate index of discrete source */
//        for (i = 0; i < NSRCEXT; i++) {
//            ixvec1_src[NSRCEXT-i-1] = floor(indx_xs_1) - i;
//            ixvec1_src[NSRCEXT+i  ] =  ceil(indx_xs_1) + i;
//            izvec1_src[NSRCEXT-i-1] = floor(indx_zs_1) - i;
//            izvec1_src[NSRCEXT+i  ] =  ceil(indx_zs_1) + i;
//
//            ixvec2_src[NSRCEXT-i-1] = floor(indx_xs_2) - i;
//            ixvec2_src[NSRCEXT+i  ] =  ceil(indx_xs_2) + i;
//            izvec2_src[NSRCEXT-i-1] = floor(indx_zs_2) - i;
//            izvec2_src[NSRCEXT+i  ] =  ceil(indx_zs_2) + i;
//        }
//
//        /* calculate discrete delta function */
//        for (i = 0; i < 2*NSRCEXT; i++) {
//            DNx1[i] =  calculate_windowed_impulse(ixvec1_src[i]-indx_xs_1, KAISERB, NSRCEXT);
//            DNz1[i] =  calculate_windowed_impulse(izvec1_src[i]-indx_zs_1, KAISERB, NSRCEXT);
//            DNx2[i] =  calculate_windowed_impulse(ixvec2_src[i]-indx_xs_2, KAISERB, NSRCEXT);
//            DNz2[i] =  calculate_windowed_impulse(izvec2_src[i]-indx_zs_2, KAISERB, NSRCEXT);
//        }
//
//        /* add source */
//        for (i = 0; i < 2*NSRCEXT; i++) {
//            for (k = 0; k < 2*NSRCEXT; k++ ) {
//
//                V = dx * dz;
//
//                //-- for 00
//                indx = izvec1_src[k] * nx + ixvec1_src[i];
//                d = DNx1[i] * DNz1[k];        /* d: Hicks, 2002 */
//                moment_src = d*moment_src/V * dt;
//                Txx_1[indx] = Txx_1[indx] - Mxx * moment_src;
//                Tzz_1[indx] = Tzz_1[indx] - Mzz * moment_src;
//                Txz_2[indx] = Txz_2[indx] - Mxz * moment_src;
//
//                //-- for 11
//                indx = izvec2_src[k] * nx + ixvec2_src[i];
//                d = DNx2[i] * DNz2[k];
//                moment_src = d*moment_src/V * dt;
//                Txx_2[indx] = Txx_2[indx] - Mxx * moment_src;
//                Tzz_2[indx] = Tzz_2[indx] - Mzz * moment_src;
//                Txz_1[indx] = Txz_1[indx] - Mxz * moment_src;
//
//            }
//        }
//
//    }
//    free(ixvec1_src); free(ixvec2_src);
//    free(izvec1_src); free(izvec2_src);
//    free(DNx1); free(DNx2); free(DNz1); free(DNz2);
//    return 0;
//}
//
///* Force source */
//int src_tti_force(float *Vx_1, float *Vz_1,
//                  float *Vx_2, float *Vz_2,
//                  float *xvec1, float *zvec1,  /* for the integer grids  */
//                  float *xvec2, float *zvec2,  /* for the half grids     */
//                  int source_impulse_method, struct Src src,
//                  float current_time, float dt,
//                  float dx, float dz, int nx,
//                  float *B01, float *B10)
//{
//    int   i_source, i, k, indx;  // i-th source
//    // index of the source location in  grid
//    float indx_xs_1, indx_zs_1, indx_xs_2, indx_zs_2;
//    int   *ixvec1_src, *ixvec2_src, *izvec1_src, *izvec2_src;
//    float *DNx1, *DNx2, *DNz1, *DNz2;
//    float force_src, V, Fx, Fz, Mxz, d;
//
//    ixvec1_src = (int*)malloc(2*NSRCEXT*sizeof(int));
//    ixvec2_src = (int*)malloc(2*NSRCEXT*sizeof(int));
//    izvec1_src = (int*)malloc(2*NSRCEXT*sizeof(int));
//    izvec2_src = (int*)malloc(2*NSRCEXT*sizeof(int));
//
//    DNx1 = (float*) malloc(2*NSRCEXT*sizeof(float));
//    DNz1 = (float*) malloc(2*NSRCEXT*sizeof(float));
//    DNx2 = (float*) malloc(2*NSRCEXT*sizeof(float));
//    DNz2 = (float*) malloc(2*NSRCEXT*sizeof(float));
//
//
//    /* If there i_source no moment source */
//    if (src.number_of_src <= 0)
//        return;
//
//    /* Loop each moment source */
//    for (i_source = 0; i_source < src.number_of_src; i_source++) {
//
//        /* the factor of the i_source-th source */
//        force_src = cal_source_time_function(src.stf_type_id[i_source], \
//                       current_time, src.stf_timefactor[i_source], src.stf_freqfactor[i_source]);
//        Fx = src.Fx[i_source];
//        Fz = src.Fz[i_source];
//
//        indx_xs_1 = (src.xs[i_source] - xvec1[0])/dx;
//        indx_zs_1 = (src.zs[i_source] - zvec1[0])/dz;
//
//        indx_xs_2 = (src.xs[i_source] - xvec2[0])/dx;
//        indx_zs_2 = (src.zs[i_source] - zvec2[0])/dz;
//
//        /* calculate index of discrete source */
//        for (i = 0; i < NSRCEXT; i++) {
//            ixvec1_src[NSRCEXT-i-1] = floor(indx_xs_1) - i;
//            ixvec1_src[NSRCEXT+i  ] =  ceil(indx_xs_1) + i;
//            izvec1_src[NSRCEXT-i-1] = floor(indx_zs_1) - i;
//            izvec1_src[NSRCEXT+i  ] =  ceil(indx_zs_1) + i;
//
//            ixvec2_src[NSRCEXT-i-1] = floor(indx_xs_2) - i;
//            ixvec2_src[NSRCEXT+i  ] =  ceil(indx_xs_2) + i;
//            izvec2_src[NSRCEXT-i-1] = floor(indx_zs_2) - i;
//            izvec2_src[NSRCEXT+i  ] =  ceil(indx_zs_2) + i;
//        }
//
//        /* calculate discrete delta function */
//        for (i = 0; i < 2*NSRCEXT; i++) {
//            DNx1[i] =  calculate_windowed_impulse(ixvec1_src[i]-indx_xs_1, KAISERB, NSRCEXT);
//            DNz1[i] =  calculate_windowed_impulse(izvec1_src[i]-indx_zs_1, KAISERB, NSRCEXT);
//            DNx2[i] =  calculate_windowed_impulse(ixvec2_src[i]-indx_xs_2, KAISERB, NSRCEXT);
//            DNz2[i] =  calculate_windowed_impulse(izvec2_src[i]-indx_zs_2, KAISERB, NSRCEXT);
//        }
//
//        /* add source */
//        for (i = 0; i < 2*NSRCEXT; i++) {
//            for (k = 0; k < 2*NSRCEXT; k++ ) {
//
//                V = dx * dz;
//
//                //-- for points 01
//                indx = izvec1_src[k] * nx + ixvec2_src[i];
//                d = DNx2[i] * DNz1[k];        /* d: Hicks, 2002 */
//                force_src = d*force_src/V * B01[indx] * dt;
//                Vx_1[indx] = Vx_1[indx] - Fx * force_src;
//                Vz_2[indx] = Vz_2[indx] - Fz * force_src;
//
//                //-- for points 10
//                indx = izvec2_src[k] * nx + ixvec1_src[i];
//                d = DNx1[i] * DNz2[k];
//                force_src = d*force_src/V * B10[indx] * dt;
//                Vx_2[indx] = Vx_2[indx] - Fx * force_src;
//                Vz_1[indx] = Vz_1[indx] - Fz * force_src;
//
//            }
//        }
//    }
//
//    free(ixvec1_src); free(ixvec2_src);
//    free(izvec1_src); free(izvec2_src);
//    free(DNx1); free(DNx2); free(DNz1); free(DNz2);
//
//    return 0;
//
//}

