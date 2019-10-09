/***************************************************************************
 *
 * This function is used for selective filter
 *
 * Authors: Luqian Jiang   Email: jlq0608@gmail.com
 *          Wei Zhang      Email: zhangwei@sustc.edu.cn
 *
 * Copyright (c) 2018 - 2019 zwlab
 *
 * History: 08/2019: Original version created by Luqian Jiang
 *
 ***************************************************************************/

#include "elastic2d_filter.h"

int filter_velocity(int filter_method, int half_fd_stencil, int nx,
                    int ni1, int ni2, int nk1, int nk2,
                    float *Vx_1, float *Vx_2, float *Vz_1, float *Vz_2)
{
    float *f_Coe=NULL;
    int  half_filter_stencil;

    f_Coe = (float*)malloc(half_fd_stencil*sizeof(float));
    filter_coe(filter_method, half_fd_stencil, f_Coe, &half_filter_stencil);

    /* Vx */
    sf01_45(half_filter_stencil, f_Coe, Vx_1, Vx_2, nx, ni1, ni2, nk1, nk2);
    sf10_45(half_filter_stencil, f_Coe, Vx_1, Vx_2, nx, ni1, ni2, nk1, nk2);
    //sf01_135(half_filter_stencil, f_Coe, Vx_1, Vx_2, nx, ni1, ni2, nk1, nk2);
    //sf10_135(half_filter_stencil, f_Coe, Vx_1, Vx_2, nx, ni1, ni2, nk1, nk2);

    /* Vz */
    sf10_45(half_filter_stencil, f_Coe, Vz_2, Vz_1, nx, ni1, ni2, nk1, nk2);
    sf01_45(half_filter_stencil, f_Coe, Vz_2, Vz_1, nx, ni1, ni2, nk1, nk2);
    //sf10_135(half_filter_stencil, f_Coe, Vz_2, Vz_1, nx, ni1, ni2, nk1, nk2);
    //sf01_135(half_filter_stencil, f_Coe, Vz_2, Vz_1, nx, ni1, ni2, nk1, nk2);

    free(f_Coe);
    return 0;
}

int filter_stresses(int filter_method, int half_fd_stencil, int nx,
                    int ni1, int ni2, int nk1, int nk2,
                    float *Txx_1, float *Tzz_1, float *Txz_1,
                    float *Txx_2, float *Tzz_2, float *Txz_2)
{
    float *f_Coe=NULL;
    int half_filter_stencil;


    f_Coe = (float*)malloc(half_fd_stencil*sizeof(float));
    filter_coe(filter_method, half_fd_stencil, f_Coe, &half_filter_stencil);

    /* Txx */
//    sf00_45(half_filter_stencil, f_Coe, Txx_1, Txx_2, nx, ni1, ni2, nk1, nk2);
//    sf11_45(half_filter_stencil, f_Coe, Txx_1, Txx_2, nx, ni1, ni2, nk1, nk2);
//    sf00_135(half_filter_stencil, f_Coe, Txx_1, Txx_2, nx, ni1, ni2, nk1, nk2);
//    sf11_135(half_filter_stencil, f_Coe, Txx_1, Txx_2, nx, ni1, ni2, nk1, nk2);
//
//    /* Tzz */
//    sf00_45(half_filter_stencil, f_Coe, Tzz_1, Tzz_2, nx, ni1, ni2, nk1, nk2);
//    sf11_45(half_filter_stencil, f_Coe, Tzz_1, Tzz_2, nx, ni1, ni2, nk1, nk2);
//    sf00_135(half_filter_stencil, f_Coe, Tzz_1, Tzz_2, nx, ni1, ni2, nk1, nk2);
//    sf11_135(half_filter_stencil, f_Coe, Tzz_1, Tzz_2, nx, ni1, ni2, nk1, nk2);
//
//    /* Txz */
//    sf00_45(half_filter_stencil, f_Coe, Txz_2, Txz_1, nx, ni1, ni2, nk1, nk2);
//    sf11_45(half_filter_stencil, f_Coe, Txz_2, Txz_1, nx, ni1, ni2, nk1, nk2);
//    sf00_135(half_filter_stencil, f_Coe, Txz_2, Txz_1, nx, ni1, ni2, nk1, nk2);
//    sf11_135(half_filter_stencil, f_Coe, Txz_2, Txz_1, nx, ni1, ni2, nk1, nk2);

    free(f_Coe);
    return 0;
}


/* 1D selective filtering for (i, j) points */
int sf00_45(int half_sf_stencil, float *sf, float *u00, float *u11,
        int nx, int ni1, int ni2, int nk1, int nk2)
{
    int ix, iz, ip, index;
    float Du;

    for (iz = nk1+half_sf_stencil; iz < nk2-half_sf_stencil-1; iz++) {
        for (ix = ni1+half_sf_stencil; ix < ni2-half_sf_stencil-1; ix++) {

            Du = sf[0] * u00[iz*nx+ix];
            for (ip = 1; ip <= half_sf_stencil; ip = ip+2)
            {
                index = (ip+1)/2;
                /*-- 45 deg --*/
                Du += sf[ip] * ( u11[ ( iz+(index-1) )*nx + ix+(index-1) ]
                               + u11[ ( iz-(index  ) )*nx + ix-(index  ) ] );
                if (ip < half_sf_stencil) {
                    Du += sf[ip+1] * ( u00[ (iz+(index))*nx + ix+(index) ]
                                     + u00[ (iz-(index))*nx + ix-(index) ] );
                }
            }

            u00[iz*nx + ix] -= SIGMAD * Du;
        }
    }
    return 0;
}

/* 2D selective filtering for (i, j) points */
int sf00_135(int half_sf_stencil, float *sf, float *u00, float *u11,
        int nx, int ni1, int ni2, int nk1, int nk2)
{
    int ix, iz, ip, index;
    float Du;

    for (iz = nk1+half_sf_stencil; iz < nk2-half_sf_stencil-1; iz++) {
        for (ix = ni1+half_sf_stencil; ix < ni2-half_sf_stencil-1; ix++) {

            Du = sf[0] * u00[iz*nx+ix];
            for (ip = 1; ip <= half_sf_stencil; ip = ip+2)
            {
                index = (ip+1)/2;

                /*-- 135 deg --*/
                Du += sf[ip] * ( u11[ (iz-(index  ))*nx + ix+(index-1) ]
                               + u11[ (iz+(index-1))*nx + ix-(index  ) ]);
                if (ip < half_sf_stencil) {
                    Du += sf[ip+1] * ( u00[ (iz+(index))*nx + ix-(index) ]
                                     + u00[ (iz-(index))*nx + ix+(index) ] );
                }

            }

            u00[iz*nx + ix] -= SIGMAD * Du;
        }
    }
    return 0;
}


/* 1D selective filtering for (i+1/2, j+1/2) points */
int sf11_45(int half_sf_stencil, float *sf, float *u00, float *u11,
         int nx, int ni1, int ni2, int nk1, int nk2)
{
    int ix, iz, ip, index;
    float Du;

    for (iz = nk1+half_sf_stencil; iz < nk2-half_sf_stencil-1; iz++) {
        for (ix = ni1+half_sf_stencil; ix < ni2-half_sf_stencil-1; ix++) {

            Du = sf[0] * u11[iz*nx+ix];
            for (ip = 1; ip <= half_sf_stencil; ip = ip+2)
            {
                index = (ip+1)/2;
                /*-- 45 deg --*/
                Du += sf[ip] * ( u00[ ( iz+(index  ) )*nx + ix+(index  ) ]
                               + u00[ ( iz-(index-1) )*nx + ix-(index-1) ] );
                if (ip < half_sf_stencil) {
                    Du += sf[ip+1] * ( u11[ (iz+(index))*nx + ix+(index) ]
                                     + u11[ (iz-(index))*nx + ix-(index) ] );
                }
            }
            u11[iz*nx + ix] -= SIGMAD * Du;
        }
    }

    return 0;
}

/* 1D selective filtering for (i+1/2, j+1/2) points */
int sf11_135(int half_sf_stencil, float *sf, float *u00, float *u11,
         int nx, int ni1, int ni2, int nk1, int nk2)
{
    int ix, iz, ip, index;
    float Du;

    for (iz = nk1+half_sf_stencil; iz < nk2-half_sf_stencil-1; iz++) {
        for (ix = ni1+half_sf_stencil; ix < ni2-half_sf_stencil-1; ix++) {

            Du = sf[0] * u11[iz*nx+ix];
            for (ip = 1; ip <= half_sf_stencil; ip = ip+2)
            {
                index = (ip+1)/2;

                /*-- 135 deg --*/
                Du += sf[ip] * ( u00[ (iz-(index-1))*nx + ix+(index  ) ]
                               + u00[ (iz+(index  ))*nx + ix-(index-1) ]);
                if (ip < half_sf_stencil) {
                    Du += sf[ip+1] * ( u11[ (iz+(index))*nx + ix-(index) ]
                                     + u11[ (iz-(index))*nx + ix+(index) ] );
                }
            }

            u11[iz*nx + ix] -= SIGMAD * Du;
        }
    }
    return 0;
}

/* 1D selective filtering for (i+1/2, j) points */
int sf01_45(int half_sf_stencil, float *sf, float *u01, float *u10,
         int nx, int ni1, int ni2, int nk1, int nk2)
{
    int ix, iz, ip, index;
    float Du;

    for (iz = nk1+half_sf_stencil; iz < nk2-half_sf_stencil-1; iz++) {
        for (ix = ni1+half_sf_stencil; ix < ni2-half_sf_stencil-1; ix++) {

            Du = sf[0] * u01[iz*nx+ix];
            for (ip = 1; ip <= half_sf_stencil; ip = ip+2)
            {
                index = (ip+1)/2;
                /*-- 45 deg --*/
                Du += sf[ip] * ( u10[ ( iz+(index-1) )*nx + ix+(index  ) ]
                               + u10[ ( iz-(index  ) )*nx + ix-(index-1) ] );
                if (ip < half_sf_stencil) {
                    Du += sf[ip+1] * ( u01[ (iz+(index))*nx + ix+(index) ]
                                     + u01[ (iz-(index))*nx + ix-(index) ] );
                }

            }

            u01[iz*nx + ix] -= SIGMAD * Du;
        }
    }

    return 0;
}

/* 1D selective filtering for (i+1/2, j) points */
int sf01_135(int half_sf_stencil, float *sf, float *u01, float *u10,
         int nx, int ni1, int ni2, int nk1, int nk2)
{
    int ix, iz, ip, index;
    float Du;

    for (iz = nk1+half_sf_stencil; iz < nk2-half_sf_stencil-1; iz++) {
        for (ix = ni1+half_sf_stencil; ix < ni2-half_sf_stencil-1; ix++) {

            //printf("%d %d %f", ix, iz, sf[0]);

            Du = sf[0] * u01[iz*nx+ix];
            for (ip = 1; ip <= half_sf_stencil; ip = ip+2)
            {
                index = (ip+1)/2;

                /*-- 135 deg --*/
                Du += sf[ip] * ( u10[ (iz-(index  ))*nx + ix+(index  ) ]
                               + u10[ (iz+(index-1))*nx + ix-(index-1) ]);
                if (ip < half_sf_stencil) {
                    Du += sf[ip+1] * ( u01[ (iz+(index))*nx + ix-(index) ]
                                     + u01[ (iz-(index))*nx + ix+(index) ] );
                }

            }

            u01[iz*nx + ix] -= SIGMAD * Du;
        }
    }
    return 0;
}


/* 1D selective filtering for (i, j+1/2) points */
int sf10_45(int half_sf_stencil, float *sf, float *u01, float *u10,
         int nx, int ni1, int ni2, int nk1, int nk2)
{
    int ix, iz, ip, index;
    float Du;

    for (iz = nk1+half_sf_stencil; iz < nk2-half_sf_stencil-1; iz++) {
        for (ix = ni1+half_sf_stencil; ix < ni2-half_sf_stencil-1; ix++) {

            Du = sf[0] * u10[iz*nx + ix];
            for (ip = 1; ip <= half_sf_stencil; ip = ip+2)
            {
                index = (ip+1)/2;
                /*-- 45 deg --*/
                Du += sf[ip] * ( u01[ ( iz+(index  ) )*nx + ix+(index-1) ]
                               + u01[ ( iz-(index-1) )*nx + ix-(index  ) ] );
                if (ip < half_sf_stencil) {
                    Du += sf[ip+1] * ( u10[ (iz+(index))*nx + ix+(index) ]
                                     + u10[ (iz-(index))*nx + ix-(index) ] );
                }

            }
            u10[iz*nx + ix] -= SIGMAD * Du;
        }
    }
    return 0;
}

/* 1D selective filtering for (i, j+1/2) points */
int sf10_135(int half_sf_stencil, float *sf, float *u01, float *u10,
         int nx, int ni1, int ni2, int nk1, int nk2)
{
    int ix, iz, ip, index;
    float Du;

    for (iz = nk1+half_sf_stencil; iz < nk2-half_sf_stencil-1; iz++) {
        for (ix = ni1+half_sf_stencil; ix < ni2-half_sf_stencil-1; ix++) {

            Du = sf[0] * u10[iz*nx + ix];
            for (ip = 1; ip <= half_sf_stencil; ip = ip+2)
            {
                index = (ip+1)/2;

                /*-- 135 deg --*/
                Du += sf[ip] * ( u01[ (iz-(index-1))*nx + ix+(index-1) ]
                               + u01[ (iz+(index  ))*nx + ix-(index  ) ]);
                if (ip < half_sf_stencil) {
                    Du += sf[ip+1] * ( u10[ (iz+(index))*nx + ix-(index) ]
                                     + u10[ (iz-(index))*nx + ix+(index) ] );
                }

            }
            u10[iz*nx + ix] -= SIGMAD * Du;
        }
    }
    return 0;
}

/* Calculate the filter coefficient, easy to change other filters */
int filter_coe(int filter_method, int half_fd_stencil, float *f_Coe, int *half_sf_stencil)
{

    switch (filter_method) {
        case FILTER_SFs:
            SFsCoe(half_fd_stencil, f_Coe, half_sf_stencil);
            break;
        case FILTER_SFo:
            SFsCoe(half_fd_stencil, f_Coe, half_sf_stencil);
            break;
        default:
            fprintf(stderr, "Please check filter method!\n" );
            return;
    }
    return 0;
}

/* standard high-order filter */
int SFsCoe(int half_fd_stencil, float *SFs, int *half_sf_stencil)
{
    if (half_fd_stencil <= 4)
    {
        *half_sf_stencil = 4;
        SFs = (float*)realloc(SFs, (*half_sf_stencil+1)*sizeof(float));
        SFs[0] = 35.0/128.0;
        SFs[1] = -7.0/32.0;
        SFs[2] = 7.0/64.0;
        SFs[3] = -1.0/32.0;
        SFs[4] = 1.0/256.0;
    }
    else if (half_fd_stencil <= 5)
    {
        *half_sf_stencil = 5;
        SFs = (float*)realloc(SFs, (*half_sf_stencil+1)*sizeof(float));
        SFs[0] =   63.0/256.0;
        SFs[1] = -105.0/512.0;
        SFs[2] =   15.0/128.0;
        SFs[3] =  -45.0/1024.0;
        SFs[4] =    5.0/512.0;
        SFs[5] =  - 1.0/1024.0;

    }
    else if (half_fd_stencil >= 6)
    {
        *half_sf_stencil = 6;
        SFs = (float*)realloc(SFs, (*half_sf_stencil+1)*sizeof(float));
        SFs[0] = 231.0/1024.0;
        SFs[1] = -99.0/ 512.0;
        SFs[2] = 495.0/4096.0;
        SFs[3] = -55.0/1024.0;
        SFs[4] =  33.0/2048.0;
        SFs[5] = - 3.0/1024.0;
        SFs[6] =   1.0/4096.0;

    }
    return 0;
}

/* optimized selective filter */
int SFoCoe(int half_fd_stencil, float *SFo, int *half_sf_stencil)
{
    if (half_fd_stencil <= 4)
    {
        *half_sf_stencil = 4;
        SFo = (float*)realloc(SFo, (*half_sf_stencil+1)*sizeof(float));
        SFo[0] =  0.243527493120;
        SFo[1] = -0.201788880640;
        SFo[2] =  0.120007591680;
        SFo[3] = -0.045211119360;
        SFo[4] =  0.008228661760;

    }
    else if (half_fd_stencil <= 5)
    {
        *half_sf_stencil = 5;
        SFo = (float*)realloc(SFo, (*half_sf_stencil+1)*sizeof(float));
        SFo[0] =  0.215044884112;
        SFo[1] = -0.187772883589;
        SFo[2] =  0.123755948787;
        SFo[3] = -0.059227575576;
        SFo[4] =  0.018721609157;
        SFo[5] = -0.002999540835;

    }
    else if (half_fd_stencil >= 6)
    {
        *half_sf_stencil = 6;
        SFo = (float*)realloc(SFo, (*half_sf_stencil+1)*sizeof(float));
        SFo[0] =  0.190899511506;
        SFo[1] = -0.171503832236;
        SFo[2] =  0.123632891797;
        SFo[3] = -0.069975429105;
        SFo[4] =  0.029662754736;
        SFo[5] = -0.008520738659;
        SFo[6] =  0.001254597714;

    }

    return 0;

}


///* 2D selective filtering for (i, j) points */
//int sf00(int half_sf_stencil, float *sf, float *u00, float *u11,
//        int nx, int ni1, int ni2, int nk1, int nk2)
//{
//    int ix, iz, ip, index;
//    float Du;
//
//    for (iz = nk1+half_sf_stencil; iz < nk2-half_sf_stencil-1; iz++) {
//        for (ix = ni1+half_sf_stencil; ix < ni2-half_sf_stencil-1; ix++) {
//
//            Du = sf[0] * u00[iz*nx+ix];
//            for (ip = 1; ip <= half_sf_stencil; ip = ip+2)
//            {
//                index = (ip+1)/2;
//                /*-- 45 deg --*/
//                Du += sf[ip] * ( u11[ ( iz+(index-1) )*nx + ix+(index-1) ]
//                               + u11[ ( iz-(index  ) )*nx + ix-(index  ) ] );
//                if (ip < half_sf_stencil) {
//                    Du += sf[ip+1] * ( u00[ (iz+(index))*nx + ix+(index) ]
//                                     + u00[ (iz-(index))*nx + ix-(index) ] );
//                }
//                /*-- 135 deg --*/
//                Du += sf[ip] * ( u11[ (iz-(index  ))*nx + ix+(index-1) ]
//                               + u11[ (iz+(index-1))*nx + ix-(index  ) ]);
//                if (ip < half_sf_stencil) {
//                    Du += sf[ip+1] * ( u00[ (iz+(index))*nx + ix-(index) ]
//                                     + u00[ (iz-(index))*nx + ix+(index) ] );
//                }
//
//            }
//
//            u00[iz*nx + ix] -= SIGMAD * Du;
//        }
//    }
//    return 0;
//}
//
///* 2D selective filtering for (i+1/2, j+1/2) points */
//int sf11(int half_sf_stencil, float *sf, float *u00, float *u11,
//         int nx, int ni1, int ni2, int nk1, int nk2)
//{
//    int ix, iz, ip, index;
//    float Du;
//
//    for (iz = nk1+half_sf_stencil; iz < nk2-half_sf_stencil-1; iz++) {
//        for (ix = ni1+half_sf_stencil; ix < ni2-half_sf_stencil-1; ix++) {
//
//            Du = sf[0] * u11[iz*nx+ix];
//            for (ip = 1; ip <= half_sf_stencil; ip = ip+2)
//            {
//                index = (ip+1)/2;
//                /*-- 45 deg --*/
//                Du = Du + sf[ip] * ( u00[ ( iz+(index  ) )*nx + ix+(index  ) ]
//                                   + u00[ ( iz-(index-1) )*nx + ix-(index-1) ] );
//                if (ip < half_sf_stencil) {
//                    Du = Du + sf[ip+1] * ( u11[ (iz+(index))*nx + ix+(index) ]
//                                         + u11[ (iz-(index))*nx + ix-(index) ] );
//                }
//
//                /*-- 135 deg --*/
//                Du = Du + sf[ip] * ( u00[ (iz-(index-1))*nx + ix+(index  ) ]
//                                   + u00[ (iz+(index  ))*nx + ix-(index-1) ]);
//                if (ip < half_sf_stencil) {
//                    Du = Du + sf[ip+1] * ( u11[ (iz+(index))*nx + ix-(index) ]
//                                         + u11[ (iz-(index))*nx + ix+(index) ] );
//                }
//
//            }
//
//            u11[iz*nx + ix] = u11[iz*nx + ix] - SIGMAD * Du;
//        }
//    }
//    return 0;
//}
//
///* 2D selective filtering for (i+1/2, j) points */
//int sf01(int half_sf_stencil, float *sf, float *u01, float *u10,
//         int nx, int ni1, int ni2, int nk1, int nk2)
//{
//    int ix, iz, ip, index;
//    float Du;
//
//    for (iz = nk1+half_sf_stencil; iz < nk2-half_sf_stencil-1; iz++) {
//        for (ix = ni1+half_sf_stencil; ix < ni2-half_sf_stencil-1; ix++) {
//
//            //printf("%d %d %f", ix, iz, sf[0]);
//
//            Du = sf[0] * u01[iz*nx+ix];
//            for (ip = 1; ip <= half_sf_stencil; ip = ip+2)
//            {
//                index = (ip+1)/2;
//                /*-- 45 deg --*/
//                Du = Du + sf[ip] * ( u10[ ( iz+(index-1) )*nx + ix+(index  ) ]
//                                   + u10[ ( iz-(index  ) )*nx + ix-(index-1) ] );
//                if (ip < half_sf_stencil) {
//                    Du = Du + sf[ip+1] * ( u01[ (iz+(index))*nx + ix+(index) ]
//                                         + u01[ (iz-(index))*nx + ix-(index) ] );
//                }
//
//                /*-- 135 deg --*/
//                Du = Du + sf[ip] * ( u10[ (iz-(index  ))*nx + ix+(index  ) ]
//                                   + u10[ (iz+(index-1))*nx + ix-(index-1) ]);
//                if (ip < half_sf_stencil) {
//                    Du = Du + sf[ip+1] * ( u01[ (iz+(index))*nx + ix-(index) ]
//                                         + u01[ (iz-(index))*nx + ix+(index) ] );
//                }
//
//            }
//
//            u01[iz*nx + ix] = u01[iz*nx + ix] - SIGMAD * Du;
//        }
//    }
//    return 0;
//}
//
///* 2D selective filtering for (i, j+1/2) points */
//int sf10(int half_sf_stencil, float *sf, float *u01, float *u10,
//         int nx, int ni1, int ni2, int nk1, int nk2)
//{
//    int ix, iz, ip, index;
//    float Du;
//
//    for (iz = nk1+half_sf_stencil; iz < nk2-half_sf_stencil-1; iz++) {
//        for (ix = ni1+half_sf_stencil; ix < ni2-half_sf_stencil-1; ix++) {
//
//            Du = sf[0] * u10[iz*nx + ix];
//            for (ip = 1; ip <= half_sf_stencil; ip = ip+2)
//            {
//                index = (ip+1)/2;
//                /*-- 45 deg --*/
//                Du = Du + sf[ip] * ( u01[ ( iz+(index  ) )*nx + ix+(index-1) ]
//                                   + u01[ ( iz-(index-1) )*nx + ix-(index  ) ] );
//                if (ip < half_sf_stencil) {
//                    Du = Du + sf[ip+1] * ( u10[ (iz+(index))*nx + ix+(index) ]
//                                         + u10[ (iz-(index))*nx + ix-(index) ] );
//                }
//
//                /*-- 135 deg --*/
//                Du = Du + sf[ip] * ( u01[ (iz-(index-1))*nx + ix+(index-1) ]
//                                   + u01[ (iz+(index  ))*nx + ix-(index  ) ]);
//                if (ip < half_sf_stencil) {
//                    Du = Du + sf[ip+1] * ( u10[ (iz+(index))*nx + ix-(index) ]
//                                         + u10[ (iz-(index))*nx + ix+(index) ] );
//                }
//
//            }
//            u10[iz*nx + ix] = u10[iz*nx + ix] - SIGMAD * Du;
//        }
//    }
//    return 0;
//}
