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


int filter_velocity(int half_fd_stencil, int nx, int nz,
                    float *Vx_1, float *Vx_2, float *Vz_1, float *Vz_2)
{
    float *f_Coe, half_filter_stencil

    half_filter_stencil = filter_coe(half_fd_stencil, f_Coe);

    /* Vx */
    sf01(half_filter_stencil, f_Coe, Vx_1, Vx_2, nx, nz);
    sf10(half_filter_stencil, f_Coe, Vx_1, Vx_2, nx, nz);

    /* Vz */
    sf10(half_filter_stencil, f_Coe, Vz_2, Vz_1, nx, nz);
    sf01(half_filter_stencil, f_Coe, Vz_2, Vz_1, nx, nz);

    return 0;
}

int filter_stresses(int half_fd_stencil, int nx, int nz,
                    float *Txx_1, float *Tzz_1, float *Txz_1,
                    float *Txx_2, float *Tzz_2, float *Txz_2)
{
    float *f_Coe, half_filter_stencil

    half_filter_stencil = filter_coe(half_fd_stencil, f_Coe);

    /* Txx */
    sf00(half_filter_stencil, f_Coe, Txx_1, Txx_2, nx, nz);
    sf11(half_filter_stencil, f_Coe, Txx_1, Txx_2, nx, nz);

    /* Tzz */
    sf00(half_filter_stencil, f_Coe, Tzz_1, Tzz_2, nx, nz);
    sf11(half_filter_stencil, f_Coe, Tzz_1, Tzz_2, nx, nz);

    /* Txz */
    sf00(half_filter_stencil, f_Coe, Txz_2, Txz_1, nx, nz);
    sf11(half_filter_stencil, f_Coe, Txz_2, Txz_1, nx, nz);

    return 0;
}

/* Calculate the filter coefficient, easy to change other filter */
int filter_coe(int half_fd_stencil, float *f_Coe)
{
    int half_sf_stencil;

    half_sf_stencil = SFoCoe(half_fd_stencil, f_Coe);

    return half_sf_stencil;
}

/* 2D selective filtering for (i, j) points */
int sf00(int half_sf_stencil, float *sf, float *u00, float *u11, int nx, int nz)
{
    int ix, iz, ip, index;
    float Du;

    for (iz = half_sf_stencil; iz < nz-half_sf_stencil; iz++) {
        for (ix = half_sf_stencil; ix < nx-half_sf_stencil; ix++) {

            Du = sf[0] * u00[iz*nx + ix];
            for (ip = 1; ip <= half_sf_stencil; ip = ip+2)
            {
                index = (ip+1)/2;
                /*-- 45 deg --*/
                Du = Du + sf[ip] * ( u11[ ( iz+(index-1) )*nx + ix+(index-1) ]
                                   + u11[ ( iz-(index  ) )*nx + ix-(index  ) ] );
                if (ip < half_sf_stencil) {
                    Du = Du + sf[ip+1] * ( u00[ (iz+(index))*nx + ix+(index) ]
                                         + u00[ (iz-(index))*nx + ix-(index) ] );
                }
                /*-- 135 deg --*/
                Du = Du + sf[ip] * ( u11[ (iz-(index  ))*nx + ix+(index-1) ]
                                   + u11[ (iz+(index-1))*nx + ix-(index  ) ]);
                if (ip < half_sf_stencil) {
                    Du = Du + sf[ip+1] * ( u00[ (iz+(index))*nx + ix-(index) ]
                                         + u00[ (iz-(index))*nx + ix+(index) ] );
                }

            }

            u00[iz*nx + ix] = u00[iz*nx + ix] - SIGMAD * Du;
        }
    }
    return 0;
}

/* 2D selective filtering for (i+1/2, j+1/2) points */
int sf11(int half_sf_stencil, float *sf, float *u00, float *u11, int nx, int nz)
{
    int ix, iz, ip, index;
    float Du;

    for (iz = half_sf_stencil; iz < nz-half_sf_stencil; iz++) {
        for (ix = half_sf_stencil; ix < nx-half_sf_stencil; ix++) {

            Du = sf[0] * u11[iz*nx + ix];
            for (ip = 1; ip <= half_sf_stencil; ip = ip+2)
            {
                index = (ip+1)/2;
                /*-- 45 deg --*/
                Du = Du + sf[ip] * ( u00[ ( iz+(index  ) )*nx + ix+(index  ) ]
                                   + u00[ ( iz-(index-1) )*nx + ix-(index-1) ] );
                if (ip < half_sf_stencil) {
                    Du = Du + sf[ip+1] * ( u11[ (iz+(index))*nx + ix+(index) ]
                                         + u11[ (iz-(index))*nx + ix-(index) ] );
                }

                /*-- 135 deg --*/
                Du = Du + sf[ip] * ( u00[ (iz-(index-1))*nx + ix+(index  ) ]
                                   + u00[ (iz+(index  ))*nx + ix-(index-1) ]);
                if (ip < half_sf_stencil) {
                    Du = Du + sf[ip+1] * ( u11[ (iz+(index))*nx + ix-(index) ]
                                         + u11[ (iz-(index))*nx + ix+(index) ] );
                }

            }

            u11[iz*nx + ix] = u11[iz*nx + ix] - SIGMAD * Du;
        }
    }
    return 0;
}

/* 2D selective filtering for (i+1/2, j) points */
int sf01(int half_sf_stencil, float *sf, float *u01, float *u10, int nx, int nz)
{
    int ix, iz, ip, index;
    float Du;

    for (iz = half_sf_stencil; iz < nz-half_sf_stencil; iz++) {
        for (ix = half_sf_stencil; ix < nx-half_sf_stencil; ix++) {

            Du = sf[0] * u01[iz*nx + ix];
            for (ip = 1; ip <= half_sf_stencil; ip = ip+2)
            {
                index = (ip+1)/2;
                /*-- 45 deg --*/
                Du = Du + sf[ip]   * ( u10[ ( iz+(index-1) )*nx + ix+(index  ) ]
                                     + u10[ ( iz-(index  ) )*nx + ix-(index-1) ] );
                if (ip < half_sf_stencil) {
                    Du = Du + sf[ip+1] * ( u01[ (iz+(index))*nx + ix+(index) ]
                                         + u01[ (iz-(index))*nx + ix-(index) ] );
                }

                /*-- 135 deg --*/
                Du = Du + sf[ip]   * ( u10[ (iz-(index  ))*nx + ix+(index  ) ]
                                     + u10[ (iz+(index-1))*nx + ix-(index-1) ]);
                if (ip < half_sf_stencil) {
                    Du = Du + sf[ip+1] * ( u01[ (iz+(index))*nx + ix-(index) ]
                                         + u01[ (iz-(index))*nx + ix+(index) ] );
                }

            }

            u01[iz*nx + ix] = u01[iz*nx + ix] - SIGMAD * Du;
        }
    }
    return 0;
}

/* 2D selective filtering for (i, j+1/2) points */
int sf10(int half_sf_stencil, float *sf, float *u01, float *u10, int nx, int nz)
{
    int ix, iz, ip, index;
    float Du;

    for (iz = half_sf_stencil; iz < nz-half_sf_stencil; iz++) {
        for (ix = half_sf_stencil; ix < nx-half_sf_stencil; ix++) {

            Du = sf[0] * u10[iz*nx + ix];
            for (ip = 1; ip <= half_sf_stencil; ip = ip+2)
            {
                index = (ip+1)/2;
                /*-- 45 deg --*/
                Du = Du + sf[ip]   * ( u01[ ( iz+(index  ) )*nx + ix+(index-1) ]
                                     + u01[ ( iz-(index-1) )*nx + ix-(index  ) ] );
                if (ip < half_sf_stencil) {
                    Du = Du + sf[ip+1] * ( u10[ (iz+(index))*nx + ix+(index) ]
                                         + u10[ (iz-(index))*nx + ix-(index) ] );
                }

                /*-- 135 deg --*/
                Du = Du + sf[ip]   * ( u01[ (iz-(index-1))*nx + ix+(index-1) ]
                                     + u01[ (iz+(index  ))*nx + ix-(index  ) ]);
                if (ip < half_sf_stencil) {
                    Du = Du + sf[ip+1] * ( u10[ (iz+(index))*nx + ix-(index) ]
                                         + u10[ (iz-(index))*nx + ix+(index) ] );
                }

            }

            u10[iz*nx + ix] = u10[iz*nx + ix] - SIGMAD * Du;
        }
    }
    return 0;
}

/* standard high-order filter */
int SFsCoe(int half_fd_stencil, float *SFs)
{
    int half_sf_stencil;
    if (half_fd_stencil <= 4)
    {
        half_sf_stencil = 4;
        SFs = (float*) malloc((half_sf_stencil+1) * sizeof(float));
        SFs = {35/128, -7/32, 7/64, -1/32, 1/256};

    }
    else if (half_fd_stencil <= 5)
    {
        half_sf_stencil = 5;
        SFs = (float*) malloc((half_sf_stencil+1) * sizeof(float));
        SFs = {63/256, -105/512, 15/128, -45/1024, 5/512, -1/1024};

    }
    else if (half_fd_stencil >= 6)
    {
        half_sf_stencil = 6;
        SFs = (float*) malloc((half_sf_stencil+1) * sizeof(float));
        SFs = {231/1024, -99/512, 495/4096, -55/1024, 33/2048, -3/1024, 1/4096};

    }
    return half_sf_stencil;
}

/* optimized selective filter */
int SFoCoe(int half_fd_stencil, float *SFo)
{
    int half_sf_stencil;
    if (half_fd_stencil <= 4)
    {
        half_sf_stencil = 4;
        SFo = (float*) malloc((half_sf_stencil+1) * sizeof(float));
        SFo = { 0.243527493120,
               -0.201788880640,
                0.120007591680,
               -0.045211119360,
                0.008228661760};

    }
    else if (half_fd_stencil <= 5)
    {
        half_sf_stencil = 5;
        SFo = (float*) malloc((half_sf_stencil+1) * sizeof(float));
        SFo = { 0.215044884112,
               -0.187772883589,
                0.123755948787,
               -0.059227575576,
                0.018721609157,
               -0.002999540835};

    }
    else if (half_fd_stencil >= 6)
    {
        half_sf_stencil = 6;
        SFo = (float*) malloc((half_sf_stencil+1) * sizeof(float));
        SFo = { 0.190899511506,
               -0.171503832236,
                0.123632891797,
               -0.069975429105,
                0.029662754736,
               -0.008520738659,
                0.001254597714};

    }

    return half_sf_stencil;

}
