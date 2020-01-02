/***************************************************************************
 *
 * This function is used for the absorbing outgoing waves based on
 *     unsplit-field ADE CFS-PML (auxiliary differential equation complex
 *     frequency shifted PML)
 * (Zhang and Shen. (2010), Geophysics, 75(4):T141-T154)
 *
 * Authors: Luqian Jiang   Email: jlq0608@gmail.com
 *          Wei Zhang      Email: zhangwei@sustc.edu.cn
 *
 * Copyright (c) 2018 - 2020 zwlab
 *
 * History: 12/2019: Original version created by Luqian Jiang
 *
 ***************************************************************************/
#include "elastic2d_math.h"
#include "elastic2d_abs_npml.h"

int abs_npml_velocity_tti(int ni1, int ni2, int nk1, int nk2, int nx, int nz,
                          int *abs_layer_number, float dt,
                          float *hVx01,float *hVx10,float *hVz10,float *hVz01,
                          float *Ax_regu, float *Bx_regu, float *Dx_regu,
                          float *Ax_half, float *Bx_half, float *Dx_half,
                          float *Az_regu, float *Bz_regu, float *Dz_regu,
                          float *Az_half, float *Bz_half, float *Dz_half,
                          float *Txx01_x1, float *Txx01_x2,
                          float *Tzz01_z1, float *Tzz01_z2,
                          float *Txz01_x1, float *Txz01_x2,
                          float *Txz01_z1, float *Txz01_z2,
                          float *Txx10_x1, float *Txx10_x2,
                          float *Tzz10_z1, float *Tzz10_z2,
                          float *Txz10_x1, float *Txz10_x2,
                          float *Txz10_z1, float *Txz10_z2,
                          float *DxTxx01, float *DzTxz01, float *DxTxz01, float *DzTzz01,
                          float *DxTxx10, float *DzTxz10, float *DxTxz10, float *DzTzz10)
{
    int ierr = 0, iz, ix, indx, indx_abs;
    float Txx01_tmp, Txx10_tmp, Tzz10_tmp, Tzz01_tmp, Txz01_tmp, Txz10_tmp;

    /* left: x1 */
    if (abs_layer_number[0] > 0) {
        for (iz = nk1; iz < nk2; iz++) {
            for (ix = ni1; ix <= ni1+abs_layer_number[0]; ix++) {
                indx = iz*nx + ix;
                indx_abs = iz * abs_layer_number[0] + ix-ni1;

                /* for v_x: 01-half; 10-integer */
                Txx01_tmp =  (2.0*Txx01_x1[indx_abs] + dt*Dx_half[ix]/Bx_half[ix]*DxTxx01[indx])
                            /(2.0 + dt*(Ax_half[ix]+Dx_half[ix]/Bx_half[ix]));
                hVx01[indx] += (1.0/Bx_half[ix]-1)*DxTxx01[indx] - Txx01_tmp/Bx_half[ix];
                Txx01_x1[indx_abs] = 2.0*Txx01_tmp - Txx01_x1[indx_abs];
                //printf("hVx01[%d,%d]: %f, Txx01_x1[%d,%d]: %f\n", iz, ix, hVx01[indx], iz, ix, Txx01_tmp);

                Txx10_tmp =  (2.0*Txx10_x1[indx_abs] + dt*Dx_regu[ix]/Bx_regu[ix]*DxTxx10[indx]) \
                            /(2.0 + dt*(Ax_regu[ix]+Dx_regu[ix]/Bx_regu[ix]));
                hVx10[indx] += (1.0/Bx_regu[ix]-1)*DxTxx10[indx] - Txx10_tmp/Bx_regu[ix];
                Txx10_x1[indx_abs] = 2.0*Txx10_tmp - Txx10_x1[indx_abs];

                /* for v_z: 01-half; 10-integer */
                Txz01_tmp =  (2.0*Txz01_x1[indx_abs] + dt*Dx_half[ix]/Bx_half[ix]*DxTxz01[indx]) \
                            /(2.0 + dt*(Ax_half[ix]+Dx_half[ix]/Bx_half[ix]));
                hVz01[indx] += (1.0/Bx_half[ix]-1)*DxTxz01[indx] - Txz01_tmp/Bx_half[ix];
                Txz01_x1[indx_abs] = 2.0*Txz01_tmp - Txz01_x1[indx_abs];

                Txz10_tmp =  (2.0*Txz10_x1[indx_abs] + dt*Dx_regu[ix]/Bx_regu[ix]*DxTxz10[indx]) \
                            /(2.0 + dt*(Ax_regu[ix]+Dx_regu[ix]/Bx_regu[ix]));
                hVz10[indx] += (1.0/Bx_regu[ix]-1)*DxTxz10[indx] - Txz10_tmp/Bx_regu[ix];
                Txz10_x1[indx_abs] = 2.0*Txz10_tmp - Txz10_x1[indx_abs];

            }
        }
    }

    /* right: x2 */
    if (abs_layer_number[1] > 0) {
        for (iz = nk1; iz < nk2; iz++) {
            for (ix = ni2-abs_layer_number[1]-1; ix < ni2; ix++) {
                indx = iz*nx + ix;
                indx_abs = iz * abs_layer_number[1] + ix-(ni2-abs_layer_number[1]-1);

                /* for v_x: 01-half; 10-integer */
                Txx01_tmp =  (2.0*Txx01_x2[indx_abs] + dt*Dx_half[ix]/Bx_half[ix]*DxTxx01[indx]) \
                            /(2.0 + dt*(Ax_half[ix]+Dx_half[ix]/Bx_half[ix]));
                hVx01[indx] += (1.0/Bx_half[ix]-1)*DxTxx01[indx] - Txx01_tmp/Bx_half[ix];
                Txx01_x2[indx_abs] = 2.0*Txx01_tmp - Txx01_x2[indx_abs];

                Txx10_tmp =  (2.0*Txx10_x2[indx_abs] + dt*Dx_regu[ix]/Bx_regu[ix]*DxTxx10[indx]) \
                            /(2.0 + dt*(Ax_regu[ix]+Dx_regu[ix]/Bx_regu[ix]));
                hVx10[indx] += (1.0/Bx_regu[ix]-1)*DxTxx10[indx] - Txx10_tmp/Bx_regu[ix];
                Txx10_x2[indx_abs] = 2.0*Txx10_tmp - Txx10_x2[indx_abs];

                /* for v_z: 01-half; 10-integer */
                Txz01_tmp =  (2.0*Txz01_x2[indx_abs] + dt*Dx_half[ix]/Bx_half[ix]*DxTxz01[indx]) \
                            /(2.0 + dt*(Ax_half[ix]+Dx_half[ix]/Bx_half[ix]));
                hVz01[indx] += (1.0/Bx_half[ix]-1)*DxTxz01[indx] - Txz01_tmp/Bx_half[ix];
                Txz01_x2[indx_abs] = 2.0*Txz01_tmp - Txz01_x2[indx_abs];

                Txz10_tmp =  (2.0*Txz10_x2[indx_abs] + dt*Dx_regu[ix]/Bx_regu[ix]*DxTxz10[indx]) \
                            /(2.0 + dt*(Ax_regu[ix]+Dx_regu[ix]/Bx_regu[ix]));
                hVz10[indx] += (1.0/Bx_regu[ix]-1)*DxTxz10[indx] - Txz10_tmp/Bx_regu[ix];
                Txz10_x2[indx_abs] = 2.0*Txz10_tmp - Txz10_x2[indx_abs];
            }
        }
    }

    /* top: z1 */
    if (abs_layer_number[2] > 0) {
        for (iz = nk1; iz <= nk1+abs_layer_number[2]; iz++) {
            for (ix = ni1; ix < ni2; ix++) {
                indx = iz*nx + ix;
                indx_abs = ix * abs_layer_number[2] + iz-nk1;

                /* for v_x: 01-integer; 10-half */
                Txz01_tmp =  (2.0*Txz01_z1[indx_abs] + dt*Dz_regu[iz]/Bz_regu[iz]*DzTxz01[indx]) \
                            /(2.0 + dt*(Az_regu[iz]+Dz_regu[iz]/Bz_regu[iz]));
                hVx01[indx] += (1.0/Bz_regu[iz]-1)*DzTxz01[indx] - Txz01_tmp/Bz_regu[iz];
                Txz01_z1[indx_abs] = 2.0*Txz01_tmp - Txz01_z1[indx_abs];

                Txz10_tmp =  (2.0*Txz10_z1[indx_abs] + dt*Dz_half[iz]/Bz_half[iz]*DzTxz10[indx]) \
                            /(2.0 + dt*(Az_half[iz]+Dz_half[iz]/Bz_half[iz]));
                hVx10[indx] = hVx10[indx] + (1.0/Bz_half[iz]-1)*DzTxz10[indx] - Txz10_tmp/Bz_half[iz];
                Txz10_z1[indx_abs] = 2.0*Txz10_tmp - Txz10_z1[indx_abs];

                /* for v_z: 01-integer; 10-half */
                Tzz01_tmp =  (2.0*Tzz01_z1[indx_abs] + dt*Dz_regu[iz]/Bz_regu[iz]*DzTzz01[indx]) \
                            /(2.0 + dt*(Az_regu[iz]+Dz_regu[iz]/Bz_regu[iz]));
                hVz01[indx] += (1.0/Bz_regu[iz]-1)*DzTzz01[indx] - Tzz01_tmp/Bz_regu[iz];
                Tzz01_z1[indx_abs] = 2.0*Tzz01_tmp - Tzz01_z1[indx_abs];

                Tzz10_tmp =  (2.0*Tzz10_z1[indx_abs] + dt*Dz_half[iz]/Bz_half[iz]*DzTzz10[indx]) \
                            /(2.0 + dt*(Az_half[iz]+Dz_half[iz]/Bz_half[iz]));
                hVz10[indx] += (1.0/Bz_half[iz]-1)*DzTzz10[indx] - Tzz10_tmp/Bz_half[iz];
                Tzz10_z1[indx_abs] = 2.0*Tzz10_tmp - Tzz10_z1[indx_abs];

            }
        }
    }

    /* bottom: z2 */
    if (abs_layer_number[3] > 0) {
        for (iz = nk2-abs_layer_number[3]-1; iz < nk2; iz++) {
            for (ix = ni1; ix < ni2; ix++) {
                indx = iz*nx + ix;
                indx_abs = ix * abs_layer_number[3] + iz-(nk2-abs_layer_number[3]-1);

                /* for v_x: 01-integer; 10-half */
                Txz01_tmp =  (2.0*Txz01_z2[indx_abs] + dt*Dz_regu[iz]/Bz_regu[iz]*DzTxz01[indx]) \
                            /(2.0 + dt*(Az_regu[iz]+Dz_regu[iz]/Bz_regu[iz]));
                hVx01[indx] += (1.0/Bz_regu[iz]-1)*DzTxz01[indx] - Txz01_tmp/Bz_regu[iz];
                Txz01_z2[indx_abs] = 2.0*Txz01_tmp - Txz01_z2[indx_abs];

                Txz10_tmp =  (2.0*Txz10_z2[indx_abs] + dt*Dz_half[iz]/Bz_half[iz]*DzTxz10[indx]) \
                            /(2.0 + dt*(Az_half[iz]+Dz_half[iz]/Bz_half[iz]));
                hVx10[indx] += (1.0/Bz_half[iz]-1)*DzTxz01[indx] - Txz01_tmp/Bz_half[iz];
                Txz10_z2[indx_abs] = 2.0*Txz10_tmp - Txz10_z2[indx_abs];

                /* for v_z: 01-integer; 10-half */
                Tzz01_tmp =  (2.0*Tzz01_z2[indx_abs] + dt*Dz_regu[iz]/Bz_regu[iz]*DzTzz01[indx]) \
                            /(2.0 + dt*(Az_regu[iz]+Dz_regu[iz]/Bz_regu[iz]));
                hVz01[indx] += (1.0/Bz_regu[iz]-1)*DzTzz01[indx] - Tzz01_tmp/Bz_regu[iz];
                Tzz01_z2[indx_abs] = 2.0*Tzz01_tmp - Tzz01_z2[indx_abs];

                Tzz10_tmp =  (2.0*Tzz10_z2[indx_abs] + dt*Dz_half[iz]/Bz_half[iz]*DzTzz10[indx]) \
                            /(2.0 + dt*(Az_half[iz]+Dz_half[iz]/Bz_half[iz]));
                hVz10[indx] += (1.0/Bz_half[iz]-1)*DzTzz10[indx] - Tzz10_tmp/Bz_half[iz];
                Tzz10_z2[indx_abs] = 2.0*Tzz10_tmp - Tzz10_z2[indx_abs];

            }
        }
    }

    return  ierr;
}

int abs_npml_stress_tti(int ni1, int ni2, int nk1, int nk2, int nx, int nz,
                        int *abs_layer_number, float dt,
                        float *hTxx00, float *hTxx11,
                        float *hTxz00, float *hTxz11,
                        float *hTzz00, float *hTzz11,
                        float *DxVx00, float *DzVz00, float *DzVx00, float *DxVz00,
                        float *DxVx11, float *DzVz11, float *DzVx11, float *DxVz11,
                        float *c11_1 , float *c13_1 , float *c15_1 , float *c33_1 , float *c35_1, float *c55_1,
                        float *c11_2 , float *c13_2 , float *c15_2 , float *c33_2 , float *c35_2, float *c55_2,
                        float *Ax_regu, float *Bx_regu, float *Dx_regu,
                        float *Ax_half, float *Bx_half, float *Dx_half,
                        float *Az_regu, float *Bz_regu, float *Dz_regu,
                        float *Az_half, float *Bz_half, float *Dz_half,
                        float *Vx00_x1, float *Vx00_x2, float *Vx11_x1, float *Vx11_x2,
                        float *Vx00_z1, float *Vx00_z2, float *Vx11_z1, float *Vx11_z2,
                        float *Vz00_x1, float *Vz00_x2, float *Vz11_x1, float *Vz11_x2,
                        float *Vz00_z1, float *Vz00_z2, float *Vz11_z1, float *Vz11_z2)
{
    int ierr = 0, iz, ix, indx, indx_abs;
    float Vx00_x1_tmp, Vz00_x1_tmp, Vx00_x2_tmp, Vz00_x2_tmp;
    float Vx00_z1_tmp, Vz00_z1_tmp, Vx00_z2_tmp, Vz00_z2_tmp;
    float Vx11_x1_tmp, Vz11_x1_tmp, Vx11_x2_tmp, Vz11_x2_tmp;
    float Vx11_z1_tmp, Vz11_z1_tmp, Vx11_z2_tmp, Vz11_z2_tmp;

    /* left: x1 */
    if (abs_layer_number[0] > 0) {
        for (iz = nk1; iz < nk2; iz++) {
            for (ix = ni1; ix <= ni1+abs_layer_number[0]; ix++) {
                indx = iz*nx + ix;
                indx_abs = iz * abs_layer_number[0] + ix-ni1;

                /* for Txx, Tzz, Txz: 00-integer, 11-half */
                Vx00_x1_tmp = (2.0*Vx00_x1[indx_abs] + dt*Dx_regu[ix]/Bx_regu[ix]*DxVx00[indx]) \
                             /(2.0 + dt*(Ax_regu[ix]+Dx_regu[ix]/Bx_regu[ix]));
                Vz00_x1_tmp = (2.0*Vz00_x1[indx_abs] + dt*Dx_regu[ix]/Bx_regu[ix]*DxVz00[indx]) \
                             /(2.0 + dt*(Ax_regu[ix]+Dx_regu[ix]/Bx_regu[ix]));
                hTxx00[indx] += c11_1[indx]*( (1.0/Bx_regu[ix]-1)*DxVx00[indx] - Vx00_x1_tmp/Bx_regu[ix] ) \
                              + c15_1[indx]*( (1.0/Bx_regu[ix]-1)*DxVz00[indx] - Vz00_x1_tmp/Bx_regu[ix] );
                hTzz00[indx] += c13_1[indx]*( (1.0/Bx_regu[ix]-1)*DxVx00[indx] - Vx00_x1_tmp/Bx_regu[ix] ) \
                              + c35_1[indx]*( (1.0/Bx_regu[ix]-1)*DxVz00[indx] - Vz00_x1_tmp/Bx_regu[ix] );
                hTxz00[indx] += c15_1[indx]*( (1.0/Bx_regu[ix]-1)*DxVx00[indx] - Vx00_x1_tmp/Bx_regu[ix] ) \
                              + c55_2[indx]*( (1.0/Bx_regu[ix]-1)*DxVz00[indx] - Vz00_x1_tmp/Bx_regu[ix] );
                Vx00_x1[indx_abs] = 2.0*Vx00_x1_tmp - Vx00_x1[indx_abs];
                Vz00_x1[indx_abs] = 2.0*Vz00_x1_tmp - Vz00_x1[indx_abs];

                Vx11_x1_tmp = (2.0*Vx11_x1[indx_abs] + dt*Dx_half[ix]/Bx_half[ix]*DxVx11[indx]) \
                             /(2.0 + dt*(Ax_half[ix]+Dx_half[ix]/Bx_half[ix]));
                Vz11_x1_tmp = (2.0*Vz11_x1[indx_abs] + dt*Dx_half[ix]/Bx_half[ix]*DxVz11[indx]) \
                             /(2.0 + dt*(Ax_half[ix]+Dx_half[ix]/Bx_half[ix]));
                hTxx11[indx] += c11_2[indx]*( (1.0/Bx_half[ix]-1)*DxVx11[indx] - Vx11_x1_tmp/Bx_half[ix] ) \
                              + c15_2[indx]*( (1.0/Bx_half[ix]-1)*DxVz11[indx] - Vz11_x1_tmp/Bx_half[ix] );
                hTzz11[indx] += c13_2[indx]*( (1.0/Bx_half[ix]-1)*DxVx11[indx] - Vx11_x1_tmp/Bx_half[ix] ) \
                              + c35_2[indx]*( (1.0/Bx_half[ix]-1)*DxVz11[indx] - Vz11_x1_tmp/Bx_half[ix] );
                hTxz11[indx] += c15_2[indx]*( (1.0/Bx_half[ix]-1)*DxVx11[indx] - Vx11_x1_tmp/Bx_half[ix] ) \
                              + c55_1[indx]*( (1.0/Bx_half[ix]-1)*DxVz11[indx] - Vz11_x1_tmp/Bx_half[ix] );
                Vx11_x1[indx_abs] = 2.0*Vx11_x1_tmp - Vx11_x1[indx_abs];
                Vz11_x1[indx_abs] = 2.0*Vz11_x1_tmp - Vz11_x1[indx_abs];
            }
        }
    }


    /* right: x2 */
    if (abs_layer_number[1] > 0) {
        for (iz = nk1; iz < nk2; iz++) {
            for (ix = ni2-abs_layer_number[1]-1; ix < ni2; ix++) {
                indx = iz*nx + ix;
                indx_abs = iz * abs_layer_number[1] + ix-(ni2-abs_layer_number[1]-1);

                /* for Txx, Tzz, Txz: 00-integer, 11-half */
                Vx00_x2_tmp = (2.0*Vx00_x2[indx_abs] + dt*Dx_regu[ix]/Bx_regu[ix]*DxVx00[indx]) \
                             /(2.0 + dt*(Ax_regu[ix]+Dx_regu[ix]/Bx_regu[ix]));
                Vz00_x2_tmp = (2.0*Vz00_x2[indx_abs] + dt*Dx_regu[ix]/Bx_regu[ix]*DxVz00[indx]) \
                             /(2.0 + dt*(Ax_regu[ix]+Dx_regu[ix]/Bx_regu[ix]));
                hTxx00[indx] += c11_1[indx]*( (1.0/Bx_regu[ix]-1)*DxVx00[indx] - Vx00_x2_tmp/Bx_regu[ix] ) \
                              + c15_1[indx]*( (1.0/Bx_regu[ix]-1)*DxVz00[indx] - Vz00_x2_tmp/Bx_regu[ix] );
                hTzz00[indx] += c13_1[indx]*( (1.0/Bx_regu[ix]-1)*DxVx00[indx] - Vx00_x2_tmp/Bx_regu[ix] ) \
                              + c35_1[indx]*( (1.0/Bx_regu[ix]-1)*DxVz00[indx] - Vz00_x2_tmp/Bx_regu[ix] );
                hTxz00[indx] += c15_1[indx]*( (1.0/Bx_regu[ix]-1)*DxVx00[indx] - Vx00_x2_tmp/Bx_regu[ix] ) \
                              + c55_2[indx]*( (1.0/Bx_regu[ix]-1)*DxVz00[indx] - Vz00_x2_tmp/Bx_regu[ix] );
                Vx00_x2[indx_abs] = 2.0*Vx00_x2_tmp - Vx00_x2[indx_abs];
                Vz00_x2[indx_abs] = 2.0*Vz00_x2_tmp - Vz00_x2[indx_abs];

                Vx11_x2_tmp = (2.0*Vx11_x2[indx_abs] + dt*Dx_half[ix]/Bx_half[ix]*DxVx11[indx]) \
                             /(2.0 + dt*(Ax_half[ix]+Dx_half[ix]/Bx_half[ix]));
                Vz11_x2_tmp = (2.0*Vz11_x2[indx_abs] + dt*Dx_half[ix]/Bx_half[ix]*DxVz11[indx]) \
                             /(2.0 + dt*(Ax_half[ix]+Dx_half[ix]/Bx_half[ix]));
                hTxx11[indx] += c11_2[indx]*( (1.0/Bx_half[ix]-1)*DxVx11[indx] - Vx11_x2_tmp/Bx_half[ix] ) \
                              + c15_2[indx]*( (1.0/Bx_half[ix]-1)*DxVz11[indx] - Vz11_x2_tmp/Bx_half[ix] );
                hTzz11[indx] += c13_2[indx]*( (1.0/Bx_half[ix]-1)*DxVx11[indx] - Vx11_x2_tmp/Bx_half[ix] ) \
                              + c35_2[indx]*( (1.0/Bx_half[ix]-1)*DxVz11[indx] - Vz11_x2_tmp/Bx_half[ix] );
                hTxz11[indx] += c15_2[indx]*( (1.0/Bx_half[ix]-1)*DxVx11[indx] - Vx11_x2_tmp/Bx_half[ix] ) \
                              + c55_1[indx]*( (1.0/Bx_half[ix]-1)*DxVz11[indx] - Vz11_x2_tmp/Bx_half[ix] );
                Vx11_x2[indx_abs] = 2.0*Vx11_x2_tmp - Vx11_x2[indx_abs];
                Vz11_x2[indx_abs] = 2.0*Vz11_x2_tmp - Vz11_x2[indx_abs];
            }
        }
    }

    /* top: z1 */
    if (abs_layer_number[2] > 0) {
        for (iz = nk1; iz <= nk1+abs_layer_number[2]; iz++) {
            for (ix = ni1; ix < ni2; ix++) {
                indx = iz*nx + ix;
                indx_abs = ix * abs_layer_number[2] + iz-nk1;

                /* for Txx, Tzz, Txz: 00-integer, 11-half */
                Vx00_z1_tmp = (2.0*Vx00_z1[indx_abs] + dt*Dz_regu[iz]/Bz_regu[iz]*DzVx00[indx]) \
                             /(2.0 + dt*(Az_regu[iz]+Dz_regu[iz]/Bz_regu[iz]));
                Vz00_z1_tmp = (2.0*Vz00_z1[indx_abs] + dt*Dz_regu[iz]/Bz_regu[iz]*DzVz00[indx]) \
                              /(2.0 + dt*(Az_regu[iz]+Dz_regu[iz]/Bz_regu[iz]));
                hTxx00[indx] += c13_1[indx] * ( (1.0/Bz_regu[iz]-1)*DzVz00[indx] - Vz00_z1_tmp/Bz_regu[iz] ) \
                              + c15_1[indx] * ( (1.0/Bz_regu[iz]-1)*DzVx00[indx] - Vx00_z1_tmp/Bz_regu[iz] );
                hTzz00[indx] += c33_1[indx] * ( (1.0/Bz_regu[iz]-1)*DzVz00[indx] - Vz00_z1_tmp/Bz_regu[iz] ) \
                              + c35_1[indx] * ( (1.0/Bz_regu[iz]-1)*DzVx00[indx] - Vx00_z1_tmp/Bz_regu[iz] );
                hTxz00[indx] += c35_1[indx] * ( (1.0/Bz_regu[iz]-1)*DzVz00[indx] - Vz00_z1_tmp/Bz_regu[iz] ) \
                              + c55_2[indx] * ( (1.0/Bz_regu[iz]-1)*DzVx00[indx] - Vx00_z1_tmp/Bz_regu[iz] );
                Vx00_z1[indx_abs] = 2.0*Vx00_z1_tmp - Vx00_z1[indx_abs];
                Vz00_z1[indx_abs] = 2.0*Vz00_z1_tmp - Vz00_z1[indx_abs];

                Vx11_z1_tmp = (2.0*Vx11_z1[indx_abs] + dt*Dz_half[iz]/Bz_half[iz]*DzVx11[indx]) \
                             /(2.0 + dt*(Az_half[iz]+Dz_half[iz]/Bz_half[iz]));
                Vz11_z1_tmp = (2.0*Vz11_z1[indx_abs] + dt*Dz_half[iz]/Bz_half[iz]*DzVz11[indx]) \
                             /(2.0 + dt*(Az_half[iz]+Dz_half[iz]/Bz_half[iz]));
                hTxx11[indx] += c13_1[indx] * ( (1.0/Bz_half[iz]-1)*DzVz11[indx] - Vz11_z1_tmp/Bz_half[iz] ) \
                              + c15_1[indx] * ( (1.0/Bz_half[iz]-1)*DzVx11[indx] - Vx11_z1_tmp/Bz_half[iz] );
                hTzz11[indx] += c33_1[indx] * ( (1.0/Bz_half[iz]-1)*DzVz11[indx] - Vz11_z1_tmp/Bz_half[iz] ) \
                              + c35_1[indx] * ( (1.0/Bz_half[iz]-1)*DzVx11[indx] - Vx11_z1_tmp/Bz_half[iz] );
                hTxz11[indx] += c35_2[indx] * ( (1.0/Bz_half[iz]-1)*DzVz11[indx] - Vz11_z1_tmp/Bz_half[iz] )
                              + c55_1[indx] * ( (1.0/Bz_half[iz]-1)*DzVx11[indx] - Vx11_z1_tmp/Bz_half[iz] );
                Vx11_z1[indx_abs] = 2.0*Vx11_z1_tmp - Vx11_z1[indx_abs];
                Vz11_z1[indx_abs] = 2.0*Vz11_z1_tmp - Vz11_z1[indx_abs];
            }
        }
    }

    /* bottom: z2 */
    if (abs_layer_number[3] > 0) {
        for (iz = nk2-abs_layer_number[3]-1; iz < nk2; iz++) {
            for (ix = ni1; ix < ni2; ix++) {
                indx = iz*nx + ix;
                indx_abs = ix * abs_layer_number[3] + iz-(nk2-abs_layer_number[3]-1);

                /* for Txx, Tzz, Txz: 00-integer, 11-half */
                Vx00_z2_tmp = (2.0*Vx00_z2[indx_abs] + dt*Dz_regu[iz]/Bz_regu[iz]*DzVx00[indx]) \
                             /(2.0 + dt*(Az_regu[iz]+Dz_regu[iz]/Bz_regu[iz]));
                Vz00_z2_tmp = (2.0*Vz00_z2[indx_abs] + dt*Dz_regu[iz]/Bz_regu[iz]*DzVz00[indx]) \
                             /(2.0 + dt*(Az_regu[iz]+Dz_regu[iz]/Bz_regu[iz]));
                hTxx00[indx] += c13_1[indx] * ( (1.0/Bz_regu[iz]-1)*DzVz00[indx] - Vz00_z2_tmp/Bz_regu[iz] ) \
                              + c15_1[indx] * ( (1.0/Bz_regu[iz]-1)*DzVx00[indx] - Vx00_z2_tmp/Bz_regu[iz] );
                hTzz00[indx] += c33_1[indx] * ( (1.0/Bz_regu[iz]-1)*DzVz00[indx] - Vz00_z2_tmp/Bz_regu[iz] ) \
                              + c35_1[indx] * ( (1.0/Bz_regu[iz]-1)*DzVx00[indx] - Vx00_z2_tmp/Bz_regu[iz] );
                hTxz00[indx] += c35_1[indx] * ( (1.0/Bz_regu[iz]-1)*DzVz00[indx] - Vz00_z2_tmp/Bz_regu[iz] ) \
                              + c55_2[indx] * ( (1.0/Bz_regu[iz]-1)*DzVx00[indx] - Vx00_z2_tmp/Bz_regu[iz] );
                Vx00_z2[indx_abs] = 2.0*Vx00_z2_tmp - Vx00_z2[indx_abs];
                Vz00_z2[indx_abs] = 2.0*Vz00_z2_tmp - Vz00_z2[indx_abs];

                Vx11_z2_tmp = (2.0*Vx11_z2[indx_abs] + dt*Dz_half[iz]/Bz_half[iz]*DzVx11[indx]) \
                             /(2.0 + dt*(Az_half[iz]+Dz_half[iz]/Bz_half[iz]));
                Vz11_z2_tmp = (2.0*Vz11_z2[indx_abs] + dt*Dz_half[iz]/Bz_half[iz]*DzVz11[indx]) \
                             /(2.0 + dt*(Az_half[iz]+Dz_half[iz]/Bz_half[iz]));
                hTxx11[indx] += c13_1[indx] * ( (1.0/Bz_half[iz]-1)*DzVz11[indx] - Vz11_z2_tmp/Bz_half[iz] ) \
                              + c15_1[indx] * ( (1.0/Bz_half[iz]-1)*DzVx11[indx] - Vx11_z2_tmp/Bz_half[iz] );
                hTzz11[indx] += c33_1[indx] * ( (1.0/Bz_half[iz]-1)*DzVz11[indx] - Vz11_z2_tmp/Bz_half[iz] ) \
                              + c35_1[indx] * ( (1.0/Bz_half[iz]-1)*DzVx11[indx] - Vx11_z2_tmp/Bz_half[iz] );
                hTxz11[indx] += c35_2[indx] * ( (1.0/Bz_half[iz]-1)*DzVz11[indx] - Vz11_z2_tmp/Bz_half[iz] ) \
                              + c55_1[indx] * ( (1.0/Bz_half[iz]-1)*DzVx11[indx] - Vx11_z2_tmp/Bz_half[iz] );
                Vx11_z2[indx_abs] = 2.0*Vx11_z2_tmp - Vx11_z2[indx_abs];
                Vz11_z2[indx_abs] = 2.0*Vz11_z2_tmp - Vz11_z2[indx_abs];
            }
        }
    }

    return ierr;
}

int abs_npml_init(int *abs_layer_number,
                  float fc, float *c11, float *c55, float *rho_recip,
                  int half_fd_stencil,
                  float *Ax_regu, float *Bx_regu, float *Dx_regu,
                  float *Ax_half, float *Bx_half, float *Dx_half,
                  float *Az_regu, float *Bz_regu, float *Dz_regu,
                  float *Az_half, float *Bz_half, float *Dz_half,
                  float *xvec1, float *xvec2, float *zvec1, float *zvec2,
                  int ni1, int ni2, int nk1, int nk2, int nx, float dx, float dz)
{
    int ierr = 0, i;
    float L0, x0, d0, alpha0, beta0[4], Rpp, Lx, Lz, vp_pml[4], vs_pml[4], PPW_FD;

    /* alpha max */
    alpha0 = cal_pml_alpha0(fc);

    /* for beta0: PPW_S/PPW_FD */
    ierr = cal_vel_pml(c55, rho_recip, vs_pml, abs_layer_number,
                       ni1, ni2, nk1, nk2, nx);

    if ( half_fd_stencil <= 2 )
        PPW_FD = 7.0;
    else if (half_fd_stencil <= 4)
        PPW_FD = 5.5;
    else if (half_fd_stencil <= 7)
        PPW_FD = 4.0;
    else
        PPW_FD = 3.0;

    printf("vs: %f\n", vs_pml[0]);
    beta0[0] = 2 * vs_pml[0]/dx/fc/PPW_FD;
    beta0[1] = 2 * vs_pml[1]/dx/fc/PPW_FD;
    beta0[2] = 2 * vs_pml[2]/dz/fc/PPW_FD;
    beta0[3] = 2 * vs_pml[3]/dz/fc/PPW_FD;

    /* vp_pml, used for calculate d0 */
    ierr = cal_vel_pml(c11, rho_recip, vp_pml, abs_layer_number,
                       ni1, ni2, nk1, nk2, nx);

    fprintf(stdout, " -------- For ADE CFS-PML -------- \n");

    // left
    if (abs_layer_number[0] > 0) {
        L0  = xvec2[ abs_layer_number[0]+ni1 ] - xvec1[ni1];
        Rpp = cal_pml_R(abs_layer_number[0]);
        d0  = cal_pml_d0(L0, vp_pml[0], Rpp);

        fprintf(stdout, "  left:   alpha0=%f, beta0=%f, d0=%f\n", \
                alpha0, beta0[0], d0);

        for (i = ni1; i <= abs_layer_number[0]+ni1; i++) {
            /* integer grid points */
            Lx = xvec2[abs_layer_number[0]+ni1] - xvec1[i];
            Dx_regu[i] = cal_pml_d(d0, Lx, L0);
            Ax_regu[i] = cal_pml_alpha(alpha0, Lx, L0);
            Bx_regu[i] = cal_pml_beta(beta0[0], Lx, L0);
            //printf("Ax_regu[%d]: %f, Bx_regu: %f, Dx_regu: %f\n", i, Ax_regu[i], Bx_regu[i], Dx_regu[i]);

            /* half grid points */
            Lx = xvec2[abs_layer_number[0]+ni1] - xvec2[i];
            Dx_half[i] = cal_pml_d(d0, Lx, L0);
            Ax_half[i] = cal_pml_alpha(alpha0, Lx, L0);
            Bx_half[i] = cal_pml_beta(beta0[0], Lx, L0);
           //printf("Ax_half[%d]: %f, Bx_half: %f, Dx_half: %f\n", i, Ax_half[i], Bx_half[i], Dx_half[i]);

        }

    }

    // right
    if (abs_layer_number[1] > 0) {
        L0  = xvec2[ni2-1] - xvec1[ni2-abs_layer_number[1]-1];
        Rpp = cal_pml_R(abs_layer_number[1]);
        d0  = cal_pml_d0(L0, vp_pml[1], Rpp);

        fprintf(stdout, "  right:  alpha0=%f, beta0=%f, d0=%f\n",
                alpha0, beta0[1], d0);

        for (i = ni2-abs_layer_number[1]-1; i < ni2; i++) {
            /* integer grid points */
            Lx = xvec1[i] - xvec1[ni2-abs_layer_number[1]-1];
            Dx_regu[i] = cal_pml_d(d0, Lx, L0);
            Ax_regu[i] = cal_pml_alpha(alpha0, Lx, L0);
            Bx_regu[i] = cal_pml_beta(beta0[1], Lx, L0);

            /* half grid points */
            Lx = xvec2[i] - xvec1[ni2-abs_layer_number[1]-1];
            Dx_half[i] = cal_pml_d(d0, Lx, L0);
            Ax_half[i] = cal_pml_alpha(alpha0, Lx, L0);
            Bx_half[i] = cal_pml_beta(beta0[1], Lx, L0);
        }
    }

    // top
    if (abs_layer_number[2] > 0) {
        L0  = zvec2[ abs_layer_number[2]+nk1 ] - zvec1[nk1];
        Rpp = cal_pml_R(abs_layer_number[2]);
        d0  = cal_pml_d0(L0, vp_pml[2], Rpp);

        fprintf(stdout, "  top:    alpha0=%f, beta0=%f, d0=%f\n",
                alpha0, beta0[2], d0);

        for (i = nk1; i <= abs_layer_number[2]+nk1; i++) {
            /* integer grid points */
            Lz = zvec2[abs_layer_number[2]+nk1] - zvec1[i];
            Dz_regu[i] = cal_pml_d(d0, Lz, L0);
            Az_regu[i] = cal_pml_alpha(alpha0, Lz, L0);
            Bz_regu[i] = cal_pml_beta(beta0[2], Lz, L0);
//            printf("Az_regu[%d]: %f, Bz_regu: %f, Dz_regu: %f\n", i, Az_regu[i], Bz_regu[i], Dz_regu[i]);

            /* half grid points */
            Lz = zvec2[abs_layer_number[2]+nk1] - zvec2[i];
            Dz_half[i] = cal_pml_d(d0, Lz, L0);
            Az_half[i] = cal_pml_alpha(alpha0, Lz, L0);
            Bz_half[i] = cal_pml_beta(beta0[2], Lz, L0);
//            printf("Az_half[%d]: %f, Bz_half: %f, Dz_half: %f\n", i, Az_half[i], Bz_half[i], Dz_half[i]);

        }

    }

    // bottom
    if (abs_layer_number[3] > 0) {
        L0  = zvec2[nk2-1] - zvec1[nk2-abs_layer_number[3]-1];
        Rpp = cal_pml_R(abs_layer_number[3]);
        d0  = cal_pml_d0(L0, vp_pml[3], Rpp);

        fprintf(stdout, "  bottom: alpha0=%f, beta0=%f, d0=%f\n",
                alpha0, beta0[3], d0);

        for (i = nk2-abs_layer_number[3]-1; i < nk2; i++) {
            /* integer grid points */
            Lz = zvec1[i] - zvec1[nk2-abs_layer_number[3]-1];
            Dz_regu[i] = cal_pml_d(d0, Lz, L0);
            Az_regu[i] = cal_pml_alpha(alpha0, Lz, L0);
            Bz_regu[i] = cal_pml_beta(beta0[3], Lz, L0);
//            printf("Az_regu[%d]: %f, Bz_regu: %f, Dz_regu: %f\n", i, Az_regu[i], Bz_regu[i], Dz_regu[i]);


            /* half grid points */
            Lz = zvec2[i] - zvec1[nk2-abs_layer_number[3]-1];
            Dz_half[i] = cal_pml_d(d0, Lz, L0);
            Az_half[i] = cal_pml_alpha(alpha0, Lz, L0);
            Bz_half[i] = cal_pml_beta(beta0[3], Lz, L0);
//            printf("Az_half[%d]: %f, Bz_half: %f, Dz_half: %f\n", i, Az_half[i], Bz_half[i], Dz_half[i]);

        }
    }

    return ierr;
}



float cal_pml_R(int N)
{
    float R;

    R = pow(10, -(log10(N)-1)/log10(2)-3 );

    return R;
}

float cal_pml_d0(float L, float vp, float R)
{
    float d0;

    d0 = -3.0*vp/2.0/L*log(R);

    return d0;
}

float cal_pml_alpha0(float fc)
{
    float alpha0;

    alpha0 = PI*fc;

    return alpha0;
}

float cal_pml_beta0(float PPW_s, float PPW_FD)
{
    float beta0;

    beta0 = 2 * PPW_s/PPW_FD;

    return beta0;
}

float cal_pml_d(float d0, float x, float L)
{
    float d = 0.0;

    if (x<0) {
        d = 0.0;
    } else {
        d = d0*(x/L)*(x/L);
    }

    return d;
}

float cal_pml_alpha(float alpha0, float x, float L)
{
    float alpha = 0.0;

    if(x<0) {
        alpha = 0.0;
    } else {
        alpha = alpha0*(1-x/L);
    }

    return alpha;
}

float cal_pml_beta(float beta0, float x, float L)
{
    float  beta = 0;

    if(x<0) {
        beta = 1.0;
    } else {
        beta = 1+(beta0-1)*(x/L)*(x/L);
    }

    return beta;
}

int cal_vel_pml(float *c11, float *rho_recip, float *vp_pml,
               int *abs_layer_number,
               int ni1, int ni2, int nk1, int nk2, int nx)
{
    int ix, iz, indx;

    /* left */
    indx = nk1*nx + ni1;
    vp_pml[0] = sqrt(c11[indx]*rho_recip[indx]);
    for (ix = ni1; ix <= ni1 + abs_layer_number[0]; ix ++) {
        for (iz = nk1; iz < nk2; iz++) {
            indx = ix + iz*nx;
            if ( vp_pml[0] < sqrt(c11[indx]*rho_recip[indx]) )
                vp_pml[0] = sqrt(c11[indx]*rho_recip[indx]);
        }
    }

    /* right */
    indx = nk1*nx + (ni2-abs_layer_number[1]-1);
    vp_pml[1] = sqrt(c11[indx]*rho_recip[indx]);
    for (ix = ni2-abs_layer_number[1]-1; ix < ni2; ix ++) {
        for (iz = nk1; iz < nk2; iz++) {
            indx = ix + iz*nx;
            if ( vp_pml[1] < sqrt(c11[indx]*rho_recip[indx]) )
                vp_pml[1] = sqrt(c11[indx]*rho_recip[indx]);
        }
    }

    /* top */
    indx = nk1*nx + ni1;
    vp_pml[2] = sqrt(c11[indx]*rho_recip[indx]);
    for (ix = ni1; ix < ni2 ; ix ++) {
        for (iz = nk1; iz <= nk1 + abs_layer_number[2]; iz++) {
            indx = ix + iz*nx;
            if ( vp_pml[2] < sqrt(c11[indx]*rho_recip[indx]) )
                vp_pml[2] = sqrt(c11[indx]*rho_recip[indx]);
        }
    }

    /* bottom */
    indx = (nk2-abs_layer_number[3]-1)*nx + ni1;
    vp_pml[3] = sqrt(c11[indx]*rho_recip[indx]);
    for (ix = ni1; ix < ni2; ix ++) {
        for (iz = nk2-abs_layer_number[3]-1; iz < nk2; iz++) {
            indx = ix + iz*nx;
            if ( vp_pml[3] < sqrt(c11[indx]*rho_recip[indx]) )
                vp_pml[3] = sqrt(c11[indx]*rho_recip[indx]);
        }
    }

    return 0;
}


int abs_npml_velocity_staggered(int ni1, int ni2, int nk1, int nk2, int nx, int nz,
                                int *abs_layer_number, float dt,
                                float *hVx, float *hVz,
                                float *Ax_regu, float *Bx_regu, float *Dx_regu,
                                float *Ax_half, float *Bx_half, float *Dx_half,
                                float *Az_regu, float *Bz_regu, float *Dz_regu,
                                float *Az_half, float *Bz_half, float *Dz_half,
                                float *Txx_x1, float *Txx_x2,
                                float *Tzz_z1, float *Tzz_z2,
                                float *Txz_x1, float *Txz_x2,
                                float *Txz_z1, float *Txz_z2,
                                float *DxTxx01, float *DzTxz01,
                                float *DxTxz10, float *DzTzz10)
{
    int ierr = 0, iz, ix, indx, indx_abs;
    float Txx01_tmp, Tzz10_tmp, Txz01_tmp, Txz10_tmp;

    /* left: x1 */
    if (abs_layer_number[0] > 0) {
        for (iz = nk1; iz < nk2; iz++) {
            for (ix = ni1; ix <= ni1+abs_layer_number[0]; ix++) {

                indx = iz*nx + ix;
                indx_abs = iz * abs_layer_number[0] + ix-ni1;

                /* for v_x: 01-half; */
                Txx01_tmp =  (2.0*Txx_x1[indx_abs] + dt*Dx_half[ix]/Bx_half[ix]*DxTxx01[indx])
                            /(2.0 + dt*(Ax_half[ix]+Dx_half[ix]/Bx_half[ix]));
                hVx[indx] += (1.0/Bx_half[ix]-1)*DxTxx01[indx] - Txx01_tmp/Bx_half[ix];
                Txx_x1[indx_abs] = 2.0*Txx01_tmp - Txx_x1[indx_abs];

                /* for v_z: 10-integer */
                Txz10_tmp =  (2.0*Txz_x1[indx_abs] + dt*Dx_half[ix]/Bx_half[ix]*DxTxz10[indx]) \
                            /(2.0 + dt*(Ax_half[ix]+Dx_half[ix]/Bx_half[ix]));
                hVz[indx] += (1.0/Bx_half[ix]-1)*DxTxz10[indx] - Txz01_tmp/Bx_half[ix];
                Txz_x1[indx_abs] = 2.0*Txz10_tmp - Txz_x1[indx_abs];

            }
        }
    }

    /* right: x2 */
    if (abs_layer_number[1] > 0) {
        for (iz = nk1; iz < nk2; iz++) {
            for (ix = ni2-abs_layer_number[1]-1; ix < ni2; ix++) {

                indx = iz*nx + ix;
                indx_abs = iz * abs_layer_number[1] + ix-(ni2-abs_layer_number[1]-1);

                /* for v_x: 01-half;  */
                Txx01_tmp =  (2.0*Txx_x2[indx_abs] + dt*Dx_half[ix]/Bx_half[ix]*DxTxx01[indx]) \
                            /(2.0 + dt*(Ax_half[ix]+Dx_half[ix]/Bx_half[ix]));
                hVx[indx] += (1.0/Bx_half[ix]-1)*DxTxx01[indx] - Txx01_tmp/Bx_half[ix];
                Txx_x2[indx_abs] = 2.0*Txx01_tmp - Txx_x2[indx_abs];

                /* for v_z: 10-integer */
                Txz10_tmp =  (2.0*Txz_x2[indx_abs] + dt*Dx_regu[ix]/Bx_regu[ix]*DxTxz10[indx]) \
                            /(2.0 + dt*(Ax_regu[ix]+Dx_regu[ix]/Bx_regu[ix]));
                hVz[indx] += (1.0/Bx_regu[ix]-1)*DxTxz10[indx] - Txz10_tmp/Bx_regu[ix];
                Txz_x2[indx_abs] = 2.0*Txz10_tmp - Txz_x2[indx_abs];

            }
        }
    }

    /* top: z1 */
    if (abs_layer_number[2] > 0) {
        for (iz = nk1; iz <= nk1+abs_layer_number[2]; iz++) {
            for (ix = ni1; ix < ni2; ix++) {
                indx = iz*nx + ix;
                indx_abs = ix * abs_layer_number[2] + iz-nk1;

                /* for v_x: 01-integer; */
                Txz01_tmp =  (2.0*Txz_z1[indx_abs] + dt*Dz_regu[iz]/Bz_regu[iz]*DzTxz01[indx]) \
                            /(2.0 + dt*(Az_regu[iz]+Dz_regu[iz]/Bz_regu[iz]));
                hVx[indx] += (1.0/Bz_regu[iz]-1)*DzTxz01[indx] - Txz01_tmp/Bz_regu[iz];
                Txz_z1[indx_abs] = 2.0*Txz01_tmp - Txz_z1[indx_abs];

                /* for v_z:  10-half */
                Tzz10_tmp =  (2.0*Tzz_z1[indx_abs] + dt*Dz_half[iz]/Bz_half[iz]*DzTzz10[indx]) \
                            /(2.0 + dt*(Az_half[iz]+Dz_half[iz]/Bz_half[iz]));
                hVz[indx] += (1.0/Bz_half[iz]-1)*DzTzz10[indx] - Tzz10_tmp/Bz_half[iz];
                Tzz_z1[indx_abs] = 2.0*Tzz10_tmp - Tzz_z1[indx_abs];

            }
        }
    }

    /* bottom: z2 */
    if (abs_layer_number[3] > 0) {
        for (iz = nk2-abs_layer_number[3]-1; iz < nk2; iz++) {
            for (ix = ni1; ix < ni2; ix++) {
                indx = iz*nx + ix;
                indx_abs = ix * abs_layer_number[3] + iz-(nk2-abs_layer_number[3]-1);

                /* for v_x: 01-integer; */
                Txz01_tmp =  (2.0*Txz_z2[indx_abs] + dt*Dz_regu[iz]/Bz_regu[iz]*DzTxz01[indx]) \
                            /(2.0 + dt*(Az_regu[iz]+Dz_regu[iz]/Bz_regu[iz]));
                hVx[indx] += (1.0/Bz_regu[iz]-1)*DzTxz01[indx] - Txz01_tmp/Bz_regu[iz];
                Txz_z2[indx_abs] = 2.0*Txz01_tmp - Txz_z2[indx_abs];

                /* for v_z:  10-half */
                Tzz10_tmp =  (2.0*Tzz_z2[indx_abs] + dt*Dz_half[iz]/Bz_half[iz]*DzTzz10[indx]) \
                            /(2.0 + dt*(Az_half[iz]+Dz_half[iz]/Bz_half[iz]));
                hVz[indx] += (1.0/Bz_half[iz]-1)*DzTzz10[indx] - Tzz10_tmp/Bz_half[iz];
                Tzz_z2[indx_abs] = 2.0*Tzz10_tmp - Tzz_z2[indx_abs];

            }
        }
    }

    return  ierr;
}

int abs_npml_stress_staggered(int ni1, int ni2, int nk1, int nk2, int nx, int nz,
                              int *abs_layer_number, float dt,
                              float *hTxx, float *hTxz, float *hTzz,
                              float *DxVx00, float *DzVz00, float *DzVx11, float *DxVz11,
                              float *c11, float *c13, float *c33, float *c55,
                              float *Ax_regu, float *Bx_regu, float *Dx_regu,
                              float *Ax_half, float *Bx_half, float *Dx_half,
                              float *Az_regu, float *Bz_regu, float *Dz_regu,
                              float *Az_half, float *Bz_half, float *Dz_half,
                              float *Vx_x1, float *Vx_x2, float *Vx_z1, float *Vx_z2,
                              float *Vz_x1, float *Vz_x2, float *Vz_z1, float *Vz_z2)
{
    int ierr = 0, iz, ix, indx, indx_abs;
    float Vx_x1_tmp, Vz_x1_tmp, Vx_x2_tmp, Vz_x2_tmp;
    float Vx_z1_tmp, Vz_z1_tmp, Vx_z2_tmp, Vz_z2_tmp;

    /* left: x1 */
    if (abs_layer_number[0] > 0) {
        for (iz = nk1; iz < nk2; iz++) {
            for (ix = ni1; ix <= ni1+abs_layer_number[0]; ix++) {
                indx = iz*nx + ix;
                indx_abs = iz * abs_layer_number[0] + ix-ni1;

                /* for Txx, Tzz, Txz: 00-integer, 11-half */
                Vx_x1_tmp = (2.0*Vx_x1[indx_abs] + dt*Dx_regu[ix]/Bx_regu[ix]*DxVx00[indx]) \
                           /(2.0 + dt*(Ax_regu[ix]+Dx_regu[ix]/Bx_regu[ix]));
                hTxx[indx] += c11[indx]*( (1.0/Bx_regu[ix]-1)*DxVx00[indx] - Vx_x1_tmp/Bx_regu[ix] );
                hTzz[indx] += c13[indx]*( (1.0/Bx_regu[ix]-1)*DxVx00[indx] - Vx_x1_tmp/Bx_regu[ix] );
                Vx_x1[indx_abs] = 2.0*Vx_x1_tmp - Vx_x1[indx_abs];

                /* for Txz, Tzz: 11-half */
                Vz_x1_tmp = (2.0*Vz_x1[indx_abs] + dt*Dx_half[ix]/Bx_half[ix]*DxVz11[indx]) \
                           /(2.0 + dt*(Ax_half[ix]+Dx_half[ix]/Bx_half[ix]));
                hTxz[indx] += c55[indx]*( (1.0/Bx_half[ix]-1)*DxVz11[indx] - Vz_x1_tmp/Bx_half[ix] );
                Vz_x1[indx_abs] = 2.0*Vz_x1_tmp - Vz_x1[indx_abs];
            }
        }
    }


    /* right: x2 */
    if (abs_layer_number[1] > 0) {
        for (iz = nk1; iz < nk2; iz++) {
            for (ix = ni2-abs_layer_number[1]-1; ix < ni2; ix++) {
                indx = iz*nx + ix;
                indx_abs = iz * abs_layer_number[1] + ix-(ni2-abs_layer_number[1]-1);

                /* for Txx, Tzz: 00-integer */
                Vx_x2_tmp = (2.0*Vx_x2[indx_abs] + dt*Dx_regu[ix]/Bx_regu[ix]*DxVx00[indx]) \
                           /(2.0 + dt*(Ax_regu[ix]+Dx_regu[ix]/Bx_regu[ix]));
                hTxx[indx] += c11[indx]*( (1.0/Bx_regu[ix]-1)*DxVx00[indx] - Vx_x2_tmp/Bx_regu[ix] );
                hTzz[indx] += c13[indx]*( (1.0/Bx_regu[ix]-1)*DxVx00[indx] - Vx_x2_tmp/Bx_regu[ix] );
                Vx_x2[indx_abs] = 2.0*Vx_x2_tmp - Vx_x2[indx_abs];

                /* for Txz, Tzz: 11-half */
                Vz_x2_tmp = (2.0*Vz_x2[indx_abs] + dt*Dx_half[ix]/Bx_half[ix]*DxVz11[indx]) \
                           /(2.0 + dt*(Ax_half[ix]+Dx_half[ix]/Bx_half[ix]));
                hTxz[indx] += c55[indx]*( (1.0/Bx_half[ix]-1)*DxVz11[indx] - Vz_x2_tmp/Bx_half[ix] );
                Vz_x2[indx_abs] = 2.0*Vz_x2_tmp - Vz_x2[indx_abs];
            }
        }
    }

    /* top: z1 */
    if (abs_layer_number[2] > 0) {
        for (iz = nk1; iz <= nk1+abs_layer_number[2]; iz++) {
            for (ix = ni1; ix < ni2; ix++) {
                indx = iz*nx + ix;
                indx_abs = ix * abs_layer_number[2] + iz-nk1;

                /* for Txx, Tzz, 00-integer */
                Vz_z1_tmp =  (2.0*Vz_z1[indx_abs] + dt*Dz_regu[iz]/Bz_regu[iz]*DzVz00[indx]) \
                            /(2.0 + dt*(Az_regu[iz]+Dz_regu[iz]/Bz_regu[iz]));
                hTxx[indx] += c13[indx] * ((1.0/Bz_regu[iz]-1)*DzVz00[indx] - Vz_z1_tmp/Bz_regu[iz]);
                hTzz[indx] += c33[indx] * ((1.0/Bz_regu[iz]-1)*DzVz00[indx] - Vz_z1_tmp/Bz_regu[iz]);
                Vz_z1[indx_abs] = 2.0*Vz_z1_tmp - Vz_z1[indx_abs];

                /* for Txz: 11-half */
                Vx_z1_tmp =  (2.0*Vx_z1[indx_abs] + dt*Dz_half[iz]/Bz_half[iz]*DzVx11[indx]) \
                            /(2.0 + dt*(Az_half[iz]+Dz_half[iz]/Bz_half[iz]));
                hTxz[indx] += c55[indx] * ((1.0/Bz_half[iz]-1.0)*DzVx11[indx] - Vx_z1_tmp/Bz_half[iz]);
                Vx_z1[indx_abs] = 2.0*Vx_z1_tmp - Vx_z1[indx_abs];
            }
        }
    }

    /* bottom: z2 */
    if (abs_layer_number[3] > 0) {
        for (iz = nk2-abs_layer_number[3]-1; iz < nk2; iz++) {
            for (ix = ni1; ix < ni2; ix++) {
                indx = iz*nx + ix;
                indx_abs = ix * abs_layer_number[3] + iz-(nk2-abs_layer_number[3]-1);

                /* for Txx, Tzz: 00-integer */
                Vz_z2_tmp =  (2.0 * Vz_z2[indx_abs] + dt*Dz_regu[iz]/Bz_regu[iz]*DzVz00[indx])
                            /(2.0 + dt*(Az_regu[iz]+Dz_regu[iz]/Bz_regu[iz]));
                hTxx[indx] += c13[indx] * ((1.0/Bz_regu[iz]-1.0)*DzVz00[indx] - Vz_z2_tmp/Bz_regu[iz]);
                hTzz[indx] += c33[indx] * ((1.0/Bz_regu[iz]-1.0)*DzVz00[indx] - Vz_z2_tmp/Bz_regu[iz]);
                Vz_z2[indx_abs] = 2.0*Vz_z2_tmp - Vz_z2[indx_abs];

                /* for Txz, Tzz: 11-half */
                Vx_z2_tmp = (2.0*Vx_z2[indx_abs] + dt*Dz_half[iz]/Bz_half[iz]*DzVx11[indx]) \
                           /(2.0 + dt*(Az_half[iz]+Dz_half[iz]/Bz_half[iz]));
                hTxz[indx] += c55[indx] * ( (1.0/Bz_half[iz]-1.0)*DzVx11[indx] - Vx_z2_tmp/Bz_half[iz] );
                Vx_z2[indx_abs] = 2.0*Vx_z2_tmp - Vx_z2[indx_abs];
            }
        }
    }

    return ierr;
}
