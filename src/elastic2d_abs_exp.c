/***************************************************************************
 *
 * This function is used for the absorbing outgoing waves based on
 *   Cerjan and Kosloff's nonreflecting boundary condition
 *   (Cerjan C., et al. (1985), Geophysics, 50(4):705-708)
 *
 *
 * Authors: Luqian Jiang   Email: jlq0608@gmail.com
 *          Wei Zhang      Email: zhangwei@sustc.edu.cn
 *
 * Copyright (c) 2018 - 2020 zwlab
 *
 * History: 08/2019: Original version created by Luqian Jiang
 *
 ***************************************************************************/

#include "elastic2d_abs_exp.h"

int abs_exp_velocity(int ni1, int ni2, int nk1, int nk2, int nx, int nz, int half_fd_stencil,
                     int *boundary_layer_number,
                     float *Vx_1, float *Vx_2, float *Vz_1, float *Vz_2)
{
    int ierr = 0;

    if( boundary_layer_number[0] > 0 ) {
        ierr = abs_exp_x1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[0], Vx_1);
        ierr = abs_exp_x1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[0], Vx_2);
        ierr = abs_exp_x1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[0], Vz_1);
        ierr = abs_exp_x1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[0], Vz_2);
    }

    if( boundary_layer_number[1] > 0) {
        ierr = abs_exp_x2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[1], Vx_1);
        ierr = abs_exp_x2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[1], Vx_2);
        ierr = abs_exp_x2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[1], Vz_1);
        ierr = abs_exp_x2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[1], Vz_2);
    }

    if( boundary_layer_number[2] > 0 ) {
        ierr = abs_exp_z1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[2], Vx_1);
        ierr = abs_exp_z1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[2], Vx_2);
        ierr = abs_exp_z1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[2], Vz_1);
        ierr = abs_exp_z1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[2], Vz_2);
    }

    if( boundary_layer_number[30] > 0 ) {
        ierr = abs_exp_z2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[3], Vx_1);
        ierr = abs_exp_z2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[3], Vx_2);
        ierr = abs_exp_z2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[3], Vz_1);
        ierr = abs_exp_z2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[3], Vz_2);
    }
    return 0;
}

int abs_exp_stresses(int ni1, int ni2, int nk1, int nk2, int nx, int nz, int half_fd_stencil,
                     int *boundary_layer_number,
                     float *Txx_1, float *Txx_2, float *Tzz_1, float *Tzz_2, float *Txz_1, float *Txz_2)
{
    int ierr = 0;

    if(boundary_layer_number[0] > 0) {
        ierr = abs_exp_x1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[0], Txx_1);
        ierr = abs_exp_x1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[0], Txx_2);
        ierr = abs_exp_x1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[0], Tzz_1);
        ierr = abs_exp_x1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[0], Tzz_2);
        ierr = abs_exp_x1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[0], Txz_1);
        ierr = abs_exp_x1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[0], Txz_2);
    }

    if(boundary_layer_number[1] > 0) {
        ierr = abs_exp_x2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[1], Txx_1);
        ierr = abs_exp_x2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[1], Txx_2);
        ierr = abs_exp_x2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[1], Tzz_1);
        ierr = abs_exp_x2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[1], Tzz_2);
        ierr = abs_exp_x2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[1], Txz_1);
        ierr = abs_exp_x2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[1], Txz_2);
    }

    if(boundary_layer_number[2] > 0) {
        ierr = abs_exp_z1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[2], Txx_1);
        ierr = abs_exp_z1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[2], Txx_2);
        ierr = abs_exp_z1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[2], Tzz_1);
        ierr = abs_exp_z1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[2], Tzz_2);
        ierr = abs_exp_z1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[2], Txz_1);
        ierr = abs_exp_z1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[2], Txz_2);
    }

    if(boundary_layer_number[3] > 0) {
        ierr = abs_exp_z2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[3], Txx_1);
        ierr = abs_exp_z2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[3], Txx_2);
        ierr = abs_exp_z2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[3], Tzz_1);
        ierr = abs_exp_z2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[3], Tzz_2);
        ierr = abs_exp_z2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[3], Txz_1);
        ierr = abs_exp_z2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[3], Txz_2);
    }


    return 0;
}

// TODO: four corner may be wrong!
int abs_exp_x1(int ni1, int ni2, int nk1, int nk2, int nx, int nz, int half_fd_stencil,
               int boundary_layer_number, float *U)
{
    int indx, ix, iz;
    int bni1, bni2, bnk1, bnk2;
    float damp;

    bni1 = ni1;
    bni2 = ni1 - 1 + boundary_layer_number;
    bnk1 = nk1;
    bnk2 = nk2;

    for (ix = bni1; ix < bni2; ix++) {
        for (iz = bnk1; iz < bnk2; iz++) {
            damp = cal_exp( boundary_layer_number-(ix-bni1) );
            U[iz*nx+ix] = U[iz*nx+ix] * damp;
        }
    }

    return 0;
}

int abs_exp_x2(int ni1, int ni2, int nk1, int nk2, int nx, int nz, int half_fd_stencil,
               int boundary_layer_number, float *U)
{
    int indx, ix, iz;
    int bni1, bni2, bnk1, bnk2;
    float damp;


    bni1 = ni2 - boundary_layer_number + 1;
    bni2 = ni2;
    bnk1 = nk1;
    bnk2 = nk2;

    for (ix = bni1; ix < bni2; ix++) {
        for (iz = bnk1; iz < bnk2; iz++) {
            damp = cal_exp( ix-(bni2-boundary_layer_number)+1 );
            U[iz*nx+ix] = U[iz*nx+ix] * damp;
        }
    }

    return 0;
}

int abs_exp_z1(int ni1, int ni2, int nk1, int nk2, int nx, int nz, int half_fd_stencil,
               int boundary_layer_number, float *U)
{
    int indx, ix, iz;
    int bni1, bni2, bnk1, bnk2;
    float damp;

    bni1 = ni1;
    bni2 = ni2;
    bnk1 = nk1;
    bnk2 = nk1 - 1 + boundary_layer_number;

    for (ix = bni1; ix < bni2; ix++) {
        for (iz = bnk1; iz < bnk2; iz++) {
            damp = cal_exp( boundary_layer_number-(iz-bnk1) );
            U[iz*nx+ix] = U[iz*nx+ix] * damp;
        }
    }

    return 0;
}

int abs_exp_z2(int ni1, int ni2, int nk1, int nk2, int nx, int nz, int half_fd_stencil,
               int boundary_layer_number, float *U)
{
    int indx, ix, iz;
    int bni1, bni2, bnk1, bnk2;
    float damp;

    bni1 = ni1;
    bni2 = ni2;
    bnk1 = nk2- boundary_layer_number + 1;
    bnk2 = nk2;

    for (ix = bni1; ix < bni2; ix++) {
        for (iz = bnk1; iz < bnk2; iz++) {
            damp = cal_exp( iz-(bnk2-boundary_layer_number)+1 );
            U[iz*nx+ix] = U[iz*nx+ix] * damp;
        }
    }

    return 0;
}

float cal_exp(int i)
{
    float e;

    e =  exp( -(ALPHA*i) * (ALPHA*i) );

    return e;
}

int abs_exp_velocity_ssg(int ni1, int ni2, int nk1, int nk2, int nx, int nz, int half_fd_stencil,
                         int *boundary_layer_number,
                         float *Vx, float *Vz)
{
    int ierr = 0;

    if(boundary_layer_number[0] > 0) {
        ierr = abs_exp_x1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[0], Vx);
        ierr = abs_exp_x1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[0], Vz);;
    }

    if(boundary_layer_number[1] > 0) {
        ierr = abs_exp_x2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[1], Vx);
        ierr = abs_exp_x2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[1], Vz);
    }

    if(boundary_layer_number[2] >0 ) {
        ierr = abs_exp_z1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[2], Vx);
        ierr = abs_exp_z1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[2], Vz);
    }

    if(boundary_layer_number[3] > 0) {
        ierr = abs_exp_z2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[3], Vx);
        ierr = abs_exp_z2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[3], Vz);
    }
    return 0;
}

int abs_exp_stresses_ssg(int ni1, int ni2, int nk1, int nk2, int nx, int nz, int half_fd_stencil,
                     int *boundary_layer_number,
                     float *Txx, float *Tzz, float *Txz)
{
    int ierr = 0;

    if(boundary_layer_number[0] > 0) {
        ierr = abs_exp_x1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[0], Txx);
        ierr = abs_exp_x1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[0], Tzz);
        ierr = abs_exp_x1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[0], Txz);
    }

    if(boundary_layer_number[1] > 0) {
        ierr = abs_exp_x2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[1], Txx);
        ierr = abs_exp_x2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[1], Tzz);
        ierr = abs_exp_x2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[1], Txz);
    }

    if(boundary_layer_number[2] > 0) {
        ierr = abs_exp_z1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[2], Txx);
        ierr = abs_exp_z1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[2], Tzz);
        ierr = abs_exp_z1(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[2], Txz);
    }

    if(boundary_layer_number[3] > 0) {
        ierr = abs_exp_z2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[3], Txx);
        ierr = abs_exp_z2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[3], Tzz);
        ierr = abs_exp_z2(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                          boundary_layer_number[3], Txz);
    }

    return 0;
}
