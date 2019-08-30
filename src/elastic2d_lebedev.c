/***************************************************************************
 *
 * This function is used to simulate 2D elastic waveform in the TTI media
 *   using the high-order Lebedev scheme
 *
 * Authors: Luqian Jiang   Email: jlq0608@gmail.com
 *          Wei Zhang      Email: zhangwei@sustc.edu.cn
 *
 * Copyright (c) 2018 - 2019 zwlab
 *
 * History: 08/2019: Original version created by Luqian Jiang
 *
 ***************************************************************************/

#include "elastic2d_lebedev.h"
#include "elastic2d_stf.c"
#include "elastic2d_coef.c"
#include "elastic2d_filter.h"


int elastic2d_lebedev(float dx, float dz, int nx, int nz, int nt, int half_fd_stencil,
    float *c11_1, float *c13_1, float *c15_1, float *c33_1, float *c35_1, float *c55_1,
    float *c11_2, float *c13_2, float *c15_2, float *c33_2, float *c35_2, float *c55_2,
    float *B01, float *B10,
    struct Src, struct Rec,
    int *boundary_type, int *boundary_layer_number)
{
    int ni1, ni2, nk1, nk2, it, i, ierr;
    float *fdx = NULL, *fdz = NULL;

    /*---------------------------------------------------------*/
    /* staggered grid differential coefficient                 */
    /*---------------------------------------------------------*/
    fdx = (int*) malloc(half_fd_stencil*sizeof(int));
    fdz = (int*) malloc(half_fd_stencil*sizeof(int));

    for (i = 0; i < half_fd_stencil; i++) {
        fdx[i] = coef(half_fd_stencil,i+1)/dx;
        fdz[i] = coef(half_fd_stencil,i+1)/dz;
    }

    /* TODO! check stability */

    /* Initial */

    /*==========================================================*/
    /* Time loop to simulate wavefield propagation              */
    /*==========================================================*/
    for (it = 0; it < nt; it++) {
        fprintf(stdout, "\n it = %d \n", it);

        /*------------------------------------------------------*/
        /* update velocity (equation of the motion)             */
        /*------------------------------------------------------*/

        /* Free surface */

        /* Moment source */
        ierr = src_tti_moment();

        /* RHS of equation of the motion */
        ierr = cal_momentum(fdx, fdz, half_fd_stencil, ni1, ni2, nk1, nk2, nx,
                            hVx_1, hVz_1, Txx_1, Txz_1, Tzz_1,
                            hVx_2, hVz_2, Txx_2, Txz_2, Tzz_2);

        /* Absorbing boundary condition */
        ierr = abs_exp_velocity(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                            boundary_type, boundary_layer_number,
                            Vx_1, Vx_2, Vz_1, Vz_2);

        /* moment equation to update velocity */
        ierr = update_velocity(nx, nz, Vx_1, Vz_1, Vx_2, Vz_2,
                            hVx_1, hVz_1, hVx_2, hVz_2, B10, B01);

        /* Filtering velocity */
        ierr = filter_velocity(half_fd_stencil, nx, nz,
                            Vx_1, Vx_2, Vz_1, Vz_2);

        /* Receiver and snapshot */


        /*------------------------------------------------------*/
        /* update stress (Hoke law)                             */
        /*------------------------------------------------------*/

        /* Free surface */

        /* Force source */
        ierr = src_tti_force();

        /* RHS of Hook law euqation (stress-strain equation) */
        ierr = cal_hook(fdx, fdz, half_fd_stencil, ni1, ni2, nk1, nk2, nx,
                        hTxx_1, hTzz_1, hTxz_1, Vx_1 , Vz_1,
                        hTxx_2, hTzz_2, hTxz_2, Vx_2 , Vz_2,
                        c11_1 , c13_1 , c15_1 , c33_1, c35_1, c55_1,
                        c11_2 , c13_2 , c15_2 , c33_2, c35_2, c55_2 );

        /* update stress */
        ierr = update_stress(nx, nz, Txx_1, Tzz_1, Txz_1, Txx_2, Tzz_2, Txz_2,
                        hTxx_1, hTzz_1, hTxz_1, hTxx_2, hTzz_2, hTxz_2 );

        /* Absorbing boundary condition */
        ierr = abs_exp_stresses(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                        boundary_type,  boundary_layer_number,
                        Txx_1, Txx_2, Tzz_1, Tzz_2, Txz_1, Txz_2);

        /* Filtering stresses */
        ierr = filter_stresses(half_fd_stencil, nx, nz,
                        Txx_1, Tzz_1, *Txz_1,
                        Txx_2, Tzz_2, *Txz_2);

        /* receiver and snapshot */
        ierr


    }



}

/* RHS of equation of motion */
int cal_momentum(float *fdx, float *fdz, int half_fd_stencil, int ni1, int ni2, int nk1, int nk2, int nx,
                 float *hVx_1, float *hVz_1, float *Txx_1, float *Txz_1, float *Tzz_1,
                 float *hVx_2, float *hVz_2, float *Txx_2, float *Txz_2, float *Tzz_2)
{
    int ix, iz, i;
    float Txx_x, Txz_z, Txz_x, Tzz_z;

    for (iz = nk1; iz < nk2; iz++) {
        for (ix = ni1; ix < ni2; ix++) {

            /*-- for the 1st set of grids --*/
            Txx_x = 0.0; Txz_z = 0.0; Txz_x = 0.0; Tzz_z = 0.0;
            for (i = 1; i <= half_fd_stencil; i++) {
                Txx_x = Txx_x + Dx01(Txx_1,fdx,i,ix,iz);
                Txz_z = Txz_z + Dz01(Txz_1,fdz,i,ix,iz);
                Txz_x = Txz_x + Dx10(Txz_1,fdx,i,ix,iz);
                Tzz_z = Tzz_z + Dz10(Tzz_1,fdz,i,ix,iz);
            }

            hVx_1[iz*nx+ix] = Txx_x + Txz_z;
            hVz_1[iz*nx+iz] = Txz_x + Tzz_z;


            /*-- for the 2nd set of grids --*/
            Txx_x = 0.0; Txz_z = 0.0; Txz_x = 0.0; Tzz_z = 0.0;
            for (i = 1; i <= half_fd_stencil; i++) {
                Txx_x = Txx_x + Dx10(Txx_2,fdx,i,ix,iz);
                Txz_z = Txz_z + Dz10(Txz_2,fdz,i,ix,iz);
                Txz_x = Txz_x + Dx01(Txz_2,fdx,i,ix,iz);
                Tzz_z = Tzz_z + Dz01(Tzz_2,fdz,i,ix,iz);
            }

            hVx_2[iz*nx+ix] = Txx_x + Txz_z;
            hVz_2[iz*nx+iz] = Txz_x + Tzz_z;


        }
    }

    return 0;
}


/* RHS of stress-strain equation (Hook law) */
int cal_hook(float *fdx, float *fdz, int half_fd_stencil, int ni1, int ni2, int nk1, int nk2, int nx,
             float *hTxx_1, float *hTzz_1, float *hTxz_1, float *Vx_1 , float *Vz_1,
             float *hTxx_2, float *hTzz_2, float *hTxz_2, float *Vx_2 , float *Vz_2,
             float *c11_1 , float *c13_1 , float *c15_1 , float *c33_1, float *c35_1, float *c55_1,
             float *c11_2 , float *c13_2 , float *c15_2 , float *c33_2, float *c35_2, float *c55_2)
{
    int ix, iz, i, index;
    float Vx_x00, Vz_z00, Vx_z00, Vz_x00; /* partial difference centered of (i,j) */
    float Vx_x11, Vz_z11, Vx_z11, Vz_x11; /* partial difference centered of (i+1/2,j+1/2) */

    for (iz = nk1; iz < nk2; iz++) {
        for (ix = ni1; ix < ni2; ix++) {

            Vx_x00 = 0.0; Vz_z00 = 0.0; Vx_z00 = 0.0; Vz_x00 = 0.0;
            Vx_x11 = 0.0; Vz_z11 = 0.0; Vx_z11 = 0.0; Vz_x11 = 0.0;

            /*-- for the 1st set of grids --*/
            for (i = 1; i <= half_fd_stencil; i++) {

                Vx_x00 = Vx_x00 + Dx00(Vx_1,fdx,i,ix,iz);
                Vz_z00 = Vz_z00 + Dz00(Vz_1,fdz,i,ix,iz);
                Vx_z11 = Vx_z11 + Dz11(Vx_1,fdz,i,ix,iz);
                Vz_x11 = Vz_x11 + Dx11(Vz_1,fdx,i,ix,iz);

                Vx_x11 = Vx_x11 + Dx11(Vx_2,fdx,i,ix,iz);
                Vz_z11 = Vz_z11 + Dz11(Vz_2,fdz,i,ix,iz);
                Vx_z00 = Vx_z00 + Dz00(Vx_2,fdz,i,ix,iz);
                Vz_x00 = Vz_x00 + Dx00(Vz_2,fdx,i,ix,iz);

            }

            index = iz*nx + ix;

            hTxx_1[index] = c11_1[index] * Vx_x00 + c13_1[index] * Vz_z00 \
                          + c15_1[index] * Vx_z00 + c15_1[index] * Vz_x00;

            hTzz_1[index] = c13_1[index] * Vx_x00 + c33_1[index] * Vz_z00 \
                          + c35_1[index] * Vx_z00 + c35_1[index] * Vz_x00;

            hTxz_1[index] = c15_2[index] * Vx_x11 + c15_2[index] * Vz_z11 \
                          + c55_1[index] * Vx_z11 + c55_1[index] * Vz_x11;

            hTxx_2[index] = c11_2[index] * Vx_x11 + c13_2[index] * Vz_z11 \
                          + c15_2[index] * Vx_z11 + c15_2[index] * Vz_x11;

            hTzz_2[index] = c13_2[index] * Vx_x11 + c33_2[index] * Vz_z11 \
                          + c35_2[index] * Vx_z11 + c35_2[index] * Vz_x11;

            hTxz_2[index] = c15_1[index] * Vx_x00 + c15_1[index] * Vz_z00 \
                          + c55_2[index] * Vx_z00 + c55_2[index] * Vz_x00;

        }
    }

    return 0;
}

int update_velocity(int nx, int nz, float *Vx_1, float *Vz_1, float *Vx_2, float *Vz_2,
        float *hVx_1, float *hVz_1,  float *hVx_2, float *hVz_2, float *B10, float *B01)
{
    int ix, iz, index;
    for (iz = 0; iz < nz; iz++) {
        for (ix = 0; ix < nx; ix++) {
            index = iz*nx + ix;

            Vx_1[index] = Vx_1[index] + hVx_1[index] * B01[index];
            Vz_1[index] = Vz_1[index] + hVz_1[index] * B10[index];

            Vx_2[index] = Vx_2[index] + hVx_2[index] * B10[index];
            Vz_2[index] = Vz_2[index] + hVz_2[index] * B01[index];
        }
    }

    return 0;
}

int update_stress(int nx, int nz, float *Txx_1, float *Tzz_1, float *Txz_1,
                  float *Txx_2 , float *Tzz_2 , float *Txz_2 ,
                  float *hTxx_1, float *hTzz_1, float *hTxz_1,
                  float *hTxx_2, float *hTzz_2, float *hTxz_2 )
{
    int ix, iz, index;
    for (iz = 0; iz < nz; iz++) {
        for (ix = 0; ix < nx; ix++) {
            index = iz*nx + ix;

            Txx_1[index] = Txx_1[index] + hTxx_1[index];
            Tzz_1[index] = Tzz_1[index] + hTzz_1[index];
            Txz_1[index] = Txz_1[index] + hTxz_1[index];

            Txx_2[index] = Txx_2[index] + hTxx_2[index];
            Tzz_2[index] = Tzz_2[index] + hTzz_2[index];
            Txz_2[index] = Txz_2[index] + hTxz_2[index];

        }
    }

    return 0;
}
